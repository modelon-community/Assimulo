/*
 * TR-BDF2 stepper (see trbdf2.h for the reference and the license).
 *
 * One step from (t, y) with step h, gamma = 2 - sqrt(2), d = gamma/2, w = sqrt(2)/4:
 *   TR stage    y1 = y + d h (f(t, y) + f(t + gamma h, y1))
 *   BDF2 stage  y2 = a1 y1 - a0 y + d h f(t + h, y2),  a1 = 1/(gamma(2-gamma)), a0 = (1-gamma)^2/(gamma(2-gamma))
 * Both stages solve with the same iteration matrix M = I - d h J (simplified Newton), so one LU
 * per step size. Butcher form: c = (0, gamma, 1), b = (w, w, d); embedded third-order weights
 * bhat = ((1-w)/3, (3w+1)/3, d/3); local error e = h sum (b_i - bhat_i) f_i, filtered with M^-1
 * (Hosea & Shampine's stiff modification) before the norm is taken.
 *
 * The stage derivatives f_i are taken from the stage equations, f_i = (y_i - c_i) / (d h), not from
 * an rhs call at the converged iterate (Hosea & Shampine, section 5): the Newton iterate is within
 * newton_tol of the stage value, and for a stiff component an rhs call there is off by h*lambda
 * times that, which pollutes the FSAL derivative, the TR stage of the next step and the dense
 * output. The algebraic value is off by (1/d) times the Newton error only. This also saves two
 * rhs calls per step. Predictors: stage 1 extrapolates the previous step's dense output, stage 2
 * the quadratic through y0, y1 with slope f1 (both measured on the reference FMUs; explicit Euler
 * for stage 1 fails on the data-center model).
 *
 * Jacobian policy as in CVode: the Jacobian is recomputed only on the first step, after a Newton
 * failure with a stale Jacobian, or after max_steps_between_jac accepted steps; the LU is redone
 * whenever h changes (h is kept when the proposed change is within [keep_lo, keep_hi]).
 * Recoverable callback failures (a model refusing a trial point) shrink h by step_factor and are
 * retried up to max_consecutive times, as CVode does.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "trbdf2.h"

static int trace_on(void) { static int v = -1; if (v < 0) { const char *e = getenv("TRBDF2_TRACE"); v = e ? atoi(e) : 0; } return v; }

#define SQRT2 1.4142135623730950488
#define NEWTON_DIV_RATE 2.0   /* Newton is declared diverging when |delta_k| >= NEWTON_DIV_RATE |delta_{k-1}| (CVode's RDIV) */
static const double GAMMA_ = 2.0 - SQRT2;         /* 0.5857864376269049 */
static const double D_     = 1.0 - SQRT2 / 2.0;    /* gamma/2 = 0.2928932188134524 */
static const double W_     = SQRT2 / 4.0;          /* 0.3535533905932738 */

struct trbdf2_mem {
    int n;
    double *rtol, *atol;                 /* n */
    double *y0, *f0, *y1, *f1, *y2, *f2; /* stage values of the current step; f0 is FSAL */
    double *fd0, *yd1, *fd1;             /* dense output: derivative at the start, state and derivative at the end of the
                                            last accepted step (y0 holds its start; y2/f2 may belong to a rejected attempt) */
    double *fbase;                       /* rhs at the base point of a finite-difference Jacobian */
    double *ytmp, *ftmp, *delta, *err, *werr, *scal;
    double *jac, *lu;                    /* n*n column-major */
    int *piv;
    double t0, h;                        /* start of the last accepted step and its size */
    double h_next;                       /* step size proposed for the next step */
    double h_lu;                         /* step size the LU was formed with (0: no LU) */
    double err_prev;                     /* error norm of the previous accepted step (PI controller), 0: none */
    int first_step;                      /* no step accepted yet in this segment: growth cap 1e4 (CVode's ETAMX1) */
    int have_lu, have_jac, have_f0;
    long steps_since_jac;
    /* options */
    double hmax, safety, fac_min, fac_max, keep_lo, keep_hi, newton_tol, fail_factor;
    int newton_max, max_between_jac, use_user_jac, fail_max;
    long max_steps;
    trbdf2_stats stats;
    char err_msg[256];
};

/* ------------------------------------------------------------------ dense LU (column-major) */
static int dec(int n, double *a, int *piv) {
    /* LU with partial pivoting, a[i + j*n]; returns 0 or the (1-based) index of a zero pivot */
    int i, j, k, p;
    for (k = 0; k < n; k++) {
        p = k;
        for (i = k + 1; i < n; i++) if (fabs(a[i + k * n]) > fabs(a[p + k * n])) p = i;
        piv[k] = p;
        if (a[p + k * n] == 0.0) return k + 1;
        if (p != k) for (j = 0; j < n; j++) { double s = a[k + j * n]; a[k + j * n] = a[p + j * n]; a[p + j * n] = s; }
        for (i = k + 1; i < n; i++) a[i + k * n] /= a[k + k * n];
        for (j = k + 1; j < n; j++) {
            double akj = a[k + j * n];
            if (akj != 0.0) for (i = k + 1; i < n; i++) a[i + j * n] -= a[i + k * n] * akj;
        }
    }
    return 0;
}

static void sol(int n, const double *a, const int *piv, double *b) {
    /* dec swaps whole rows (L included), so all interchanges are applied to b first, as LAPACK's dgetrs does */
    int i, k;
    for (k = 0; k < n; k++) if (piv[k] != k) { double s = b[k]; b[k] = b[piv[k]]; b[piv[k]] = s; }
    for (k = 0; k < n; k++) for (i = k + 1; i < n; i++) b[i] -= a[i + k * n] * b[k];
    for (k = n - 1; k >= 0; k--) {
        b[k] /= a[k + k * n];
        for (i = 0; i < k; i++) b[i] -= a[i + k * n] * b[k];
    }
}

/* ------------------------------------------------------------------ helpers */
static void set_scale(trbdf2_mem *m, const double *y) {
    int i;
    for (i = 0; i < m->n; i++) m->scal[i] = m->rtol[i] * fabs(y[i]) + m->atol[i];
}

static double wrms(trbdf2_mem *m, const double *v) {
    int i; double s = 0.0, q;
    for (i = 0; i < m->n; i++) { q = v[i] / m->scal[i]; s += q * q; }
    return sqrt(s / m->n);
}

static void fail(trbdf2_mem *m, const char *msg, double t) {
    snprintf(m->err_msg, sizeof(m->err_msg), "%s (at t = %g)", msg, t);
}

/* ------------------------------------------------------------------ API: memory and options */
int trbdf2_create(int n, trbdf2_mem **mem_out) {
    trbdf2_mem *m;
    if (n < 1 || mem_out == NULL) return TRBDF2_ERROR_INPUT;
    m = (trbdf2_mem *)calloc(1, sizeof(trbdf2_mem));
    if (m == NULL) return TRBDF2_ERROR_MEM;
    m->n = n;
    m->rtol = (double *)malloc(n * sizeof(double)); m->atol = (double *)malloc(n * sizeof(double));
    m->y0 = (double *)malloc(n * sizeof(double)); m->f0 = (double *)malloc(n * sizeof(double));
    m->y1 = (double *)malloc(n * sizeof(double)); m->f1 = (double *)malloc(n * sizeof(double));
    m->y2 = (double *)malloc(n * sizeof(double)); m->f2 = (double *)malloc(n * sizeof(double));
    m->fd0 = (double *)malloc(n * sizeof(double)); m->yd1 = (double *)malloc(n * sizeof(double)); m->fd1 = (double *)malloc(n * sizeof(double));
    m->fbase = (double *)malloc(n * sizeof(double));
    m->ytmp = (double *)malloc(n * sizeof(double)); m->ftmp = (double *)malloc(n * sizeof(double));
    m->delta = (double *)malloc(n * sizeof(double)); m->err = (double *)malloc(n * sizeof(double));
    m->werr = (double *)malloc(n * sizeof(double)); m->scal = (double *)malloc(n * sizeof(double));
    m->jac = (double *)malloc((size_t)n * n * sizeof(double)); m->lu = (double *)malloc((size_t)n * n * sizeof(double));
    m->piv = (int *)malloc(n * sizeof(int));
    if (!m->rtol || !m->atol || !m->y0 || !m->f0 || !m->y1 || !m->f1 || !m->y2 || !m->f2 || !m->fd0 || !m->yd1 || !m->fd1 || !m->fbase || !m->ytmp || !m->ftmp ||
        !m->delta || !m->err || !m->werr || !m->scal || !m->jac || !m->lu || !m->piv) {
        trbdf2_free(&m); return TRBDF2_ERROR_MEM;
    }
    { int i; for (i = 0; i < n; i++) { m->rtol[i] = 1e-6; m->atol[i] = 1e-8; } }
    m->hmax = 0.0; m->safety = 0.9; m->fac_min = 0.2; m->fac_max = 5.0; m->keep_lo = 0.8; m->keep_hi = 1.25;
    m->newton_tol = 0.1; m->newton_max = 6; m->max_between_jac = 50; m->use_user_jac = 1;
    m->fail_factor = 0.25; m->fail_max = 40; m->max_steps = 100000;
    *mem_out = m;
    return trbdf2_reinit(m);
}

void trbdf2_free(trbdf2_mem **mem) {
    trbdf2_mem *m;
    if (mem == NULL || *mem == NULL) return;
    m = *mem;
    free(m->rtol); free(m->atol); free(m->y0); free(m->f0); free(m->y1); free(m->f1); free(m->y2); free(m->f2); free(m->fd0); free(m->yd1); free(m->fd1); free(m->fbase);
    free(m->ytmp); free(m->ftmp); free(m->delta); free(m->err); free(m->werr); free(m->scal);
    free(m->jac); free(m->lu); free(m->piv);
    free(m); *mem = NULL;
}

int trbdf2_set_tolerances(trbdf2_mem *m, const double *rtol, const double *atol) {
    int i;
    if (!m) return TRBDF2_ERROR_MEM;
    for (i = 0; i < m->n; i++) {
        if (rtol[i] <= 0.0 || atol[i] < 0.0) return TRBDF2_ERROR_INPUT;
        m->rtol[i] = rtol[i]; m->atol[i] = atol[i];
    }
    return TRBDF2_OK;
}
int trbdf2_set_hmax(trbdf2_mem *m, double hmax) { if (!m) return TRBDF2_ERROR_MEM; m->hmax = hmax > 0.0 ? hmax : 0.0; return TRBDF2_OK; }
int trbdf2_set_max_steps(trbdf2_mem *m, long nmax) { if (!m || nmax < 1) return TRBDF2_ERROR_INPUT; m->max_steps = nmax; return TRBDF2_OK; }
int trbdf2_set_newton(trbdf2_mem *m, int max_iter, double tol) {
    if (!m || max_iter < 1 || tol <= 0.0) return TRBDF2_ERROR_INPUT;
    m->newton_max = max_iter; m->newton_tol = tol; return TRBDF2_OK;
}
int trbdf2_set_step_control(trbdf2_mem *m, double safety, double fac_min, double fac_max, double keep_lo, double keep_hi) {
    if (!m || safety <= 0.0 || safety > 1.0 || fac_min <= 0.0 || fac_min >= 1.0 || fac_max <= 1.0 || keep_lo > 1.0 || keep_lo <= 0.0 || keep_hi < 1.0)
        return TRBDF2_ERROR_INPUT;
    m->safety = safety; m->fac_min = fac_min; m->fac_max = fac_max; m->keep_lo = keep_lo; m->keep_hi = keep_hi;
    return TRBDF2_OK;
}
int trbdf2_set_jac_policy(trbdf2_mem *m, int max_steps_between_jac, int use_user_jac) {
    if (!m || max_steps_between_jac < 1) return TRBDF2_ERROR_INPUT;
    m->max_between_jac = max_steps_between_jac; m->use_user_jac = use_user_jac ? 1 : 0; return TRBDF2_OK;
}
int trbdf2_set_failure_policy(trbdf2_mem *m, double step_factor, int max_consecutive) {
    if (!m || step_factor <= 0.0 || step_factor >= 1.0 || max_consecutive < 1) return TRBDF2_ERROR_INPUT;
    m->fail_factor = step_factor; m->fail_max = max_consecutive; return TRBDF2_OK;
}
int trbdf2_reinit(trbdf2_mem *m) {
    if (!m) return TRBDF2_ERROR_MEM;
    m->h = 0.0; m->h_next = 0.0; m->h_lu = 0.0; m->have_lu = 0; m->have_jac = 0; m->have_f0 = 0; m->err_prev = 0.0;
    m->steps_since_jac = 0; m->first_step = 1; m->err_msg[0] = '\0';
    memset(&m->stats, 0, sizeof(m->stats));
    return TRBDF2_OK;
}
int trbdf2_reset_stats(trbdf2_mem *m) { if (!m) return TRBDF2_ERROR_MEM; memset(&m->stats, 0, sizeof(m->stats)); return TRBDF2_OK; }
double trbdf2_get_h(trbdf2_mem *m) { return m ? m->h_next : 0.0; }
void trbdf2_get_stats(trbdf2_mem *m, trbdf2_stats *s) { if (m && s) *s = m->stats; }
const char *trbdf2_get_err_msg(trbdf2_mem *m) { return m ? m->err_msg : "no memory"; }

/* ------------------------------------------------------------------ Jacobian */
static int form_jacobian(trbdf2_mem *m, trbdf2_rhs_fn rhs, trbdf2_jac_fn jac, void *user, double t, const double *y) {
    int n = m->n, i, j, ret;
    const double *f = m->fbase;
    if (m->use_user_jac && jac != NULL) {
        ret = jac(n, t, y, m->jac, user);
        if (ret != 0) return ret;
    } else {
        /* forward differences, one column per component, from an rhs call at the base point
           (f0 is the algebraic stage derivative, not accurate enough for a difference quotient) */
        ret = rhs(n, t, y, m->fbase, user);
        m->stats.nfcnjac++;
        if (ret != 0) return ret;
        for (j = 0; j < n; j++) {
            double dy = sqrt(2.2e-16) * fmax(fabs(y[j]), m->atol[j] / m->rtol[j] > 0.0 ? m->atol[j] / m->rtol[j] : 1e-5);
            if (dy == 0.0) dy = 1e-8;
            memcpy(m->ytmp, y, n * sizeof(double));
            m->ytmp[j] += dy;
            ret = rhs(n, t, m->ytmp, m->ftmp, user);
            m->stats.nfcnjac++;
            if (ret != 0) return ret;
            for (i = 0; i < n; i++) m->jac[i + j * n] = (m->ftmp[i] - f[i]) / dy;
        }
    }
    m->stats.njac++;
    m->have_jac = 1; m->steps_since_jac = 0; m->have_lu = 0;
    return 0;
}

static int form_lu(trbdf2_mem *m, double h) {
    int n = m->n, i, j;
    double dh = D_ * h;
    for (j = 0; j < n; j++) for (i = 0; i < n; i++) m->lu[i + j * n] = -dh * m->jac[i + j * n];
    for (i = 0; i < n; i++) m->lu[i + i * n] += 1.0;
    m->stats.nlu++;
    if (dec(n, m->lu, m->piv) != 0) { m->have_lu = 0; return TRBDF2_ERROR_SINGULAR; }
    m->have_lu = 1; m->h_lu = h;
    return 0;
}

/* ------------------------------------------------------------------ Newton for one stage
 * Solves  y - c_y - dh f(ts, y) = rhs_const  with the fixed matrix M = I - dh J.
 * y holds the predictor on entry and the solution on exit; fy the stage derivative (y - cvec) / dh.
 * Returns 0 converged, 1 not converged (diverged / too many iterations), > 1 recoverable rhs failure,
 * < 0 unrecoverable. */
static int newton_stage(trbdf2_mem *m, trbdf2_rhs_fn rhs, void *user, double ts, double dh,
                        const double *cvec, double *y, double *fy, int stage) {
    int n = m->n, i, k, ret;
    double dnorm, dnorm_old = 0.0, rate = 0.0;
    for (k = 0; k < m->newton_max; k++) {
        ret = rhs(n, ts, y, fy, user);
        m->stats.nfcn++;
        if (ret != 0) return ret > 0 ? 2 : ret;
        /* residual r = y - cvec - dh*fy ; delta = -M^-1 r */
        for (i = 0; i < n; i++) m->delta[i] = -(y[i] - cvec[i] - dh * fy[i]);
        sol(n, m->lu, m->piv, m->delta);
        m->stats.nsolve++;
        m->stats.nnewton++;
        for (i = 0; i < n; i++) y[i] += m->delta[i];
        dnorm = wrms(m, m->delta);
        if (k > 0) rate = dnorm / dnorm_old;
        if (trace_on() > 1) fprintf(stderr, "      stage %d newton k=%d |delta|=%.3e rate=%.3f\n", stage, k, dnorm, rate);
        /* converged when the correction is small, or when the estimated remaining error
           |delta| * rate/(1-rate) is; the divergence test only applies to corrections that are
           not small already (a converged iterate's next correction is roundoff with a random ratio) */
        if (dnorm <= m->newton_tol || (k > 0 && rate < 1.0 && dnorm * rate / (1.0 - rate) <= m->newton_tol)) {
            for (i = 0; i < n; i++) fy[i] = (y[i] - cvec[i]) / dh;    /* stage derivative from the stage equation (see the file comment) */
            return 0;
        }
        if (k > 0 && rate >= NEWTON_DIV_RATE) { m->stats.nnfail_div++; return 1; }      /* diverging */
        dnorm_old = dnorm;
    }
    m->stats.nnfail_iter++;
    return 1;
}

/* ------------------------------------------------------------------ initial step (Hairer's hinit) */
static int initial_step(trbdf2_mem *m, trbdf2_rhs_fn rhs, void *user, double t, const double *y, double tend, double *h0) {
    int n = m->n, i, ret;
    double dnf, dny, h, der2, h1, dir = tend > t ? 1.0 : -1.0, span = fabs(tend - t);
    set_scale(m, y);
    dnf = wrms(m, m->f0); dny = wrms(m, y);
    h = (dnf <= 1e-10 || dny <= 1e-10) ? 1e-6 : 0.01 * sqrt(dny / dnf);   /* Hairer's hinit */
    h = fmin(h, span);
    if (m->hmax > 0.0) h = fmin(h, m->hmax);
    for (i = 0; i < n; i++) m->ytmp[i] = y[i] + dir * h * m->f0[i];
    ret = rhs(n, t + dir * h, m->ytmp, m->ftmp, user);
    m->stats.nfcn++;
    if (ret != 0) { *h0 = h * 0.1; return ret > 0 ? 0 : ret; }
    for (i = 0; i < n; i++) m->ftmp[i] = (m->ftmp[i] - m->f0[i]) / h;
    der2 = wrms(m, m->ftmp);
    h1 = (fmax(der2, dnf) <= 1e-15) ? fmax(1e-6, h * 1e-3) : pow(0.01 / fmax(der2, dnf), 1.0 / 3.0);
    h = fmin(100.0 * h, h1);
    h = fmax(h, 1e-8 * span);                  /* a too-small guess is corrected by the first step's growth (see fac_first) */
    h = fmin(h, span);
    if (m->hmax > 0.0) h = fmin(h, m->hmax);
    *h0 = h;
    return 0;
}

/* ------------------------------------------------------------------ dense output */
int trbdf2_interpolate(trbdf2_mem *m, double t, double *y_out) {
    int n, i; double s, a, b, c, dd, h;
    if (!m) return TRBDF2_ERROR_MEM;
    if (m->h == 0.0) return TRBDF2_ERROR_INPUT;
    n = m->n; h = m->h;
    s = (t - m->t0) / h;
    /* cubic Hermite over the whole step through the two stiffly accurate end points of the last accepted step */
    a = (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s); b = s * (1.0 - s) * (1.0 - s); c = s * s * (3.0 - 2.0 * s); dd = s * s * (s - 1.0);
    for (i = 0; i < n; i++) y_out[i] = a * m->y0[i] + b * h * m->fd0[i] + c * m->yd1[i] + dd * h * m->fd1[i];
    return TRBDF2_OK;
}

/* ------------------------------------------------------------------ the integration loop */
int trbdf2_solve(trbdf2_mem *m, trbdf2_rhs_fn rhs, trbdf2_jac_fn jac, trbdf2_solout_fn solout, void *user,
                 double *t, double *y, double tend, double h0) {
    int n, i, ret, nfail_consec = 0, last, jac_fresh_this_attempt;
    double h, dir, tn, hnew, errnorm, fac, a1, a0;
    const double c1 = (4.0 * W_ - 1.0) / 3.0, c2 = -1.0 / 3.0, c3 = 2.0 * D_ / 3.0;   /* b - bhat */
    if (!m || !rhs || !t || !y) return TRBDF2_ERROR_INPUT;
    n = m->n;
    a1 = 1.0 / (GAMMA_ * (2.0 - GAMMA_));
    a0 = (1.0 - GAMMA_) * (1.0 - GAMMA_) / (GAMMA_ * (2.0 - GAMMA_));
    dir = tend >= *t ? 1.0 : -1.0;
    if (fabs(tend - *t) <= 1e-13 * fmax(fabs(tend), 1.0)) { *t = tend; return TRBDF2_OK; }   /* nothing to integrate */

    /* rhs at the start (FSAL after an accepted step) */
    if (!m->have_f0) {
        ret = rhs(n, *t, y, m->f0, user);
        m->stats.nfcn++;
        if (ret != 0) { fail(m, "rhs failed at the start of the segment", *t); return ret > 0 ? TRBDF2_ERROR_REPEATED_FAILURE : TRBDF2_ERROR_CALLBACK; }
        m->have_f0 = 1;
    }
    h = m->h_next;
    if (h <= 0.0) {
        if (h0 > 0.0) h = h0;
        else { ret = initial_step(m, rhs, user, *t, y, tend, &h); if (ret < 0) { fail(m, "rhs failed during the initial step estimate", *t); return TRBDF2_ERROR_CALLBACK; } }
    }
    if (m->hmax > 0.0) h = fmin(h, m->hmax);
    if (trace_on()) fprintf(stderr, "solve: t=%.10g tend=%.10g h0=%.3e h_next=%.3e -> h=%.3e hmax=%.3e\n", *t, tend, h0, m->h_next, h, m->hmax);

    for (;;) {
        if (m->stats.nsteps >= m->max_steps) { fail(m, "maximum number of steps reached", *t); return TRBDF2_ERROR_MAX_STEPS; }
        /* land exactly on tend */
        last = 0;
        if (dir * (*t + dir * h - tend) >= -1e-10 * fmax(fabs(tend), 1.0)) { h = fabs(tend - *t); last = 1; }
        if (h < 1e-14 * fmax(fabs(*t), 1.0)) {
            if (last) { *t = tend; return TRBDF2_OK; }      /* the remaining segment is roundoff */
            fail(m, "step size too small", *t); return TRBDF2_ERROR_STEP_TOO_SMALL;
        }
        m->stats.nsteps++;
        set_scale(m, y);

        /* Jacobian / LU */
        jac_fresh_this_attempt = 0;
        if (!m->have_jac || m->steps_since_jac >= m->max_between_jac) {
            ret = form_jacobian(m, rhs, jac, user, *t, y);
            if (ret > 0) { m->stats.nrhsfail++; goto recoverable; }
            if (ret < 0) { fail(m, "Jacobian evaluation failed", *t); return TRBDF2_ERROR_CALLBACK; }
            jac_fresh_this_attempt = 1;
        }
        /* the LU is reused while h stays within [keep_lo, keep_hi] of the h it was formed with
           (an inexact iteration matrix only slows Newton down; CVode does the same with gamma) */
        if (!m->have_lu || dir * h / m->h_lu < m->keep_lo || dir * h / m->h_lu > m->keep_hi) {
            ret = form_lu(m, dir * h);
            if (ret != 0) {
                if (!jac_fresh_this_attempt) { m->have_jac = 0; continue; }   /* retry with a fresh Jacobian */
                fail(m, "iteration matrix is singular", *t); return TRBDF2_ERROR_SINGULAR;
            }
        }

        /* ---- stage 1 (TR): y1 - y - dh f0 - dh f(y1) = 0  ->  cvec = y + dh f0
           predictor: the previous step's dense output extrapolated (explicit Euler on the first step) */
        for (i = 0; i < n; i++) m->ytmp[i] = y[i] + D_ * dir * h * m->f0[i];
        if (m->h != 0.0 && m->t0 + m->h == *t) trbdf2_interpolate(m, *t + dir * GAMMA_ * h, m->y1);
        else for (i = 0; i < n; i++) m->y1[i] = y[i] + GAMMA_ * dir * h * m->f0[i];
        ret = newton_stage(m, rhs, user, *t + dir * GAMMA_ * h, D_ * dir * h, m->ytmp, m->y1, m->f1, 1);
        if (ret < 0) { fail(m, "rhs failed (unrecoverable) in stage 1", *t); return TRBDF2_ERROR_CALLBACK; }
        if (ret == 2) { m->stats.nrhsfail++; goto recoverable; }
        if (ret == 1) goto newton_failure;

        /* ---- stage 2 (BDF2): y2 - (a1 y1 - a0 y) - dh f(y2) = 0
           predictor: the quadratic through y (s = 0) and y1 (s = gamma) with slope f1 at y1, at s = 1 */
        for (i = 0; i < n; i++) m->ytmp[i] = a1 * m->y1[i] - a0 * y[i];
        { double sg = 1.0 - GAMMA_, q = sg * sg / (GAMMA_ * GAMMA_);
          for (i = 0; i < n; i++) m->y2[i] = m->y1[i] + sg * dir * h * m->f1[i] + q * (y[i] - m->y1[i] + GAMMA_ * dir * h * m->f1[i]); }
        ret = newton_stage(m, rhs, user, *t + dir * h, D_ * dir * h, m->ytmp, m->y2, m->f2, 2);
        if (ret < 0) { fail(m, "rhs failed (unrecoverable) in stage 2", *t); return TRBDF2_ERROR_CALLBACK; }
        if (ret == 2) { m->stats.nrhsfail++; goto recoverable; }
        if (ret == 1) goto newton_failure;

        /* ---- error estimate, filtered with M^-1 (Hosea & Shampine) */
        for (i = 0; i < n; i++) m->err[i] = dir * h * (c1 * m->f0[i] + c2 * m->f1[i] + c3 * m->f2[i]);
        sol(n, m->lu, m->piv, m->err);
        m->stats.nsolve++;
        /* weights from the new solution too, so that growing components are not under-weighted */
        for (i = 0; i < n; i++) m->scal[i] = m->rtol[i] * fmax(fabs(y[i]), fabs(m->y2[i])) + m->atol[i];
        errnorm = wrms(m, m->err);
        for (i = 0; i < n; i++) m->werr[i] = m->err[i] / m->scal[i];
        /* PI step control (Gustafsson): h_new = h safety err^-0.7/3 err_prev^0.4/3 after an accepted
           step, the plain err^-1/3 formula after a rejection or on the first step */
        if (errnorm <= 1.0 && m->err_prev > 0.0)
            fac = m->safety * pow(fmax(errnorm, 1e-10), -0.7 / 3.0) * pow(m->err_prev, 0.4 / 3.0);
        else
            fac = m->safety * pow(fmax(errnorm, 1e-10), -1.0 / 3.0);
        fac = fmin(m->first_step ? 1e4 : m->fac_max, fmax(m->fac_min, fac));

        if (trace_on()) fprintf(stderr, "  t=%.7f h=%.3e err=%.3g %s newton_iters_so_far=%ld\n", *t, h, errnorm, errnorm > 1.0 ? "REJ" : "acc", m->stats.nnewton);
        if (errnorm > 1.0) {
            /* ---- rejected */
            m->stats.nreject++;
            nfail_consec = 0;
            m->err_prev = 0.0;
            h *= fmin(fac, 0.9);
            continue;
        }

        /* ---- accepted */
        m->stats.naccpt++;
        nfail_consec = 0;
        m->first_step = 0;
        m->err_prev = fmax(errnorm, 1e-10);
        m->steps_since_jac++;
        m->t0 = *t; m->h = dir * h;
        memcpy(m->y0, y, n * sizeof(double));
        memcpy(m->fd0, m->f0, n * sizeof(double));
        memcpy(m->yd1, m->y2, n * sizeof(double));
        memcpy(m->fd1, m->f2, n * sizeof(double));
        tn = last ? tend : *t + dir * h;
        memcpy(y, m->y2, n * sizeof(double));
        memcpy(m->f0, m->f2, n * sizeof(double));            /* FSAL */
        hnew = h * fac;
        if (m->hmax > 0.0) hnew = fmin(hnew, m->hmax);
        m->h_next = hnew;
        *t = tn;
        if (solout != NULL) {
            ret = solout((int)m->stats.naccpt, m->t0, *t, y, m->werr, user);
            if (ret < 0) { fail(m, "solout callback failed", *t); return TRBDF2_ERROR_CALLBACK; }
            if (ret > 0) return TRBDF2_STOP;
        }
        if (last) return TRBDF2_OK;
        h = hnew;
        continue;

    newton_failure:
        m->stats.nnfail++;
        if (trace_on()) fprintf(stderr, "  t=%.7f h=%.3e NEWTON FAIL (jac fresh %d, steps_since_jac %ld, h_lu %.3e)\n", *t, h, jac_fresh_this_attempt, m->steps_since_jac, m->h_lu);
        if (!jac_fresh_this_attempt && m->steps_since_jac > 0) { m->have_jac = 0; continue; }   /* stale Jacobian: refresh, same h */
        h *= m->fail_factor;                                                                    /* fresh Jacobian: smaller step */
        m->have_lu = 0;
        m->stats.nreject++;
        if (h < 1e-14 * fmax(fabs(*t), 1.0)) { fail(m, "Newton iteration does not converge", *t); return TRBDF2_ERROR_NEWTON; }
        continue;

    recoverable:
        /* the model refused a trial point: shrink and retry, CVode-like, with a budget */
        if (++nfail_consec >= m->fail_max) { fail(m, "repeated recoverable failures of the model at trial points", *t); return TRBDF2_ERROR_REPEATED_FAILURE; }
        h *= m->fail_factor;
        m->stats.nreject++;
        continue;
    }
}
