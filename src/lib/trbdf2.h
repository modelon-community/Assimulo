/*
 * TR-BDF2: a second-order, L-stable, stiffly accurate ESDIRK method with a third-order
 * embedded error estimate and a piecewise cubic Hermite dense output.
 *
 * Hosea, M. E. and Shampine, L. F., "Analysis and implementation of TR-BDF2",
 * Applied Numerical Mathematics 20 (1996) 21-37.
 *
 * Copyright (C) 2026 Modelon AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef TRBDF2_H
#define TRBDF2_H

/* return codes of trbdf2_solve */
#define TRBDF2_OK                        0   /* reached tend */
#define TRBDF2_STOP                      1   /* solout asked to stop (event) */
#define TRBDF2_ERROR_MEM               -10
#define TRBDF2_ERROR_INPUT             -11
#define TRBDF2_ERROR_MAX_STEPS         -12
#define TRBDF2_ERROR_STEP_TOO_SMALL    -13
#define TRBDF2_ERROR_SINGULAR          -14
#define TRBDF2_ERROR_REPEATED_FAILURE  -15   /* recoverable callback failures exceeded the budget */
#define TRBDF2_ERROR_NEWTON            -16   /* Newton did not converge even with a fresh Jacobian and a tiny step */
#define TRBDF2_ERROR_CALLBACK          -20   /* unrecoverable error from a callback (rhs, jac, solout) */

/* callbacks: return 0 on success, > 0 for a recoverable failure (the trial point is refused),
   < 0 for an unrecoverable one. jac fills the column-major n x n matrix J[i + j*n] = df_i/dy_j.
   solout is called after every accepted step; 0 = continue, 1 = stop, < 0 = unrecoverable. */
typedef int (*trbdf2_rhs_fn)(int n, double t, const double *y, double *ydot, void *user);
typedef int (*trbdf2_jac_fn)(int n, double t, const double *y, double *jac, void *user);
typedef int (*trbdf2_solout_fn)(int naccpt, double told, double t, const double *y, const double *werr, void *user);

typedef struct trbdf2_stats {
    long nsteps, naccpt, nreject, nfcn, nfcnjac, njac, nlu, nsolve, nnewton, nnfail, nrhsfail;
    long nnfail_div, nnfail_iter;        /* Newton failures by divergence / by the iteration limit */
} trbdf2_stats;

typedef struct trbdf2_mem trbdf2_mem;

int    trbdf2_create(int n, trbdf2_mem **mem_out);
void   trbdf2_free(trbdf2_mem **mem);
int    trbdf2_set_tolerances(trbdf2_mem *mem, const double *rtol, const double *atol);   /* vectors of length n */
int    trbdf2_set_hmax(trbdf2_mem *mem, double hmax);                 /* <= 0: no limit */
int    trbdf2_set_max_steps(trbdf2_mem *mem, long nmax);
int    trbdf2_set_newton(trbdf2_mem *mem, int max_iter, double tol);  /* tol in the weighted norm, default 0.05 */
int    trbdf2_set_step_control(trbdf2_mem *mem, double safety, double fac_min, double fac_max,
                               double keep_lo, double keep_hi);       /* reuse the LU while keep_lo <= h/h_lu <= keep_hi */
int    trbdf2_set_jac_policy(trbdf2_mem *mem, int max_steps_between_jac, int use_user_jac);
int    trbdf2_set_failure_policy(trbdf2_mem *mem, double step_factor, int max_consecutive);
int    trbdf2_reinit(trbdf2_mem *mem);                                /* forget step size, Jacobian and stats of the segment */
int    trbdf2_reset_stats(trbdf2_mem *mem);                           /* zero the counters only */
int    trbdf2_solve(trbdf2_mem *mem, trbdf2_rhs_fn rhs, trbdf2_jac_fn jac, trbdf2_solout_fn solout, void *user,
                    double *t, double *y, double tend, double h0);    /* h0 <= 0: estimate */
int    trbdf2_interpolate(trbdf2_mem *mem, double t, double *y_out);  /* on the last accepted step */
double trbdf2_get_h(trbdf2_mem *mem);                                 /* step size proposed for the next step */
void   trbdf2_get_stats(trbdf2_mem *mem, trbdf2_stats *stats);
const char *trbdf2_get_err_msg(trbdf2_mem *mem);

#endif
