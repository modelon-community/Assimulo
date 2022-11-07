/* based on f2c (version 20100827) translation of radau_decsol.f. */
/* Note: Due to this, matrices (double*) are stored in Fortran-style column major format */
/* Also: This is the source of some odd looking pointer increments */

#include <stdio.h>
#include <math.h>
#include "radau5_c.h"

#ifdef __RADAU5_WITH_SUPERLU
#include "superlu_double.h"
#include "superlu_complex.h"
#include "superlu_util.h"
#endif /*__RADAU5_WITH_SUPERLU*/

#define TRUE_ (1)
#define FALSE_ (0)
#define radau5_abs(x) ((x) >= 0 ? (x) : -(x))
#define radau_min(a,b) ((a) <= (b) ? (a) : (b))
#define radau_max(a,b) ((a) >= (b) ? (a) : (b))
#define copysign(a,b) (((a < 0 && b > 0) || (a > 0 && b < 0)) ? (-a) : (a))

/* Butcher array related constants */
#define t11 .091232394870892942792
#define t12 -.14125529502095420843
#define t13 -.030029194105147424492
#define t21 .24171793270710701896
#define t22 .20412935229379993199
#define t23 .38294211275726193779
#define t31 .96604818261509293619

#define ti11 4.325579890063155351
#define ti12 .33919925181580986954
#define ti13 .54177053993587487119
#define ti21 -4.1787185915519047273
#define ti22 -.32768282076106238708
#define ti23 .47662355450055045196
#define ti31 -.50287263494578687595
#define ti32 2.5719269498556054292
#define ti33 -.59603920482822492497

/* other magical numbers */
#define expm .66666666666666663

/* forward declarations of private functions */
static int _radcor(radau_mem_t *rmem, int n, FP_CB_f fcn, void *fcn_EXT,
			double *x, double *y, double *xend, double *h__,
			double *rtol, double *atol,
			FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void *jac_EXT, int *ijac,
			FP_CB_solout solout, void *solout_EXT, int *iout, int *idid,
			double *z1, double *z2, double *z3,
			double *y0, double *scal,
			double *f1, double *f2, double *f3,
			double *fjac, double *e1, double *e2r, double *e2i,
			double *cont, double *werr);

static int _dec(int n, double *a, int *ip, int *ier);
static int _sol(int n, double *a, double *b, int *ip);

static int _decc(int n, double *ar, double *ai, int *ip, int *ier);
static int _solc(int n, double *ar, double *ai, double *br, double *bi, int *ip);

static int _decomr(radau_linsol_mem_t *mem, int n, double *fjac, double fac1, double *e1, int *ier);
static int _decomc(radau_linsol_mem_t *mem, int n, double *fjac, double alphn, double betan,
	double *e2r, double *e2i, int *ier);

static int _slvrad(radau_mem_t *rmem, int n, double fac1, double alphn, double betan, 
	double *e1, double *e2r, double *e2i, 
	double *z1, double *z2, double *z3,
	double *f1, double *f2, double *f3);

static int _estrad(radau_mem_t *rmem, int n, double h,
	double dd1, double dd2, double dd3,
	FP_CB_f fcn, void *fcn_EXT,
	double *y0, double *y,
	double x,double *e1,
	double *z1, double *z2, double *z3,
	double *cont, double *werr,
	double *f1, double *f2,
	double *err, int first, int reject);

int radau5_solve(void *radau_mem, FP_CB_f fcn, void *fcn_EXT,
	double *x, double *y, double *xend, double *h__,
	double *rtol, double *atol,
	FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void *jac_EXT, int *ijac,
	FP_CB_solout solout, void *solout_EXT, int *iout, int *idid)
{
    int i;
    double quot;
    double tolst;
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!radau_mem){ return RADAU_PARA_RADAU_MEM_NULL; }
/* ---------------------------------------------------------- */
/*     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC) */
/*     SYSTEM OF FIRST ORDER ORDINARY DIFFERENTIAL EQUATIONS */
/*                     Y'=F(X,Y). */
/*     THE METHOD USED IS AN IMPLICIT RUNGE-KUTTA METHOD (RADAU IIA) */
/*     OF ORDER 5 WITH STEP SIZE CONTROL AND CONTINUOUS OUTPUT. */
/*     CF. SECTION IV.8 */

/*     AUTHORS: E. HAIRER AND G. WANNER */
/*              UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES */
/*              CH-1211 GENEVE 24, SWITZERLAND */
/*              E-MAIL:  Ernst.Hairer@math.unige.ch */
/*                       Gerhard.Wanner@math.unige.ch */

/*     THIS CODE IS PART OF THE BOOK: */
/*         E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL */
/*         EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS. */
/*         SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14, */
/*         SPRINGER-VERLAG 1991, SECOND EDITION 1996. */

/*     VERSION OF JULY 9, 1996 */
/*     (latest small correction: January 18, 2002) */

/*     INPUT PARAMETERS */
/*     ---------------- */
/*     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE */
/*                 VALUE OF F(X,Y): */
/*                    SUBROUTINE FCN(N,X,Y,F) */
/*                    DOUBLE PRECISION X,Y(N),F(N) */
/*                    F(1)=...   ETC. */

/*     X           INITIAL X-VALUE */

/*     Y(N)        INITIAL VALUES FOR Y */

/*     XEND        FINAL X-VALUE (XEND-X MAY BE POSITIVE OR NEGATIVE) */

/*     H           INITIAL STEP SIZE GUESS; */
/*                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT, */
/*                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD. */
/*                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS */
/*                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6). */

/*     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. VECTORS OF LENGTH N. */

/*     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES */
/*                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y */
/*                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1; SUPPLY */
/*                 A DUMMY SUBROUTINE IN THE CASE IJAC=0). */
/*                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM */
/*                    SUBROUTINE JAC(N,X,Y,DFY) */
/*                    DOUBLE PRECISION X,Y(N),DFY(N,N) */
/*                    DFY(1,1)= ... */

/*     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN: */
/*                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE */
/*                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED. */
/*                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC. */

/*     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE */
/*                 NUMERICAL SOLUTION DURING INTEGRATION. */
/*                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP. */
/*                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. */
/*                 IT MUST HAVE THE FORM */
/*                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N, */
/*                                       IRTRN) */
/*                    DOUBLE PRECISION X,Y(N) */
/*                    .... */
/*                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH */
/*                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS */
/*                    THE FIRST GRID-POINT). */
/*                 "XOLD" IS THE PRECEEDING GRID-POINT. */
/*                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN */
/*                    IS SET <0, RADAU5 RETURNS TO THE CALLING PROGRAM. */

/*     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT: */
/*                    IOUT=0: SUBROUTINE IS NEVER CALLED */
/*                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT. */

/* ----------------------------------------------------------------------- */

/*     OUTPUT PARAMETERS */
/*     RETURN FLAGS, SEE HEADER FILE*/

/*     ----------------- */
/*     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED */
/*                 (AFTER SUCCESSFUL RETURN X=XEND). */

/*     Y(N)        NUMERICAL SOLUTION AT X */

/*     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP */

/*     IDID        IF RETURN == RADAU_SUCCESS_SOLOUT_INTERRUPT, IDID = RETURN FLAG FROM SOLOUT CALLBACK */

/* ----------------------------------------------------------------------- */

	/* -------- CHECK AND CHANGE THE TOLERANCES */
	for (i = 0; i < rmem->n; ++i) {
		if (atol[i] <= 0. || rtol[i] <= rmem->input->uround * 10.) {
			printf("TOLERANCES (%i) ARE TOO SMALL \n", i);
			return RADAU_ERROR_INCONSISTENT_INPUT;
		} else {
			quot = atol[i] / rtol[i];
			rmem->rtol[i] = pow(rtol[i], expm) * .1;
			rmem->atol[i] = rmem->rtol[i] * quot;
		}
	}

	tolst = rmem->rtol[0]; /* internal representation of rtol */
	/* check inputs for consistency */
	if (!rmem->input->_checked){
		/* sanity check on fnewt */
		if (rmem->input->fnewt <= rmem->input->uround / tolst){
			return RADAU_ERROR_INCONSISTENT_INPUT;
		}
		rmem->input->_checked = TRUE_;
	}
	/* automatically set input parameters that have not been set directly */
	if (!rmem->input->hmax_set){
		rmem->input->hmax = *xend - *x;
	}
	if (!rmem->input->fnewt_set){
		
		/* Note: This is using the transformed tolerances */
		rmem->input->fnewt = radau_max(rmem->input->uround * 10 / tolst, radau_min(.03, pow(tolst, .5)));
	}

	if (rmem->lin_sol->sparseLU && !(*ijac)){
		printf("CURIOUS INPUT; ANALYTICAL JACOBIAN DISABLED, IJAC = %i, WHICH IS REQUIRED FOR SPARSE SOLVER\n", *ijac);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}
	/* POSSIBLE ADDITIONAL RETURN FLAG */
	*idid = 0;
	/* -------- CALL TO CORE INTEGRATOR ------------ */
    return _radcor(rmem, rmem->n, (FP_CB_f)fcn, fcn_EXT, x, y, xend, h__, rmem->rtol, rmem->atol, 
					(FP_CB_jac)jac, (FP_CB_jac_sparse) jac_sparse, jac_EXT, ijac,
					(FP_CB_solout)solout, solout_EXT, iout, idid,
					rmem->z1, rmem->z2, rmem->z3, rmem->y0,
					rmem->scal, rmem->f1, rmem->f2, rmem->f3,
					rmem->lin_sol->jac, rmem->lin_sol->e1, rmem->lin_sol->e2r, rmem->lin_sol->e2i, 
					rmem->cont, rmem->werr);
} /* radau5_solve */


static int _radcor(radau_mem_t *rmem, int n, FP_CB_f fcn, void *fcn_EXT,
	double *x, double *y, double *xend, double *h__,
	double *rtol, double *atol,
	FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void *jac_EXT, int *ijac, 
	FP_CB_solout solout, void *solout_EXT, int *iout, int *idid,
	double *z1, double *z2, double *z3,
	double *y0, double *scal,
	double *f1, double *f2, double *f3,
	double *fjac, double *e1, double *e2r, double *e2i,
	double *cont, double *werr)
{
	int ret = RADAU_OK;
    int i, j;
    double a1, a2, a3; /* auxiliary variables */
	double ak, ak1, ak2, ak3, acont3; /* auxiliary variables */
	double qt, f1i, f2i, f3i, c1q, c2q, c3q, z1i, z2i, z3i; /* auxiliary variables */

    int ier;
	/* Newton related parameters */
	int newt; /* number of newton iterations taken */
    double xph; /* x + h, used for newton rhs */
	double thq, cfac;
	double dyno, dyth, theta;
	double qnewt;
	double dynold = 0;
	double thqold = 0;

	/* Linear system related parameters */
	double fac1 = 0; /* diagonal factor for real linear system */
	double betan = 0; /* real diagonal factor for complex linear system */
	double alphn = 0; /* complex diagonal factor for complex linear system */

	/* stepsize control related parameters */
	/* TODO: Separate into predive controller and normal one */
	int first; /* switch if first timestep */
    int last; /* switch if current step is last one */
	int reject; /* switch if step is to be rejected */
    double hopt; /* backup stepsize for output */
	double xold; /* previous time-point */
    double hold; /* previos stepsize */
	double hnew; /* new stepsize */
	double err; /* local error norm, return from estrad */
	double fac; /* stepsize adjustment factor */
	double quot; /* auxiliary; a quotient */ 
	double hhfac;
	double hmaxn; /* maximum stepsize in current step */
	double facgus;
	double posneg; /* sign */

	/* misc */
    double delt; /* perturbation in jacobian finite diffs */
	double ysafe; /* auxiliary value in jacobians*/

	/* solout callback related */
	int irtrn; /*additional solout return TODO: use return value instead */

	/* tracking for various failure states */
	int nunexpect; /* number of unexpected step failures, max = 10 */
	int nsing; /* number of singular jacobians */

    /* Parameter adjustments */
    --cont;
    --f3;
    --f2;
    --f1;
    --scal;
    --y0;
    --z3;
    --z2;
    --z1;
    --y;
    --rtol;
    --atol;
    --werr;
    fjac -= 1 + n;
    e2i -= 1 + n;
    e2r -= 1 + n;
    e1 -= 1 + n;

    posneg = copysign(1., *xend - *x);
    hmaxn = radau_min(radau5_abs(rmem->input->hmax), radau5_abs(*xend - *x));
    if (radau5_abs(*h__) <= rmem->input->uround * 10.) {
		*h__ = 1e-6;
    }
    *h__ = radau_min(radau5_abs(*h__), hmaxn);
    *h__ = copysign(*h__, posneg);
    hold = *h__;
    reject = FALSE_;
    first = TRUE_;
    last = FALSE_;
    if ((*x + *h__ * 1.0001 - *xend) * posneg >= 0.) {
		*h__ = *xend - *x;
		last = TRUE_;
    }
    hopt = *h__;
    cfac = rmem->input->step_size_safety * ((rmem->input->nmax_newton << 1) + 1);
    nsing = 0;
    nunexpect = 0;
    xold = *x;

	for (i = 1; i <= n; ++i) {
		scal[i] = atol[i] + rtol[i] * radau5_abs(y[i]);
	}
    hhfac = *h__;
    ier = (*fcn)(n, *x, &y[1], &y0[1], fcn_EXT);
	rmem->stats->nfcn++;
	if (ier < 0){
		goto L79;
	}
/* --- BASIC INTEGRATION STEP/REPEAT STEP WITH FRESH JACOBIAN */
L10:
/* *** *** *** *** *** *** *** */
/*  COMPUTATION OF THE JACOBIAN */
/* *** *** *** *** *** *** *** */
	if (!rmem->new_jac_req){
		goto L30; /* no new jacobian required; reuse old one + LU factorization */
	}
    rmem->stats->njac++;
    if (*ijac == 0) {
		/* --- COMPUTE JACOBIAN MATRIX NUMERICALLY */
		/* --- JACOBIAN IS FULL */
		for (i = 1; i <= n; ++i) {
			ysafe = y[i];
			delt = sqrt(rmem->input->uround * radau_max(1e-5, radau5_abs(ysafe)));
			y[i] = ysafe + delt;
			ier = (*fcn)(n, *x, &y[1], &cont[1], fcn_EXT);
			if (ier != RADAU_CALLBACK_OK) {
				y[i] = ysafe - delt;
				ier =  (*fcn)(n, *x, &y[1], &cont[1], fcn_EXT);
				if (ier != RADAU_CALLBACK_OK) {
					y[i] = ysafe;
					goto L79;
				}
				for (j = 1; j <= n; ++j) {
					fjac[j + i * n] = (y0[j] - cont[j]) / delt;
				}
			} else {
				for (j =  1; j <= n; ++j) {
					fjac[j + i * n] = (cont[j] - y0[j]) / delt;
				}
			}
			y[i] = ysafe;
		}
    } else {
		/* --- COMPUTE JACOBIAN MATRIX ANALYTICALLY */
		if (rmem->lin_sol->sparseLU){
			#ifdef __RADAU5_WITH_SUPERLU
			rmem->lin_sol->nnz_actual = rmem->lin_sol->nnz;
			ier = (*jac_sparse)(n, *x, &y[1], &(rmem->lin_sol->nnz_actual), rmem->lin_sol->jac, rmem->lin_sol->jac_indices, rmem->lin_sol->jac_indptr, jac_EXT);
			if (ier != RADAU_CALLBACK_OK){
				goto L79;
			}
			ier = sparse_csc_add_diagonal(n, &(rmem->lin_sol->nnz_actual), rmem->lin_sol->jac, rmem->lin_sol->jac_indices, rmem->lin_sol->jac_indptr);
			if (ier != 0){
				goto L183;
			}
			rmem->lin_sol->LU_with_fresh_jac = TRUE_;
			#endif /*__RADAU5_WITH_SUPERLU*/
		} else {
			ier = (*jac)(n, *x, &y[1], &fjac[1 + n], jac_EXT);
			if (ier != RADAU_CALLBACK_OK){
				goto L79;
			}
		}
    }
    rmem->jac_is_fresh = TRUE_; /* flag that Jacobian is freshly computed this timestep */
	rmem->new_jac_req = FALSE_; /* no new Jacobian required */
/* --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS */
L20:
    fac1 = rmem->mconst->u1 / *h__;
    _decomr(rmem->lin_sol, n, fjac, fac1, e1, &ier);
    if (ier != 0) {
		goto L185;
    }
	alphn = rmem->mconst->alph / *h__;
    betan = rmem->mconst->beta / *h__;
    _decomc(rmem->lin_sol, n, fjac, alphn, betan, e2r, e2i, &ier);
    if (ier != 0) {
		goto L185;
    }
	if (rmem->lin_sol->sparseLU){
		#ifdef __RADAU5_WITH_SUPERLU
			rmem->lin_sol->LU_with_fresh_jac = FALSE_; /* has once been used to create a decomposition now */
		#endif /*__RADAU5_WITH_SUPERLU*/
	}
    rmem->stats->ludecomps++; /* increment LU decompositions counter */
/* --- NEXT STEP, NO NEW JAC/LU */
L30:
    rmem->stats->nsteps++;
    if (rmem->stats->nsteps > rmem->input->nmax) {
		goto L178;
    }
    if (radau5_abs(*h__) * .1 <= radau5_abs(*x) * rmem->input->uround) {
		goto L177;
    }
    xph = *x + *h__;
	/* *** *** *** *** *** *** *** */
	/*  STARTING VALUES FOR NEWTON ITERATION */
	/* *** *** *** *** *** *** *** */
    if (first || rmem->input->newton_start_zero) {
		for (i = 1; i <= n; ++i) {
			z1[i] = 0.;
			z2[i] = 0.;
			z3[i] = 0.;
			f1[i] = 0.;
			f2[i] = 0.;
			f3[i] = 0.;
		}
    } else {
		c3q = *h__ / hold;
		c1q = rmem->mconst->c1 * c3q;
		c2q = rmem->mconst->c2 * c3q;
		for (i = 1; i <= n; ++i) {
			ak1 = cont[i + n];
			ak2 = cont[i + 2*n];
			ak3 = cont[i + 3*n];
			z1i = c1q * (ak1 + (c1q - rmem->mconst->c2m1) * (ak2 + (c1q - rmem->mconst->c1m1) * ak3));
			z2i = c2q * (ak1 + (c2q - rmem->mconst->c2m1) * (ak2 + (c2q - rmem->mconst->c1m1) * ak3));
			z3i = c3q * (ak1 + (c3q - rmem->mconst->c2m1) * (ak2 + (c3q - rmem->mconst->c1m1) * ak3));
			z1[i] = z1i;
			z2[i] = z2i;
			z3[i] = z3i;
			f1[i] = ti11 * z1i + ti12 * z2i + ti13 * z3i;
			f2[i] = ti21 * z1i + ti22 * z2i + ti23 * z3i;
			f3[i] = ti31 * z1i + ti32 * z2i + ti33 * z3i;
		}
    }
	/* *** *** *** *** *** *** *** */
	/*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
	/* *** *** *** *** *** *** *** */
    newt = 0;
    rmem->faccon = pow(radau_max(rmem->faccon,rmem->input->uround), .8);
    theta = radau5_abs(rmem->input->theta_jac_recomp);
	/* --- NEWTON ITERATION */
L40:
    if (newt >= rmem->input->nmax_newton) {
		goto L78;
    }
	/* ---     COMPUTE THE RIGHT-HAND SIDE */
    for (i = 1; i <= n; ++i) {
		cont[i] = y[i] + z1[i];
    }
    ier = (*fcn)(n, *x + rmem->mconst->c1 * *h__, &cont[1], &z1[1], fcn_EXT);
    rmem->stats->nfcn++;
    if (ier != RADAU_CALLBACK_OK) {
		goto L79;
    }
    for (i = 1; i <= n; ++i) {
		cont[i] = y[i] + z2[i];
    }
    ier = (*fcn)(n, *x + rmem->mconst->c2 * *h__, &cont[1], &z2[1], fcn_EXT);
    rmem->stats->nfcn++;
    if (ier != RADAU_CALLBACK_OK) {
		goto L79;
    }
    for (i = 1; i <= n; ++i) {
		cont[i] = y[i] + z3[i];
    }
    ier = (*fcn)(n, xph, &cont[1], &z3[1], fcn_EXT);
    rmem->stats->nfcn++;
    if (ier != RADAU_CALLBACK_OK) {
		goto L79;
    }
	/* ---     SOLVE THE LINEAR SYSTEMS */
    for (i = 1; i <= n; ++i) {
		a1 = z1[i];
		a2 = z2[i];
		a3 = z3[i];
		z1[i] = ti11 * a1 + ti12 * a2 + ti13 * a3;
		z2[i] = ti21 * a1 + ti22 * a2 + ti23 * a3;
		z3[i] = ti31 * a1 + ti32 * a2 + ti33 * a3;
    }
    ier = _slvrad(rmem, n, fac1, alphn, betan, &e1[1 + n],
			&e2r[1 + n], &e2i[1 + n],
			&z1[1], &z2[1], &z3[1], &f1[1], &f2[1], &f3[1]);
	if (ier != 0){
		goto L184;
	}
    ++newt;
    dyno = 0.;
    for (i = 1; i <= n; ++i) {
		dyno = dyno + (z1[i] / scal[i] * z1[i] / scal[i]) 
		            + (z2[i] / scal[i] * z2[i] / scal[i]) 
					+ (z3[i] / scal[i] * z3[i] / scal[i]);
    }
    dyno = sqrt(dyno / (3*n));
	/* ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE */
    if (newt > 1 && newt < rmem->input->nmax_newton) {
		thq = dyno / dynold;
		if (newt == 2) {
			theta = thq;
		} else {
			theta = sqrt(thq * thqold);
		}
		thqold = thq;
		if (theta < .99) {
			rmem->faccon = theta / (1. - theta);
			dyth = rmem->faccon * dyno * pow(theta, rmem->input->nmax_newton - 1 - newt) / rmem->input->fnewt;
			if (dyth >= 1.) {
				qnewt = radau_max(1e-4, radau_min(20., dyth));
				hhfac = pow(qnewt, -1. / (rmem->input->nmax_newton + 4. - 1 - newt)) * .8;
				*h__ = hhfac * *h__;
				reject = TRUE_;
				last = FALSE_;
				if (rmem->jac_is_fresh) {
					goto L20;
				}
				rmem->new_jac_req = TRUE_;
				goto L10;
			}
		} else {
			goto L78;
		}
    }
    dynold = radau_max(dyno, rmem->input->uround);
    for (i = 1; i <= n; ++i) { /* update Newton RHS ? */
		f1i = f1[i] + z1[i];
		f2i = f2[i] + z2[i];
		f3i = f3[i] + z3[i];
		f1[i] = f1i;
		f2[i] = f2i;
		f3[i] = f3i;
		z1[i] = t11 * f1i + t12 * f2i + t13 * f3i;
		z2[i] = t21 * f1i + t22 * f2i + t23 * f3i;
		z3[i] = t31 * f1i + f2i;
    }
    if (rmem->faccon * dyno > rmem->input->fnewt) {
		goto L40;
    }
	/* --- ERROR ESTIMATION */
    ret = _estrad(rmem, n, *h__, rmem->mconst->dd1, rmem->mconst->dd2, rmem->mconst->dd3, (FP_CB_f) fcn, fcn_EXT,
				   &y0[1], &y[1], *x, &e1[1 + n],
				   &z1[1], &z2[1], &z3[1], &cont[1], &werr[1], &f1[1], &f2[1],
				   &err, first, reject);
	if (ret < 0){
		goto L184;
	}
	/* --- COMPUTATION OF HNEW */
	/* --- WE REQUIRE .2<=HNEW/H<=8. */
    fac = radau_min(rmem->input->step_size_safety, cfac / (newt + (rmem->input->nmax_newton << 1)));
    quot = radau_max(rmem->input->fac_upper ,radau_min(rmem->input->fac_lower, pow(err, .25) / fac));
    hnew = *h__ / quot;
	/* *** *** *** *** *** *** *** */
	/*  IS THE ERROR SMALL ENOUGH ? */
	/* *** *** *** *** *** *** *** */
    if (err < 1.) {
	/* --- STEP IS ACCEPTED */
		first = FALSE_;
		rmem->stats->naccpt++;
		if (rmem->input->pred_step_control) {
			/*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
			if (rmem->stats->naccpt > 1) {
				facgus = rmem->h_old / *h__ * pow(err * err / rmem->erracc, .25) / rmem->input->step_size_safety;
				facgus = radau_max(rmem->input->fac_upper, radau_min(rmem->input->fac_lower, facgus));
				quot = radau_max(quot,facgus);
				hnew = *h__ / quot;
			}
			rmem->h_old = *h__;
			rmem->erracc = radau_max(.01, err);
		}
		xold = *x;
		hold = *h__;
		*x = xph;
		for (i = 1; i <= n; ++i) {
			y[i] += z3[i];
			z2i = z2[i];
			z1i = z1[i];
			cont[i + n] = (z2i - z3[i]) / rmem->mconst->c2m1;
			ak = (z1i - z2i) / rmem->mconst->c1mc2;
			acont3 = z1i / rmem->mconst->c1;
			acont3 = (ak - acont3) / rmem->mconst->c2;
			cont[i + 2*n] = (ak - cont[i + n]) / rmem->mconst->c1m1;
			cont[i + 3*n] = cont[i + 2*n] - acont3;
		}
		for (i = 1; i <= n; ++i) {
			scal[i] = atol[i] + rtol[i] * radau5_abs(y[i]);
		}
		if (*iout != 0) {
			rmem->xsol = *x;
			for (i = 1; i <= n; ++i) {
				cont[i] = y[i];
			}
			rmem->hsol = hold;

			irtrn = 0;
			rmem->_dense_output_valid = TRUE_;
			(*solout)(rmem->stats->naccpt + 1, xold, rmem->xsol, &y[1], &werr[1], n, &irtrn, solout_EXT);
			if (irtrn < 0) {
				goto L179;
			}
			rmem->_dense_output_valid = FALSE_;
		}
		rmem->jac_is_fresh = FALSE_;
		if (last) {
			*h__ = hopt;
			return RADAU_SUCCESS;
		}
		(*fcn)(n, *x, &y[1], &y0[1], fcn_EXT);
		rmem->stats->nfcn++;
		hnew = posneg * radau_min(radau5_abs(hnew), hmaxn);
		hopt = hnew;
		hopt = radau_min(*h__,hnew);
		if (reject) {
			hnew = posneg * radau_min(radau5_abs(hnew), radau5_abs(*h__));
		}
		reject = FALSE_;
		if ((*x + hnew / rmem->input->quot1 - *xend) * posneg >= 0.) {
			*h__ = *xend - *x;
			last = TRUE_;
		} else {
			qt = hnew / *h__;
			hhfac = *h__;
			if (theta <= rmem->input->theta_jac_recomp && qt >= rmem->input->quot1 && qt <= rmem->input->quot2) {
				goto L30;
			}
			*h__ = hnew;
		}
		hhfac = *h__;
		if (theta <= rmem->input->theta_jac_recomp) {
			goto L20;
		}
		rmem->new_jac_req = TRUE_;
		goto L10;
    } else {
		/* --- STEP IS REJECTED */
		reject = TRUE_;
		last = FALSE_;
		if (first) {
			*h__ *= .1;
			hhfac = .1;
		} else {
			hhfac = hnew / *h__;
			*h__ = hnew;
		}
		if (rmem->stats->naccpt >= 1) {
			rmem->stats->nreject++;
		}
		if (rmem->jac_is_fresh) {
			goto L20;
		}
		rmem->new_jac_req = TRUE_;
		goto L10;
    }
/* --- UNEXPECTED STEP-REJECTION OR SINGULAR JACOBIAN */
L78:
    if (ier != 0) {
		++nsing;
		if (nsing >= 5) {
			goto L176;
		}
    }
    *h__ *= .5;
    hhfac = .5;
    reject = TRUE_;
    last = FALSE_;
    if (rmem->jac_is_fresh) {
		goto L20;
    }
	rmem->new_jac_req = TRUE_;
    goto L10;
/* --- UNEXPECTED STEP-REJECTION */
L79:
	if (ier < RADAU_CALLBACK_ERROR_INVALID_NNZ){
		#ifdef __RADAU5_WITH_SUPERLU
			printf("FAILURE IN JACOBIAN EVALUATIONS, NNZ TOO SMALL, SPECIFIED NNZ: %i, ACTUAL: %i \n", rmem->lin_sol->nnz, -(ier - RADAU_CALLBACK_ERROR_INVALID_NNZ));
		#endif /*__RADAU5_WITH_SUPERLU*/
		return RADAU_ERROR_NNZ_TOO_SMALL;
	}

	switch(ier){
		case RADAU_CALLBACK_ERROR_RECOVERABLE:
			++nunexpect;
			if (nunexpect >= 10) {
				goto L175;
			}
			*h__ *= .5;
			hhfac = .5;
			reject = TRUE_;
			last = FALSE_;
			if (rmem->jac_is_fresh) {
				goto L20;
			}
			rmem->new_jac_req = TRUE_;
			goto L10;
			break;

		case RADAU_CALLBACK_ERROR_NONRECOVERABLE:
			goto L186;
			break;

		case RADAU_CALLBACK_ERROR_INVALID_JAC_FORMAT:
			goto L182;
			break;
	}
/* --- FAIL EXIT */
L175:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("REPEATEDLY UNEXPECTED STEP REJECTIONS\n");
	return RADAU_ERROR_REP_STEP_REJECT;
/* --- FAIL EXIT, REPEATED SINGULAR JACOBIAN*/
L176:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("MATRIX IS REPEATEDLY SINGULAR IER= %i\n", ier);
	return RADAU_ERROR_JAC_SINGULAR;
/* --- FAIL EXIT, STEP SIZE TOO SMALL*/
L177:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("STEP SIZE TOO SMALL, H= %e\n", *h__);
	return RADAU_ERROR_STEPSIZE_TOO_SMALL;
/* --- FAIL EXIT, NMAX EXCEEDED*/
L178:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("MORE THAN NMAX = %i STEPS ARE NEEDED\n", rmem->input->nmax);
	return RADAU_ERROR_NMAX_TOO_SMALL;
/* --- FAILURE EXIT, ERROR IN SPARSE JACOBIAN CALLBACK*/
L182:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("JACOBIAN GIVEN IN WRONG FORMAT, REQUIRED SPARSE FORMAT: CSC\n");
	return RADAU_ERROR_WRONG_SPARSE_JAC_FORMAT;
#ifdef __RADAU5_WITH_SUPERLU /* TODO */
	L183:
		printf("EXIT OF RADAU5 AT X = %e \n", *x);
		printf("UNEXPECTED MALLOC FAILURE\n");
		return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;
#endif /*__RADAU5_WITH_SUPERLU*/
/* --- FAIL EXIT, UNEXPECTED SUPERLU FAILURE*/
L184:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("UNEXPECTED FAILURE OF SUPERLU FUNCTION CALL, error code = %i \n", ret);
	return RADAU_ERROR_UNEXPECTED_SUPERLU_FAILURE;
/* --- FAIL EXIT, OTHER SUPERLU FAILURE*/
L185:
	if (ier < 0){ /* incorrect input to function call */
		goto L184;
	}
	if (ier <= n){ /* factorization singular */
		goto L78;
	}else{ /* memory allocation failure */
		printf("EXIT OF RADAU5 AT X = %e \n", *x);
		printf("SUPERLU MEMORY ALLOCATION FAILURE, NUMBER OF BYTES ALLOCATED AT POINT OF FAILURE: %i.\n", ier - n);
		return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;
	}
/* --- FAIL EXIT, UNRECOVERABLE EXCEPTION IN FCN OR JAC CALLBACK*/
L186:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("UNRECOVERABLE EXCEPTION ENCOUNTERED IN PROBLEM CALLBACK\n");
	return RADAU_ERROR_UNRECOVERABLE_CALLBACK_ERROR;
/* --- EXIT CAUSED BY SOLOUT */
L179:
    *idid = irtrn; /* return from solout callback */
	return RADAU_SUCCESS_SOLOUT_INTERRUPT;
} /* _radcor */

int radau_get_cont_output_single(void *radau_mem, int i, double x, double *out){
	/*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN */
	/*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
	/*     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR */
	/*     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU5). */
	double s;
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}
	if (!rmem->_dense_output_valid){return -100;}

	s = (x - rmem->xsol) / rmem->hsol;
    *out = rmem->cont[i] + s * (rmem->cont[i + rmem->n] + (s - rmem->mconst->c2m1)
	       * (rmem->cont[i + 2*rmem->n] + (s - rmem->mconst->c1m1) * rmem->cont[i + 3*rmem->n]));
	return RADAU_OK;
}

int radau_get_cont_output(void *radau_mem, double x, double *out){
	/* see radau_get_cont_output_single; outputs solution in all components */
	int i;
	double s;
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}
	if (!rmem->_dense_output_valid){return -100;}
	
	s  = (x - rmem->xsol) / rmem->hsol;
	for(i = 0; i < rmem->n; i++){
		out[i] = rmem->cont[i] + s * (rmem->cont[i + rmem->n] + (s - rmem->mconst->c2m1)
	       		 * (rmem->cont[i + 2*rmem->n] + (s - rmem->mconst->c1m1) * rmem->cont[i + 3*rmem->n]));
	}
	return RADAU_OK;
}

static int _dec(int n, double *a, int *ip, int *ier)
{
    int i, j, k, m;
    double t;
    int kp1;

/* VERSION REAL DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION. */
/*  INPUT.. */
/*     N = ORDER OF MATRIX. */
/*     A = MATRIX TO BE TRIANGULARIZED. */
/*  OUTPUT.. */
/*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
/*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
/*           SINGULAR AT STAGE K. */
/*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N). */
/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;

    /* Function Body */
    *ier = 0;
    ip[n] = 1;
    if (n == 1) {
		goto L70;
    }
    for (k = 1; k <= n - 1; ++k) {
		kp1 = k + 1;
		m = k;
		for (i = kp1; i <= n; ++i) {
			if (radau5_abs(a[i + k * n]) > radau5_abs(a[m + k * n])) {
			m = i;
			}
		}
		ip[k] = m;
		t = a[m + k * n];
		if (m == k) {
			goto L20;
		}
		ip[n] = -ip[n];
		a[m + k * n] = a[k + k * n];
		a[k + k * n] = t;
L20:
		if (t == 0.) {
			goto L80;
		}
		t = 1. / t;
		for (i = kp1; i <= n; ++i) {
			a[i + k * n] = -a[i + k * n] * t;
		}
		for (j = kp1; j <= n; ++j) {
			t = a[m + j * n];
			a[m + j * n] = a[k + j * n];
			a[k + j * n] = t;
			if (t == 0.) {
				goto L45;
			}
			for (i = kp1; i <= n; ++i) {
				a[i + j * n] += a[i + k * n] * t;
			}
L45:
			;
		}
    }
L70:
    k = n;
    if (a[n + n * n] == 0.) {
		goto L80;
    }
    return RADAU_OK;
L80:
    *ier = k;
    ip[n] = 0;
    return RADAU_OK;
} /* _dec */


static int _sol(int n, double *a, double *b, int *ip)
{
    int i, k, m;
    double t;
    int kb, km1;

/* VERSION REAL DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX. */
/*    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
/*    B = RIGHT HAND SIDE VECTOR. */
/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
/*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    B = SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    --b;
    a -= 1 + n;

    if (n == 1) {
		goto L50;
    }
    for (k = 1; k <= n - 1; ++k) {
		m = ip[k];
		t = b[m];
		b[m] = b[k];
		b[k] = t;
		for (i = k + 1; i <= n; ++i) {
			b[i] += a[i + k * n] * t;
		}
    }
    for (kb = 1; kb <= n - 1; ++kb) {
		km1 = n - kb;
		k = km1 + 1;
		b[k] /= a[k + k * n];
		t = -b[k];
		for (i = 1; i <= km1; ++i) {
			b[i] += a[i + k * n] * t;
		}
    }
L50:
    b[1] /= a[n + 1];
    return RADAU_OK;
} /* _sol */


static int _decc(int n, double *ar, double *ai, int *ip, int *ier)
{
    int i, j, k, m;
    double ti, tr;
    int kp1;
    double den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
/*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
/*  INPUT.. */
/*     N = ORDER OF MATRIX. */
/*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
/*  OUTPUT.. */
/*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
/*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
/*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*                                                    REAL PART. */
/*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*                                                    IMAGINARY PART. */
/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
/*           SINGULAR AT STAGE K. */
/*  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;

    /* Function Body */
    *ier = 0;
    ip[n] = 1;
    if (n == 1) {
		goto L70;
    }
    for (k = 1; k <= n - 1; ++k) {
		kp1 = k + 1;
		m = k;
		for (i = kp1; i <= n; ++i) {
			if (radau5_abs(ar[i + k * n]) + radau5_abs(ai[i + k * n]) > radau5_abs(ar[m + k * n]) + radau5_abs(ai[m + k * n])) {
				m = i;
			}
		}
		ip[k] = m;
		tr = ar[m + k * n];
		ti = ai[m + k * n];
		if (m == k) {
			goto L20;
		}
		ip[n] = -ip[n];
		ar[m + k * n] = ar[k + k * n];
		ai[m + k * n] = ai[k + k * n];
		ar[k + k * n] = tr;
		ai[k + k * n] = ti;
L20:
		if (radau5_abs(tr) + radau5_abs(ti) == 0.) {
			goto L80;
		}
		den = tr * tr + ti * ti;
		tr /= den;
		ti = -ti / den;
		for (i = kp1; i <= n; ++i) {
			prodr = ar[i + k * n] * tr - ai[i + k * n] * ti;
			prodi = ai[i + k * n] * tr + ar[i + k * n] * ti;
			ar[i + k * n] = -prodr;
			ai[i + k * n] = -prodi;
		}
		for (j = kp1; j <= n; ++j) {
			tr = ar[m + j * n];
			ti = ai[m + j * n];
			ar[m + j * n] = ar[k + j * n];
			ai[m + j * n] = ai[k + j * n];
			ar[k + j * n] = tr;
			ai[k + j * n] = ti;
			if (radau5_abs(tr) + radau5_abs(ti) == 0.) {
				goto L48;
			}
			if (ti == 0.) {
				for (i = kp1; i <= n; ++i) {
					prodr = ar[i + k * n] * tr;
					prodi = ai[i + k * n] * tr;
					ar[i + j * n] += prodr;
					ai[i + j * n] += prodi;
				}
				goto L48;
			}
			if (tr == 0.) {
				for (i = kp1; i <= n; ++i) {
					prodr = -ai[i + k * n] * ti;
					prodi = ar[i + k * n] * ti;
					ar[i + j * n] += prodr;
					ai[i + j * n] += prodi;
				}
				goto L48;
			}
			for (i = kp1; i <= n; ++i) {
				prodr = ar[i + k * n] * tr - ai[i + k * n] * ti;
				prodi = ai[i + k * n] * tr + ar[i + k * n] * ti;
				ar[i + j * n] += prodr;
				ai[i + j * n] += prodi;
			}
L48:
			;
		}
    }
L70:
    k = n;
    if (radau5_abs(ar[n + n * n]) + radau5_abs(ai[n + n * n]) == 0.) {
		goto L80;
    }
    return RADAU_OK;
L80:
    *ier = k;
    ip[n] = 0;
    return RADAU_OK;
} /* _decc */


static int _solc(int n, double *ar, double *ai, double *br, double *bi, int *ip)
{
    int i, k, m, kb;
    double ti, tr;
    int km1, kp1;
    double den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX. */
/*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
/*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
/*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    (BR,BI) = SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    --bi;
    --br;
    ai -= 1 + n;
    ar -= 1 + n;

    /* Function Body */
    if (n == 1) {
		goto L50;
    }
    for (k = 1; k <= n - 1; ++k) {
		kp1 = k + 1;
		m = ip[k];
		tr = br[m];
		ti = bi[m];
		br[m] = br[k];
		bi[m] = bi[k];
		br[k] = tr;
		bi[k] = ti;
		for (i = kp1; i <= n; ++i) {
			prodr = ar[i + k * n] * tr - ai[i + k * n] * ti;
			prodi = ai[i + k * n] * tr + ar[i + k * n] * ti;
			br[i] += prodr;
			bi[i] += prodi;
		}
    }
    for (kb = 1; kb <= n - 1; ++kb) {
		km1 = n - kb;
		k = km1 + 1;
		den = ar[k + k * n] * ar[k + k * n] + ai[k + k * n] 
			* ai[k + k * n];
		prodr = br[k] * ar[k + k * n] + bi[k] * ai[k + k * n];
		prodi = bi[k] * ar[k + k * n] - br[k] * ai[k + k * n];
		br[k] = prodr / den;
		bi[k] = prodi / den;
		tr = -br[k];
		ti = -bi[k];
		for (i = 1; i <= km1; ++i) {
			prodr = ar[i + k * n] * tr - ai[i + k * n] * ti;
			prodi = ai[i + k * n] * tr + ar[i + k * n] * ti;
			br[i] += prodr;
			bi[i] += prodi;
		}
    }
L50:
    den = ar[n + 1] * ar[n + 1] + ai[n + 1] * ai[n + 1];
    prodr = br[1] * ar[n + 1] + bi[1] * ai[n + 1];
    prodi = bi[1] * ar[n + 1] - br[1] * ai[n + 1];
    br[1] = prodr / den;
    bi[1] = prodi / den;
    return RADAU_OK;
} /* _solc */


static int _decomr(radau_linsol_mem_t *lmem, int n, double *fjac, double fac1, double *e1, int *ier)
{
    int i, j;

	if(lmem->sparseLU){
		#ifdef __RADAU5_WITH_SUPERLU
			superlu_setup_d((SuperLU_aux_d*)lmem->slu_aux_d, fac1, lmem->jac, lmem->jac_indices, lmem->jac_indptr, lmem->LU_with_fresh_jac, lmem->nnz_actual);
			*ier = superlu_factorize_d((SuperLU_aux_d*)lmem->slu_aux_d);
		#endif /*__RADAU5_WITH_SUPERLU*/
	}else{
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= n; ++i) {
				e1[i + j * n] = -fjac[i + j * n];
			}
			e1[j + j * n] += fac1;
		}
		_dec(n, e1, lmem->ip1, ier);
	}
	return RADAU_OK;
} /* _decomr */


int _decomc(radau_linsol_mem_t *lmem, int n, double *fjac, double alphn, double betan,
	double *e2r, double *e2i, int *ier)
{
    int i, j;

	if (lmem->sparseLU){
		#ifdef __RADAU5_WITH_SUPERLU
		superlu_setup_z((SuperLU_aux_z*)lmem->slu_aux_z, alphn, betan, lmem->jac, lmem->jac_indices, lmem->jac_indptr, lmem->LU_with_fresh_jac, lmem->nnz_actual);
		*ier = superlu_factorize_z((SuperLU_aux_z*)lmem->slu_aux_z);
		#endif /*__RADAU5_WITH_SUPERLU*/
	}else{
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= n; ++i) {
				e2r[i + j * n] = -fjac[i + j * n];
				e2i[i + j * n] = 0.;
			}
			e2r[j + j * n] += alphn;
			e2i[j + j * n] = betan;
		}
		_decc(n, e2r, e2i, lmem->ip2, ier);
	}
	return RADAU_OK;
} /* _decomc */


static int _slvrad(radau_mem_t *rmem, int n, double fac1, double alphn, double betan, 
	double *e1, double *e2r, double *e2i, 
	double *z1, double *z2, double *z3,
	double *f1, double *f2, double *f3)
{
    int i;
    double s2, s3;
	int ret = RADAU_OK;

    for (i = 0; i < n; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * fac1;
		z2[i] = z2[i] + s2 * alphn - s3 * betan;
		z3[i] = z3[i] + s3 * alphn + s2 * betan;
    }

	if (rmem->lin_sol->sparseLU){
		#ifdef __RADAU5_WITH_SUPERLU
			ret = superlu_solve_d((SuperLU_aux_d*)rmem->lin_sol->slu_aux_d, z1);
			if (ret != 0) { return ret; }
			ret = superlu_solve_z((SuperLU_aux_z*)rmem->lin_sol->slu_aux_z, z2, z3);
		#endif /*__RADAU5_WITH_SUPERLU*/
	}else{
		_sol(n, e1, z1, rmem->lin_sol->ip1);
    	_solc(n, e2r, e2i, z2, z3, rmem->lin_sol->ip2);
	}
	rmem->stats->lusolves++; /* increment factorization counter */
	return ret;
} /* _slvrad */


static int _estrad(radau_mem_t *rmem, int n, double h,
	double dd1, double dd2, double dd3,
	FP_CB_f fcn, void *fcn_EXT,
	double *y0, double *y,
	double x, double *e1,
	double *z1, double *z2, double *z3,
	double *cont, double *werr,
	double *f1, double *f2,
	double *err, int first, int reject)
{
    int i;
	int ret = RADAU_OK;
    double hee1 = dd1 / h;
    double hee2 = dd2 / h;
    double hee3 = dd3 / h;

	for (i = 0; i < n; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
	}

	if (rmem->lin_sol->sparseLU){
		#ifdef __RADAU5_WITH_SUPERLU
			ret = superlu_solve_d((SuperLU_aux_d*)rmem->lin_sol->slu_aux_d, cont);
			if (ret < 0){
				return ret;
			}
		#endif /*__RADAU5_WITH_SUPERLU*/
	}else{
		_sol(n, e1, cont, rmem->lin_sol->ip1);
	}

    *err = 0.;
    for (i = 0; i < n; ++i) {
		werr[i] = cont[i] / rmem->scal[i];
		*err += werr[i] * werr[i];
    }
    *err = radau_max(sqrt(*err / n), 1e-10);

    if (*err < 1.) {
		return RADAU_OK;
    }
    if (first || reject) {
		for (i = 0; i < n; ++i) {
			cont[i] = y[i] + cont[i];
		}
		(*fcn)(n, x, cont, f1, fcn_EXT);
		rmem->stats->nfcn++;
		for (i = 0; i < n; ++i) {
			cont[i] = f1[i] + f2[i];
		}

		if (rmem->lin_sol->sparseLU){
			#ifdef __RADAU5_WITH_SUPERLU
				ret = superlu_solve_d((SuperLU_aux_d*)rmem->lin_sol->slu_aux_d, cont);
				if (ret < 0){
					return ret;
				}
			#endif /*__RADAU5_WITH_SUPERLU*/
		}else{
			_sol(n, e1, cont, rmem->lin_sol->ip1);
		}

		*err = 0.;
		for (i = 0; i < n; ++i) {
			werr[i] = cont[i] / rmem->scal[i];
			*err += werr[i] * werr[i];
		}
		*err = radau_max(sqrt(*err / n), 1e-10);
    }
    return ret;
} /* _estrad */

