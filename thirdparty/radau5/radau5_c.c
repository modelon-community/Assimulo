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

int radau5_c(void *radau_mem, FP_CB_f fcn, void *fcn_EXT,
	double *x, double *y, double *xend, double *h__,
	double *rtol, double *atol,
	FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void *jac_EXT, int *ijac,
	FP_CB_solout solout, void *solout_EXT, int *iout,
	double *work, int *idid)
{
    int i;
    double facl;
    double facr;
    double hmax;
    double thet;
    double quot;
    double quot1, quot2;
    double fnewt;
    double tolst;
    double uround;
	int radcor_ret = 0;
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
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

/*     WORK        ARRAY OF WORKING SPACE OF LENGTH 20. */
/*                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS */
/*                 FOR THE CODE. FOR STANDARD USE OF THE CODE */
/*                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE */
/*                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE. */

/* ---------------------------------------------------------------------- */

/*     SOPHISTICATED SETTING OF PARAMETERS */
/*    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16. */

/*    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION, */
/*              DEFAULT 0.9D0. */

/*    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
/*              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS */
/*              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER */
/*              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO */
/*              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP. */
/*              DEFAULT 0.001D0. */

/*    WORK(4)   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. */
/*              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER. */
/*              DEFAULT MIN(0.03D0,RTOL(1)**0.5D0) */

/*    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE */
/*              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A */
/*              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR */
/*              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE */
/*              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS */
/*              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD. */
/*              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 . */

/*    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X. */

/*    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION */
/*              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION */
/*                 WORK(8) <= HNEW/HOLD <= WORK(9) */
/*              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0 */

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
	if (!radau_mem){
		return RADAU_PARA_RADAU_MEM_NULL;
	}

    /* Parameter adjustments */
    --work;

	/* -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0 */
    if (work[1] == 0.) {
		uround = 1e-16;
    } else {
		uround = work[1];
		if (uround <= 1e-19 || uround >= 1.) {
			printf(" COEFFICIENTS HAVE 20 DIGITS, UROUND= \t %e \n", uround);
			return RADAU_ERROR_INCONSISTENT_INPUT;
		}
    }
	/* -------- CHECK AND CHANGE THE TOLERANCES */
	for (i = 0; i < rmem->n; ++i) {
		if (atol[i] <= 0. || rtol[i] <= uround * 10.) {
			printf("TOLERANCES (%i) ARE TOO SMALL \n", i);
			return RADAU_ERROR_INCONSISTENT_INPUT;
		} else {
			quot = atol[i] / rtol[i];
			rmem->rtol[i] = pow(rtol[i], expm) * .1;
			rmem->atol[i] = rmem->rtol[i] * quot;
		}
	}
	/* ------ THET, DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
    if (work[3] == 0.) {
		thet = .001;
    } else {
		thet = work[3];
		if (thet >= 1.) {
			printf("CURIOUS INPUT FOR WORK(3)= %f \n", thet);
			return RADAU_ERROR_INCONSISTENT_INPUT;
		}
    }
	/* --- FNEWT, STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. */
    tolst = rmem->rtol[0];
    if (work[4] == 0.) {
		fnewt = radau_max(uround * 10 / tolst, radau_min(.03, pow(tolst, .5)));
    } else {
		fnewt = work[4];
		if (fnewt <= uround / tolst) {
			printf("CURIOUS INPUT FOR WORK(4)= %f \n", fnewt);
			return RADAU_ERROR_INCONSISTENT_INPUT;
		}
    }
	/* --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST. */
    if (work[5] == 0.) {
		quot1 = 1.;
    } else {
		quot1 = work[5];
    }
    if (work[6] == 0.) {
		quot2 = 1.2;
    } else {
		quot2 = work[6];
    }
    if (quot1 > 1. || quot2 < 1.) {
		printf("CURIOUS INPUT FOR WORK(5, 6)= %f \t %f \n", quot1, quot2);
		return RADAU_ERROR_INCONSISTENT_INPUT;
    }
	/* -------- MAXIMAL STEP SIZE */
    if (work[7] == 0.) {
		hmax = *xend - *x;
    } else {
		hmax = work[7];
    }
	/* -------  FACL,FACR, PARAMETERS FOR STEP SIZE SELECTION */
    if (work[8] == 0.) {
		facl = 5.;
    } else {
		facl = 1. / work[8];
    }
    if (work[9] == 0.) {
		facr = .125;
    } else {
		facr = 1. / work[9];
    }
    if (facl < 1. || facr > 1.) {
		printf("CURIOUS INPUT FOR WORK(8, 9)= %f \t %f \n", facl, facr);
		return RADAU_ERROR_INCONSISTENT_INPUT;
    }
	if (rmem->lin_sol->sparseLU && !(*ijac)){
		printf("CURIOUS INPUT; ANALYTICAL JACOBIAN DISABLED, IJAC = %i, WHICH IS REQUIRED FOR SPARSE SOLVER\n", *ijac);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}
	/* POSSIBLE ADDITIONAL RETURN FLAG */
	*idid = 0;
	/* -------- CALL TO CORE INTEGRATOR ------------ */
    radcor_ret = radcor_(rmem, rmem->n, (FP_CB_f)fcn, fcn_EXT, x, y, xend, &hmax, h__, rmem->rtol, rmem->atol, 
						 (FP_CB_jac)jac, (FP_CB_jac_sparse) jac_sparse, jac_EXT, ijac,
						 (FP_CB_solout)solout, solout_EXT, iout, idid, &uround, &thet, &fnewt,
						 &quot1, &quot2, &facl, &facr,
						 rmem->z1, rmem->z2, rmem->z3, rmem->y0,
						 rmem->scal, rmem->f1, rmem->f2, rmem->f3,
						 rmem->lin_sol->jac, rmem->lin_sol->e1, rmem->lin_sol->e2r, rmem->lin_sol->e2i, 
						 rmem->cont, rmem->werr);
	return radcor_ret;
} /* radau5_ */


int radcor_(radau_mem_t *rmem, int n, FP_CB_f fcn, void *fcn_EXT,
	double *x, double *y, double *xend, double *hmax, double *h__,
	double *rtol, double *atol,
	FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void *jac_EXT, int *ijac, 
	FP_CB_solout solout, void *solout_EXT, int *iout, int *idid,
	double *uround, double *thet,
	double *fnewt, double *quot1, double *quot2,
	double *facl, double *facr,
	double *z1, double *z2, double *z3,
	double *y0, double *scal,
	double *f1, double *f2, double *f3,
	double *fjac, double *e1, double *e2r, double *e2i,
	double *cont, double *werr)
{
	int ret = RADAU_OK;
    int i, j;
    double a1, a2, a3;
    int nunexpect;
    double ak;
    double qt, ak1, ak2, ak3, f1i, f2i, f3i, c1q, 
	    c2q, c3q, z1i, z2i, z3i, fac;
    int ier;
    double xph, thq, err, fac1, cfac;
    double hold;
    double delt, hnew;
    int last;
    double hopt, xold;
    int newt;
    double dyno, dyth, quot, hhfac, betan, alphn, theta, ysafe, hmaxn;
    int nsing;
    int first;
    int irtrn, nrsol, nsolu;
    double qnewt, xosol, acont3;
    double faccon;
    int reject;
    double facgus;
    double posneg;

	fac1 = 0; alphn = 0; betan = 0;  /* avoid warnings */

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
    hmaxn = radau_min(radau5_abs(*hmax), radau5_abs(*xend - *x));
    if (radau5_abs(*h__) <= *uround * 10.) {
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
    faccon = 1.;
    cfac = rmem->para->step_size_safety * ((rmem->para->nmax_newton << 1) + 1);
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
	if (!rmem->para->new_jac_req){
		goto L30; /* no new jacobian required; reuse old one + LU factorization */
	}
    rmem->stats->njac++;
    if (*ijac == 0) {
		/* --- COMPUTE JACOBIAN MATRIX NUMERICALLY */
		/* --- JACOBIAN IS FULL */
		for (i = 1; i <= n; ++i) {
			ysafe = y[i];
			delt = sqrt(*uround * radau_max(1e-5, radau5_abs(ysafe)));
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
    rmem->para->jac_is_fresh = TRUE_; /* flag that Jacobian is freshly computed this timestep */
	rmem->para->new_jac_req = FALSE_; /* no new Jacobian required */
/* --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS */
L20:
    fac1 = rmem->mconst->u1 / *h__;
    decomr_(rmem->lin_sol, n, fjac, fac1, e1, &ier);
    if (ier != 0) {
		goto L185;
    }
	alphn = rmem->mconst->alph / *h__;
    betan = rmem->mconst->beta / *h__;
    decomc_(rmem->lin_sol, n, fjac, alphn, betan, e2r, e2i, &ier);
    if (ier != 0) {
		goto L185;
    }
	if (rmem->lin_sol->sparseLU){
		#ifdef __RADAU5_WITH_SUPERLU
			rmem->lin_sol->LU_with_fresh_jac = FALSE_; /* has once been used to create a decomposition now */
		#endif /*__RADAU5_WITH_SUPERLU*/
	}
    rmem->stats->ludecomps++; /* increment LU decompositions counter */
/* --- COMPUTE STEPSIZE; step accepted*/
L30:
    rmem->stats->nsteps++;
    if (rmem->stats->nsteps > rmem->para->nmax) {
		goto L178;
    }
    if (radau5_abs(*h__) * .1 <= radau5_abs(*x) * *uround) {
		goto L177;
    }
    xph = *x + *h__;
	/* *** *** *** *** *** *** *** */
	/*  STARTING VALUES FOR NEWTON ITERATION */
	/* *** *** *** *** *** *** *** */
    if (first || rmem->para->newton_start_zero) {
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
    faccon = pow(radau_max(faccon,*uround), .8);
    theta = radau5_abs(*thet);
	/* --- NEWTON ITERATION */
L40:
    if (newt >= rmem->para->nmax_newton) {
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
    ier = slvrad_(rmem, n, fac1, alphn, betan, &e1[1 + n],
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
    if (newt > 1 && newt < rmem->para->nmax_newton) {
		thq = dyno / rmem->para->dynold;
		if (newt == 2) {
			theta = thq;
		} else {
			theta = sqrt(thq * rmem->para->thqold);
		}
		rmem->para->thqold = thq;
		if (theta < .99) {
			faccon = theta / (1. - theta);
			dyth = faccon * dyno * pow(theta, rmem->para->nmax_newton - 1 - newt) / *fnewt;
			if (dyth >= 1.) {
				qnewt = radau_max(1e-4, radau_min(20., dyth));
				hhfac = pow(qnewt, -1. / (rmem->para->nmax_newton + 4. - 1 - newt)) * .8;
				*h__ = hhfac * *h__;
				reject = TRUE_;
				last = FALSE_;
				if (rmem->para->jac_is_fresh) {
					goto L20;
				}
				rmem->para->new_jac_req = TRUE_;
				goto L10;
			}
		} else {
			goto L78;
		}
    }
    rmem->para->dynold = radau_max(dyno, *uround);
    for (i = 1; i <= n; ++i) {
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
    if (faccon * dyno > *fnewt) {
		goto L40;
    }
	/* --- ERROR ESTIMATION */
    ret = estrad_(rmem, n, *h__, rmem->mconst->dd1, rmem->mconst->dd2, rmem->mconst->dd3, (FP_CB_f) fcn, fcn_EXT,
				   &y0[1], &y[1], *x, &e1[1 + n],
				   &z1[1], &z2[1], &z3[1], &cont[1], &werr[1], &f1[1], &f2[1],
				   &err, first, reject);
	if (ret < 0){
		goto L184;
	}
	/* --- COMPUTATION OF HNEW */
	/* --- WE REQUIRE facl<=HNEW/H<=facr */
    fac = radau_min(rmem->para->step_size_safety, cfac / (newt + (rmem->para->nmax_newton << 1)));
    quot = radau_max(*facr ,radau_min(*facl, pow(err, .25) / fac));
    hnew = *h__ / quot;
	/* *** *** *** *** *** *** *** */
	/*  IS THE ERROR SMALL ENOUGH ? */
	/* *** *** *** *** *** *** *** */
    if (err < 1.) {
	/* --- STEP IS ACCEPTED */
		first = FALSE_;
		rmem->stats->naccpt++;
		if (rmem->para->pred_step_control) {
			/*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
			if (rmem->stats->naccpt > 1) {
				facgus = rmem->para->h_old / *h__ * pow(err * err / rmem->para->erracc, .25) / rmem->para->step_size_safety;
				facgus = radau_max(*facr, radau_min(*facl, facgus));
				quot = radau_max(quot,facgus);
				hnew = *h__ / quot;
			}
			rmem->para->h_old = *h__;
			rmem->para->erracc = radau_max(.01, err);
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
			nrsol = rmem->stats->naccpt + 1;
			rmem->xsol = *x;
			xosol = xold;
			for (i = 1; i <= n; ++i) {
				cont[i] = y[i];
			}
			nsolu = n;
			rmem->hsol = hold;

			irtrn = 0;
			(*solout)(nrsol, xosol, rmem->xsol, &y[1], &werr[1], &nsolu, &irtrn, solout_EXT);
			if (irtrn < 0) {
				goto L179;
			}
		}
		rmem->para->jac_is_fresh = FALSE_;
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
		if ((*x + hnew / *quot1 - *xend) * posneg >= 0.) {
			*h__ = *xend - *x;
			last = TRUE_;
		} else {
			qt = hnew / *h__;
			hhfac = *h__;
			if (theta <= *thet && qt >= *quot1 && qt <= *quot2) {
				goto L30;
			}
			*h__ = hnew;
		}
		hhfac = *h__;
		if (theta <= *thet) {
			goto L20;
		}
		rmem->para->new_jac_req = TRUE_;
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
		if (rmem->para->jac_is_fresh) {
			goto L20;
		}
		rmem->para->new_jac_req = TRUE_;
		goto L10;
    }
/* --- UNEXPECTED STEP-REJECTION, SINGULAR JACOBIAN */
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
    if (rmem->para->jac_is_fresh) {
		goto L20;
    }
	rmem->para->new_jac_req = TRUE_;
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
			if (rmem->para->jac_is_fresh) {
				goto L20;
			}
			rmem->para->new_jac_req = TRUE_;
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
	printf("MORE THAN NMAX = %i STEPS ARE NEEDED\n", rmem->para->nmax);
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
} /* radcor_ */

int radau_get_cont_output_single(void *radau_mem, int i, double x, double *out){
	/*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN */
	/*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
	/*     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR */
	/*     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU5). */
	double s;
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}
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
	
	s  = (x - rmem->xsol) / rmem->hsol;
	for(i = 0; i < rmem->n; i++){
		out[i] = rmem->cont[i] + s * (rmem->cont[i + rmem->n] + (s - rmem->mconst->c2m1)
	       		 * (rmem->cont[i + 2*rmem->n] + (s - rmem->mconst->c1m1) * rmem->cont[i + 3*rmem->n]));
	}
	return RADAU_OK;
}

int dec_(int n, double *a, int *ip, int *ier)
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
} /* dec_ */


int sol_(int n, double *a, double *b, int *ip)
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
} /* sol_ */


int decc_(int n, double *ar, double *ai, int *ip, int *ier)
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
} /* decc_ */


int solc_(int n, double *ar, double *ai, double *br, double *bi, int *ip)
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
} /* solc_ */


int decomr_(radau_linsol_mem_t *lmem, int n, double *fjac, double fac1, double *e1, int *ier)
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
		dec_(n, e1, lmem->ip1, ier);
	}
	return RADAU_OK;
} /* decomr_ */


int decomc_(radau_linsol_mem_t *lmem, int n, double *fjac, double alphn, double betan,
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
		decc_(n, e2r, e2i, lmem->ip2, ier);
	}
	return RADAU_OK;
} /* decomc_ */


int slvrad_(radau_mem_t *rmem, int n, double fac1, double alphn, double betan, 
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
		sol_(n, e1, z1, rmem->lin_sol->ip1);
    	solc_(n, e2r, e2i, z2, z3, rmem->lin_sol->ip2);
	}
	rmem->stats->lusolves++; /* increment factorization counter */
	return ret;
} /* slvrad_ */


int estrad_(radau_mem_t *rmem, int n, double h,
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
		sol_(n, e1, cont, rmem->lin_sol->ip1);
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
			sol_(n, e1, cont, rmem->lin_sol->ip1);
		}

		*err = 0.;
		for (i = 0; i < n; ++i) {
			werr[i] = cont[i] / rmem->scal[i];
			*err += werr[i] * werr[i];
		}
		*err = radau_max(sqrt(*err / n), 1e-10);
    }
    return ret;
} /* estrad_ */

