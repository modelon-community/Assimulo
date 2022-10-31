/* based on f2c (version 20100827) translation of radau_decsol.f. */
/* Note: Due to this, matrices (double*) are stored in Fortran-style column major format */
/* Also: This is the source of some odd looking pointer increments */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "radau5_c.h"

/* Common Block Declarations */
struct conra5_{
    int nn, nn2, nn3, nn4;
    double xsol, hsol, c2m1, c1m1;
} conra5_1;

/* Table of constant values */
#define c_b54 .5
#define c_b91 81.
#define c_b92 .33333333333333331
#define c_b93 9.
#define c_b114 .8
#define c_b116 .25

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

int radau5_c(void* radau_mem, int n, FP_CB_f fcn, void* fcn_PY,
	double *x, double *y, double *xend, double *h__,
	double *rtol, double *atol, int *itol,
	FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, int *ijac, int sparse_LU,
	FP_CB_solout solout, void* solout_PY, int *iout,
	double *work, int *lwork, int *iwork, int *liwork, int *idid,
	Radau_SuperLU_aux* radau_slu_aux)
{
    static int i, nit;
    static double facl;
    static int ndec, njac;
    static double facr, safe;
    static int nfcn;
    static int pred;
    static double hmax;
    static int nmax;
    static double thet, expm;
    static int nsol;
    static double quot;
    static double quot1, quot2;
    static double fnewt;
    static int nstep;
    static double tolst;
    static int naccpt;
    static int nrejct;
    static int startn;
    static double uround;
	int radcor_ret = 0;
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
/* ---------------------------------------------------------- */
/*     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC) */
/*     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS */
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
/*     N           DIMENSION OF THE SYSTEM */

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

/*     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY */
/*                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N. */

/*     ITOL        SWITCH FOR RTOL AND ATOL: */
/*                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS. */
/*                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF */
/*                     Y(I) BELOW RTOL*ABS(Y(I))+ATOL */
/*                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS. */
/*                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) BELOW */
/*                     RTOL(I)*ABS(Y(I))+ATOL(I). */

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
/*                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N, */
/*                                       IRTRN) */
/*                    DOUBLE PRECISION X,Y(N),CONT(LRC) */
/*                    .... */
/*                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH */
/*                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS */
/*                    THE FIRST GRID-POINT). */
/*                 "XOLD" IS THE PRECEEDING GRID-POINT. */
/*                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN */
/*                    IS SET <0, RADAU5 RETURNS TO THE CALLING PROGRAM. */

/*          -----  CONTINUOUS OUTPUT: ----- */
/*                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION */
/*                 FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH */
/*                 THE FUNCTION */
/*                        >>>   CONTR5(I,S,CONT,LRC)   <<< */
/*                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH */
/*                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE */
/*                 S SHOULD LIE IN THE INTERVAL [XOLD,X]. */
/*                 DO NOT CHANGE THE ENTRIES OF CONT(LRC), IF THE */
/*                 DENSE OUTPUT FUNCTION IS USED. */

/*     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT: */
/*                    IOUT=0: SUBROUTINE IS NEVER CALLED */
/*                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT. */

/*     WORK        ARRAY OF WORKING SPACE OF LENGTH "LWORK". */
/*                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS */
/*                 FOR THE CODE. FOR STANDARD USE OF THE CODE */
/*                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE */
/*                 CALLING. SEE BELOW FOR A MORE SOPHISTICATED USE. */
/*                 WORK(21),..,WORK(LWORK) SERVE AS WORKING SPACE */
/*                 FOR ALL VECTORS AND MATRICES. */
/*                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL */
/*				   THE MINIMUM STORAGE REQUIREMENT IS */
/*                             LWORK = 4*N*N+13*N+20. */

/*     LWORK       DECLARED LENGTH OF ARRAY "WORK". */

/*     IWORK       int WORKING SPACE OF LENGTH "LIWORK". */
/*                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS */
/*                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),.., */
/*                 IWORK(20) TO ZERO BEFORE CALLING. */
/*                 IWORK(21),...,IWORK(LIWORK) SERVE AS WORKING AREA. */
/*                 "LIWORK" MUST BE AT LEAST 3*N+20. */

/*     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". */

/* ---------------------------------------------------------------------- */

/*     SOPHISTICATED SETTING OF PARAMETERS */
/*     ----------------------------------- */
/*              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK */
/*              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),... */
/*              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO. */
/*              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES: */

/*    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. */
/*              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000. */

/*    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE */
/*              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP. */
/*              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7. */

/*    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION */
/*              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD. */
/*              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED. */
/*              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS */
/*              DIFFICULTIES WITH CONVERGENCE (THIS IS THE CASE WHEN */
/*              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.). */
/*              DEFAULT IS IWORK(4)=0. */

/*    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY */
/*              IF IWORK(8).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON) */
/*              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL */
/*              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1. */
/*              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS; */
/*              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES */
/*              OFTEN SLIGHTLY FASTER RUNS */
/* ---------- */

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

/*   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL */
/*                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED) */
/*   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY */
/*                      OR NUMERICALLY) */
/*   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS */
/*   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS */
/*   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST), */
/*                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED) */
/*   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF JACOBIAN MATRICES */
/*   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS */
/*                      THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS, */
/*                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED */
/* ----------------------------------------------------------------------- */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*          DECLARATIONS */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/* *** *** *** *** *** *** *** */
/*        SETTING THE PARAMETERS */
/* *** *** *** *** *** *** *** */
	if (!radau_mem){
		printf("RADAU_MEM not setup properly.\n");
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}

	if (sparse_LU){
		#ifndef __RADAU5_WITH_SUPERLU
			printf("SUPERLU support not enabled.\n");
			return RADAU_ERROR_SUPERLU_NOT_ENABLED;
		#else /*__RADAU5_WITH_SUPERLU*/
			;
		#endif /*__RADAU5_WITH_SUPERLU*/
	}

    /* Parameter adjustments */
    --y;
    --rtol;
    --atol;
    --work;
    --iwork;

    /* Function Body */
    nfcn = 0;
    njac = 0;
    nstep = 0;
    naccpt = 0;
    nrejct = 0;
    ndec = 0;
    nsol = 0;
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
    expm = .66666666666666663;
    if (*itol == 0) {
		if (atol[1] <= 0. || rtol[1] <= uround * 10.) {
			printf(" TOLERANCES ARE TOO SMALL \n");
			return RADAU_ERROR_INCONSISTENT_INPUT;
		} else {
			quot = atol[1] / rtol[1];
			rtol[1] = pow(rtol[1], expm) * .1;
			atol[1] = rtol[1] * quot;
		}
    } else {
		for (i = 1; i <= n; ++i) {
			if (atol[i] <= 0. || rtol[i] <= uround * 10.) {
				printf("TOLERANCES (%i) ARE TOO SMALL \n", i);
				return RADAU_ERROR_INCONSISTENT_INPUT;
			} else {
				quot = atol[i] / rtol[i];
				rtol[i] = pow(rtol[i], expm) * .1;
				atol[i] = rtol[i] * quot;
			}
		}
    }
	/* -------- NMAX, THE MAXIMAL NUMBER OF STEPS ----- */
    if (iwork[2] == 0) {
		nmax = 100000;
    } else {
		nmax = iwork[2];
		if (nmax <= 0) {
			printf("WRONG INPUT IWORK(2)= %i \n", nmax);
			return RADAU_ERROR_INCONSISTENT_INPUT;
		}
    }
	/* -------- NIT, MAXIMAL NUMBER OF NEWTON ITERATIONS */
    if (iwork[3] == 0) {
		nit = 7;
    } else {
		nit = iwork[3];
		if (nit <= 0) {
			printf("CURIOUS INPUT IWORK(3)= %i \n", nit);
			return RADAU_ERROR_INCONSISTENT_INPUT;
		}
    }
	/* -------- STARTN SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS */
    if (iwork[4] == 0) {
		startn = FALSE_;
    } else {
		startn = TRUE_;
    }
	/* -------- PRED STEP SIZE CONTROL */
    if (iwork[8] <= 1) {
		pred = TRUE_;
    } else {
		pred = FALSE_;
    }
	/* --------- SAFE, SAFETY FACTOR IN STEP SIZE PREDICTION */
    if (work[2] == 0.) {
		safe = .9;
    } else {
		safe = work[2];
		if (safe <= .001 || safe >= 1.) {
			printf("CURIOUS INPUT FOR WORK(2)= %f \n", safe);
			return RADAU_ERROR_INCONSISTENT_INPUT;
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
    tolst = rtol[1];
    if (work[4] == 0.) {
		fnewt = radau_max(uround * 10 / tolst, radau_min(.03, pow(tolst, c_b54)));
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
	if (sparse_LU && !(*ijac)){
		printf("CURIOUS INPUT; ANALYTICAL JACOBIAN DISABLED, IJAC = %i, WHICH IS REQUIRED FOR SPARSE SOLVER\n", *ijac);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}
	/* POSSIBLE ADDITIONAL RETURN FLAG */
	*idid = 0;
	/* -------- CALL TO CORE INTEGRATOR ------------ */
    radcor_ret = radcor_(rmem, n, (FP_CB_f)fcn, fcn_PY, x, &y[1], xend, &hmax, h__, &rtol[1], &atol[1], 
						 itol, (FP_CB_jac)jac, (FP_CB_jac_sparse) jac_sparse, jac_PY, ijac,
						 (FP_CB_solout)solout, solout_PY, iout, idid, &nmax, &uround, &safe, &thet, &fnewt,
						 &quot1, &quot2, &nit, sparse_LU, &startn,
						 &pred, &facl, &facr,
						 rmem->z1, rmem->z2, rmem->z3, rmem->y0,
						 rmem->scal, rmem->f1, rmem->f2, rmem->f3,
						 rmem->jac, rmem->e1, rmem->e2r, rmem->e2i, 
						 rmem->ip1, rmem->ip2, rmem->con,
						 &nfcn, &njac, &nstep, &naccpt, &nrejct, &ndec, &nsol,
						 rmem->werr,
						 radau_slu_aux);
    iwork[14] = nfcn;
    iwork[15] = njac;
    iwork[16] = nstep;
    iwork[17] = naccpt;
    iwork[18] = nrejct;
    iwork[19] = ndec;
    iwork[20] = nsol;
	/* -------- RESTORE TOLERANCES */
    expm = 1. / expm;
    if (*itol == 0) {
		quot = atol[1] / rtol[1];
		rtol[1] = pow(rtol[1] * 10., expm);
		atol[1] = rtol[1] * quot;
    } else {
		for (i = 1; i <= n; ++i) {
			quot = atol[i] / rtol[i];
			rtol[i] = pow(rtol[i] * 10., expm);
			atol[i] = rtol[i] * quot;
		}
    }
	return radcor_ret;
} /* radau5_ */


int radcor_(radau_mem_t *rmem, int n, FP_CB_f fcn, void* fcn_PY,
	double *x, double *y, double *xend, double *hmax, double *h__,
	double *rtol, double *atol, int *itol,
	FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, int *ijac, 
	FP_CB_solout solout, void* solout_PY, int *iout, int *idid,
	int *nmax, double *uround, double *safe, double *thet,
	double *fnewt, double *quot1, double *quot2, int *nit,
	int sparse_LU, int *startn, int *pred,
	double *facl, double *facr,
	double *z1, double *z2, double *z3,
	double *y0, double *scal,
	double *f1, double *f2, double *f3,
	double *fjac, double *e1, double *e2r, double *e2i,
	int *ip1, int *ip2, double *cont,
	int *nfcn, int *njac, int *nstep, int *naccpt,
	int *nrejct, int *ndec, int *nsol,
	double *werr, Radau_SuperLU_aux *radau_slu_aux)
{
    double d__1;

    static int i, j;
    static double a1, a2, c1, c2, a3;
    static int n2, n3;
    static double u1;
    static int nunexpect;
    static double ak;
    static double qt, dd1, dd2, dd3, ak1, ak2, ak3, f1i, f2i, f3i, c1q, 
	    c2q, c3q, z1i, z2i, z3i, sq6, fac, cno;
    static int lrc;
    static int ier;
    static double xph, thq, err, fac1, cfac, hacc, c1mc2, beta;
    static double alph, hold;
    static double delt, hnew;
    static int last;
    static double hopt, xold;
    static int newt;
    static double dyno, dyth, quot, hhfac, betan, alphn, theta, 
	    ysafe, hmaxn;
    static int nsing;
    static int first;
    static int irtrn, nrsol, nsolu;
    static double qnewt, xosol, acont3;
	static int caljac;
    static double faccon;
    static double erracc;
    static int reject;
    static double facgus;
    static double dynold, posneg;
    static double thqold;
	/* ---------------------------------------------------------- */
	/*     CORE INTEGRATOR FOR RADAU5 */
	/*     PARAMETERS SAME AS IN RADAU5 WITH WORKSPACE ADDED */
	/* ---------------------------------------------------------- */
	/*         DECLARATIONS */
	/* ---------------------------------------------------------- */
	/* *** *** *** *** *** *** *** */
	/*  INITIALISATIONS */
	/* *** *** *** *** *** *** *** */
	/* --------- DUPLIFY N FOR COMMON BLOCK CONT ----- */

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

    /* Function Body */
    conra5_1.nn = n;
    conra5_1.nn2 = n << 1;
    conra5_1.nn3 = n * 3;
    lrc = n << 2;
	/* ---------- CONSTANTS --------- */
    sq6 = sqrt(6.);
    c1 = (4. - sq6) / 10.;
    c2 = (sq6 + 4.) / 10.;
    conra5_1.c1m1 = c1 - 1.;
    conra5_1.c2m1 = c2 - 1.;
    c1mc2 = c1 - c2;
    dd1 = -(sq6 * 7. + 13.) / 3.;
    dd2 = (sq6 * 7. - 13.) / 3.;
    dd3 = -.33333333333333331;
    u1 = 1./((pow(c_b91, c_b92) + 6. - pow(c_b93, c_b92)) / 30.);
    alph = (12. - pow(c_b91, c_b92) + pow(c_b93, c_b92)) / 60.;
    beta = (pow(c_b91, c_b92) + pow(c_b93, c_b92)) * sqrt(3.) / 60.;
    cno = alph * alph + beta * beta;
    alph /= cno;
    beta /= cno;
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
    cfac = *safe * ((*nit << 1) + 1);
    nsing = 0;
    nunexpect = 0;
    xold = *x;
    if (*iout != 0) {
		nrsol = 1;
		xosol = xold;
		conra5_1.xsol = *x;
		for (i = 1; i <= n; ++i) {
			werr[i] = 0.;
			cont[i] = y[i];
		}
		nsolu = n;
		conra5_1.hsol = hold;

		irtrn = 0;
		(*solout)(&nrsol, &xosol, &conra5_1.xsol, &y[1], &cont[1], &werr[1], &
				  lrc, &nsolu, &irtrn, solout_PY);
		if (irtrn < 0) {
			goto L179;
		}
    }
    n2 = n << 1;
    n3 = n * 3;
    if (*itol == 0) {
		for (i = 1; i <= n; ++i) {
			scal[i] = atol[1] + rtol[1] * radau5_abs(y[i]);
		}
    } else {
		for (i = 1; i <= n; ++i) {
			scal[i] = atol[i] + rtol[i] * radau5_abs(y[i]);
		}
    }
    hhfac = *h__;
    ier = (*fcn)(n, x, &y[1], &y0[1], fcn_PY);
	if (ier < 0){
		goto L79;
	}
    ++(*nfcn);
/* --- BASIC INTEGRATION STEP/REPEAT STEP WITH FRESH JACOBIAN */
L10:
/* *** *** *** *** *** *** *** */
/*  COMPUTATION OF THE JACOBIAN */
/* *** *** *** *** *** *** *** */
    ++(*njac);
    if (*ijac == 0) {
		/* --- COMPUTE JACOBIAN MATRIX NUMERICALLY */
		/* --- JACOBIAN IS FULL */
		for (i = 1; i <= n; ++i) {
			ysafe = y[i];
			delt = sqrt(*uround * radau_max(1e-5, radau5_abs(ysafe)));
			y[i] = ysafe + delt;
			ier = (*fcn)(n, x, &y[1], &cont[1], fcn_PY);
			if (ier != RADAU_CALLBACK_OK) {
				y[i] = ysafe - delt;
				ier =  (*fcn)(n, x, &y[1], &cont[1], fcn_PY);
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
		if (sparse_LU){
			#ifdef __RADAU5_WITH_SUPERLU
			radau_slu_aux->nnz_actual = radau_slu_aux->nnz;
			ier = (*jac_sparse)(n, x, &y[1], &(radau_slu_aux->nnz_actual), radau_slu_aux->jac_data, radau_slu_aux->jac_indices, radau_slu_aux->jac_indptr, jac_PY);
			if (ier != RADAU_CALLBACK_OK){
				goto L79;
			}
			ier = sparse_csc_add_diagonal(n, &(radau_slu_aux->nnz_actual), radau_slu_aux->jac_data, radau_slu_aux->jac_indices, radau_slu_aux->jac_indptr);
			if (ier != 0){
				goto L183;
			}
			radau_slu_aux->fresh_jacobian = 1;
			#endif /*__RADAU5_WITH_SUPERLU*/
		} else {
			ier = (*jac)(n, x, &y[1], &fjac[1 + n], jac_PY);
			if (ier != RADAU_CALLBACK_OK){
				goto L79;
			}
		}
    }
    caljac = TRUE_;
/* --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS */
L20:
    fac1 = u1 / *h__;
    alphn = alph / *h__;
    betan = beta / *h__;
    decomr_(n, fjac,
			&fac1, e1, ip1, &ier, 
			sparse_LU, radau_slu_aux);
    if (ier != 0) {
		goto L185;
    }
    decomc_(n, fjac,
			&alphn, &betan, e2r, e2i,
			ip2, &ier, sparse_LU, radau_slu_aux);
    if (ier != 0) {
		goto L185;
    }
	if (sparse_LU){
		#ifdef __RADAU5_WITH_SUPERLU
			radau_slu_aux->fresh_jacobian = 0; /* has once been used to create a decomposition now */
		#endif /*__RADAU5_WITH_SUPERLU*/
	}
    ++(*ndec);
/* --- COMPUTE STEPSIZE */
L30:
    ++(*nstep);
    if (*nstep > *nmax) {
		goto L178;
    }
    if (radau5_abs(*h__) * .1 <= radau5_abs(*x) * *uround) {
		goto L177;
    }
    xph = *x + *h__;
	/* *** *** *** *** *** *** *** */
	/*  STARTING VALUES FOR NEWTON ITERATION */
	/* *** *** *** *** *** *** *** */
    if (first || *startn) {
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
		c1q = c1 * c3q;
		c2q = c2 * c3q;
		for (i = 1; i <= n; ++i) {
			ak1 = cont[i + n];
			ak2 = cont[i + n2];
			ak3 = cont[i + n3];
			z1i = c1q * (ak1 + (c1q - conra5_1.c2m1) * (ak2 + (c1q - conra5_1.c1m1) * ak3));
			z2i = c2q * (ak1 + (c2q - conra5_1.c2m1) * (ak2 + (c2q - conra5_1.c1m1) * ak3));
			z3i = c3q * (ak1 + (c3q - conra5_1.c2m1) * (ak2 + (c3q - conra5_1.c1m1) * ak3));
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
    faccon = pow(radau_max(faccon,*uround), c_b114);
    theta = radau5_abs(*thet);
	/* --- NEWTON */
L40:
    if (newt >= *nit) {
		goto L78;
    }
	/* ---     COMPUTE THE RIGHT-HAND SIDE */
    for (i = 1; i <= n; ++i) {
		cont[i] = y[i] + z1[i];
    }
    d__1 = *x + c1 * *h__;
    ier = (*fcn)(n, &d__1, &cont[1], &z1[1], fcn_PY);
    ++(*nfcn);
    if (ier != RADAU_CALLBACK_OK) {
		goto L79;
    }
    for (i = 1; i <= n; ++i) {
		cont[i] = y[i] + z2[i];
    }
    d__1 = *x + c2 * *h__;
    ier = (*fcn)(n, &d__1, &cont[1], &z2[1], fcn_PY);
    ++(*nfcn);
    if (ier != RADAU_CALLBACK_OK) {
		goto L79;
    }
    for (i = 1; i <= n; ++i) {
		cont[i] = y[i] + z3[i];
    }
    ier = (*fcn)(n, &xph, &cont[1], &z3[1], fcn_PY);
    ++(*nfcn);
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
    ier = slvrad_(n, &fac1, &alphn, &betan, &e1[1 + n],
			&e2r[1 + n], &e2i[1 + n],
			&z1[1], &z2[1], &z3[1], &f1[1], &f2[1], &f3[1], ip1, ip2,
			sparse_LU, radau_slu_aux);
	if (ier != 0){
		goto L184;
	}
    ++(*nsol);
    ++newt;
    dyno = 0.;
    for (i = 1; i <= n; ++i) {
		dyno = dyno + (z1[i] / scal[i] * z1[i] / scal[i]) 
		            + (z2[i] / scal[i] * z2[i] / scal[i]) 
					+ (z3[i] / scal[i] * z3[i] / scal[i]);
    }
    dyno = sqrt(dyno / n3);
	/* ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE */
    if (newt > 1 && newt < *nit) {
		thq = dyno / dynold;
		if (newt == 2) {
			theta = thq;
		} else {
			theta = sqrt(thq * thqold);
		}
		thqold = thq;
		if (theta < .99) {
			faccon = theta / (1. - theta);
			dyth = faccon * dyno * pow(theta, *nit - 1 - newt) / *fnewt;
			if (dyth >= 1.) {
				qnewt = radau_max(1e-4, radau_min(20., dyth));
				hhfac = pow(qnewt, -1. / (*nit + 4. - 1 - newt)) * .8;
				*h__ = hhfac * *h__;
				reject = TRUE_;
				last = FALSE_;
				if (caljac) {
					goto L20;
				}
				goto L10;
			}
		} else {
			goto L78;
		}
    }
    dynold = radau_max(dyno, *uround);
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
    estrad_(n, h__, &dd1, &dd2, &dd3, (FP_CB_f) fcn, fcn_PY, nfcn,
			&y0[1], &y[1], sparse_LU, x, &e1[1 + n],
			&z1[1], &z2[1], &z3[1], &cont[1], &werr[1], &f1[1], &f2[1], ip1, 
			&scal[1], &err, &first, &reject,
			radau_slu_aux, &ier);
	if (ier){
		goto L184;
	}
	/* --- COMPUTATION OF HNEW */
	/* --- WE REQUIRE .2<=HNEW/H<=8. */
    fac = radau_min(*safe, cfac / (newt + (*nit << 1)));
    quot = radau_max(*facr ,radau_min(*facl, pow(err, c_b116) / fac));
    hnew = *h__ / quot;
	/* *** *** *** *** *** *** *** */
	/*  IS THE ERROR SMALL ENOUGH ? */
	/* *** *** *** *** *** *** *** */
    if (err < 1.) {
	/* --- STEP IS ACCEPTED */
		first = FALSE_;
		++(*naccpt);
		if (*pred) {
			/*       --- PREDICTIVE CONTROLLER OF GUSTAFSSON */
			if (*naccpt > 1) {
				facgus = hacc / *h__ * pow(err * err / erracc, c_b116) / *safe;
				facgus = radau_max(*facr, radau_min(*facl, facgus));
				quot = radau_max(quot,facgus);
				hnew = *h__ / quot;
			}
			hacc = *h__;
			erracc = radau_max(.01, err);
		}
		xold = *x;
		hold = *h__;
		*x = xph;
		for (i = 1; i <= n; ++i) {
			y[i] += z3[i];
			z2i = z2[i];
			z1i = z1[i];
			cont[i + n] = (z2i - z3[i]) / conra5_1.c2m1;
			ak = (z1i - z2i) / c1mc2;
			acont3 = z1i / c1;
			acont3 = (ak - acont3) / c2;
			cont[i + n2] = (ak - cont[i + n]) / conra5_1.c1m1;
			cont[i + n3] = cont[i + n2] - acont3;
		}
		if (*itol == 0) {
			for (i = 1; i <= n; ++i) {
				scal[i] = atol[1] + rtol[1] * radau5_abs(y[i]);
			}
		} else {
			for (i = 1; i <= n; ++i) {
				scal[i] = atol[i] + rtol[i] * radau5_abs(y[i]);
			}
		}
		if (*iout != 0) {
			nrsol = *naccpt + 1;
			conra5_1.xsol = *x;
			xosol = xold;
			for (i = 1; i <= n; ++i) {
				cont[i] = y[i];
			}
			nsolu = n;
			conra5_1.hsol = hold;

			irtrn = 0;
			(*solout)(&nrsol, &xosol, &conra5_1.xsol, &y[1], &cont[1], &werr[1], &
								lrc, &nsolu, &irtrn, solout_PY);
			if (irtrn < 0) {
				goto L179;
			}
		}
		caljac = FALSE_;
		if (last) {
			*h__ = hopt;
			return RADAU_SUCCESS;
		}
		(*fcn)(n, x, &y[1], &y0[1], fcn_PY);
		++(*nfcn);
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
		if (*naccpt >= 1) {
			++(*nrejct);
		}
		if (caljac) {
			goto L20;
		}
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
    if (caljac) {
		goto L20;
    }
    goto L10;
/* --- UNEXPECTED STEP-REJECTION */
L79:
	if (ier < RADAU_CALLBACK_ERROR_INVALID_NNZ){
		#ifdef __RADAU5_WITH_SUPERLU
			printf("FAILURE IN JACOBIAN EVALUATIONS, NNZ TOO SMALL, SPECIFIED NNZ: %i, ACTUAL: %i \n", radau_slu_aux->nnz, -(ier - RADAU_CALLBACK_ERROR_INVALID_NNZ));
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
			if (caljac) {
				goto L20;
			}
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
	printf("MORE THAN NMAX = %i STEPS ARE NEEDED\n", *nmax);
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
	printf("UNEXPECTED FAILURE OF SUPERLU FUNCTION CALL, ier = %i \n", ier);
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


double contr5_c(int *i, double *x, double *cont, int *lrc)
{
    double ret_val;
    static double s;
	if(*lrc){;}; /* TODO: Remove; only here to remove "unused" warning */

/* ---------------------------------------------------------- */
/*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN */
/*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
/*     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR */
/*     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU5). */
/* ---------------------------------------------------------- */
    /* Parameter adjustments */
    --cont;

    s = (*x - conra5_1.xsol) / conra5_1.hsol;
    ret_val = cont[*i] + s * (cont[*i + conra5_1.nn] + (s - conra5_1.c2m1)
	     * (cont[*i + conra5_1.nn2] + (s - conra5_1.c1m1) * cont[*i + conra5_1.nn3]));
    return ret_val;
} /* contr5_ */


int dec_(int n, double *a, int *ip, int *ier)
{
    static int i, j, k, m;
    static double t;
    static int kp1;

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
    return 0;
L80:
    *ier = k;
    ip[n] = 0;
    return 0;
} /* dec_ */


int sol_(int n, double *a, double *b, int *ip)
{
    static int i, k, m;
    static double t;
    static int kb, km1, kp1;

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
		kp1 = k + 1;
		m = ip[k];
		t = b[m];
		b[m] = b[k];
		b[k] = t;
		for (i = kp1; i <= n; ++i) {
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
    return 0;
} /* sol_ */


int decc_(int n, double *ar, double *ai, int *ip, int *ier)
{
    static int i, j, k, m;
    static double ti, tr;
    static int kp1;
    static double den, prodi, prodr;

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
    return 0;
L80:
    *ier = k;
    ip[n] = 0;
    return 0;
} /* decc_ */


int solc_(int n, double *ar, double *ai, double *br, double *bi, int *ip)
{
    static int i, k, m, kb;
    static double ti, tr;
    static int km1, kp1;
    static double den, prodi, prodr;

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
    return 0;
} /* solc_ */


int decomr_(int n, double *fjac,double *fac1, double *e1,
	int *ip1, int *ier, int sparse_LU, Radau_SuperLU_aux *radau_slu_aux)
{
    static int i, j;

	if(sparse_LU){
		#ifdef __RADAU5_WITH_SUPERLU
			superlu_setup_d(radau_slu_aux->slu_aux_d, *fac1, radau_slu_aux->jac_data, radau_slu_aux->jac_indices, radau_slu_aux->jac_indptr, radau_slu_aux->fresh_jacobian, radau_slu_aux->nnz_actual);
			*ier = superlu_factorize_d(radau_slu_aux->slu_aux_d);
		#endif /*__RADAU5_WITH_SUPERLU*/
	}else{
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= n; ++i) {
				e1[i + j * n] = -fjac[i + j * n];
			}
			e1[j + j * n] += *fac1;
		}
		dec_(n, e1, ip1, ier);
	}
	return 0;
} /* decomr_ */


int decomc_(int n, double *fjac, double *alphn, double *betan,
	double *e2r, double *e2i, int *ip2,
	int *ier, int sparse_LU, Radau_SuperLU_aux *radau_slu_aux)
{
    static int i, j;

	if (sparse_LU){
		#ifdef __RADAU5_WITH_SUPERLU
		superlu_setup_z(radau_slu_aux->slu_aux_z, *alphn, *betan, radau_slu_aux->jac_data, radau_slu_aux->jac_indices, radau_slu_aux->jac_indptr, radau_slu_aux->fresh_jacobian, radau_slu_aux->nnz_actual);
		*ier = superlu_factorize_z(radau_slu_aux->slu_aux_z);
		#endif /*__RADAU5_WITH_SUPERLU*/
	}else{
		for (j = 1; j <= n; ++j) {
			for (i = 1; i <= n; ++i) {
				e2r[i + j * n] = -fjac[i + j * n];
				e2i[i + j * n] = 0.;
			}
			e2r[j + j * n] += *alphn;
			e2i[j + j * n] = *betan;
		}
		decc_(n, e2r, e2i, ip2, ier);
	}
	return 0;
} /* decomc_ */


int slvrad_(int n, double *fac1, double *alphn, double *betan, 
	double *e1, double *e2r, double *e2i, 
	double *z1, double *z2, double *z3,
	double *f1, double *f2, double *f3,
	int *ip1, int *ip2, int sparse_LU,
	Radau_SuperLU_aux* radau_slu_aux)
{
    static int i;
    static double s2, s3;
	int ier = 0;

    for (i = 0; i < n; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }

	if (sparse_LU){
		#ifdef __RADAU5_WITH_SUPERLU
			ier = superlu_solve_d(radau_slu_aux->slu_aux_d, z1);
			if (ier != 0) {
				return ier;
			}
			ier = superlu_solve_z(radau_slu_aux->slu_aux_z, z2, z3);
		#endif /*__RADAU5_WITH_SUPERLU*/
	}else{
		sol_(n, e1, z1, ip1);
    	solc_(n, e2r, e2i, z2, z3, ip2);
	}
	return ier;
} /* slvrad_ */


int estrad_(int n, double *h__,
	double *dd1, double *dd2, double *dd3,
	FP_CB_f fcn, void* fcn_PY, int *nfcn,
	double *y0, double *y, int sparse_LU,
	double *x,double *e1,
	double *z1, double *z2, double *z3,
	double *cont, double *werr,
	double *f1, double *f2, int *ip1,
	double *scal, double *err, int *first, int *reject, 
	Radau_SuperLU_aux* radau_slu_aux, int *ier)
{
    static int i;
    static double hee1, hee2, hee3;

    hee1 = *dd1 / *h__;
    hee2 = *dd2 / *h__;
    hee3 = *dd3 / *h__;

	for (i = 0; i < n; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
	}

	if (sparse_LU){
		#ifdef __RADAU5_WITH_SUPERLU
			*ier = superlu_solve_d(radau_slu_aux->slu_aux_d, cont);
			if (*ier){
				return 0;
			}
		#endif /*__RADAU5_WITH_SUPERLU*/
	}else{
		sol_(n, e1, cont, ip1);
	}

    *err = 0.;
    for (i = 0; i < n; ++i) {
		werr[i] = cont[i] / scal[i];
		*err += werr[i] * werr[i];
    }
    *err = radau_max(sqrt(*err / n), 1e-10);

    if (*err < 1.) {
		return 0;
    }
    if (*first || *reject) {
		for (i = 0; i < n; ++i) {
			cont[i] = y[i] + cont[i];
		}
		(*fcn)(n, x, cont, f1, fcn_PY);
		++(*nfcn);
		for (i = 0; i < n; ++i) {
			cont[i] = f1[i] + f2[i];
		}

		if (sparse_LU){
			#ifdef __RADAU5_WITH_SUPERLU
				*ier = superlu_solve_d(radau_slu_aux->slu_aux_d, cont);
				if (*ier){
					return 0;
				}
			#endif /*__RADAU5_WITH_SUPERLU*/
		}else{
			sol_(n, e1, cont, ip1);
		}

		*err = 0.;
		for (i = 0; i < n; ++i) {
			werr[i] = cont[i] / scal[i];
			*err += werr[i] * werr[i];
		}
		*err = radau_max(sqrt(*err / n), 1e-10);
    }
    return 0;
} /* estrad_ */

#ifdef __RADAU5_WITH_SUPERLU
Radau_SuperLU_aux* radau_superlu_aux_setup(int n, int nnz, int nprocs, int* info){
	Radau_SuperLU_aux *radau_slu_aux = malloc(sizeof(Radau_SuperLU_aux));
	*info = 0;

	if (n < 1)     { *info = RADAU_SUPERLU_INVALID_INPUT_N;}
	if (nnz < 0)   { *info = RADAU_SUPERLU_INVALID_INPUT_NNZ;}
	if (nnz > n*n) { *info = RADAU_SUPERLU_INVALID_INPUT_NNZ_TOO_LARGE;}
	if (nprocs < 1){ *info = RADAU_SUPERLU_INVALID_INPUT_NPROC;}
	if (*info < 0) { return radau_slu_aux;}

	radau_slu_aux->n = n;
	radau_slu_aux->nnz = nnz;
	radau_slu_aux->nprocs = nprocs;

	radau_slu_aux->nnz_actual = nnz;
	radau_slu_aux->fresh_jacobian = 0;

	/* allocate memory for sparse Jacobian structure */
	radau_slu_aux->jac_data = (double*)malloc((nnz + n)*sizeof(double));
	radau_slu_aux->jac_indices = (int*)malloc((nnz + n)*sizeof(int));
	radau_slu_aux->jac_indptr = (int*)malloc((n + 1)*sizeof(int));

	/* Create auxiliary superLU structures */
	radau_slu_aux->slu_aux_d = superlu_init_d(nprocs, n, nnz);
	radau_slu_aux->slu_aux_z = superlu_init_z(nprocs, n, nnz);

	return radau_slu_aux;
}
#else
Radau_SuperLU_aux* radau_superlu_aux_setup(int n, int nnz, int nprocs, int* info){
	*info = RADAU_ERROR_SUPERLU_NOT_ENABLED;
	return NULL;
}
#endif /*__RADAU5_WITH_SUPERLU*/

#ifdef __RADAU5_WITH_SUPERLU
int radau_superlu_aux_finalize(Radau_SuperLU_aux *radau_slu_aux){
	free(radau_slu_aux->jac_data);
	free(radau_slu_aux->jac_indices);
	free(radau_slu_aux->jac_indptr);

	superlu_finalize_d(radau_slu_aux->slu_aux_d);
	superlu_finalize_z(radau_slu_aux->slu_aux_z);

	free(radau_slu_aux);
	return 0;
}
#else
int radau_superlu_aux_finalize(Radau_SuperLU_aux *radau_slu_aux){
	return 0;
}
#endif /*__RADAU5_WITH_SUPERLU*/

void *setup_radau_mem(int n){
	radau_mem_t *rmem;
	int n_sq = n*n;
	
	if (n < 1){return NULL;} /* sanity check */

	rmem = (radau_mem_t*)malloc(sizeof(radau_mem_t));
	if (!rmem){return NULL;}

	rmem->work = (double*)malloc(20*sizeof(double)); /* TODO */
	rmem->werr = (double*)malloc(n*sizeof(double));
	rmem->z1 = (double*)malloc(n*sizeof(double));
	rmem->z2 = (double*)malloc(n*sizeof(double));
	rmem->z3 = (double*)malloc(n*sizeof(double));
	rmem->y0 = (double*)malloc(n*sizeof(double));
	rmem->scal = (double*)malloc(n*sizeof(double));
	rmem->f1 = (double*)malloc(n*sizeof(double));
	rmem->f2 = (double*)malloc(n*sizeof(double));
	rmem->f3 = (double*)malloc(n*sizeof(double));
	rmem->con = (double*)malloc(4*n*sizeof(double));

	if(!rmem->work || !rmem->werr || !rmem->z1 || !rmem->z3 || !rmem->y0 || !rmem->scal || !rmem->f1 || !rmem->f2 || !rmem->f3 || !rmem->con){
		return NULL;
	}

	rmem->jac = (double*)malloc(n_sq*sizeof(double));
	rmem->e1 = (double*)malloc(n_sq*sizeof(double));
	rmem->e2r = (double*)malloc(n_sq*sizeof(double));
	rmem->e2i = (double*)malloc(n_sq*sizeof(double));

	if(!rmem->jac || !rmem->e1 || !rmem->e2r || !rmem->e2i ){
		return NULL;
	}

	rmem->iwork = (int*)malloc(20*sizeof(int)); /* TODO */
	rmem->ip1 = (int*)malloc(n*sizeof(int));
	rmem->ip2 = (int*)malloc(n*sizeof(int));

	if(!rmem->iwork || !rmem->ip1 || !rmem->ip2){
		return NULL;
	}

	return (void*)rmem;
}

void free_radau_mem(void **radau_mem){
	radau_mem_t *rmem = (radau_mem_t*) *radau_mem;
	if(!rmem){
		return;
	}

	free(rmem->work);
	free(rmem->werr);
	free(rmem->z1);
	free(rmem->z2);
	free(rmem->z3);
	free(rmem->y0);
	free(rmem->scal);
	free(rmem->f1);
	free(rmem->f2);
	free(rmem->f3);
	free(rmem->con);

	free(rmem->jac);
	free(rmem->e1);
	free(rmem->e2r);
	free(rmem->e2i);

	free(rmem->iwork);
	free(rmem->ip1);
	free(rmem->ip2);

	free(rmem);
}
