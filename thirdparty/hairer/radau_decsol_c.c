// based on f2c (version 20100827) translation of radau_decsol.f.
// Note: Due to this, matrices (doublereal*) are stored in Fortran-style column major format

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "radau_decsol_c.h"
#include <inttypes.h>

/* Common Block Declarations */

struct conra5_{
    integer nn, nn2, nn3, nn4;
    doublereal xsol, hsol, c2m1, c1m1;
} conra5_1;

struct linal_{
    integer mle, mue, mbjac, mbb, mdiag, mdiff, mbdiag;
} linal_1;

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b54 = .5;
static doublereal c_b91 = 81.;
static doublereal c_b92 = .33333333333333331;
static doublereal c_b93 = 9.;
static doublereal c_b114 = .8;
static doublereal c_b116 = .25;

/* Subroutine */ int radau5_c(integer n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *
	atol, integer *itol, FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, integer *ijac, integer *mljac, integer 
	*mujac, integer *imas, integer *mlmas, integer *mumas, FP_CB_solout 
	solout, void* solout_PY, integer *iout, doublereal *work, integer *lwork, integer *
	iwork, integer *liwork, doublereal *rpar, integer *ipar, integer *idid,
	Radau_SuperLU_aux* radau_slu_aux)
{
    /* Local variables */
	static integer iewerr;
    static integer i, m1, m2, nm1, nit, iee1, ief1, lde1, ief2, ief3, iey0, 
	    iez1, iez2, iez3;
    static doublereal facl;
    static integer ndec, njac;
    static doublereal facr, safe;
    static integer ijob, nfcn; // ijob is identifier for type of LU decomposition used
    static logical pred;
    static doublereal hmax;
    static integer nmax;
    static doublereal thet, expm;
    static integer nsol;
    static doublereal quot;
    static integer iee2i, iee2r, ieip1, ieip2, nind1, nind2, nind3;
    static doublereal quot1, quot2;
    static integer iejac, ldjac;
    static logical jband;
    static integer iecon, iemas, ldmas, ieiph;
    static logical arret;
    static doublereal fnewt;
    static integer nstep;
    static doublereal tolst;
    static integer ldmas2, iescal, naccpt;
    extern /* Subroutine */ int radcor_(integer, FP_CB_f, void*, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, FP_CB_jac, FP_CB_jac_sparse, void*, integer *, integer *,
		integer *, integer *, integer *, FP_CB_solout, void*, integer *, integer *,
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, integer *, integer *, logical *, 
	    integer *, integer *, integer *, logical *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, logical *, logical *,
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer *,
		integer *, integer *, integer *, doublereal *, integer *, doublereal *, 
		Radau_SuperLU_aux*);
    static integer nrejct;
    static logical implct;
    static integer istore;
    static logical startn;
    static doublereal uround;
/* ---------------------------------------------------------- */
/*     NUMERICAL SOLUTION OF A STIFF (OR DIFFERENTIAL ALGEBRAIC) */
/*     SYSTEM OF FIRST 0RDER ORDINARY DIFFERENTIAL EQUATIONS */
/*                     M*Y'=F(X,Y). */
/*     THE SYSTEM CAN BE (LINEARLY) IMPLICIT (MASS-MATRIX M .NE. I) */
/*     OR EXPLICIT (M=I). */
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
/*                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR) */
/*                    DOUBLE PRECISION X,Y(N),F(N) */
/*                    F(1)=...   ETC. */
/*                 RPAR, IPAR (SEE BELOW) */

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
/*                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR) */
/*                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N) */
/*                    DFY(1,1)= ... */
/*                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS */
/*                 FURNISHED BY THE CALLING PROGRAM. */
/*                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO */
/*                    BE FULL AND THE PARTIAL DERIVATIVES ARE */
/*                    STORED IN DFY AS */
/*                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J) */
/*                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND */
/*                    THE PARTIAL DERIVATIVES ARE STORED */
/*                    DIAGONAL-WISE AS */
/*                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J). */

/*     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN: */
/*                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE */
/*                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED. */
/*                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC. */

/*     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN: */
/*                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR */
/*                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. */
/*                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN */
/*                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW */
/*                       THE MAIN DIAGONAL). */

/*     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON- */
/*                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). */
/*                 NEED NOT BE DEFINED IF MLJAC=N. */

/*     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      ----- */
/*     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): - */

/*     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS- */
/*                 MATRIX M. */
/*                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY */
/*                 MATRIX AND NEEDS NOT TO BE DEFINED; */
/*                 SUPPLY A DUMMY SUBROUTINE IN THIS CASE. */
/*                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM */
/*                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR) */
/*                    DOUBLE PRECISION AM(LMAS,N) */
/*                    AM(1,1)= .... */
/*                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED */
/*                    AS FULL MATRIX LIKE */
/*                         AM(I,J) = M(I,J) */
/*                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED */
/*                    DIAGONAL-WISE AS */
/*                         AM(I-J+MUMAS+1,J) = M(I,J). */

/*     IMAS       GIVES INFORMATION ON THE MASS-MATRIX: */
/*                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY */
/*                       MATRIX, MAS IS NEVER CALLED. */
/*                    IMAS=1: MASS-MATRIX  IS SUPPLIED. */

/*     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX: */
/*                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR */
/*                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. */
/*                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE */
/*                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW */
/*                       THE MAIN DIAGONAL). */
/*                 MLMAS IS SUPPOSED TO BE .LE. MLJAC. */

/*     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON- */
/*                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). */
/*                 NEED NOT BE DEFINED IF MLMAS=N. */
/*                 MUMAS IS SUPPOSED TO BE .LE. MUJAC. */

/*     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE */
/*                 NUMERICAL SOLUTION DURING INTEGRATION. */
/*                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP. */
/*                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. */
/*                 IT MUST HAVE THE FORM */
/*                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,CONT,LRC,N, */
/*                                       RPAR,IPAR,IRTRN) */
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
/*                 "LWORK" MUST BE AT LEAST */
/*                             N*(LJAC+LMAS+3*LE+13)+20 */
/*                 WHERE */
/*                    LJAC=N              IF MLJAC=N (FULL JACOBIAN) */
/*                    LJAC=MLJAC+MUJAC+1  IF MLJAC<N (BANDED JAC.) */
/*                 AND */
/*                    LMAS=0              IF IMAS=0 */
/*                    LMAS=N              IF IMAS=1 AND MLMAS=N (FULL) */
/*                    LMAS=MLMAS+MUMAS+1  IF MLMAS<N (BANDED MASS-M.) */
/*                 AND */
/*                    LE=N               IF MLJAC=N (FULL JACOBIAN) */
/*                    LE=2*MLJAC+MUJAC+1 IF MLJAC<N (BANDED JAC.) */

/*                 IN THE USUAL CASE WHERE THE JACOBIAN IS FULL AND THE */
/*                 MASS-MATRIX IS THE INDENTITY (IMAS=0), THE MINIMUM */
/*                 STORAGE REQUIREMENT IS */
/*                             LWORK = 4*N*N+13*N+20. */
/*                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST */
/*                          N*(LJAC+13)+(N-M1)*(LMAS+3*LE)+20 */
/*                 WHERE IN THE DEFINITIONS OF LJAC, LMAS AND LE THE */
/*                 NUMBER N CAN BE REPLACED BY N-M1. */

/*     LWORK       DECLARED LENGTH OF ARRAY "WORK". */

/*     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK". */
/*                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS */
/*                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),.., */
/*                 IWORK(20) TO ZERO BEFORE CALLING. */
/*                 IWORK(21),...,IWORK(LIWORK) SERVE AS WORKING AREA. */
/*                 "LIWORK" MUST BE AT LEAST 3*N+20. */

/*     LIWORK      DECLARED LENGTH OF ARRAY "IWORK". */

/*     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH */
/*                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING */
/*                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES. */

/* ---------------------------------------------------------------------- */

/*     SOPHISTICATED SETTING OF PARAMETERS */
/*     ----------------------------------- */
/*              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK */
/*              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),... */
/*              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO. */
/*              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES: */

/*    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN */
/*              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY */
/*              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN. */
/*              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N) */
/*              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1). */

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

/*       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR */
/*       DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1. */
/*       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT */
/*       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER. */
/*       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE */
/*       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2. */

/*    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR */
/*              ODE'S THIS EQUALS THE DIMENSION OF THE SYSTEM. */
/*              DEFAULT IWORK(5)=N. */

/*    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0. */

/*    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0. */

/*    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY */
/*              IF IWORK(8).EQ.1  MOD. PREDICTIVE CONTROLLER (GUSTAFSSON) */
/*              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL */
/*              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1. */
/*              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS; */
/*              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES */
/*              OFTEN SLIGHTLY FASTER RUNS */

/*       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT */
/*            Y(I)' = Y(I+M2)   FOR  I=1,...,M1, */
/*       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME */
/*       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10). */
/*       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE */
/*       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2. */
/*       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS: */
/*       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE */
/*              JACOBIAN HAVE TO BE STORED */
/*              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL */
/*                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J) */
/*                FOR I=1,N-M1 AND J=1,N. */
/*              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM ) */
/*                 DFY(I-J+MUJAC+1,J+K*M2) = PARTIAL F(I+M1) / PARTIAL Y(J+K*M2) */
/*                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM. */
/*       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL */
/*                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM) */
/*                     PARTIAL F(I+M1) / PARTIAL Y(J+K*M2),  I,J=1,M2 */
/*                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH */
/*                    OF THESE MM+1 SUBMATRICES */
/*       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES */
/*                NEED NOT BE DEFINED IF MLJAC=N-M1 */
/*       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND */
/*              NEED NOT BE DEFINED. SUPPLY A DUMMY SUBROUTINE IN THIS CASE. */
/*              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF */
/*              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX. */
/*              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL */
/*                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1. */
/*              ELSE, THE MASS MATRIX IS BANDED */
/*                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1) */
/*       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL */
/*                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX */
/*       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX */
/*                NEED NOT BE DEFINED IF MLMAS=N-M1 */

/*    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0. */

/*    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1. */

/*    IWORK(11) SWITCH FOR USAGE OF SPARSE LINEAR SOLVER (SUPERLU),
				CURRENTLY ONLY COMPATIBLE WITH B = I OR DIAGONAL B, I.E. 
				MLMAS = MUMAS = 0 */
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
/*     ----------------- */
/*     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED */
/*                 (AFTER SUCCESSFUL RETURN X=XEND). */

/*     Y(N)        NUMERICAL SOLUTION AT X */

/*     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP */

/*     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN: */
/*                   IDID= 1  COMPUTATION SUCCESSFUL, */
/*                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT) */
/*                   IDID=-1  INPUT IS NOT CONSISTENT, */
/*                   IDID=-2  LARGER NMAX IS NEEDED, */
/*                   IDID=-3  STEP SIZE BECOMES TOO SMALL, */
/*                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR. */

/*   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL */
/*                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED) */
/*   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY */
/*                      OR NUMERICALLY) */
/*   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS */
/*   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS */
/*   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST), */
/*                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED) */
/*   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF BOTH MATRICES */
/*   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF BOTH */
/*                      SYSTEMS; THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS, */
/*                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED */
/* ----------------------------------------------------------------------- */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*          DECLARATIONS */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/* *** *** *** *** *** *** *** */
/*        SETTING THE PARAMETERS */
/* *** *** *** *** *** *** *** */
    /* Parameter adjustments */
    --y;
    --rtol;
    --atol;
    --work;
    --iwork;
    --rpar;
    --ipar;

    /* Function Body */
    nfcn = 0;
    njac = 0;
    nstep = 0;
    naccpt = 0;
    nrejct = 0;
    ndec = 0;
    nsol = 0;
    arret = FALSE_;
	/* -------- UROUND   SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0 */
    if (work[1] == 0.) {
		uround = 1e-16;
    } else {
		uround = work[1];
		if (uround <= 1e-19 || uround >= 1.) {
			printf(" COEFFICIENTS HAVE 20 DIGITS, UROUND= \t %e \n", uround);
			arret = TRUE_;
		}
    }
	/* -------- CHECK AND CHANGE THE TOLERANCES */
    expm = .66666666666666663;
    if (*itol == 0) {
		if (atol[1] <= 0. || rtol[1] <= uround * 10.) {
			printf(" TOLERANCES ARE TOO SMALL \n");
			arret = TRUE_;
		} else {
			quot = atol[1] / rtol[1];
			rtol[1] = pow(rtol[1], expm) * .1;
			atol[1] = rtol[1] * quot;
		}
    } else {
		for (i = 1; i <= n; ++i) {
			if (atol[i] <= 0. || rtol[i] <= uround * 10.) {
				printf("TOLERANCES (%i) ARE TOO SMALL \n", i);
				arret = TRUE_;
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
			arret = TRUE_;
		}
    }
	/* -------- NIT, MAXIMAL NUMBER OF NEWTON ITERATIONS */
    if (iwork[3] == 0) {
		nit = 7;
    } else {
		nit = iwork[3];
		if (nit <= 0) {
			printf("CURIOUS INPUT IWORK(3)= %i \n", nit);
			arret = TRUE_;
		}
    }
	/* -------- STARTN SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS */
    if (iwork[4] == 0) {
		startn = FALSE_;
    } else {
		startn = TRUE_;
    }
	/* -------- PARAMETER FOR DIFFERENTIAL-ALGEBRAIC COMPONENTS */
    nind1 = iwork[5];
    nind2 = iwork[6];
    nind3 = iwork[7];
    if (nind1 == 0) {
		nind1 = n;
    }
    if (nind1 + nind2 + nind3 != n) {
		printf("CURIOUS INPUT FOR IWORK(5,6,7)= \t %i\t %i\t %i\n", nind1, nind2, nind3);
		arret = TRUE_;
    }
	/* -------- PRED STEP SIZE CONTROL */
    if (iwork[8] <= 1) {
		pred = TRUE_;
    } else {
		pred = FALSE_;
    }
	/* -------- PARAMETER FOR SECOND ORDER EQUATIONS */
    m1 = iwork[9];
    m2 = iwork[10];
    nm1 = n - m1;
    if (m1 == 0) {
		m2 = n;
    }
    if (m2 == 0) {
		m2 = m1;
    }
    if (m1 < 0 || m2 < 0 || m1 + m2 > n) {
		printf("CURIOUS INPUT FOR IWORK(9,10)= \t %i\t %i\n", m1, m2);
		arret = TRUE_;
    }
	/* --------- SAFE, SAFETY FACTOR IN STEP SIZE PREDICTION */
    if (work[2] == 0.) {
		safe = .9;
    } else {
		safe = work[2];
		if (safe <= .001 || safe >= 1.) {
			printf("CURIOUS INPUT FOR WORK(2)= %f \n", safe);
			arret = TRUE_;
		}
    }
	/* ------ THET, DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
    if (work[3] == 0.) {
		thet = .001;
    } else {
		thet = work[3];
		if (thet >= 1.) {
			printf("CURIOUS INPUT FOR WORK(3)= %f \n", thet);
			arret = TRUE_;
		}
    }
	/* --- FNEWT, STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. */
    tolst = rtol[1];
    if (work[4] == 0.) {
		fnewt = max(uround * 10 / tolst, min(.03, pow(tolst, c_b54)));
    } else {
		fnewt = work[4];
		if (fnewt <= uround / tolst) {
			printf("CURIOUS INPUT FOR WORK(4)= %f \n", fnewt);
			arret = TRUE_;
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
		arret = TRUE_;
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
		arret = TRUE_;
    }
	/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
	/*         COMPUTATION OF ARRAY ENTRIES */
	/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
	/* ---- IMPLICIT, BANDED OR NOT ? */
    implct = *imas != 0;
    jband = *mljac < nm1;
	/* -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- */
	/* -- JACOBIAN AND MATRICES E1, E2 */
    if (jband) {
		ldjac = *mljac + *mujac + 1;
		lde1 = *mljac + ldjac;
    } else {
		*mljac = nm1;
		*mujac = nm1;
		ldjac = nm1;
		lde1 = nm1;
    }
	/* -- MASS MATRIX */
    if (implct) {
		if (*mlmas != nm1) {
			ldmas = *mlmas + *mumas + 1;
			if (jband) {
				ijob = 4;
			} else {
				ijob = 3;
			}
		} else {
			*mumas = nm1;
			ldmas = nm1;
			ijob = 5;
		}
	/* ------ BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC" */
		if (*mlmas > *mljac || *mumas > *mujac) {
			printf("BANDWITH OF \"MAS\" NOT SMALLER THAN BANDWITH OF \"JAC\"\n");
			arret = TRUE_;
		}
    } else {
		ldmas = 0;
		if (jband) {
			ijob = 2;
		} else {
			ijob = 1;
			if (n > 2 && iwork[1] != 0) {
				ijob = 7;
			}
		}
    }
    ldmas2 = max(1,ldmas);
	/* ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN */
    if ((implct || jband) && ijob == 7) {
		printf(" HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS WITH FULL JACOBIAN\n");
		arret = TRUE_;
    }
	/* -------- SPARSE LU DECOMPOSITION OPTIONS */
	if (iwork[11]){ // SPARSE LU
		ijob = 8;
		if (implct){
			printf("SPARSE LU DECOMPOSITION IS NOT SUPPORTED FOR IMPLICIT ODES\n");
			arret = TRUE_;
		}
		if (!(*ijac)){
			printf("CURIOUS INPUT; ANALYTICAL JACOBIAN DISABLED, IJAC = %i, WHICH IS REQUIRED FOR SPARSE SOLVER\n", *ijac);
			arret = TRUE_;
		}
	}
	/* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
	iewerr = 21;
    iez1 = iewerr + n;
    iez2 = iez1 + n;
    iez3 = iez2 + n;
    iey0 = iez3 + n;
    iescal = iey0 + n;
    ief1 = iescal + n;
    ief2 = ief1 + n;
    ief3 = ief2 + n;
    iecon = ief3 + n;
    iejac = iecon + (n << 2);
    iemas = iejac + n * ldjac;
    iee1 = iemas + nm1 * ldmas;
    iee2r = iee1 + nm1 * lde1;
    iee2i = iee2r + nm1 * lde1;
	/* ------ TOTAL STORAGE REQUIREMENT ----------- */
    istore = iee2i + nm1 * lde1 - 1;
    if (istore > *lwork) {
		printf("INSUFFICIENT STORAGE FOR WORK, MIN. LWORK= %i\n", istore);
		arret = TRUE_;
    }
	/* ------- ENTRY POINTS FOR INTEGER WORKSPACE ----- */
    ieip1 = 21;
    ieip2 = ieip1 + nm1;
    ieiph = ieip2 + nm1;
	/* --------- TOTAL REQUIREMENT --------------- */
    istore = ieiph + nm1 - 1;
    if (istore > *liwork) {
		printf("INSUFF. STORAGE FOR IWORK, MIN. LIWORK= %i\n", istore);
		arret = TRUE_;
    }
	/* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
    if (arret) {
		*idid = -1;
		return 0;
    }
	/* -------- CALL TO CORE INTEGRATOR ------------ */
	radcor_(n, (FP_CB_f)fcn, fcn_PY, x, &y[1], xend, &hmax, h__, &rtol[1], &atol[1], 
		itol, (FP_CB_jac)jac, (FP_CB_jac_sparse) jac_sparse, jac_PY, ijac, mljac, mujac, mlmas, mumas,
		(FP_CB_solout)solout, solout_PY, iout, idid, &nmax, &uround, &safe, &thet, &fnewt,
		&quot1, &quot2, &nit, &ijob, &startn, &nind1, &nind2, &nind3,
		&pred, &facl, &facr, &m1, &m2, &nm1, &implct, &jband, &ldjac,
		&lde1, &ldmas2, &work[iez1], &work[iez2], &work[iez3], &work[iey0],
		&work[iescal], &work[ief1], &work[ief2], &work[ief3], &work[iejac],
		&work[iee1], &work[iee2r], &work[iee2i], &work[iemas],
		&iwork[ieip1], &iwork[ieip2], &iwork[ieiph], &work[iecon], &nfcn,
		&njac, &nstep, &naccpt, &nrejct, &ndec, &nsol, &rpar[1], &ipar[1],
		&work[iewerr], radau_slu_aux);
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
	/* ----------- RETURN ----------- */
    return 0;
} /* radau5_ */


/*     END OF SUBROUTINE RADAU5 */

/* *********************************************************** */

/* Subroutine */ int radcor_(integer n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *hmax, doublereal *h__, doublereal *
	rtol, doublereal *atol, integer *itol, FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, integer *ijac, 
	integer *mljac, integer *mujac, integer *mlmas, integer *
	mumas, FP_CB_solout solout, void* solout_PY, integer *iout, integer *idid, integer *nmax, 
	doublereal *uround, doublereal *safe, doublereal *thet, doublereal *
	fnewt, doublereal *quot1, doublereal *quot2, integer *nit, integer *
	ijob, logical *startn, integer *nind1, integer *nind2, integer *nind3,
	 logical *pred, doublereal *facl, doublereal *facr, integer *m1, 
	integer *m2, integer *nm1, logical *implct, logical *banded, integer *
	ldjac, integer *lde1, integer *ldmas, doublereal *z1, doublereal *z2, 
	doublereal *z3, doublereal *y0, doublereal *scal, doublereal *f1, 
	doublereal *f2, doublereal *f3, doublereal *fjac, doublereal *e1, 
	doublereal *e2r, doublereal *e2i, doublereal *fmas, integer *ip1, 
	integer *ip2, integer *iphes, doublereal *cont, integer *nfcn, 
	integer *njac, integer *nstep, integer *naccpt, integer *nrejct, 
	integer *ndec, integer *nsol, doublereal *rpar, integer *ipar,
	doublereal *werr, Radau_SuperLU_aux* radau_slu_aux)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, e2r_dim1, e2r_offset, e2i_dim1, e2i_offset;
    doublereal d__1;

    /* Local variables */
    static integer i, j, k, l;
    static doublereal a1, a2, c1, c2, a3;
    static integer j1, n2, n3;
    static doublereal u1;
    static integer nunexpect;
    static doublereal ak;
    static integer md;
    static doublereal t11, t12, t13, t21, t22, t23, t31;
    static integer mm;
    static doublereal qt, dd1, dd2, dd3, ak1, ak2, ak3, f1i, f2i, f3i, c1q, 
	    c2q, c3q, z1i, z2i, z3i, sq6, fac, ti11, cno;
    static integer lrc;
    static doublereal ti12, ti13, ti21, ti22, ti23, ti31, ti32, ti33;
    static integer ier;
    static doublereal xph, thq, err, fac1, cfac, hacc, c1mc2, beta;
    static integer lbeg;
    static doublereal alph, hold;
    static integer lend;
    static doublereal delt, hnew;
    static logical last;
    static doublereal hopt, xold;
    static integer newt;
    static doublereal dyno, dyth, quot, hhfac, betan, alphn, theta, 
	    ysafe, hmaxn;
    static integer nsing;
    static logical first;
    static integer irtrn, nrsol, nsolu;
    static doublereal qnewt, xosol, acont3;
	static logical index2, index3, caljac;
    static doublereal faccon;
    extern /* Subroutine */ int decomc_(integer, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer *,
	    integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *,
		Radau_SuperLU_aux*);
    static logical calhes;
    static doublereal erracc;
    static integer mujacj;
    extern /* Subroutine */ int decomr_(integer, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer *,
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, logical *, integer *,
		Radau_SuperLU_aux*);
    static logical reject;
    static doublereal facgus;
    static integer mujacp;
    extern /* Subroutine */ int estrad_(integer, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer *,
	    doublereal *, doublereal *, doublereal *, doublereal *, FP_CB_f, void*, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, logical *, logical *, doublereal *, 
	    doublereal *, integer *, Radau_SuperLU_aux*, integer *);
    static doublereal dynold, posneg;
    extern /* Subroutine */ int slvrad_(integer, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer *,
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, Radau_SuperLU_aux*);
    static doublereal thqold;

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
    --iphes;
    --ip2;
    --ip1;
    --werr;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    e2i_dim1 = *lde1;
    e2i_offset = 1 + e2i_dim1;
    e2i -= e2i_offset;
    e2r_dim1 = *lde1;
    e2r_offset = 1 + e2r_dim1;
    e2r -= e2r_offset;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    --rpar;
    --ipar;

    /* Function Body */
    conra5_1.nn = n;
    conra5_1.nn2 = n << 1;
    conra5_1.nn3 = n * 3;
    lrc = n << 2;
	/* -------- CHECK THE INDEX OF THE PROBLEM ----- */
    // index1 = *nind1 != 0;
    index2 = *nind2 != 0;
    index3 = *nind3 != 0;
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
    u1 = (pow(c_b91, c_b92) + 6. - pow(c_b93, c_b92)) / 30.;
    alph = (12. - pow(c_b91, c_b92) + pow(c_b93, c_b92)) / 60.;
    beta = (pow(c_b91, c_b92) + pow(c_b93, c_b92)) * sqrt(3.) / 60.;
    cno = alph * alph + beta * beta;
    u1 = 1. / u1;
    alph /= cno;
    beta /= cno;
    t11 = .091232394870892942792;
    t12 = -.14125529502095420843;
    t13 = -.030029194105147424492;
    t21 = .24171793270710701896;
    t22 = .20412935229379993199;
    t23 = .38294211275726193779;
    t31 = .96604818261509293619;
    ti11 = 4.325579890063155351;
    ti12 = .33919925181580986954;
    ti13 = .54177053993587487119;
    ti21 = -4.1787185915519047273;
    ti22 = -.32768282076106238708;
    ti23 = .47662355450055045196;
    ti31 = -.50287263494578687595;
    ti32 = 2.5719269498556054292;
    ti33 = -.59603920482822492497;
    if (*m1 > 0) {
		*ijob += 10;
    }
    posneg = copysign(1., *xend - *x);
    hmaxn = min(abs(*hmax), abs(*xend - *x));
    if (abs(*h__) <= *uround * 10.) {
		*h__ = 1e-6;
    }
    *h__ = min(abs(*h__), hmaxn);
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
		irtrn = 1;
		nrsol = 1;
		xosol = xold;
		conra5_1.xsol = *x;
		for (i = 1; i <= n; ++i) {
			werr[i] = 0.;
			cont[i] = y[i];
		}
		nsolu = n;
		conra5_1.hsol = hold;
		(*solout)(&nrsol, &xosol, &conra5_1.xsol, &y[1], &cont[1], &werr[1], &
			lrc, &nsolu, &rpar[1], &ipar[1], &irtrn, solout_PY);
		if (irtrn < 0) {
			goto L179;
		}
    }
    linal_1.mle = *mljac;
    linal_1.mue = *mujac;
    linal_1.mbjac = *mljac + *mujac + 1;
    linal_1.mbb = *mlmas + *mumas + 1;
    linal_1.mdiag = linal_1.mle + linal_1.mue + 1;
    linal_1.mdiff = linal_1.mle + linal_1.mue - *mumas;
    linal_1.mbdiag = *mumas + 1;
    n2 = n << 1;
    n3 = n * 3;
    if (*itol == 0) {
		for (i = 1; i <= n; ++i) {
			scal[i] = atol[1] + rtol[1] * abs(y[i]);
		}
    } else {
		for (i = 1; i <= n; ++i) {
			scal[i] = atol[i] + rtol[i] * abs(y[i]);
		}
    }
    hhfac = *h__;
    (*fcn)(n, x, &y[1], &y0[1], &rpar[1], &ipar[1], fcn_PY);
	if (ipar[1] < 0){
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
		if (*banded) {
		/* --- JACOBIAN IS BANDED */
			mujacp = *mujac + 1;
			md = min(linal_1.mbjac,*m2);
			for (mm = 1; mm <= *m1 / *m2 + 1; ++mm) {
				for (k = 1; k <= md; ++k) {
					j = k + (mm - 1) * *m2;
L12:
					f1[j] = y[j];
					f2[j] = sqrt(*uround * max(1e-5, abs(y[j])));
					y[j] += f2[j];
					j += md;
					if (j <= mm * *m2) {
						goto L12;
					}
					(*fcn)(n, x, &y[1], &cont[1], &rpar[1], &ipar[1], fcn_PY);
					j = k + (mm - 1) * *m2;
					j1 = k;
					lbeg = max(1, j1 - *mujac) + *m1;
L14:
					lend = min(*m2, j1 + *mljac) + *m1;
					y[j] = f1[j];
					mujacj = mujacp - j1 - *m1;
					for (l = lbeg; l <= lend; ++l) {
						fjac[l + mujacj + j * fjac_dim1] = (cont[l] - y0[l]) / f2[j];
					}
					j += md;
					j1 += md;
					lbeg = lend + 1;
					if (j <= mm * *m2) {
						goto L14;
					}
				}
			}
		} else {
			/* --- JACOBIAN IS FULL */
			for (i = 1; i <= n; ++i) {
				ysafe = y[i];
				delt = sqrt(*uround * max(1e-5, abs(ysafe)));
				y[i] = ysafe + delt;
				(*fcn)(n, x, &y[1], &cont[1], &rpar[1], &ipar[1], fcn_PY);
				if (ipar[1] < 0) {
					y[i] = ysafe - delt;
					(*fcn)(n, x, &y[1], &cont[1], &rpar[1], &ipar[1], fcn_PY);
					if (ipar[1] < 0) {
						y[i] = ysafe;
						goto L79;
					}
					for (j = *m1 + 1; j <= n; ++j) {
						fjac[j - *m1 + i * fjac_dim1] = (y0[j] - cont[j]) / delt;
					}
				} else {
					for (j = *m1 + 1; j <= n; ++j) {
						fjac[j - *m1 + i * fjac_dim1] = (cont[j] - y0[j]) / delt;
					}
				}
				y[i] = ysafe;
			}
		}
    } else {
		/* --- COMPUTE JACOBIAN MATRIX ANALYTICALLY */
		if (*ijob == 8){ // Sparse LU
			radau_slu_aux->nnz_actual = radau_slu_aux->nnz;
			(*jac_sparse)(n, x, &y[1], &(radau_slu_aux->nnz_actual), radau_slu_aux->jac_data, radau_slu_aux->jac_indices, radau_slu_aux->jac_indptr, &rpar[1], &ipar[1], jac_PY);
			if (ipar[1] < 0){
				goto L79;
			}
			ier = sparse_csc_add_diagonal(n, &(radau_slu_aux->nnz_actual), radau_slu_aux->jac_data, radau_slu_aux->jac_indices, radau_slu_aux->jac_indptr);
			if (ier != 0){
				goto L183;
			}
			radau_slu_aux->fresh_jacobian = 1;
		} else { // dense jacobian
			(*jac)(n, x, &y[1], &fjac[fjac_offset], ldjac, &rpar[1], &ipar[1], jac_PY);
			if (ipar[1] < 0){
				goto L79;
			}
		}
    }
    caljac = TRUE_;
    calhes = TRUE_;
/* --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS */
L20:
    fac1 = u1 / *h__;
    alphn = alph / *h__;
    betan = beta / *h__;
	decomr_(n, &fjac[fjac_offset], ldjac, &fmas[fmas_offset], ldmas, mlmas, 
			mumas, m1, m2, nm1, &fac1, &e1[e1_offset], lde1, &ip1[1], &ier, 
			ijob, &calhes, &iphes[1], radau_slu_aux);
    if (ier) {
		goto L185;
    }
	decomc_(n, &fjac[fjac_offset], ldjac, &fmas[fmas_offset], ldmas, mlmas, 
			mumas, m1, m2, nm1, &alphn, &betan, &e2r[e2r_offset], &e2i[e2i_offset],
			lde1, &ip2[1], &ier, ijob, radau_slu_aux);
    if (ier) {
		goto L185;
    }
	if (*ijob == 8){ // Sparse LU
		radau_slu_aux->fresh_jacobian = 0; // has once been used to create a decomposition now
	}
    ++(*ndec);
/* --- COMPUTE STEPSIZE */
L30:
    ++(*nstep);
    if (*nstep > *nmax) {
		goto L178;
    }
    if (abs(*h__) * .1 <= abs(*x) * *uround) {
		goto L177;
    }
    if (index2) {
		for (i = *nind1 + 1; i <= *nind1 + *nind2; ++i) {
			scal[i] /= hhfac;
		}
    }
    if (index3) {
		for (i = *nind1 + *nind2 + 1; i <= *nind1 + *nind2 + *nind3; ++i) {
			scal[i] /= hhfac * hhfac;
		}
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
    faccon = pow(max(faccon,*uround), c_b114);
    theta = abs(*thet);
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
	(*fcn)(n, &d__1, &cont[1], &z1[1], &rpar[1], &ipar[1], fcn_PY);
    ++(*nfcn);
    if (ipar[1] < 0) {
		goto L79;
    }
    for (i = 1; i <= n; ++i) {
		cont[i] = y[i] + z2[i];
    }
    d__1 = *x + c2 * *h__;
	(*fcn)(n, &d__1, &cont[1], &z2[1], &rpar[1], &ipar[1], fcn_PY);
    ++(*nfcn);
    if (ipar[1] < 0) {
		goto L79;
    }
    for (i = 1; i <= n; ++i) {
		cont[i] = y[i] + z3[i];
    }
	(*fcn)(n, &xph, &cont[1], &z3[1], &rpar[1], &ipar[1], fcn_PY);
    ++(*nfcn);
    if (ipar[1] < 0) {
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
	slvrad_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset], 
			ldmas, mlmas, mumas, m1, m2, nm1, &fac1, &alphn, &betan, &e1[e1_offset],
			&e2r[e2r_offset], &e2i[e2i_offset], lde1,
			&z1[1], &z2[1], &z3[1], &f1[1], &f2[1], &f3[1], &cont[1], &ip1[1], &ip2[1],
			&iphes[1], &ier, ijob, radau_slu_aux);
	if (ier){
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
				qnewt = max(1e-4, min(20., dyth));
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
    dynold = max(dyno, *uround);
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
	estrad_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset],  
			ldmas, mlmas, mumas, h__, &dd1, &dd2, &dd3, (FP_CB_f) fcn, fcn_PY, nfcn,
			&y0[1], &y[1], ijob, x, m1, m2, nm1, &e1[e1_offset], lde1,
			&z1[1], &z2[1], &z3[1], &cont[1], &werr[1], &f1[1], &f2[1], &ip1[1], 
			&iphes[1], &scal[1], &err, &first, &reject, &fac1, &rpar[1], &ipar[1],
			radau_slu_aux, &ier);
	if (ier){
		goto L184;
	}
	/* --- COMPUTATION OF HNEW */
	/* --- WE REQUIRE .2<=HNEW/H<=8. */
    fac = min(*safe, cfac / (newt + (*nit << 1)));
    quot = max(*facr ,min(*facl, pow(err, c_b116) / fac));
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
				facgus = max(*facr, min(*facl, facgus));
				quot = max(quot,facgus);
				hnew = *h__ / quot;
			}
			hacc = *h__;
			erracc = max(.01, err);
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
				scal[i] = atol[1] + rtol[1] * abs(y[i]);
			}
		} else {
			for (i = 1; i <= n; ++i) {
				scal[i] = atol[i] + rtol[i] * abs(y[i]);
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
			(*solout)(&nrsol, &xosol, &conra5_1.xsol, &y[1], &cont[1], &werr[1],
			 		  &lrc, &nsolu, &rpar[1], &ipar[1], &irtrn, solout_PY);
			if (irtrn < 0) {
				goto L179;
			}
		}
		caljac = FALSE_;
		if (last) {
			*h__ = hopt;
			*idid = 1;
			goto L181;
		}
		(*fcn)(n, x, &y[1], &y0[1], &rpar[1], &ipar[1], fcn_PY);
		++(*nfcn);
		hnew = posneg * min(abs(hnew), hmaxn);
		hopt = hnew;
		hopt = min(*h__,hnew);
		if (reject) {
			hnew = posneg * min(abs(hnew), abs(*h__));
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
	if (ipar[1] < RADAU_CALLBACK_ERROR_INVALID_NNZ){
		printf("FAILURE IN JACOBIAN EVALUATIONS, NNZ TOO SMALL, SPECIFIED NNZ: %i, ACTUAL: %i \n", radau_slu_aux->nnz, -(ipar[1] - RADAU_CALLBACK_ERROR_INVALID_NNZ));
		*idid = -6;
		goto L181;
	}

	switch(ipar[1]){
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
    *idid = -5;
	goto L181;
/* --- FAIL EXIT, REPEATED SINGULAR JACOBIAN*/
L176:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("MATRIX IS REPEATEDLY SINGULAR IER= %i\n", ier);
    *idid = -4;
	goto L181;
/* --- FAIL EXIT, STEP SIZE TOO SMALL*/
L177:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("STEP SIZE TOO SMALL, H= %e\n", *h__);
    *idid = -3;
	goto L181;
/* --- FAIL EXIT, NMAX EXCEEDED*/
L178:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("MORE THAN NMAX = %i STEPS ARE NEEDED\n", *nmax);
    *idid = -2;
	goto L181;
/* --- FAILURE EXIT, ERROR IN SPARSE JACOBIAN CALLBACK*/
L182:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("JACOBIAN GIVEN IN WRONG FORMAT, REQUIRED SPARSE FORMAT: CSC\n");
	*idid = -7;
	goto L181;
L183:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("UNEXPECTED MALLOC FAILURE\n");
	*idid = -9;
	goto L181;
/* --- FAIL EXIT, UNEXPECTED SUPERLU FAILURE*/
L184:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("UNEXPECTED FAILURE OF SUPERLU FUNCTION CALL, ier = %i \n", ier);
	*idid = -8;
	goto L181;
/* --- FAIL EXIT, OTHER SUPERLU FAILURE*/
L185:
	if (ier < 0){ // incorrect input to function call
		goto L184;
	}
	if (ier <= n){ // factorization singular
		goto L78;
	}else{ // memory allocation failure
		printf("EXIT OF RADAU5 AT X = %e \n", *x);
		printf("SUPERLU MEMORY ALLOCATION FAILURE, NUMBER OF BYTES ALLOCATED AT POINT OF FAILURE: %i.\n", ier - n);
		*idid = -9;
		goto L181;
	}
/* --- FAIL EXIT, UNRECOVERABLE EXCEPTION IN FCN OR JAC CALLBACK*/
L186:
	printf("EXIT OF RADAU5 AT X = %e \n", *x);
	printf("UNRECOVERABLE EXCEPTION ENCOUNTERED IN PROBLEM CALLBACK\n");
	*idid = -10;
	return 0;
/* --- EXIT CAUSED BY SOLOUT */
L179:
/*      WRITE(6,979)X */
    *idid = 2;
	goto L181;

/* --- EXIT OF RADCOR */
L181:
	return 0;
} /* radcor_ */


/*     END OF SUBROUTINE RADCOR */

/* *********************************************************** */

doublereal contr5_c(integer *i, doublereal *x, doublereal *cont, integer *lrc)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal s;

/* ---------------------------------------------------------- */
/*     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN */
/*     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X. */
/*     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR */
/*     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU5). */
/* ---------------------------------------------------------- */
    /* Parameter adjustments */
    --cont;

    /* Function Body */
    s = (*x - conra5_1.xsol) / conra5_1.hsol;
    ret_val = cont[*i] + s * (cont[*i + conra5_1.nn] + (s - conra5_1.c2m1)
	     * (cont[*i + conra5_1.nn2] + (s - conra5_1.c1m1) * cont[*i + 
	    conra5_1.nn3]));
    return ret_val;
} /* contr5_ */


/*     END OF FUNCTION CONTR5 */

/* *********************************************************** */
/* Subroutine */ int dec_(integer n, integer *ndim, doublereal *a, integer *
	ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    static integer i, j, k, m;
    static doublereal t;
    static integer nm1, kp1;

/* VERSION REAL DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION. */
/*  INPUT.. */
/*     N = ORDER OF MATRIX. */
/*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
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
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
    ip[n] = 1;
    if (n == 1) {
		goto L70;
    }
    nm1 = n - 1;
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = k;
		for (i = kp1; i <= n; ++i) {
			if (abs(a[i + k * a_dim1]) > abs(a[m + k * a_dim1])) {
			m = i;
			}
	/* L10: */
		}
		ip[k] = m;
		t = a[m + k * a_dim1];
		if (m == k) {
			goto L20;
		}
		ip[n] = -ip[n];
		a[m + k * a_dim1] = a[k + k * a_dim1];
		a[k + k * a_dim1] = t;
L20:
		if (t == 0.) {
			goto L80;
		}
		t = 1. / t;
		for (i = kp1; i <= n; ++i) {
			/* L30: */
			a[i + k * a_dim1] = -a[i + k * a_dim1] * t;
		}
		for (j = kp1; j <= n; ++j) {
			t = a[m + j * a_dim1];
			a[m + j * a_dim1] = a[k + j * a_dim1];
			a[k + j * a_dim1] = t;
			if (t == 0.) {
				goto L45;
			}
			for (i = kp1; i <= n; ++i) {
				/* L40: */
				a[i + j * a_dim1] += a[i + k * a_dim1] * t;
			}
L45:
	/* L50: */
			;
		}
	/* L60: */
    }
L70:
    k = n;
    if (a[n + n * a_dim1] == 0.) {
		goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DEC ------------------------- */
} /* dec_ */

/* Subroutine */ int sol_(integer n, integer *ndim, doublereal *a, 
	doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    static integer i, k, m;
    static doublereal t;
    static integer kb, km1, nm1, kp1;

/* VERSION REAL DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX. */
/*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
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
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (n == 1) {
		goto L50;
    }
    nm1 = n - 1;
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = ip[k];
		t = b[m];
		b[m] = b[k];
		b[k] = t;
		for (i = kp1; i <= n; ++i) {
			/* L10: */
			b[i] += a[i + k * a_dim1] * t;
		}
		/* L20: */
    }
    for (kb = 1; kb <= nm1; ++kb) {
		km1 = n - kb;
		k = km1 + 1;
		b[k] /= a[k + k * a_dim1];
		t = -b[k];
		for (i = 1; i <= km1; ++i) {
			/* L30: */
			b[i] += a[i + k * a_dim1] * t;
		}
		/* L40: */
    }
L50:
    b[1] /= a[a_dim1 + 1];
    return 0;
/* ----------------------- END OF SUBROUTINE SOL ------------------------- */
} /* sol_ */

/* Subroutine */ int dech_(integer n, integer *ndim, doublereal *a, integer *
	lb, integer *ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    static integer i, j, k, m;
    static doublereal t;
    static integer na, nm1, kp1;

/* VERSION REAL DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A HESSENBERG */
/*  MATRIX WITH LOWER BANDWIDTH LB */
/*  INPUT.. */
/*     N = ORDER OF MATRIX A. */
/*     NDIM = DECLARED DIMENSION OF ARRAY  A . */
/*     A = MATRIX TO BE TRIANGULARIZED. */
/*     LB = LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED, LB.GE.1). */
/*  OUTPUT.. */
/*     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U . */
/*     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW. */
/*     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE */
/*           SINGULAR AT STAGE K. */
/*  USE  SOLH  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N). */
/*  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     THIS IS A SLIGHT MODIFICATION OF */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
    ip[n] = 1;
    if (n == 1) {
		goto L70;
    }
    nm1 = n - 1;
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = k;
		na = min(n, *lb + k);
		for (i = kp1; i <= na; ++i) {
			if (abs(a[i + k * a_dim1]) > abs(a[m + k * a_dim1])) {
				m = i;
			}
		/* L10: */
		}
		ip[k] = m;
		t = a[m + k * a_dim1];
		if (m == k) {
			goto L20;
		}
		ip[n] = -ip[n];
		a[m + k * a_dim1] = a[k + k * a_dim1];
		a[k + k * a_dim1] = t;
L20:
		if (t == 0.) {
			goto L80;
		}
		t = 1. / t;
		for (i = kp1; i <= na; ++i) {
			/* L30: */
			a[i + k * a_dim1] = -a[i + k * a_dim1] * t;
		}
		for (j = kp1; j <= n; ++j) {
			t = a[m + j * a_dim1];
			a[m + j * a_dim1] = a[k + j * a_dim1];
			a[k + j * a_dim1] = t;
			if (t == 0.) {
				goto L45;
			}
			for (i = kp1; i <= na; ++i) {
				/* L40: */
				a[i + j * a_dim1] += a[i + k * a_dim1] * t;
			}
L45:
	/* L50: */
			;
		}
	/* L60: */
    }
L70:
    k = n;
    if (a[n + n * a_dim1] == 0.) {
		goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECH ------------------------ */
} /* dech_ */

/* Subroutine */ int solh_(integer n, integer *ndim, doublereal *a, integer *
	lb, doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    static integer i, k, m;
    static doublereal t;
    static integer kb, na, km1, nm1, kp1;

/* VERSION REAL DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX A. */
/*    NDIM = DECLARED DIMENSION OF ARRAY  A . */
/*    A = TRIANGULARIZED MATRIX OBTAINED FROM DECH. */
/*    LB = LOWER BANDWIDTH OF A. */
/*    B = RIGHT HAND SIDE VECTOR. */
/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
/*  DO NOT USE IF DECH HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    B = SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    --b;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (n == 1) {
		goto L50;
    }
    nm1 = n - 1;
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = ip[k];
		t = b[m];
		b[m] = b[k];
		b[k] = t;
		na = min(n, *lb + k);
		for (i = kp1; i <= na; ++i) {
			/* L10: */
			b[i] += a[i + k * a_dim1] * t;
		}
		/* L20: */
    }
    for (kb = 1; kb <= nm1; ++kb) {
		km1 = n - kb;
		k = km1 + 1;
		b[k] /= a[k + k * a_dim1];
		t = -b[k];
		for (i = 1; i <= km1; ++i) {
			/* L30: */
			b[i] += a[i + k * a_dim1] * t;
		}
		/* L40: */
    }
L50:
    b[1] /= a[a_dim1 + 1];
    return 0;
/* ----------------------- END OF SUBROUTINE SOLH ------------------------ */
} /* solh_ */

/* Subroutine */ int decc_(integer n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *ip, integer *ier)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset;

    /* Local variables */
    static integer i, j, k, m;
    static doublereal ti, tr;
    static integer nm1, kp1;
    static doublereal den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
/*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
/*  INPUT.. */
/*     N = ORDER OF MATRIX. */
/*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
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
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    *ier = 0;
    ip[n] = 1;
    if (n == 1) {
		goto L70;
    }
    nm1 = n - 1;
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = k;
		for (i = kp1; i <= n; ++i) {
			if (abs(ar[i + k * ar_dim1]) + abs(ai[i + k * ai_dim1]) > abs(ar[m + k * ar_dim1]) + abs(ai[m + k * ai_dim1])) {
				m = i;
			}
		/* L10: */
		}
		ip[k] = m;
		tr = ar[m + k * ar_dim1];
		ti = ai[m + k * ai_dim1];
		if (m == k) {
			goto L20;
		}
		ip[n] = -ip[n];
		ar[m + k * ar_dim1] = ar[k + k * ar_dim1];
		ai[m + k * ai_dim1] = ai[k + k * ai_dim1];
		ar[k + k * ar_dim1] = tr;
		ai[k + k * ai_dim1] = ti;
L20:
		if (abs(tr) + abs(ti) == 0.) {
			goto L80;
		}
		den = tr * tr + ti * ti;
		tr /= den;
		ti = -ti / den;
		for (i = kp1; i <= n; ++i) {
			prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
			prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
			ar[i + k * ar_dim1] = -prodr;
			ai[i + k * ai_dim1] = -prodi;
		/* L30: */
		}
		for (j = kp1; j <= n; ++j) {
			tr = ar[m + j * ar_dim1];
			ti = ai[m + j * ai_dim1];
			ar[m + j * ar_dim1] = ar[k + j * ar_dim1];
			ai[m + j * ai_dim1] = ai[k + j * ai_dim1];
			ar[k + j * ar_dim1] = tr;
			ai[k + j * ai_dim1] = ti;
			if (abs(tr) + abs(ti) == 0.) {
				goto L48;
			}
			if (ti == 0.) {
				for (i = kp1; i <= n; ++i) {
					prodr = ar[i + k * ar_dim1] * tr;
					prodi = ai[i + k * ai_dim1] * tr;
					ar[i + j * ar_dim1] += prodr;
					ai[i + j * ai_dim1] += prodi;
				/* L40: */
				}
				goto L48;
			}
			if (tr == 0.) {
				for (i = kp1; i <= n; ++i) {
					prodr = -ai[i + k * ai_dim1] * ti;
					prodi = ar[i + k * ar_dim1] * ti;
					ar[i + j * ar_dim1] += prodr;
					ai[i + j * ai_dim1] += prodi;
					/* L45: */
				}
				goto L48;
			}
			for (i = kp1; i <= n; ++i) {
				prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
				prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
				ar[i + j * ar_dim1] += prodr;
				ai[i + j * ai_dim1] += prodi;
				/* L47: */	
			}
L48:
	/* L50: */
			;
		}
	/* L60: */
    }
L70:
    k = n;
    if (abs(ar[n + n * ar_dim1]) + abs(ai[n + n * ai_dim1]) == 0.) {
		goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECC ------------------------ */
} /* decc_ */

/* Subroutine */ int solc_(integer n, integer *ndim, doublereal *ar, 
	doublereal *ai, doublereal *br, doublereal *bi, integer *ip)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset;

    /* Local variables */
    static integer i, k, m, kb;
    static doublereal ti, tr;
    static integer km1, nm1, kp1;
    static doublereal den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX. */
/*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
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
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    if (n == 1) {
		goto L50;
    }
    nm1 = n - 1;
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = ip[k];
		tr = br[m];
		ti = bi[m];
		br[m] = br[k];
		bi[m] = bi[k];
		br[k] = tr;
		bi[k] = ti;
		for (i = kp1; i <= n; ++i) {
			prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
			prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
			br[i] += prodr;
			bi[i] += prodi;
			/* L10: */
		}
		/* L20: */
    }
    for (kb = 1; kb <= nm1; ++kb) {
		km1 = n - kb;
		k = km1 + 1;
		den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] + ai[k + k * ai_dim1] 
			* ai[k + k * ai_dim1];
		prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
		prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
		br[k] = prodr / den;
		bi[k] = prodi / den;
		tr = -br[k];
		ti = -bi[k];
		for (i = 1; i <= km1; ++i) {
			prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
			prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
			br[i] += prodr;
			bi[i] += prodi;
			/* L30: */
		}
		/* L40: */
    }
L50:
    den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 + 1];
    prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
    prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
    br[1] = prodr / den;
    bi[1] = prodi / den;
    return 0;
/* ----------------------- END OF SUBROUTINE SOLC ------------------------ */
} /* solc_ */

/* Subroutine */ int dechc_(integer n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *lb, integer *ip, integer *ier)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset;

    /* Local variables */
    static integer i, j, k, m, na;
    static doublereal ti, tr;
    static integer nm1, kp1;
    static doublereal den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION */
/*  ------ MODIFICATION FOR COMPLEX MATRICES -------- */
/*  INPUT.. */
/*     N = ORDER OF MATRIX. */
/*     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI . */
/*     (AR, AI) = MATRIX TO BE TRIANGULARIZED. */
/*  OUTPUT.. */
/*     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART. */
/*     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART. */
/*     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*                                                    REAL PART. */
/*     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L. */
/*                                                    IMAGINARY PART. */
/*     LB = LOWER BANDWIDTH OF A (DIAGONAL NOT COUNTED), LB.GE.1. */
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
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    *ier = 0;
    ip[n] = 1;
    if (*lb == 0) {
		goto L70;
    }
    if (n == 1) {
		goto L70;
    }
    nm1 = n - 1;
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = k;
		na = min(n, *lb + k);
		for (i = kp1; i <= na; ++i) {
			if (abs(ar[i + k * ar_dim1]) + abs(ai[i + k * ai_dim1]) > abs(ar[m + k * ar_dim1]) + abs(ai[m + k * ai_dim1])) {
				m = i;
			}
			/* L10: */
		}
		ip[k] = m;
		tr = ar[m + k * ar_dim1];
		ti = ai[m + k * ai_dim1];
		if (m == k) {
			goto L20;
		}
		ip[n] = -ip[n];
		ar[m + k * ar_dim1] = ar[k + k * ar_dim1];
		ai[m + k * ai_dim1] = ai[k + k * ai_dim1];
		ar[k + k * ar_dim1] = tr;
		ai[k + k * ai_dim1] = ti;
L20:
		if (abs(tr) + abs(ti) == 0.) {
			goto L80;
		}
		den = tr * tr + ti * ti;
		tr /= den;
		ti = -ti / den;
		for (i = kp1; i <= na; ++i) {
			prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
			prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
			ar[i + k * ar_dim1] = -prodr;
			ai[i + k * ai_dim1] = -prodi;
			/* L30: */
		}
		for (j = kp1; j <= n; ++j) {
			tr = ar[m + j * ar_dim1];
			ti = ai[m + j * ai_dim1];
			ar[m + j * ar_dim1] = ar[k + j * ar_dim1];
			ai[m + j * ai_dim1] = ai[k + j * ai_dim1];
			ar[k + j * ar_dim1] = tr;
			ai[k + j * ai_dim1] = ti;
			if (abs(tr) + abs(ti) == 0.) {
				goto L48;
			}
			if (ti == 0.) {
				for (i = kp1; i <= na; ++i) {
					prodr = ar[i + k * ar_dim1] * tr;
					prodi = ai[i + k * ai_dim1] * tr;
					ar[i + j * ar_dim1] += prodr;
					ai[i + j * ai_dim1] += prodi;
					/* L40: */
				}
				goto L48;
			}
			if (tr == 0.) {
				for (i = kp1; i <= na; ++i) {
					prodr = -ai[i + k * ai_dim1] * ti;
					prodi = ar[i + k * ar_dim1] * ti;
					ar[i + j * ar_dim1] += prodr;
					ai[i + j * ai_dim1] += prodi;
					/* L45: */
				}
				goto L48;
			}
			for (i = kp1; i <= na; ++i) {
				prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
				prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
				ar[i + j * ar_dim1] += prodr;
				ai[i + j * ai_dim1] += prodi;
				/* L47: */
			}
L48:
	/* L50: */
			;
		}
	/* L60: */
    }
L70:
    k = n;
    if (abs(ar[n + n * ar_dim1]) + abs(ai[n + n * ai_dim1]) == 0.) {
		goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECHC ----------------------- */
} /* dechc_ */

/* Subroutine */ int solhc_(integer n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *lb, doublereal *br, doublereal *bi, integer *
	ip)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset;

    /* Local variables */
    static integer i, k, m, kb;
    static doublereal ti, tr;
    static integer km1, nm1, kp1;
    static doublereal den, prodi, prodr;

/* VERSION COMPLEX DOUBLE PRECISION */
/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N = ORDER OF MATRIX. */
/*    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI. */
/*    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC. */
/*    (BR,BI) = RIGHT HAND SIDE VECTOR. */
/*    LB = LOWER BANDWIDTH OF A. */
/*    IP = PIVOT VECTOR OBTAINED FROM DEC. */
/*  DO NOT USE IF DEC HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    (BR,BI) = SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    --bi;
    --br;
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    if (n == 1) {
		goto L50;
    }
    nm1 = n - 1;
    if (*lb == 0) {
		goto L25;
    }
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = ip[k];
		tr = br[m];
		ti = bi[m];
		br[m] = br[k];
		bi[m] = bi[k];
		br[k] = tr;
		bi[k] = ti;
		for (i = kp1; i <= min(n, *lb + k); ++i) {
			prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
			prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
			br[i] += prodr;
			bi[i] += prodi;
			/* L10: */
		}
		/* L20: */
    }
L25:
    for (kb = 1; kb <= nm1; ++kb) {
		km1 = n - kb;
		k = km1 + 1;
		den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] + ai[k + k * ai_dim1] 
			* ai[k + k * ai_dim1];
		prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
		prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
		br[k] = prodr / den;
		bi[k] = prodi / den;
		tr = -br[k];
		ti = -bi[k];
		for (i = 1; i <= km1; ++i) {
			prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
			prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
			br[i] += prodr;
			bi[i] += prodi;
			/* L30: */
		}
		/* L40: */
    }
L50:
    den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 + 1];
    prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
    prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
    br[1] = prodr / den;
    bi[1] = prodi / den;
    return 0;
/* ----------------------- END OF SUBROUTINE SOLHC ----------------------- */
} /* solhc_ */

/* Subroutine */ int decb_(integer n, integer *ndim, doublereal *a, integer *
	ml, integer *mu, integer *ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    static integer i, j, k, m;
    static doublereal t;
    static integer md, jk, mm, ju, md1, nm1, kp1, mdl, ijk;

/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED */
/*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
/*  INPUT.. */
/*     N       ORDER OF THE ORIGINAL MATRIX A. */
/*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
/*     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND */
/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
/*                ML+1 THROUGH 2*ML+MU+1 OF  A. */
/*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*  OUTPUT.. */
/*     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
/*     IP      INDEX VECTOR OF PIVOT INDICES. */
/*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
/*                SINGULAR AT STAGE K. */
/*  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1. */
/*  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     THIS IS A MODIFICATION OF */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
    ip[n] = 1;
    md = *ml + *mu + 1;
    md1 = md + 1;
    ju = 0;
    if (*ml == 0) {
		goto L70;
    }
    if (n == 1) {
		goto L70;
    }
    if (n < *mu + 2) {
		goto L7;
    }
    for (j = *mu + 2; j <= n; ++j) {
		for (i = 1; i <= *ml; ++i) {
			/* L5: */
			a[i + j * a_dim1] = 0.;
		}
    }
L7:
    nm1 = n - 1;
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = md;
		mdl = min(*ml, n - k) + md;
		for (i = md1; i <= mdl; ++i) {
			if (abs(a[i + k * a_dim1]) > abs(a[m + k * a_dim1])) {
				m = i;
			}
		/* L10: */
		}
		ip[k] = m + k - md;
		t = a[m + k * a_dim1];
		if (m == md) {
			goto L20;
		}
		ip[n] = -ip[n];
		a[m + k * a_dim1] = a[md + k * a_dim1];
		a[md + k * a_dim1] = t;
L20:
		if (t == 0.) {
			goto L80;
		}
		t = 1. / t;
		for (i = md1; i <= mdl; ++i) {
			/* L30: */
			a[i + k * a_dim1] = -a[i + k * a_dim1] * t;
		}
		ju = min(max(ju, *mu + ip[k]), n);
		mm = md;
		if (ju < kp1) {
			goto L55;
		}
		for (j = kp1; j <= ju; ++j) {
			--m;
			--mm;
			t = a[m + j * a_dim1];
			if (m == mm) {
				goto L35;
			}
			a[m + j * a_dim1] = a[mm + j * a_dim1];
			a[mm + j * a_dim1] = t;
L35:
			if (t == 0.) {
				goto L45;
			}
			jk = j - k;
			for (i = md1; i <= mdl; ++i) {
				ijk = i - jk;
				/* L40: */
				a[ijk + j * a_dim1] += a[i + k * a_dim1] * t;
			}
L45:
	/* L50: */
			;
		}
L55:
	/* L60: */
		;
    }
L70:
    k = n;
    if (a[md + n * a_dim1] == 0.) {
		goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECB ------------------------ */
} /* decb_ */

/* Subroutine */ int solb_(integer n, integer *ndim, doublereal *a, integer *
	ml, integer *mu, doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    static integer i, k, m;
    static doublereal t;
    static integer kb, md, lm, md1, nm1, imd, kmd, mdl, mdm;

/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B . */
/*  INPUT.. */
/*    N      ORDER OF MATRIX A. */
/*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
/*    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB. */
/*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*    B      RIGHT HAND SIDE VECTOR. */
/*    IP     PIVOT VECTOR OBTAINED FROM DECB. */
/*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    B      SOLUTION VECTOR, X . */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    --b;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    md = *ml + *mu + 1;
    md1 = md + 1;
    mdm = md - 1;
    nm1 = n - 1;
    if (*ml == 0) {
		goto L25;
    }
    if (n == 1) {
		goto L50;
    }
    for (k = 1; k <= nm1; ++k) {
		m = ip[k];
		t = b[m];
		b[m] = b[k];
		b[k] = t;
		mdl = min(*ml, n - k) + md;
		for (i = md1; i <= mdl; ++i) {
			imd = i + k - md;
			/* L10: */
			b[imd] += a[i + k * a_dim1] * t;
		}
		/* L20: */
    }
L25:
    for (kb = 1; kb <= nm1; ++kb) {
		k = n + 1 - kb;
		b[k] /= a[md + k * a_dim1];
		t = -b[k];
		kmd = md - k;
		lm = max(1, kmd + 1);
		for (i = lm; i <= mdm; ++i) {
			imd = i - kmd;
			/* L30: */
			b[imd] += a[i + k * a_dim1] * t;
		}
		/* L40: */
    }
L50:
    b[1] /= a[md + a_dim1];
    return 0;
/* ----------------------- END OF SUBROUTINE SOLB ------------------------ */
} /* solb_ */

/* Subroutine */ int decbc_(integer n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *ml, integer *mu, integer *ip, integer *ier)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset;

    /* Local variables */
    static integer i, j, k, m, md, jk, mm;
    static doublereal ti;
    static integer ju;
    static doublereal tr;
    static integer md1, nm1, kp1;
    static doublereal den;
    static integer mdl, ijk;
    static doublereal prodi, prodr;

/* ----------------------------------------------------------------------- */
/*  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED COMPLEX */
/*  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU */
/*  INPUT.. */
/*     N       ORDER OF THE ORIGINAL MATRIX A. */
/*     NDIM    DECLARED DIMENSION OF ARRAY  A. */
/*     AR, AI     CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS */
/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  AR (REAL */
/*                PART) AND AI (IMAGINARY PART)  AND */
/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
/*                ML+1 THROUGH 2*ML+MU+1 OF  AR AND AI. */
/*     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*  OUTPUT.. */
/*     AR, AI  AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
/*     IP      INDEX VECTOR OF PIVOT INDICES. */
/*     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O . */
/*     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE */
/*                SINGULAR AT STAGE K. */
/*  USE  SOLBC  TO OBTAIN SOLUTION OF LINEAR SYSTEM. */
/*  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1. */
/*  IF IP(N)=O, A IS SINGULAR, SOLBC WILL DIVIDE BY ZERO. */

/*  REFERENCE.. */
/*     THIS IS A MODIFICATION OF */
/*     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER, */
/*     C.A.C.M. 15 (1972), P. 274. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    *ier = 0;
    ip[n] = 1;
    md = *ml + *mu + 1;
    md1 = md + 1;
    ju = 0;
    if (*ml == 0) {
		goto L70;
    }
    if (n == 1) {
		goto L70;
    }
    if (n < *mu + 2) {
		goto L7;
    }
    for (j = *mu + 2; j <= n; ++j) {
		for (i = 1; i <= *ml; ++i) {
			ar[i + j * ar_dim1] = 0.;
			ai[i + j * ai_dim1] = 0.;
			/* L5: */
		}
    }
L7:
    nm1 = n - 1;
    for (k = 1; k <= nm1; ++k) {
		kp1 = k + 1;
		m = md;
		mdl = min(*ml, n - k) + md;
		for (i = md1; i <= mdl; ++i) {
			if (abs(ar[i + k * ar_dim1]) + abs(ai[i + k * ai_dim1]) > abs(ar[m + k * ar_dim1]) + abs(ai[m + k * ai_dim1])) {
				m = i;
			}
			/* L10: */
		}
		ip[k] = m + k - md;
		tr = ar[m + k * ar_dim1];
		ti = ai[m + k * ai_dim1];
		if (m == md) {
			goto L20;
		}
		ip[n] = -ip[n];
		ar[m + k * ar_dim1] = ar[md + k * ar_dim1];
		ai[m + k * ai_dim1] = ai[md + k * ai_dim1];
		ar[md + k * ar_dim1] = tr;
		ai[md + k * ai_dim1] = ti;
L20:
		if (abs(tr) + abs(ti) == 0.) {
			goto L80;
		}
		den = tr * tr + ti * ti;
		tr /= den;
		ti = -ti / den;
		for (i = md1; i <= mdl; ++i) {
			prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
			prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
			ar[i + k * ar_dim1] = -prodr;
			ai[i + k * ai_dim1] = -prodi;
		/* L30: */
		}
		ju = min(max(ju, *mu + ip[k]), n);
		mm = md;
		if (ju < kp1) {
			goto L55;
		}
		for (j = kp1; j <= ju; ++j) {
			--m;
			--mm;
			tr = ar[m + j * ar_dim1];
			ti = ai[m + j * ai_dim1];
			if (m == mm) {
				goto L35;
			}
			ar[m + j * ar_dim1] = ar[mm + j * ar_dim1];
			ai[m + j * ai_dim1] = ai[mm + j * ai_dim1];
			ar[mm + j * ar_dim1] = tr;
			ai[mm + j * ai_dim1] = ti;
L35:
			if (abs(tr) + abs(ti) == 0.) {
				goto L48;
			}
			jk = j - k;
			if (ti == 0.) {
				for (i = md1; i <= mdl; ++i) {
					ijk = i - jk;
					prodr = ar[i + k * ar_dim1] * tr;
					prodi = ai[i + k * ai_dim1] * tr;
					ar[ijk + j * ar_dim1] += prodr;
					ai[ijk + j * ai_dim1] += prodi;
				/* L40: */
				}
				goto L48;
			}
			if (tr == 0.) {
				for (i = md1; i <= mdl; ++i) {
					ijk = i - jk;
					prodr = -ai[i + k * ai_dim1] * ti;
					prodi = ar[i + k * ar_dim1] * ti;
					ar[ijk + j * ar_dim1] += prodr;
					ai[ijk + j * ai_dim1] += prodi;
				/* L45: */
				}
				goto L48;
			}
			for (i = md1; i <= mdl; ++i) {
				ijk = i - jk;
				prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
				prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
				ar[ijk + j * ar_dim1] += prodr;
				ai[ijk + j * ai_dim1] += prodi;
				/* L47: */
			}
L48:
	/* L50: */
			;
		}
L55:
	/* L60: */
		;
    }
L70:
    k = n;
    if (abs(ar[md + n * ar_dim1]) + abs(ai[md + n * ai_dim1]) == 0.) {
		goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECBC ------------------------ */
} /* decbc_ */

/* Subroutine */ int solbc_(integer n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *ml, integer *mu, doublereal *br, doublereal *
	bi, integer *ip)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset;

    /* Local variables */
    static integer i, k, m, kb, md, lm;
    static doublereal ti, tr;
    static integer md1, nm1;
    static doublereal den;
    static integer imd, kmd, mdl, mdm;
    static doublereal prodi, prodr;

/* ----------------------------------------------------------------------- */
/*  SOLUTION OF LINEAR SYSTEM, A*X = B , */
/*                  VERSION BANDED AND COMPLEX-DOUBLE PRECISION. */
/*  INPUT.. */
/*    N      ORDER OF MATRIX A. */
/*    NDIM   DECLARED DIMENSION OF ARRAY  A . */
/*    AR, AI TRIANGULARIZED MATRIX OBTAINED FROM DECB (REAL AND IMAG. PART). */
/*    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED). */
/*    BR, BI RIGHT HAND SIDE VECTOR (REAL AND IMAG. PART). */
/*    IP     PIVOT VECTOR OBTAINED FROM DECBC. */
/*  DO NOT USE IF DECB HAS SET IER .NE. 0. */
/*  OUTPUT.. */
/*    BR, BI SOLUTION VECTOR, X (REAL AND IMAG. PART). */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    --bi;
    --br;
    ai_dim1 = *ndim;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *ndim;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;

    /* Function Body */
    md = *ml + *mu + 1;
    md1 = md + 1;
    mdm = md - 1;
    nm1 = n - 1;
    if (*ml == 0) {
		goto L25;
    }
    if (n == 1) {
		goto L50;
    }
    for (k = 1; k <= nm1; ++k) {
		m = ip[k];
		tr = br[m];
		ti = bi[m];
		br[m] = br[k];
		bi[m] = bi[k];
		br[k] = tr;
		bi[k] = ti;
		mdl = min(*ml, n - k) + md;
		for (i = md1; i <= mdl; ++i) {
			imd = i + k - md;
			prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
			prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
			br[imd] += prodr;
			bi[imd] += prodi;
			/* L10: */
		}
		/* L20: */
    }
L25:
    for (kb = 1; kb <= nm1; ++kb) {
		k = n + 1 - kb;
		den = ar[md + k * ar_dim1] * ar[md + k * ar_dim1] + ai[md + k * 
			ai_dim1] * ai[md + k * ai_dim1];
		prodr = br[k] * ar[md + k * ar_dim1] + bi[k] * ai[md + k * ai_dim1];
		prodi = bi[k] * ar[md + k * ar_dim1] - br[k] * ai[md + k * ai_dim1];
		br[k] = prodr / den;
		bi[k] = prodi / den;
		tr = -br[k];
		ti = -bi[k];
		kmd = md - k;
		lm = max(1, kmd + 1);
		for (i = lm; i <= mdm; ++i) {
			imd = i - kmd;
			prodr = ar[i + k * ar_dim1] * tr - ai[i + k * ai_dim1] * ti;
			prodi = ai[i + k * ai_dim1] * tr + ar[i + k * ar_dim1] * ti;
			br[imd] += prodr;
			bi[imd] += prodi;
			/* L30: */
		}
		/* L40: */
    }
    den = ar[md + ar_dim1] * ar[md + ar_dim1] + ai[md + ai_dim1] * ai[md + ai_dim1];
    prodr = br[1] * ar[md + ar_dim1] + bi[1] * ai[md + ai_dim1];
    prodi = bi[1] * ar[md + ar_dim1] - br[1] * ai[md + ai_dim1];
    br[1] = prodr / den;
    bi[1] = prodi / den;
L50:
    return 0;
/* ----------------------- END OF SUBROUTINE SOLBC ------------------------ */
} /* solbc_ */

/* Subroutine */ int elmhes_(integer *nm, integer n, integer *low, integer *
	igh, doublereal *a, integer *int__)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    static integer i, j, m;
    static doublereal x, y;
    static integer la, mm1, kp1, mp1;

/*     this subroutine is a translation of the algol procedure elmhes, */
/*     num. math. 12, 349-368(1968) by martin and wilkinson. */
/*     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971). */

/*     given a real general matrix, this subroutine */
/*     reduces a submatrix situated in rows and columns */
/*     low through igh to upper hessenberg form by */
/*     stabilized elementary similarity transformations. */

/*     on input: */

/*      nm must be set to the row dimension of two-dimensional */
/*        array parameters as declared in the calling program */
/*        dimension statement; */

/*      n is the order of the matrix; */

/*      low and igh are integers determined by the balancing */
/*        subroutine  balanc.      if  balanc  has not been used, */
/*        set low=1, igh=n; */

/*      a contains the input matrix. */

/*     on output: */

/*      a contains the hessenberg matrix.  the multipliers */
/*        which were used in the reduction are stored in the */
/*        remaining triangle under the hessenberg matrix; */

/*      int contains information on the rows and columns */
/*        interchanged in the reduction. */
/*        only elements low through igh are used. */

/*     questions and comments should be directed to b. s. garbow, */
/*     applied mathematics division, argonne national laboratory */

/*     ------------------------------------------------------------------ */

    /* Parameter adjustments */
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --int__;

    /* Function Body */
    la = *igh - 1;
    kp1 = *low + 1;
    if (la < kp1) {
		goto L200;
    }

    for (m = kp1; m <= la; ++m) {
		mm1 = m - 1;
		x = 0.;
		i = m;

		for (j = m; j <= *igh; ++j) {
			if (abs(a[j + mm1 * a_dim1]) <= abs(x)) {
				goto L100;
			}
			x = a[j + mm1 * a_dim1];
			i = j;
L100:
			;
		}

		int__[m] = i;
		if (i == m) {
			goto L130;
		}
		/*    :::::::::: interchange rows and columns of a :::::::::: */
		for (j = mm1; j <= n; ++j) {
			y = a[i + j * a_dim1];
			a[i + j * a_dim1] = a[m + j * a_dim1];
			a[m + j * a_dim1] = y;
			/* L110: */
		}

		for (j = 1; j <= *igh; ++j) {
			y = a[j + i * a_dim1];
			a[j + i * a_dim1] = a[j + m * a_dim1];
			a[j + m * a_dim1] = y;
			/* L120: */
		}
		/*    :::::::::: end interchange :::::::::: */
L130:
		if (x == 0.) {
			goto L180;
		}
		mp1 = m + 1;

		for (i = mp1; i <= *igh; ++i) {
			y = a[i + mm1 * a_dim1];
			if (y == 0.) {
				goto L160;
			}
			y /= x;
			a[i + mm1 * a_dim1] = y;

			for (j = m; j <= n; ++j) {
				/* L140: */
				a[i + j * a_dim1] -= y * a[m + j * a_dim1];
			}

			for (j = 1; j <= *igh; ++j) {
				/* L150: */
				a[j + m * a_dim1] += y * a[j + i * a_dim1];
			}

L160:
			;
		}

L180:
		;
    }

L200:
    return 0;
/*    :::::::::: last card of elmhes :::::::::: */
} /* elmhes_ */

/* ****************************************** */
/*     VERSION OF SEPTEMBER 18, 1995 */
/* ****************************************** */

/* Subroutine */ int decomr_(integer n, doublereal *fjac, integer *ldjac, 
	doublereal *fmas, integer *ldmas, integer *mlmas, integer *mumas, 
	integer *m1, integer *m2, integer *nm1, doublereal *fac1, doublereal *e1,
	integer *lde1, integer *ip1, integer *ier, integer *ijob, logical *calhes,
	integer *iphes, Radau_SuperLU_aux* radau_slu_aux)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset;

    /* Local variables */
    static integer i, j, k, j1, ib, mm, jm1;
	static doublereal sum;
    extern /* Subroutine */ int dec_(integer, integer *, doublereal *, 
	    integer *, integer *);
    extern /* Subroutine */ int decb_(integer, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *);
    extern /* Subroutine */ int dech_(integer, integer *, doublereal *,
		integer *, integer *, integer *);
    extern /* Subroutine */ int elmhes_(integer *, integer, integer *,
	    integer *, doublereal *, integer *);

    /* Parameter adjustments */
    --iphes;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip1;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;

    /* Function Body */
    switch (*ijob) {
		case 1:  goto L1;
		case 2:  goto L2;
		case 3:  goto L3;
		case 4:  goto L4;
		case 5:  goto L5;
		case 6:  goto L6;
		case 7:  goto L7;
		case 8:  goto L8;
		case 10:  goto L55;
		case 11:  goto L11;
		case 12:  goto L12;
		case 13:  goto L13;
		case 14:  goto L14;
		case 15:  goto L15;
    }

/* ----------------------------------------------------------- */

L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= n; ++i) {
			e1[i + j * e1_dim1] = -fjac[i + j * fjac_dim1];
		}
		e1[j + j * e1_dim1] += *fac1;
    }
    dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= *nm1; ++i) {
			e1[i + j * e1_dim1] = -fjac[i + jm1 * fjac_dim1];
		}
		e1[j + j * e1_dim1] += *fac1;
    }
L45:
    mm = *m1 / *m2;
    for (j = 1; j <= *m2; ++j) {
		for (i = 1; i <= *nm1; ++i) {
			sum = 0.;
			for (k = 0; k <= mm - 1; ++k) {
				sum = (sum + fjac[i + (j + k * *m2) * fjac_dim1]) / *fac1;
			}
			e1[i + j * e1_dim1] -= sum;
		}
    }
    dec_(*nm1, lde1, &e1[e1_offset], &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= linal_1.mbjac; ++i) {
			e1[i + linal_1.mle + j * e1_dim1] = -fjac[i + j * fjac_dim1];
		}
		e1[linal_1.mdiag + j * e1_dim1] += *fac1;
    }
    decb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= linal_1.mbjac; ++i) {
			e1[i + linal_1.mle + j * e1_dim1] = -fjac[i + jm1 * fjac_dim1];
		}
		e1[linal_1.mdiag + j * e1_dim1] += *fac1;
    }
L46:
    mm = *m1 / *m2;
    for (j = 1; j <= *m2; ++j) {
		for (i = 1; i <= linal_1.mbjac; ++i) {
			sum = 0.;
			for (k = 0; k <= mm - 1; ++k) {
				sum = (sum + fjac[i + (j + k * *m2) * fjac_dim1]) / *fac1;
			}
			e1[i + linal_1.mle + j * e1_dim1] -= sum;
		}
    }
    decb_(*nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= n; ++i) {
			e1[i + j * e1_dim1] = -fjac[i + j * fjac_dim1];
		}
		for (i = max(1, j - *mumas); i <= min(n, j + *mlmas); ++i) {
			e1[i + j * e1_dim1] += *fac1 * fmas[i - j + linal_1.mbdiag + j * fmas_dim1];
		}
    }
    dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= *nm1; ++i) {
			e1[i + j * e1_dim1] = -fjac[i + jm1 * fjac_dim1];
		}
		for (i = max(1, j - *mumas); i <= min(*nm1, j + *mlmas); ++i) {
			e1[i + j * e1_dim1] += *fac1 * fmas[i - j + linal_1.mbdiag + j * fmas_dim1];
		}
    }
    goto L45;
/* ----------------------------------------------------------- */

L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= linal_1.mbjac; ++i) {
			e1[i + linal_1.mle + j * e1_dim1] = -fjac[i + j * fjac_dim1];
		}
		for (i = 1; i <= linal_1.mbb; ++i) {
			ib = i + linal_1.mdiff;
			e1[ib + j * e1_dim1] += *fac1 * fmas[i + j * fmas_dim1];
		}
    }
    decb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L14:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= linal_1.mbjac; ++i) {
			e1[i + linal_1.mle + j * e1_dim1] = -fjac[i + jm1 * fjac_dim1];
		}
		for (i = 1; i <= linal_1.mbb; ++i) {
			ib = i + linal_1.mdiff;
			e1[ib + j * e1_dim1] += *fac1 * fmas[i + j * fmas_dim1];
		}
    }
    goto L46;

/* ----------------------------------------------------------- */

L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= n; ++i) {
			e1[i + j * e1_dim1] = fmas[i + j * fmas_dim1] * *fac1 - fjac[i + j * fjac_dim1];
		}
    }
    dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= *nm1; ++i) {
			e1[i + j * e1_dim1] = fmas[i + j * fmas_dim1] * *fac1 - fjac[i + jm1 * fjac_dim1];
		}
    }
    goto L45;

/* ----------------------------------------------------------- */

L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
    return 0;

/* ----------------------------------------------------------- */

L7:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
    if (*calhes) {
		elmhes_(ldjac, n, &c__1, &n, &fjac[fjac_offset], &iphes[1]);
    }
    *calhes = FALSE_;
    for (j = 1; j <= n - 1; ++j) {
		j1 = j + 1;
		e1[j1 + j * e1_dim1] = -fjac[j1 + j * fjac_dim1];
    }
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= j; ++i) {
			e1[i + j * e1_dim1] = -fjac[i + j * fjac_dim1];
		}
		e1[j + j * e1_dim1] += *fac1;
    }
    dech_(n, lde1, &e1[e1_offset], &c__1, &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L8:
/* ---  B=IDENTITY, SPARSE LU */
	superlu_setup_d(radau_slu_aux->slu_aux_d, *fac1, radau_slu_aux->jac_data, radau_slu_aux->jac_indices, radau_slu_aux->jac_indptr, radau_slu_aux->fresh_jacobian, radau_slu_aux->nnz_actual);
	*ier = superlu_factorize_d(radau_slu_aux->slu_aux_d);
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* decomr_ */

/*     END OF SUBROUTINE DECOMR */

/* *********************************************************** */

/* Subroutine */ int decomc_(integer n, doublereal *fjac, integer *ldjac, 
	doublereal *fmas, integer *ldmas, integer *mlmas, integer *mumas, 
	integer *m1, integer *m2, integer *nm1, doublereal *alphn, doublereal *betan,
	doublereal *e2r, doublereal *e2i, integer *lde1, integer *ip2,
	integer *ier, integer *ijob, Radau_SuperLU_aux* radau_slu_aux)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e2r_dim1, 
	    e2r_offset, e2i_dim1, e2i_offset;

    /* Local variables */
    static integer i, j, k, j1;
    static doublereal bb;
    static integer ib, mm, jm1;
    static doublereal bet, alp;
    extern /* Subroutine */ int decc_(integer, integer *, doublereal *, 
	    doublereal *, integer *, integer *);
    static doublereal ffma, abno;
    static integer imle;
    static doublereal sumi, sumr, sums;
    extern /* Subroutine */ int decbc_(integer, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *), dechc_(
	    integer, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *);


    /* Parameter adjustments */
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip2;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e2i_dim1 = *lde1;
    e2i_offset = 1 + e2i_dim1;
    e2i -= e2i_offset;
    e2r_dim1 = *lde1;
    e2r_offset = 1 + e2r_dim1;
    e2r -= e2r_offset;

    /* Function Body */
    switch (*ijob) {
		case 1:  goto L1;
		case 2:  goto L2;
		case 3:  goto L3;
		case 4:  goto L4;
		case 5:  goto L5;
		case 6:  goto L6;
		case 7:  goto L7;
		case 8:  goto L8;
		case 10:  goto L55;
		case 11:  goto L11;
		case 12:  goto L12;
		case 13:  goto L13;
		case 14:  goto L14;
		case 15:  goto L15;
    }

/* ----------------------------------------------------------- */

L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= n; ++i) {
			e2r[i + j * e2r_dim1] = -fjac[i + j * fjac_dim1];
			e2i[i + j * e2i_dim1] = 0.;
		}
		e2r[j + j * e2r_dim1] += *alphn;
		e2i[j + j * e2i_dim1] = *betan;
    }
    decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= *nm1; ++i) {
			e2r[i + j * e2r_dim1] = -fjac[i + jm1 * fjac_dim1];
			e2i[i + j * e2i_dim1] = 0.;
		}
		e2r[j + j * e2r_dim1] += *alphn;
		e2i[j + j * e2i_dim1] = *betan;
    }
L45:
    mm = *m1 / *m2;
    abno = *alphn * *alphn + *betan * *betan;
    alp = *alphn / abno;
    bet = *betan / abno;
    for (j = 1; j <= *m2; ++j) {
		for (i = 1; i <= *nm1; ++i) {
			sumr = 0.;
			sumi = 0.;
			for (k = 0; k <= mm - 1; ++k) {
				sums = sumr + fjac[i + (j + k * *m2) * fjac_dim1];
				sumr = sums * alp + sumi * bet;
				sumi = sumi * alp - sums * bet;
			}
			e2r[i + j * e2r_dim1] -= sumr;
			e2i[i + j * e2i_dim1] -= sumi;
		}
    }
    decc_(*nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= linal_1.mbjac; ++i) {
			imle = i + linal_1.mle;
			e2r[imle + j * e2r_dim1] = -fjac[i + j * fjac_dim1];
			e2i[imle + j * e2i_dim1] = 0.;
		}
		e2r[linal_1.mdiag + j * e2r_dim1] += *alphn;
		e2i[linal_1.mdiag + j * e2i_dim1] = *betan;
    }
    decbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= linal_1.mbjac; ++i) {
			e2r[i + linal_1.mle + j * e2r_dim1] = -fjac[i + jm1 * fjac_dim1];
			e2i[i + linal_1.mle + j * e2i_dim1] = 0.;
		}
		e2r[linal_1.mdiag + j * e2r_dim1] += *alphn;
		e2i[linal_1.mdiag + j * e2i_dim1] += *betan;
    }
L46:
    mm = *m1 / *m2;
    abno = *alphn * *alphn + *betan * *betan;
    alp = *alphn / abno;
    bet = *betan / abno;
    for (j = 1; j <= *m2; ++j) {
		for (i = 1; i <= linal_1.mbjac; ++i) {
			sumr = 0.;
			sumi = 0.;
			for (k = 0; k <= mm - 1; ++k) {
				sums = sumr + fjac[i + (j + k * *m2) * fjac_dim1];
				sumr = sums * alp + sumi * bet;
				sumi = sumi * alp - sums * bet;
			}
			imle = i + linal_1.mle;
			e2r[imle + j * e2r_dim1] -= sumr;
			e2i[imle + j * e2i_dim1] -= sumi;
		}
    }
    decbc_(*nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &linal_1.mue, &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= n; ++i) {
			e2r[i + j * e2r_dim1] = -fjac[i + j * fjac_dim1];
			e2i[i + j * e2i_dim1] = 0.;
		}
    }
    for (j = 1; j <= n; ++j) {
		for (i = max(1, j - *mumas); i <= min(n, j + *mlmas); ++i) {
			bb = fmas[i - j + linal_1.mbdiag + j * fmas_dim1];
			e2r[i + j * e2r_dim1] += *alphn * bb;
			e2i[i + j * e2i_dim1] = *betan * bb;
		}
    }
    decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= *nm1; ++i) {
			e2r[i + j * e2r_dim1] = -fjac[i + jm1 * fjac_dim1];
			e2i[i + j * e2i_dim1] = 0.;
		}
		for (i = max(1, j - *mumas); i <= min(*nm1, j + *mlmas); ++i) {
			ffma = fmas[i - j + linal_1.mbdiag + j * fmas_dim1];
			e2r[i + j * e2r_dim1] += *alphn * ffma;
			e2i[i + j * e2i_dim1] += *betan * ffma;
		}
    }
    goto L45;

/* ----------------------------------------------------------- */

L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= linal_1.mbjac; ++i) {
			imle = i + linal_1.mle;
			e2r[imle + j * e2r_dim1] = -fjac[i + j * fjac_dim1];
			e2i[imle + j * e2i_dim1] = 0.;
		}
		for (i = max(1, *mumas + 2 - j); i <= min(linal_1.mbb, *mumas + 1 - j + n); ++i) {
			ib = i + linal_1.mdiff;
			bb = fmas[i + j * fmas_dim1];
			e2r[ib + j * e2r_dim1] += *alphn * bb;
			e2i[ib + j * e2i_dim1] = *betan * bb;
		}
    }
    decbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L14:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= linal_1.mbjac; ++i) {
			e2r[i + linal_1.mle + j * e2r_dim1] = -fjac[i + jm1 * fjac_dim1];
			e2i[i + linal_1.mle + j * e2i_dim1] = 0.;
		}
		for (i = 1; i <= linal_1.mbb; ++i) {
			ib = i + linal_1.mdiff;
			ffma = fmas[i + j * fmas_dim1];
			e2r[ib + j * e2r_dim1] += *alphn * ffma;
			e2i[ib + j * e2i_dim1] += *betan * ffma;
		}
    }
    goto L46;

/* ----------------------------------------------------------- */

L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= n; ++i) {
			bb = fmas[i + j * fmas_dim1];
			e2r[i + j * e2r_dim1] = bb * *alphn - fjac[i + j * fjac_dim1];
			e2i[i + j * e2i_dim1] = bb * *betan;
		}
    }
    decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (j = 1; j <= *nm1; ++j) {
		jm1 = j + *m1;
		for (i = 1; i <= *nm1; ++i) {
			e2r[i + j * e2r_dim1] = *alphn * fmas[i + j * fmas_dim1] - fjac[i + jm1 * fjac_dim1];
			e2i[i + j * e2i_dim1] = *betan * fmas[i + j * fmas_dim1];
		}
    }
    goto L45;

/* ----------------------------------------------------------- */

L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
    return 0;

/* ----------------------------------------------------------- */

L7:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
    for (j = 1; j <= n - 1; ++j) {
		j1 = j + 1;
		e2r[j1 + j * e2r_dim1] = -fjac[j1 + j * fjac_dim1];
		e2i[j1 + j * e2i_dim1] = 0.;
    }
    for (j = 1; j <= n; ++j) {
		for (i = 1; i <= j; ++i) {
			e2i[i + j * e2i_dim1] = 0.;
			e2r[i + j * e2r_dim1] = -fjac[i + j * fjac_dim1];
		}
		e2r[j + j * e2r_dim1] += *alphn;
		e2i[j + j * e2i_dim1] = *betan;
    }
    dechc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L8:
/* ---  B=IDENTITY, SPARSE LU */
	superlu_setup_z(radau_slu_aux->slu_aux_z, *alphn, *betan, radau_slu_aux->jac_data, radau_slu_aux->jac_indices, radau_slu_aux->jac_indptr, radau_slu_aux->fresh_jacobian, radau_slu_aux->nnz_actual);
	*ier = superlu_factorize_z(radau_slu_aux->slu_aux_z);
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* decomc_ */

/*     END OF SUBROUTINE DECOMC */

/* *********************************************************** */

/* Subroutine */ int slvrad_(integer n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *nm1,
	doublereal *fac1, doublereal *alphn, doublereal *betan, 
	doublereal *e1, doublereal *e2r, doublereal *e2i, integer *lde1, 
	doublereal *z1, doublereal *z2, doublereal *z3, doublereal *f1, 
	doublereal *f2, doublereal *f3, doublereal *cont, integer *ip1, 
	integer *ip2, integer *iphes, integer *ier, integer *ijob,
	Radau_SuperLU_aux* radau_slu_aux)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, e2r_dim1, e2r_offset, e2i_dim1, e2i_offset;

    /* Local variables */
    static integer i, j, k;
    static doublereal s1, s2, s3, bb;
    static integer mm, mp, im1, jm1, mp1;
    static doublereal z2i, z3i;
    static integer jkm, mpi;
    extern /* Subroutine */ int sol_(integer, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum1, sum2, sum3, ffja, abno;
    extern /* Subroutine */ int solb_(integer, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solc_(integer, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *), solh_(integer, integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal sumh, e1imp;
    extern /* Subroutine */ int solbc_(integer, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal zsafe;
    extern /* Subroutine */ int solhc_(integer, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);

    /* Parameter adjustments */
    --iphes;
    --f3;
    --f2;
    --f1;
    --z3;
    --z2;
    --z1;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip2;
    --ip1;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e2i_dim1 = *lde1;
    e2i_offset = 1 + e2i_dim1;
    e2i -= e2i_offset;
    e2r_dim1 = *lde1;
    e2r_offset = 1 + e2r_dim1;
    e2r -= e2r_offset;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;

    /* Function Body */
    switch (*ijob) {
		case 1:  goto L1;
		case 2:  goto L2;
		case 3:  goto L3;
		case 4:  goto L4;
		case 5:  goto L5;
		case 6:  goto L6;
		case 7:  goto L7;
		case 8:  goto L8;
		case 10:  goto L55;
		case 11:  goto L11;
		case 12:  goto L12;
		case 13:  goto L13;
		case 14:  goto L13;
		case 15:  goto L15;
    }

/* ----------------------------------------------------------- */

L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
    for (i = 1; i <= n; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]);
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (i = 1; i <= n; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
L48:
    abno = *alphn * *alphn + *betan * *betan;
    mm = *m1 / *m2;
    for (j = 1; j <= *m2; ++j) {
		sum1 = 0.;
		sum2 = 0.;
		sum3 = 0.;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum1 = (z1[jkm] + sum1) / *fac1;
			sumh = (z2[jkm] + sum2) / abno;
			sum3 = (z3[jkm] + sum3) / abno;
			sum2 = sumh * *alphn + sum3 * *betan;
			sum3 = sum3 * *alphn - sumh * *betan;
			for (i = 1; i <= *nm1; ++i) {
				im1 = i + *m1;
				z1[im1] += fjac[i + jkm * fjac_dim1] * sum1;
				z2[im1] += fjac[i + jkm * fjac_dim1] * sum2;
				z3[im1] += fjac[i + jkm * fjac_dim1] * sum3;
			}
		}
    }
    sol_(*nm1, lde1, &e1[e1_offset], &z1[*m1 + 1], &ip1[1]);
    solc_(*nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
L49:
    for (i = *m1; i >= 1; --i) {
		mpi = *m2 + i;
		z1[i] = (z1[i] + z1[mpi]) / *fac1;
		z2i = z2[i] + z2[mpi];
		z3i = z3[i] + z3[mpi];
		z3[i] = (z3i * *alphn - z2i * *betan) / abno;
		z2[i] = (z2i * *alphn + z3i * *betan) / abno;
    }
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    for (i = 1; i <= n; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]);
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &linal_1.mue, &z2[1], &z3[1], &ip2[1]);
    return 0;

/* ----------------------------------------------------------- */

L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    for (i = 1; i <= n; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
L45:
    abno = *alphn * *alphn + *betan * *betan;
    mm = *m1 / *m2;
    for (j = 1; j <= *m2; ++j) {
		sum1 = 0.;
		sum2 = 0.;
		sum3 = 0.;
		for (k = mm - 1; k >= 0; --k) {
			jkm = j + k * *m2;
			sum1 = (z1[jkm] + sum1) / *fac1;
			sumh = (z2[jkm] + sum2) / abno;
			sum3 = (z3[jkm] + sum3) / abno;
			sum2 = sumh * *alphn + sum3 * *betan;
			sum3 = sum3 * *alphn - sumh * *betan;
			for (i = max(1, j - *mujac); i <= min(*nm1, j + *mljac); ++i) {
				im1 = i + *m1;
				ffja = fjac[i + *mujac + 1 - j + jkm * fjac_dim1];
				z1[im1] += ffja * sum1;
				z2[im1] += ffja * sum2;
				z3[im1] += ffja * sum3;
			}
		}
    }
    solb_(*nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[*m1 + 1],&ip1[1]);
    solbc_(*nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	       linal_1.mue, &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
    goto L49;

/* ----------------------------------------------------------- */

L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    for (i = 1; i <= n; ++i) {
		s1 = 0.;
		s2 = 0.;
		s3 = 0.;
		for (j = max(1, i - *mlmas); j <= min(n, i + *mumas); ++j) {
			bb = fmas[i - j + linal_1.mbdiag + j * fmas_dim1];
			s1 -= bb * f1[j];
			s2 -= bb * f2[j];
			s3 -= bb * f3[j];
		}
		z1[i] += s1 * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]);
    return 0;

/* ----------------------------------------------------------- */

L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (i = 1; i <= *m1; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
    for (i = 1; i <= *nm1; ++i) {
		im1 = i + *m1;
		s1 = 0.;
		s2 = 0.;
		s3 = 0.;
		for (j = max(1, i - *mlmas); j <= min(*nm1, i + *mumas); ++j) {
			jm1 = j + *m1;
			bb = fmas[i - j + linal_1.mbdiag + j * fmas_dim1];
			s1 -= bb * f1[jm1];
			s2 -= bb * f2[jm1];
			s3 -= bb * f3[jm1];
		}
		z1[im1] += s1 * *fac1;
		z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
		z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
    }
    if (*ijob == 14) {
		goto L45;
    }
    goto L48;

/* ----------------------------------------------------------- */

L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    for (i = 1; i <= n; ++i) {
		s1 = 0.;
		s2 = 0.;
		s3 = 0.;
		for (j = max(1, i - *mlmas); j <= min(n, i + *mumas); ++j) {
			bb = fmas[i - j + linal_1.mbdiag + j * fmas_dim1];
			s1 -= bb * f1[j];
			s2 -= bb * f2[j];
			s3 -= bb * f3[j];
		}
		z1[i] += s1 * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]);
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[1], &z3[1], &ip2[1]);
    return 0;

/* ----------------------------------------------------------- */

L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    for (i = 1; i <= n; ++i) {
		s1 = 0.;
		s2 = 0.;
		s3 = 0.;
		for (j = 1; j <= n; ++j) {
			bb = fmas[i + j * fmas_dim1];
			s1 -= bb * f1[j];
			s2 -= bb * f2[j];
			s3 -= bb * f3[j];
		}
		z1[i] += s1 * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]);
    return 0;

/* ----------------------------------------------------------- */

L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (i = 1; i <= *m1; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
    for (i = 1; i <= *nm1; ++i) {
		im1 = i + *m1;
		s1 = 0.;
		s2 = 0.;
		s3 = 0.;
		for (j = 1; j <= *nm1; ++j) {
			jm1 = j + *m1;
			bb = fmas[i + j * fmas_dim1];
			s1 -= bb * f1[jm1];
			s2 -= bb * f2[jm1];
			s3 -= bb * f3[jm1];
		}
		z1[im1] += s1 * *fac1;
		z2[im1] = z2[im1] + s2 * *alphn - s3 * *betan;
		z3[im1] = z3[im1] + s3 * *alphn + s2 * *betan;
    }
    goto L48;

/* ----------------------------------------------------------- */

L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
    return 0;

/* ----------------------------------------------------------- */

L7:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
    for (i = 1; i <= n; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
    for (mm = n - 2; mm >= 1; --mm) {
		mp = n - mm;
		mp1 = mp - 1;
		i = iphes[mp];
		if (i == mp) {
			goto L746;
		}
		zsafe = z1[mp];
		z1[mp] = z1[i];
		z1[i] = zsafe;
		zsafe = z2[mp];
		z2[mp] = z2[i];
		z2[i] = zsafe;
		zsafe = z3[mp];
		z3[mp] = z3[i];
		z3[i] = zsafe;
L746:
		for (i = mp + 1; i <= n; ++i) {
			e1imp = fjac[i + mp1 * fjac_dim1];
			z1[i] -= e1imp * z1[mp];
			z2[i] -= e1imp * z2[mp];
			z3[i] -= e1imp * z3[mp];
		}
    }
    solh_(n, lde1, &e1[e1_offset], &c__1, &z1[1], &ip1[1]);
    solhc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &z2[1], &z3[1],&ip2[1]);
    for (mm = 1; mm <= n - 2; ++mm) {
		mp = n - mm;
		mp1 = mp - 1;
		for (i = mp + 1; i <= n; ++i) {
			e1imp = fjac[i + mp1 * fjac_dim1];
			z1[i] += e1imp * z1[mp];
			z2[i] += e1imp * z2[mp];
			z3[i] += e1imp * z3[mp];
		}
		i = iphes[mp];
		if (i == mp) {
			goto L750;
		}
		zsafe = z1[mp];
		z1[mp] = z1[i];
		z1[i] = zsafe;
		zsafe = z2[mp];
		z2[mp] = z2[i];
		z2[i] = zsafe;
		zsafe = z3[mp];
		z3[mp] = z3[i];
		z3[i] = zsafe;
L750:
		;
    }
    return 0;

/* ----------------------------------------------------------- */

L8:
/* ---  B=IDENTITY, SPARSE LU */
    for (i = 1; i <= n; ++i) {
		s2 = -f2[i];
		s3 = -f3[i];
		z1[i] -= f1[i] * *fac1;
		z2[i] = z2[i] + s2 * *alphn - s3 * *betan;
		z3[i] = z3[i] + s3 * *alphn + s2 * *betan;
    }
	*ier = superlu_solve_d(radau_slu_aux->slu_aux_d, &z1[1]);
	if (*ier) {
		return 0;
	}
	*ier = superlu_solve_z(radau_slu_aux->slu_aux_z, &z2[1], &z3[1]);
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* slvrad_ */

/*     END OF SUBROUTINE SLVRAD */

/* *********************************************************** */

/* Subroutine */ int estrad_(integer n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, doublereal *h__, doublereal *dd1, 
	doublereal *dd2, doublereal *dd3, FP_CB_f fcn, void* fcn_PY, integer *nfcn, doublereal 
	*y0, doublereal *y, integer *ijob, doublereal *x, integer *m1, 
	integer *m2, integer *nm1, doublereal *e1, integer *lde1, doublereal *
	z1, doublereal *z2, doublereal *z3, doublereal *cont, doublereal *
	werr, doublereal *f1, doublereal *f2, integer *ip1, integer *iphes, 
	doublereal *scal, doublereal *err, logical *first, logical *reject, 
	doublereal *fac1, doublereal *rpar, integer *ipar,
	Radau_SuperLU_aux* radau_slu_aux, integer *ier)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset;

    /* Local variables */
    static integer i, j, k, mm, mp, im1;
    extern /* Subroutine */ int sol_(integer, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum, hee1, hee2, hee3, sum1;
    extern /* Subroutine */ int solb_(integer, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solh_(integer, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal zsafe;

    /* Parameter adjustments */
    --scal;
    --iphes;
    --f2;
    --f1;
    --werr;
    --cont;
    --z3;
    --z2;
    --z1;
    --y;
    --y0;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip1;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;
    --rpar;
    --ipar;

    /* Function Body */
    hee1 = *dd1 / *h__;
    hee2 = *dd2 / *h__;
    hee3 = *dd3 / *h__;
    switch (*ijob) {
		case 1:  goto L1;
		case 2:  goto L2;
		case 3:  goto L3;
		case 4:  goto L4;
		case 5:  goto L5;
		case 6:  goto L6;
		case 7:  goto L7;
		case 8:  goto L8;
		case 10:  goto L55;
		case 11:  goto L11;
		case 12:  goto L12;
		case 13:  goto L13;
		case 14:  goto L14;
		case 15:  goto L15;
    }

L1:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
    for (i = 1; i <= n; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
    }
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L77;

L11:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (i = 1; i <= n; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
    }
L48:
    mm = *m1 / *m2;
    for (j = 1; j <= *m2; ++j) {
		sum1 = 0.;
		for (k = mm - 1; k >= 0; --k) {
			sum1 = (cont[j + k * *m2] + sum1) / *fac1;
			for (i = 1; i <= *nm1; ++i) {
				im1 = i + *m1;
				cont[im1] += fjac[i + (j + k * *m2) * fjac_dim1] * sum1;
			}
		}
    }
    sol_(*nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
    for (i = *m1; i >= 1; --i) {
		cont[i] = (cont[i] + cont[*m2 + i]) / *fac1;
    }
    goto L77;

L2:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    for (i = 1; i <= n; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[1]);
    goto L77;

L12:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    for (i = 1; i <= n; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
    }
L45:
    mm = *m1 / *m2;
    for (j = 1; j <= *m2; ++j) {
		sum1 = 0.;
		for (k = mm - 1; k >= 0; --k) {
			sum1 = (cont[j + k * *m2] + sum1) / *fac1;
			for (i = max(1, j - *mujac); i <= min(*nm1, j + *mljac); ++i) {
				im1 = i + *m1;
				cont[im1] += fjac[i + *mujac + 1 - j + (j + k * *m2) * fjac_dim1] * sum1;
			}
		}
    }
    solb_(*nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 1], &ip1[1]);
    for (i = *m1; i >= 1; --i) {
		cont[i] = (cont[i] + cont[*m2 + i]) / *fac1;
    }
    goto L77;

L3:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    for (i = 1; i <= n; ++i) {
		f1[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
    }
    for (i = 1; i <= n; ++i) {
		sum = 0.;
		for (j = max(1, i - *mlmas); j <= min(n, i + *mumas); ++j) {
			sum += fmas[i - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
		}
		f2[i] = sum;
		cont[i] = sum + y0[i];
    }
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L77;

L13:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (i = 1; i <= *m1; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
    }
    for (i = *m1 + 1; i <= n; ++i) {
		f1[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
    }
    for (i = 1; i <= *nm1; ++i) {
		sum = 0.;
		for (j = max(1, i - *mlmas); j <= min(*nm1, i + *mumas); ++j) {
			sum += fmas[i - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *m1];
		}
		im1 = i + *m1;
		f2[im1] = sum;
		cont[im1] = sum + y0[im1];
    }
    goto L48;

L4:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    for (i = 1; i <= n; ++i) {
		f1[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
    }
    for (i = 1; i <= n; ++i) {
		sum = 0.;
		for (j = max(1, i - *mlmas); j <= min(n, i + *mumas); ++j) {
			sum += fmas[i - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
		}
		f2[i] = sum;
		cont[i] = sum + y0[i];
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[1]);
    goto L77;

L14:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    for (i = 1; i <= *m1; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
    }
    for (i = *m1 + 1; i <= n; ++i) {
		f1[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
    }
    for (i = 1; i <= *nm1; ++i) {
		sum = 0.;
		for (j = max(1, i - *mlmas); j <= min(*nm1, i + *mumas); ++j) {
			sum += fmas[i - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *m1];
		}
		im1 = i + *m1;
		f2[im1] = sum;
		cont[im1] = sum + y0[im1];
    }
    goto L45;

L5:
/* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    for (i = 1; i <= n; ++i) {
		f1[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
    }
    for (i = 1; i <= n; ++i) {
		sum = 0.;
		for (j = 1; j <= n; ++j) {
			sum += fmas[i + j * fmas_dim1] * f1[j];
		}
		f2[i] = sum;
		cont[i] = sum + y0[i];
    }
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L77;

L15:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    for (i = 1; i <= *m1; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
    }
    for (i = *m1 + 1; i <= n; ++i) {
		f1[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
    }
    for (i = 1; i <= *nm1; ++i) {
		sum = 0.;
		for (j = 1; j <= *nm1; ++j) {
			sum += fmas[i + j * fmas_dim1] * f1[j + *m1];
		}
		im1 = i + *m1;
		f2[im1] = sum;
		cont[im1] = sum + y0[im1];
    }
    goto L48;

L6:
/* ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ------  THIS OPTION IS NOT PROVIDED */
    return 0;

L7:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
    for (i = 1; i <= n; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
    }
    for (mm = n - 2; mm >= 1; --mm) {
		mp = n - mm;
		i = iphes[mp];
		if (i == mp) {
			goto L310;
		}
		zsafe = cont[mp];
		cont[mp] = cont[i];
		cont[i] = zsafe;
L310:
		for (i = mp + 1; i <= n; ++i) {
			cont[i] -= fjac[i + (mp - 1) * fjac_dim1] * cont[mp];
		}
    }
    solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
    for (mm = 1; mm <= n - 2; ++mm) {
		mp = n - mm;
		for (i = mp + 1; i <= n; ++i) {
			cont[i] += fjac[i + (mp - 1) * fjac_dim1] * cont[mp];
		}
		i = iphes[mp];
		if (i == mp) {
			goto L440;
		}
		zsafe = cont[mp];
		cont[mp] = cont[i];
		cont[i] = zsafe;
L440:
		;
    }

/* ---  B=IDENTITY MATRIX, SPARSE LU */
L8:
    for (i = 1; i <= n; ++i) {
		f2[i] = hee1 * z1[i] + hee2 * z2[i] + hee3 * z3[i];
		cont[i] = f2[i] + y0[i];
    }
	*ier = superlu_solve_d(radau_slu_aux->slu_aux_d, &cont[1]);
	if (*ier){
		return 0;
	}
	goto L77;

/* -------------------------------------- */

L77:
    *err = 0.;
    for (i = 1; i <= n; ++i) {
		werr[i] = cont[i] / scal[i];
		*err += werr[i] * werr[i];
    }
    *err = max(sqrt(*err / n), 1e-10);

    if (*err < 1.) {
		return 0;
    }
    if (*first || *reject) {
		for (i = 1; i <= n; ++i) {
			cont[i] = y[i] + cont[i];
		}
		(*fcn)(n, x, &cont[1], &f1[1], &rpar[1], &ipar[1], fcn_PY);
		++(*nfcn);
		for (i = 1; i <= n; ++i) {
			cont[i] = f1[i] + f2[i];
		}
		switch (*ijob) {
			case 1:  goto L31;
			case 2:  goto L32;
			case 3:  goto L31;
			case 4:  goto L32;
			case 5:  goto L31;
			case 6:  goto L32;
			case 7:  goto L33;
			case 8:  goto L38;
			// case 9:  goto L39;
			case 10:  goto L55;
			case 11:  goto L41;
			case 12:  goto L42;
			case 13:  goto L41;
			case 14:  goto L42;
			case 15:  goto L41;
		}
	/* ------ FULL MATRIX OPTION */
L31:
		sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
		goto L88;
	/* ------ FULL MATRIX OPTION, SECOND ORDER */
L41:
		for (j = 1; j <= *m2; ++j) {
			sum1 = 0.;
			for (k = mm - 1; k >= 0; --k) {
				sum1 = (cont[j + k * *m2] + sum1) / *fac1;
				for (i = 1; i <= *nm1; ++i) {
					im1 = i + *m1;
					cont[im1] += fjac[i + (j + k * *m2) * fjac_dim1] * sum1;
				}
			}
		}
		sol_(*nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
		for (i = *m1; i >= 1; --i) {
			cont[i] = (cont[i] + cont[*m2 + i]) / *fac1;
		}
		goto L88;
	/* ------ BANDED MATRIX OPTION */
	L32:
		solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[1]);
		goto L88;
	/* ------ BANDED MATRIX OPTION, SECOND ORDER */
	L42:
		for (j = 1; j <= *m2; ++j) {
			sum1 = 0.;
			for (k = mm - 1; k >= 0; --k) {
				sum1 = (cont[j + k * *m2] + sum1) / *fac1;
				for (i = max(1, j - *mujac); i <= min(*nm1, j + *mljac); ++i) {
					im1 = i + *m1;
					cont[im1] += fjac[i + *mujac + 1 - j + (j + k * *m2) * 
						fjac_dim1] * sum1;
				}
			}
		}
		solb_(*nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 1], &ip1[1]);
		for (i = *m1; i >= 1; --i) {
			cont[i] = (cont[i] + cont[*m2 + i]) / *fac1;
		}
		goto L88;
		/* ------ HESSENBERG MATRIX OPTION */
L33:
		for (mm = n - 2; mm >= 1; --mm) {
			mp = n - mm;
			i = iphes[mp];
			if (i == mp) {
				goto L510;
			}
			zsafe = cont[mp];
			cont[mp] = cont[i];
			cont[i] = zsafe;
L510:
			for (i = mp + 1; i <= n; ++i) {
				cont[i] -= fjac[i + (mp - 1) * fjac_dim1] * cont[mp];
			}
		}
		solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
		for (mm = 1; mm <= n - 2; ++mm) {
			mp = n - mm;
			for (i = mp + 1; i <= n; ++i) {
				cont[i] += fjac[i + (mp - 1) * fjac_dim1] * cont[mp];
			}
			i = iphes[mp];
			if (i == mp) {
				goto L640;
			}
			zsafe = cont[mp];
			cont[mp] = cont[i];
			cont[i] = zsafe;
L640:
			;
		}
L38:
/* ---  B=IDENTITY, SPARSE LU */
		*ier = superlu_solve_d(radau_slu_aux->slu_aux_d, &cont[1]);
		if (*ier){
			return 0;
		}
    	goto L88;

/* ----------------------------------------------------------- */
L88:
		*err = 0.;
		for (i = 1; i <= n; ++i) {
			werr[i] = cont[i] / scal[i];
			*err += werr[i] * werr[i];
		}
		*err = max(sqrt(*err / n), 1e-10);
    }
    return 0;
/* ----------------------------------------------------------- */
L55:
    return 0;
} /* estrad_ */

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

	// allocate memory for sparse Jacobian structure
	radau_slu_aux->jac_data = (double*)malloc((nnz + n)*sizeof(double));
	radau_slu_aux->jac_indices = (int*)malloc((nnz + n)*sizeof(int));
	radau_slu_aux->jac_indptr = (int*)malloc((n + 1)*sizeof(int));

	// Create auxiliary superLU structures
	radau_slu_aux->slu_aux_d = superlu_init_d(nprocs, n, nnz);
	radau_slu_aux->slu_aux_z = superlu_init_z(nprocs, n, nnz);

	return radau_slu_aux;
}

int radau_superlu_aux_finalize(Radau_SuperLU_aux *radau_slu_aux){
	free(radau_slu_aux->jac_data);
	free(radau_slu_aux->jac_indices);
	free(radau_slu_aux->jac_indptr);

	superlu_finalize_d(radau_slu_aux->slu_aux_d);
	superlu_finalize_z(radau_slu_aux->slu_aux_z);

	free(radau_slu_aux);
	return 0;
}
