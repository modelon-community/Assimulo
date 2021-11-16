/* translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "radau_decsol_c.h"

/* Common Block Declarations */

struct {
    integer nn, nn2, nn3, nn4;
    doublereal xsol, hsol, c2m1, c1m1;
} conra5_;

#define conra5_1 conra5_

struct {
    integer mle, mue, mbjac, mbb, mdiag, mdiff, mbdiag;
} linal_;

#define linal_1 linal_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__3 = 3;
static doublereal c_b54 = .5;
static doublereal c_b91 = 81.;
static doublereal c_b92 = .33333333333333331;
static doublereal c_b93 = 9.;
static doublereal c_b103 = 1.;
static doublereal c_b114 = .8;
static doublereal c_b116 = .25;

/* Subroutine */ int radau5_c(integer *n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *
	atol, integer *itol, FP_CB_jac jac, void* jac_PY, integer *ijac, integer *mljac, integer 
	*mujac, FP_CB_mas mas, void* mas_PY, integer *imas, integer *mlmas, integer *mumas, FP_CB_solout 
	solout, void* solout_PY, integer *iout, doublereal *work, integer *lwork, integer *
	iwork, integer *liwork, doublereal *rpar, integer *ipar, integer *
	idid)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, m1, m2, nm1, nit, iee1, ief1, lde1, ief2, ief3, iey0, 
	    iez1, iez2, iez3;
    static doublereal facl;
    static integer ndec, njac;
    static doublereal facr, safe;
    static integer ijob, nfcn;
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
    extern /* Subroutine */ int radcor_(integer *, FP_CB_f, void*, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, FP_CB_jac, void*, integer *, integer *,
	     integer *, FP_CB_mas, void*, integer *, integer *, FP_CB_solout, void*, integer *, integer *
	    , integer *, doublereal *, doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, logical *, 
	    integer *, integer *, integer *, logical *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, logical *, logical 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *);
    static integer nrejct;
    static logical implct;
    static integer istore;
    static logical startn;
    static doublereal uround;

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };
    static cilist io___19 = { 0, 6, 0, 0, 0 };
    static cilist io___24 = { 0, 6, 0, 0, 0 };
    static cilist io___29 = { 0, 6, 0, 0, 0 };
    static cilist io___31 = { 0, 6, 0, 0, 0 };
    static cilist io___33 = { 0, 6, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___52 = { 0, 6, 0, 0, 0 };
    static cilist io___68 = { 0, 6, 0, 0, 0 };
    static cilist io___72 = { 0, 6, 0, 0, 0 };


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
/*                             N*(LJAC+LMAS+3*LE+12)+20 */
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
/*                             LWORK = 4*N*N+12*N+20. */
/*                 IF IWORK(9)=M1>0 THEN "LWORK" MUST BE AT LEAST */
/*                          N*(LJAC+12)+(N-M1)*(LMAS+3*LE)+20 */
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
	    s_wsle(&io___10);
	    do_lio(&c__9, &c__1, " COEFFICIENTS HAVE 20 DIGITS, UROUND=", (
		    ftnlen)37);
	    do_lio(&c__5, &c__1, (char *)&work[1], (ftnlen)sizeof(doublereal))
		    ;
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* -------- CHECK AND CHANGE THE TOLERANCES */
    expm = .66666666666666663;
    if (*itol == 0) {
	if (atol[1] <= 0. || rtol[1] <= uround * 10.) {
	    s_wsle(&io___12);
	    do_lio(&c__9, &c__1, " TOLERANCES ARE TOO SMALL", (ftnlen)25);
	    e_wsle();
	    arret = TRUE_;
	} else {
	    quot = atol[1] / rtol[1];
	    rtol[1] = pow_dd(&rtol[1], &expm) * .1;
	    atol[1] = rtol[1] * quot;
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (atol[i__] <= 0. || rtol[i__] <= uround * 10.) {
		s_wsle(&io___15);
		do_lio(&c__9, &c__1, " TOLERANCES(", (ftnlen)12);
		do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, ") ARE TOO SMALL", (ftnlen)15);
		e_wsle();
		arret = TRUE_;
	    } else {
		quot = atol[i__] / rtol[i__];
		rtol[i__] = pow_dd(&rtol[i__], &expm) * .1;
		atol[i__] = rtol[i__] * quot;
	    }
	}
    }
/* -------- NMAX , THE MAXIMAL NUMBER OF STEPS ----- */
    if (iwork[2] == 0) {
	nmax = 100000;
    } else {
	nmax = iwork[2];
	if (nmax <= 0) {
	    s_wsle(&io___17);
	    do_lio(&c__9, &c__1, " WRONG INPUT IWORK(2)=", (ftnlen)22);
	    do_lio(&c__3, &c__1, (char *)&iwork[2], (ftnlen)sizeof(integer));
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* -------- NIT    MAXIMAL NUMBER OF NEWTON ITERATIONS */
    if (iwork[3] == 0) {
	nit = 7;
    } else {
	nit = iwork[3];
	if (nit <= 0) {
	    s_wsle(&io___19);
	    do_lio(&c__9, &c__1, " CURIOUS INPUT IWORK(3)=", (ftnlen)24);
	    do_lio(&c__3, &c__1, (char *)&iwork[3], (ftnlen)sizeof(integer));
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* -------- STARTN  SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS */
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
	nind1 = *n;
    }
    if (nind1 + nind2 + nind3 != *n) {
	s_wsle(&io___24);
	do_lio(&c__9, &c__1, " CURIOUS INPUT FOR IWORK(5,6,7)=", (ftnlen)32);
	do_lio(&c__3, &c__1, (char *)&nind1, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&nind2, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&nind3, (ftnlen)sizeof(integer));
	e_wsle();
	arret = TRUE_;
    }
/* -------- PRED   STEP SIZE CONTROL */
    if (iwork[8] <= 1) {
	pred = TRUE_;
    } else {
	pred = FALSE_;
    }
/* -------- PARAMETER FOR SECOND ORDER EQUATIONS */
    m1 = iwork[9];
    m2 = iwork[10];
    nm1 = *n - m1;
    if (m1 == 0) {
	m2 = *n;
    }
    if (m2 == 0) {
	m2 = m1;
    }
    if (m1 < 0 || m2 < 0 || m1 + m2 > *n) {
	s_wsle(&io___29);
	do_lio(&c__9, &c__1, " CURIOUS INPUT FOR IWORK(9,10)=", (ftnlen)31);
	do_lio(&c__3, &c__1, (char *)&m1, (ftnlen)sizeof(integer));
	do_lio(&c__3, &c__1, (char *)&m2, (ftnlen)sizeof(integer));
	e_wsle();
	arret = TRUE_;
    }
/* --------- SAFE     SAFETY FACTOR IN STEP SIZE PREDICTION */
    if (work[2] == 0.) {
	safe = .9;
    } else {
	safe = work[2];
	if (safe <= .001 || safe >= 1.) {
	    s_wsle(&io___31);
	    do_lio(&c__9, &c__1, " CURIOUS INPUT FOR WORK(2)=", (ftnlen)27);
	    do_lio(&c__5, &c__1, (char *)&work[2], (ftnlen)sizeof(doublereal));
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; */
    if (work[3] == 0.) {
	thet = .001;
    } else {
	thet = work[3];
	if (thet >= 1.) {
	    s_wsle(&io___33);
	    do_lio(&c__9, &c__1, " CURIOUS INPUT FOR WORK(3)=", (ftnlen)27);
	    do_lio(&c__5, &c__1, (char *)&work[3], (ftnlen)sizeof(doublereal))
		    ;
	    e_wsle();
	    arret = TRUE_;
	}
    }
/* --- FNEWT   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. */
    tolst = rtol[1];
    if (work[4] == 0.) {
/* Computing MAX */
/* Computing MIN */
	d__3 = .03, d__4 = pow_dd(&tolst, &c_b54);
	d__1 = uround * 10 / tolst, d__2 = min(d__3,d__4);
	fnewt = max(d__1,d__2);
    } else {
	fnewt = work[4];
	if (fnewt <= uround / tolst) {
	    s_wsle(&io___36);
	    do_lio(&c__9, &c__1, " CURIOUS INPUT FOR WORK(4)=", (ftnlen)27);
	    do_lio(&c__5, &c__1, (char *)&work[4], (ftnlen)sizeof(doublereal))
		    ;
	    e_wsle();
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
	s_wsle(&io___39);
	do_lio(&c__9, &c__1, " CURIOUS INPUT FOR WORK(5,6)=", (ftnlen)29);
	do_lio(&c__5, &c__1, (char *)&quot1, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&quot2, (ftnlen)sizeof(doublereal));
	e_wsle();
	arret = TRUE_;
    }
/* -------- MAXIMAL STEP SIZE */
    if (work[7] == 0.) {
	hmax = *xend - *x;
    } else {
	hmax = work[7];
    }
/* -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION */
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
	s_wsle(&io___43);
	do_lio(&c__9, &c__1, " CURIOUS INPUT WORK(8,9)=", (ftnlen)25);
	do_lio(&c__5, &c__1, (char *)&work[8], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&work[9], (ftnlen)sizeof(doublereal));
	e_wsle();
	arret = TRUE_;
    }
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/*         COMPUTATION OF ARRAY ENTRIES */
/* *** *** *** *** *** *** *** *** *** *** *** *** *** */
/* ---- IMPLICIT, BANDED OR NOT ? */
    implct = *imas != 0;
    jband = *mljac < nm1;
/* -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- */
/* -- JACOBIAN  AND  MATRICES E1, E2 */
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
	    s_wsle(&io___50);
	    do_lio(&c__9, &c__1, "BANDWITH OF \"MAS\" NOT SMALLER THAN BANDW"
		    "ITH OF \"JAC\"", (ftnlen)52);
	    e_wsle();
	    arret = TRUE_;
	}
    } else {
	ldmas = 0;
	if (jband) {
	    ijob = 2;
	} else {
	    ijob = 1;
	    if (*n > 2 && iwork[1] != 0) {
		ijob = 7;
	    }
	}
    }
    ldmas2 = max(1,ldmas);
/* ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN */
    if ((implct || jband) && ijob == 7) {
	s_wsle(&io___52);
	do_lio(&c__9, &c__1, " HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS"
		" WITH FULL JACOBIAN", (ftnlen)65);
	e_wsle();
	arret = TRUE_;
    }
/* ------- PREPARE THE ENTRY-POINTS FOR THE ARRAYS IN WORK ----- */
    iez1 = 21;
    iez2 = iez1 + *n;
    iez3 = iez2 + *n;
    iey0 = iez3 + *n;
    iescal = iey0 + *n;
    ief1 = iescal + *n;
    ief2 = ief1 + *n;
    ief3 = ief2 + *n;
    iecon = ief3 + *n;
    iejac = iecon + (*n << 2);
    iemas = iejac + *n * ldjac;
    iee1 = iemas + nm1 * ldmas;
    iee2r = iee1 + nm1 * lde1;
    iee2i = iee2r + nm1 * lde1;
/* ------ TOTAL STORAGE REQUIREMENT ----------- */
    istore = iee2i + nm1 * lde1 - 1;
    if (istore > *lwork) {
	s_wsle(&io___68);
	do_lio(&c__9, &c__1, " INSUFFICIENT STORAGE FOR WORK, MIN. LWORK=", (
		ftnlen)43);
	do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(integer));
	e_wsle();
	arret = TRUE_;
    }
/* ------- ENTRY POINTS FOR INTEGER WORKSPACE ----- */
    ieip1 = 21;
    ieip2 = ieip1 + nm1;
    ieiph = ieip2 + nm1;
/* --------- TOTAL REQUIREMENT --------------- */
    istore = ieiph + nm1 - 1;
    if (istore > *liwork) {
	s_wsle(&io___72);
	do_lio(&c__9, &c__1, " INSUFF. STORAGE FOR IWORK, MIN. LIWORK=", (
		ftnlen)40);
	do_lio(&c__3, &c__1, (char *)&istore, (ftnlen)sizeof(integer));
	e_wsle();
	arret = TRUE_;
    }
/* ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 */
    if (arret) {
	*idid = -1;
	return 0;
    }
/* -------- CALL TO CORE INTEGRATOR ------------ */
    radcor_(n, (FP_CB_f)fcn, fcn_PY, x, &y[1], xend, &hmax, h__, &rtol[1], &atol[1], 
	    itol, (FP_CB_jac)jac, jac_PY, ijac, mljac, mujac, (FP_CB_mas)mas, mas_PY, mlmas, mumas, (
	    FP_CB_solout)solout, solout_PY, iout, idid, &nmax, &uround, &safe, &thet, &fnewt, &
	    quot1, &quot2, &nit, &ijob, &startn, &nind1, &nind2, &nind3, &
	    pred, &facl, &facr, &m1, &m2, &nm1, &implct, &jband, &ldjac, &
	    lde1, &ldmas2, &work[iez1], &work[iez2], &work[iez3], &work[iey0],
	     &work[iescal], &work[ief1], &work[ief2], &work[ief3], &work[
	    iejac], &work[iee1], &work[iee2r], &work[iee2i], &work[iemas], &
	    iwork[ieip1], &iwork[ieip2], &iwork[ieiph], &work[iecon], &nfcn, &
	    njac, &nstep, &naccpt, &nrejct, &ndec, &nsol, &rpar[1], &ipar[1]);
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
	d__1 = rtol[1] * 10.;
	rtol[1] = pow_dd(&d__1, &expm);
	atol[1] = rtol[1] * quot;
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    quot = atol[i__] / rtol[i__];
	    d__1 = rtol[i__] * 10.;
	    rtol[i__] = pow_dd(&d__1, &expm);
	    atol[i__] = rtol[i__] * quot;
	}
    }
/* ----------- RETURN ----------- */
    return 0;
} /* radau5_ */


/*     END OF SUBROUTINE RADAU5 */

/* *********************************************************** */

/* Subroutine */ int radcor_(integer *n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *hmax, doublereal *h__, doublereal *
	rtol, doublereal *atol, integer *itol, FP_CB_jac jac, void* jac_PY, integer *ijac, 
	integer *mljac, integer *mujac, FP_CB_mas mas, void* mas_PY, integer *mlmas, integer *
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
	integer *ndec, integer *nsol, doublereal *rpar, integer *ipar)
{
    /* Format strings */
    static char fmt_979[] = "(\002 EXIT OF RADAU5 AT X=\002,e18.4)";

    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, e2r_dim1, e2r_offset, e2i_dim1, e2i_offset, i__1, i__2,
	     i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), pow_dd(doublereal *, doublereal *), d_sign(
	    doublereal *, doublereal *), pow_di(doublereal *, integer *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);

    /* Local variables */
    static integer i__, j, k, l;
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
    static doublereal dyno, dyth, quot, hhfac, betan, alphn, denom, theta, 
	    ysafe, hmaxn;
    static integer nsing;
    static logical first;
    static integer irtrn, nrsol, nsolu;
    static doublereal qnewt, xosol, acont3;
    static logical index1, index2, index3, caljac;
    static doublereal faccon;
    extern /* Subroutine */ int decomc_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *);
    static logical calhes;
    static doublereal erracc;
    static integer mujacj;
    extern /* Subroutine */ int decomr_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, logical *, integer *);
    static logical reject;
    static doublereal facgus;
    static integer mujacp;
    extern /* Subroutine */ int estrad_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, doublereal *, FP_CB_f, void*, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    doublereal *, doublereal *, logical *, logical *, doublereal *, 
	    doublereal *, integer *);
    static doublereal dynold, posneg;
    extern /* Subroutine */ int slvrad_(integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, integer *);
    static doublereal thqold;

    /* Fortran I/O blocks */
    static cilist io___178 = { 0, 6, 0, fmt_979, 0 };
    static cilist io___179 = { 0, 6, 0, 0, 0 };
    static cilist io___180 = { 0, 6, 0, fmt_979, 0 };
    static cilist io___181 = { 0, 6, 0, 0, 0 };
    static cilist io___182 = { 0, 6, 0, fmt_979, 0 };
    static cilist io___183 = { 0, 6, 0, 0, 0 };
    static cilist io___184 = { 0, 6, 0, fmt_979, 0 };
    static cilist io___185 = { 0, 6, 0, 0, 0 };


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
    doublereal *werr = (doublereal*) malloc(*n * sizeof(doublereal));
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
    conra5_1.nn = *n;
    conra5_1.nn2 = *n << 1;
    conra5_1.nn3 = *n * 3;
    lrc = *n << 2;
/* -------- CHECK THE INDEX OF THE PROBLEM ----- */
    index1 = *nind1 != 0;
    index2 = *nind2 != 0;
    index3 = *nind3 != 0;
/* ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- */
    if (*implct) {
    (*mas)(nm1, &fmas[fmas_offset], ldmas, &rpar[1], &ipar[1], mas_PY);
    }
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
    u1 = (pow_dd(&c_b91, &c_b92) + 6. - pow_dd(&c_b93, &c_b92)) / 30.;
    alph = (12. - pow_dd(&c_b91, &c_b92) + pow_dd(&c_b93, &c_b92)) / 60.;
    beta = (pow_dd(&c_b91, &c_b92) + pow_dd(&c_b93, &c_b92)) * sqrt(3.) / 60.;
/* Computing 2nd power */
    d__1 = alph;
/* Computing 2nd power */
    d__2 = beta;
    cno = d__1 * d__1 + d__2 * d__2;
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
    d__1 = *xend - *x;
    posneg = d_sign(&c_b103, &d__1);
/* Computing MIN */
    d__2 = abs(*hmax), d__3 = (d__1 = *xend - *x, abs(d__1));
    hmaxn = min(d__2,d__3);
    if (abs(*h__) <= *uround * 10.) {
	*h__ = 1e-6;
    }
/* Computing MIN */
    d__1 = abs(*h__);
    *h__ = min(d__1,hmaxn);
    *h__ = d_sign(h__, &posneg);
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
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
 	    werr[i__] = 0.;
	    cont[i__] = y[i__];
	}
	nsolu = *n;
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
    n2 = *n << 1;
    n3 = *n * 3;
    if (*itol == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scal[i__] = atol[1] + rtol[1] * (d__1 = y[i__], abs(d__1));
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scal[i__] = atol[i__] + rtol[i__] * (d__1 = y[i__], abs(d__1));
	}
    }
    hhfac = *h__;
    (*fcn)(n, x, &y[1], &y0[1], &rpar[1], &ipar[1], fcn_PY);
    ++(*nfcn);
/* --- BASIC INTEGRATION STEP */
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
	    i__1 = *m1 / *m2 + 1;
	    for (mm = 1; mm <= i__1; ++mm) {
		i__2 = md;
		for (k = 1; k <= i__2; ++k) {
		    j = k + (mm - 1) * *m2;
L12:
		    f1[j] = y[j];
/* Computing MAX */
		    d__2 = 1e-5, d__3 = (d__1 = y[j], abs(d__1));
		    f2[j] = sqrt(*uround * max(d__2,d__3));
		    y[j] += f2[j];
		    j += md;
		    if (j <= mm * *m2) {
			goto L12;
		    }
            (*fcn)(n, x, &y[1], &cont[1], &rpar[1], &ipar[1], fcn_PY);
		    j = k + (mm - 1) * *m2;
		    j1 = k;
/* Computing MAX */
		    i__3 = 1, i__4 = j1 - *mujac;
		    lbeg = max(i__3,i__4) + *m1;
L14:
/* Computing MIN */
		    i__3 = *m2, i__4 = j1 + *mljac;
		    lend = min(i__3,i__4) + *m1;
		    y[j] = f1[j];
		    mujacj = mujacp - j1 - *m1;
		    i__3 = lend;
		    for (l = lbeg; l <= i__3; ++l) {
			fjac[l + mujacj + j * fjac_dim1] = (cont[l] - y0[l]) /
				 f2[j];
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
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		ysafe = y[i__];
/* Computing MAX */
		d__1 = 1e-5, d__2 = abs(ysafe);
		delt = sqrt(*uround * max(d__1,d__2));
		y[i__] = ysafe + delt;
        (*fcn)(n, x, &y[1], &cont[1], &rpar[1], &ipar[1], fcn_PY);
		if (ipar[1] < 0) {
		    y[i__] = ysafe - delt;
            (*fcn)(n, x, &y[1], &cont[1], &rpar[1], &ipar[1], fcn_PY);
		    if (ipar[1] < 0) {
			y[i__] = ysafe;
			goto L79;
		    }
		    i__2 = *n;
		    for (j = *m1 + 1; j <= i__2; ++j) {
			fjac[j - *m1 + i__ * fjac_dim1] = (y0[j] - cont[j]) / 
				delt;
		    }
		} else {
		    i__2 = *n;
		    for (j = *m1 + 1; j <= i__2; ++j) {
			fjac[j - *m1 + i__ * fjac_dim1] = (cont[j] - y0[j]) / 
				delt;
		    }
		}
		y[i__] = ysafe;
	    }
	}
    } else {
/* --- COMPUTE JACOBIAN MATRIX ANALYTICALLY */
    (*jac)(n, x, &y[1], &fjac[fjac_offset], ldjac, &rpar[1], &ipar[1], jac_PY);
    }
    caljac = TRUE_;
    calhes = TRUE_;
L20:
/* --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS */
    fac1 = u1 / *h__;
    alphn = alph / *h__;
    betan = beta / *h__;
    decomr_(n, &fjac[fjac_offset], ldjac, &fmas[fmas_offset], ldmas, mlmas, 
	    mumas, m1, m2, nm1, &fac1, &e1[e1_offset], lde1, &ip1[1], &ier, 
	    ijob, &calhes, &iphes[1]);
    if (ier != 0) {
	goto L78;
    }
    decomc_(n, &fjac[fjac_offset], ldjac, &fmas[fmas_offset], ldmas, mlmas, 
	    mumas, m1, m2, nm1, &alphn, &betan, &e2r[e2r_offset], &e2i[
	    e2i_offset], lde1, &ip2[1], &ier, ijob);
    if (ier != 0) {
	goto L78;
    }
    ++(*ndec);
L30:
    ++(*nstep);
    if (*nstep > *nmax) {
	goto L178;
    }
    if (abs(*h__) * .1 <= abs(*x) * *uround) {
	goto L177;
    }
    if (index2) {
	i__1 = *nind1 + *nind2;
	for (i__ = *nind1 + 1; i__ <= i__1; ++i__) {
	    scal[i__] /= hhfac;
	}
    }
    if (index3) {
	i__1 = *nind1 + *nind2 + *nind3;
	for (i__ = *nind1 + *nind2 + 1; i__ <= i__1; ++i__) {
	    scal[i__] /= hhfac * hhfac;
	}
    }
    xph = *x + *h__;
/* *** *** *** *** *** *** *** */
/*  STARTING VALUES FOR NEWTON ITERATION */
/* *** *** *** *** *** *** *** */
    if (first || *startn) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z1[i__] = 0.;
	    z2[i__] = 0.;
	    z3[i__] = 0.;
	    f1[i__] = 0.;
	    f2[i__] = 0.;
	    f3[i__] = 0.;
	}
    } else {
	c3q = *h__ / hold;
	c1q = c1 * c3q;
	c2q = c2 * c3q;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ak1 = cont[i__ + *n];
	    ak2 = cont[i__ + n2];
	    ak3 = cont[i__ + n3];
	    z1i = c1q * (ak1 + (c1q - conra5_1.c2m1) * (ak2 + (c1q - 
		    conra5_1.c1m1) * ak3));
	    z2i = c2q * (ak1 + (c2q - conra5_1.c2m1) * (ak2 + (c2q - 
		    conra5_1.c1m1) * ak3));
	    z3i = c3q * (ak1 + (c3q - conra5_1.c2m1) * (ak2 + (c3q - 
		    conra5_1.c1m1) * ak3));
	    z1[i__] = z1i;
	    z2[i__] = z2i;
	    z3[i__] = z3i;
	    f1[i__] = ti11 * z1i + ti12 * z2i + ti13 * z3i;
	    f2[i__] = ti21 * z1i + ti22 * z2i + ti23 * z3i;
	    f3[i__] = ti31 * z1i + ti32 * z2i + ti33 * z3i;
	}
    }
/* *** *** *** *** *** *** *** */
/*  LOOP FOR THE SIMPLIFIED NEWTON ITERATION */
/* *** *** *** *** *** *** *** */
    newt = 0;
    d__1 = max(faccon,*uround);
    faccon = pow_dd(&d__1, &c_b114);
    theta = abs(*thet);
L40:
    if (newt >= *nit) {
	goto L78;
    }
/* ---     COMPUTE THE RIGHT-HAND SIDE */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cont[i__] = y[i__] + z1[i__];
    }
    d__1 = *x + c1 * *h__;
    (*fcn)(n, &d__1, &cont[1], &z1[1], &rpar[1], &ipar[1], fcn_PY);
    ++(*nfcn);
    if (ipar[1] < 0) {
	goto L79;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cont[i__] = y[i__] + z2[i__];
    }
    d__1 = *x + c2 * *h__;
    (*fcn)(n, &d__1, &cont[1], &z2[1], &rpar[1], &ipar[1], fcn_PY);
    ++(*nfcn);
    if (ipar[1] < 0) {
	goto L79;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	cont[i__] = y[i__] + z3[i__];
    }
    (*fcn)(n, &xph, &cont[1], &z3[1], &rpar[1], &ipar[1], fcn_PY);
    ++(*nfcn);
    if (ipar[1] < 0) {
	goto L79;
    }
/* ---     SOLVE THE LINEAR SYSTEMS */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a1 = z1[i__];
	a2 = z2[i__];
	a3 = z3[i__];
	z1[i__] = ti11 * a1 + ti12 * a2 + ti13 * a3;
	z2[i__] = ti21 * a1 + ti22 * a2 + ti23 * a3;
	z3[i__] = ti31 * a1 + ti32 * a2 + ti33 * a3;
    }
    slvrad_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset], 
	    ldmas, mlmas, mumas, m1, m2, nm1, &fac1, &alphn, &betan, &e1[
	    e1_offset], &e2r[e2r_offset], &e2i[e2i_offset], lde1, &z1[1], &z2[
	    1], &z3[1], &f1[1], &f2[1], &f3[1], &cont[1], &ip1[1], &ip2[1], &
	    iphes[1], &ier, ijob);
    ++(*nsol);
    ++newt;
    dyno = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	denom = scal[i__];
/* Computing 2nd power */
	d__1 = z1[i__] / denom;
/* Computing 2nd power */
	d__2 = z2[i__] / denom;
/* Computing 2nd power */
	d__3 = z3[i__] / denom;
	dyno = dyno + d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
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
	    i__1 = *nit - 1 - newt;
	    dyth = faccon * dyno * pow(theta, i__1) / *fnewt;
	    if (dyth >= 1.) {
/* Computing MAX */
		d__1 = 1e-4, d__2 = min(20.,dyth);
		qnewt = max(d__1,d__2);
		d__1 = -1. / (*nit + 4. - 1 - newt);
		hhfac = pow_dd(&qnewt, &d__1) * .8;
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
    dynold = max(dyno,*uround);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f1i = f1[i__] + z1[i__];
	f2i = f2[i__] + z2[i__];
	f3i = f3[i__] + z3[i__];
	f1[i__] = f1i;
	f2[i__] = f2i;
	f3[i__] = f3i;
	z1[i__] = t11 * f1i + t12 * f2i + t13 * f3i;
	z2[i__] = t21 * f1i + t22 * f2i + t23 * f3i;
	z3[i__] = t31 * f1i + f2i;
    }
    if (faccon * dyno > *fnewt) {
	goto L40;
    }
/* --- ERROR ESTIMATION */
    estrad_(n, &fjac[fjac_offset], ldjac, mljac, mujac, &fmas[fmas_offset], 
	    ldmas, mlmas, mumas, h__, &dd1, &dd2, &dd3, (FP_CB_f) fcn, fcn_PY, nfcn, &y0[
	    1], &y[1], ijob, x, m1, m2, nm1, &e1[e1_offset], lde1, &z1[1], &
	    z2[1], &z3[1], &cont[1], &werr[1], &f1[1], &f2[1], &ip1[1], &
	    iphes[1], &scal[1], &err, &first, &reject, &fac1, &rpar[1], &ipar[
	    1]);
/* --- COMPUTATION OF HNEW */
/* --- WE REQUIRE .2<=HNEW/H<=8. */
/* Computing MIN */
    d__1 = *safe, d__2 = cfac / (newt + (*nit << 1));
    fac = min(d__1,d__2);
/* Computing MAX */
/* Computing MIN */
    d__3 = *facl, d__4 = pow_dd(&err, &c_b116) / fac;
    d__1 = *facr, d__2 = min(d__3,d__4);
    quot = max(d__1,d__2);
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
/* Computing 2nd power */
		d__2 = err;
		d__1 = d__2 * d__2 / erracc;
		facgus = hacc / *h__ * pow_dd(&d__1, &c_b116) / *safe;
/* Computing MAX */
		d__1 = *facr, d__2 = min(*facl,facgus);
		facgus = max(d__1,d__2);
		quot = max(quot,facgus);
		hnew = *h__ / quot;
	    }
	    hacc = *h__;
	    erracc = max(.01,err);
	}
	xold = *x;
	hold = *h__;
	*x = xph;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    y[i__] += z3[i__];
	    z2i = z2[i__];
	    z1i = z1[i__];
	    cont[i__ + *n] = (z2i - z3[i__]) / conra5_1.c2m1;
	    ak = (z1i - z2i) / c1mc2;
	    acont3 = z1i / c1;
	    acont3 = (ak - acont3) / c2;
	    cont[i__ + n2] = (ak - cont[i__ + *n]) / conra5_1.c1m1;
	    cont[i__ + n3] = cont[i__ + n2] - acont3;
	}
	if (*itol == 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		scal[i__] = atol[1] + rtol[1] * (d__1 = y[i__], abs(d__1));
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		scal[i__] = atol[i__] + rtol[i__] * (d__1 = y[i__], abs(d__1))
			;
	    }
	}
	if (*iout != 0) {
	    nrsol = *naccpt + 1;
	    conra5_1.xsol = *x;
	    xosol = xold;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		cont[i__] = y[i__];
	    }
	    nsolu = *n;
	    conra5_1.hsol = hold;
	    (*solout)(&nrsol, &xosol, &conra5_1.xsol, &y[1], &cont[1], &werr[
		    1], &lrc, &nsolu, &rpar[1], &ipar[1], &irtrn, solout_PY);
	    if (irtrn < 0) {
		goto L179;
	    }
	}
	caljac = FALSE_;
	if (last) {
	    *h__ = hopt;
	    *idid = 1;
	    return 0;
	}
	(*fcn)(n, x, &y[1], &y0[1], &rpar[1], &ipar[1], fcn_PY);
	++(*nfcn);
/* Computing MIN */
	d__1 = abs(hnew);
	hnew = posneg * min(d__1,hmaxn);
	hopt = hnew;
	hopt = min(*h__,hnew);
	if (reject) {
/* Computing MIN */
	    d__1 = abs(hnew), d__2 = abs(*h__);
	    hnew = posneg * min(d__1,d__2);
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
/* --- UNEXPECTED STEP-REJECTION */
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
L79:
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
/* --- FAIL EXIT */
L175:
    s_wsfe(&io___178);
    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsle(&io___179);
    do_lio(&c__9, &c__1, " REPEATEDLY UNEXPECTED STEP REJECTIONS", (ftnlen)38)
	    ;
    e_wsle();
    *idid = -5;
    return 0;
L176:
    s_wsfe(&io___180);
    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsle(&io___181);
    do_lio(&c__9, &c__1, " MATRIX IS REPEATEDLY SINGULAR, IER=", (ftnlen)36);
    do_lio(&c__3, &c__1, (char *)&ier, (ftnlen)sizeof(integer));
    e_wsle();
    *idid = -4;
    return 0;
L177:
    s_wsfe(&io___182);
    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsle(&io___183);
    do_lio(&c__9, &c__1, " STEP SIZE T0O SMALL, H=", (ftnlen)24);
    do_lio(&c__5, &c__1, (char *)&(*h__), (ftnlen)sizeof(doublereal));
    e_wsle();
    *idid = -3;
    return 0;
L178:
    s_wsfe(&io___184);
    do_fio(&c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsle(&io___185);
    do_lio(&c__9, &c__1, " MORE THAN NMAX =", (ftnlen)17);
    do_lio(&c__3, &c__1, (char *)&(*nmax), (ftnlen)sizeof(integer));
    do_lio(&c__9, &c__1, "STEPS ARE NEEDED", (ftnlen)16);
    e_wsle();
    *idid = -2;
    return 0;
/* --- EXIT CAUSED BY SOLOUT */
L179:
/*      WRITE(6,979)X */
    *idid = 2;
    return 0;
} /* radcor_ */


/*     END OF SUBROUTINE RADCOR */

/* *********************************************************** */

doublereal contr5_c(integer *i__, doublereal *x, doublereal *cont, integer *
	lrc)
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
    ret_val = cont[*i__] + s * (cont[*i__ + conra5_1.nn] + (s - conra5_1.c2m1)
	     * (cont[*i__ + conra5_1.nn2] + (s - conra5_1.c1m1) * cont[*i__ + 
	    conra5_1.nn3]));
    return ret_val;
} /* contr5_ */


/*     END OF FUNCTION CONTR5 */

/* *********************************************************** */
/* Subroutine */ int dec_(integer *n, integer *ndim, doublereal *a, integer *
	ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, m;
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
    ip[*n] = 1;
    if (*n == 1) {
	goto L70;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = k;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    if ((d__1 = a[i__ + k * a_dim1], abs(d__1)) > (d__2 = a[m + k * 
		    a_dim1], abs(d__2))) {
		m = i__;
	    }
/* L10: */
	}
	ip[k] = m;
	t = a[m + k * a_dim1];
	if (m == k) {
	    goto L20;
	}
	ip[*n] = -ip[*n];
	a[m + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L20:
	if (t == 0.) {
	    goto L80;
	}
	t = 1. / t;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* L30: */
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[m + j * a_dim1];
	    a[m + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
	    if (t == 0.) {
		goto L45;
	    }
	    i__3 = *n;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
/* L40: */
		a[i__ + j * a_dim1] += a[i__ + k * a_dim1] * t;
	    }
L45:
/* L50: */
	    ;
	}
/* L60: */
    }
L70:
    k = *n;
    if (a[*n + *n * a_dim1] == 0.) {
	goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[*n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DEC ------------------------- */
} /* dec_ */



/* Subroutine */ int sol_(integer *n, integer *ndim, doublereal *a, 
	doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, m;
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
    if (*n == 1) {
	goto L50;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = ip[k];
	t = b[m];
	b[m] = b[k];
	b[k] = t;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* L10: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/* L20: */
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	km1 = *n - kb;
	k = km1 + 1;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L30: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/* L40: */
    }
L50:
    b[1] /= a[a_dim1 + 1];
    return 0;
/* ----------------------- END OF SUBROUTINE SOL ------------------------- */
} /* sol_ */



/* Subroutine */ int dech_(integer *n, integer *ndim, doublereal *a, integer *
	lb, integer *ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, m;
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
    ip[*n] = 1;
    if (*n == 1) {
	goto L70;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = k;
/* Computing MIN */
	i__2 = *n, i__3 = *lb + k;
	na = min(i__2,i__3);
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    if ((d__1 = a[i__ + k * a_dim1], abs(d__1)) > (d__2 = a[m + k * 
		    a_dim1], abs(d__2))) {
		m = i__;
	    }
/* L10: */
	}
	ip[k] = m;
	t = a[m + k * a_dim1];
	if (m == k) {
	    goto L20;
	}
	ip[*n] = -ip[*n];
	a[m + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L20:
	if (t == 0.) {
	    goto L80;
	}
	t = 1. / t;
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* L30: */
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[m + j * a_dim1];
	    a[m + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
	    if (t == 0.) {
		goto L45;
	    }
	    i__3 = na;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
/* L40: */
		a[i__ + j * a_dim1] += a[i__ + k * a_dim1] * t;
	    }
L45:
/* L50: */
	    ;
	}
/* L60: */
    }
L70:
    k = *n;
    if (a[*n + *n * a_dim1] == 0.) {
	goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[*n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECH ------------------------ */
} /* dech_ */



/* Subroutine */ int solh_(integer *n, integer *ndim, doublereal *a, integer *
	lb, doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, m;
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
    if (*n == 1) {
	goto L50;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = ip[k];
	t = b[m];
	b[m] = b[k];
	b[k] = t;
/* Computing MIN */
	i__2 = *n, i__3 = *lb + k;
	na = min(i__2,i__3);
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* L10: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/* L20: */
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	km1 = *n - kb;
	k = km1 + 1;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L30: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/* L40: */
    }
L50:
    b[1] /= a[a_dim1 + 1];
    return 0;
/* ----------------------- END OF SUBROUTINE SOLH ------------------------ */
} /* solh_ */


/* Subroutine */ int decc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *ip, integer *ier)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k, m;
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
    ip[*n] = 1;
    if (*n == 1) {
	goto L70;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = k;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    if ((d__1 = ar[i__ + k * ar_dim1], abs(d__1)) + (d__2 = ai[i__ + 
		    k * ai_dim1], abs(d__2)) > (d__3 = ar[m + k * ar_dim1], 
		    abs(d__3)) + (d__4 = ai[m + k * ai_dim1], abs(d__4))) {
		m = i__;
	    }
/* L10: */
	}
	ip[k] = m;
	tr = ar[m + k * ar_dim1];
	ti = ai[m + k * ai_dim1];
	if (m == k) {
	    goto L20;
	}
	ip[*n] = -ip[*n];
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
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
	    ar[i__ + k * ar_dim1] = -prodr;
	    ai[i__ + k * ai_dim1] = -prodi;
/* L30: */
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
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
		i__3 = *n;
		for (i__ = kp1; i__ <= i__3; ++i__) {
		    prodr = ar[i__ + k * ar_dim1] * tr;
		    prodi = ai[i__ + k * ai_dim1] * tr;
		    ar[i__ + j * ar_dim1] += prodr;
		    ai[i__ + j * ai_dim1] += prodi;
/* L40: */
		}
		goto L48;
	    }
	    if (tr == 0.) {
		i__3 = *n;
		for (i__ = kp1; i__ <= i__3; ++i__) {
		    prodr = -ai[i__ + k * ai_dim1] * ti;
		    prodi = ar[i__ + k * ar_dim1] * ti;
		    ar[i__ + j * ar_dim1] += prodr;
		    ai[i__ + j * ai_dim1] += prodi;
/* L45: */
		}
		goto L48;
	    }
	    i__3 = *n;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
		prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * 
			ti;
		prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * 
			ti;
		ar[i__ + j * ar_dim1] += prodr;
		ai[i__ + j * ai_dim1] += prodi;
/* L47: */
	    }
L48:
/* L50: */
	    ;
	}
/* L60: */
    }
L70:
    k = *n;
    if ((d__1 = ar[*n + *n * ar_dim1], abs(d__1)) + (d__2 = ai[*n + *n * 
	    ai_dim1], abs(d__2)) == 0.) {
	goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[*n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECC ------------------------ */
} /* decc_ */



/* Subroutine */ int solc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, doublereal *br, doublereal *bi, integer *ip)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, m, kb;
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
    if (*n == 1) {
	goto L50;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = ip[k];
	tr = br[m];
	ti = bi[m];
	br[m] = br[k];
	bi[m] = bi[k];
	br[k] = tr;
	bi[k] = ti;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
	    br[i__] += prodr;
	    bi[i__] += prodi;
/* L10: */
	}
/* L20: */
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	km1 = *n - kb;
	k = km1 + 1;
	den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] + ai[k + k * ai_dim1] 
		* ai[k + k * ai_dim1];
	prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
	prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
	br[k] = prodr / den;
	bi[k] = prodi / den;
	tr = -br[k];
	ti = -bi[k];
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
	    br[i__] += prodr;
	    bi[i__] += prodi;
/* L30: */
	}
/* L40: */
    }
L50:
    den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 + 
	    1];
    prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
    prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
    br[1] = prodr / den;
    bi[1] = prodi / den;
    return 0;
/* ----------------------- END OF SUBROUTINE SOLC ------------------------ */
} /* solc_ */



/* Subroutine */ int dechc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *lb, integer *ip, integer *ier)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k, m, na;
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
    ip[*n] = 1;
    if (*lb == 0) {
	goto L70;
    }
    if (*n == 1) {
	goto L70;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = k;
/* Computing MIN */
	i__2 = *n, i__3 = *lb + k;
	na = min(i__2,i__3);
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    if ((d__1 = ar[i__ + k * ar_dim1], abs(d__1)) + (d__2 = ai[i__ + 
		    k * ai_dim1], abs(d__2)) > (d__3 = ar[m + k * ar_dim1], 
		    abs(d__3)) + (d__4 = ai[m + k * ai_dim1], abs(d__4))) {
		m = i__;
	    }
/* L10: */
	}
	ip[k] = m;
	tr = ar[m + k * ar_dim1];
	ti = ai[m + k * ai_dim1];
	if (m == k) {
	    goto L20;
	}
	ip[*n] = -ip[*n];
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
	i__2 = na;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
	    ar[i__ + k * ar_dim1] = -prodr;
	    ai[i__ + k * ai_dim1] = -prodi;
/* L30: */
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
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
		i__3 = na;
		for (i__ = kp1; i__ <= i__3; ++i__) {
		    prodr = ar[i__ + k * ar_dim1] * tr;
		    prodi = ai[i__ + k * ai_dim1] * tr;
		    ar[i__ + j * ar_dim1] += prodr;
		    ai[i__ + j * ai_dim1] += prodi;
/* L40: */
		}
		goto L48;
	    }
	    if (tr == 0.) {
		i__3 = na;
		for (i__ = kp1; i__ <= i__3; ++i__) {
		    prodr = -ai[i__ + k * ai_dim1] * ti;
		    prodi = ar[i__ + k * ar_dim1] * ti;
		    ar[i__ + j * ar_dim1] += prodr;
		    ai[i__ + j * ai_dim1] += prodi;
/* L45: */
		}
		goto L48;
	    }
	    i__3 = na;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
		prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * 
			ti;
		prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * 
			ti;
		ar[i__ + j * ar_dim1] += prodr;
		ai[i__ + j * ai_dim1] += prodi;
/* L47: */
	    }
L48:
/* L50: */
	    ;
	}
/* L60: */
    }
L70:
    k = *n;
    if ((d__1 = ar[*n + *n * ar_dim1], abs(d__1)) + (d__2 = ai[*n + *n * 
	    ai_dim1], abs(d__2)) == 0.) {
	goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[*n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECHC ----------------------- */
} /* dechc_ */



/* Subroutine */ int solhc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *lb, doublereal *br, doublereal *bi, integer *
	ip)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, k, m, kb;
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
    if (*n == 1) {
	goto L50;
    }
    nm1 = *n - 1;
    if (*lb == 0) {
	goto L25;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = ip[k];
	tr = br[m];
	ti = bi[m];
	br[m] = br[k];
	bi[m] = bi[k];
	br[k] = tr;
	bi[k] = ti;
/* Computing MIN */
	i__3 = *n, i__4 = *lb + k;
	i__2 = min(i__3,i__4);
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
	    br[i__] += prodr;
	    bi[i__] += prodi;
/* L10: */
	}
/* L20: */
    }
L25:
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	km1 = *n - kb;
	k = km1 + 1;
	den = ar[k + k * ar_dim1] * ar[k + k * ar_dim1] + ai[k + k * ai_dim1] 
		* ai[k + k * ai_dim1];
	prodr = br[k] * ar[k + k * ar_dim1] + bi[k] * ai[k + k * ai_dim1];
	prodi = bi[k] * ar[k + k * ar_dim1] - br[k] * ai[k + k * ai_dim1];
	br[k] = prodr / den;
	bi[k] = prodi / den;
	tr = -br[k];
	ti = -bi[k];
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
	    br[i__] += prodr;
	    bi[i__] += prodi;
/* L30: */
	}
/* L40: */
    }
L50:
    den = ar[ar_dim1 + 1] * ar[ar_dim1 + 1] + ai[ai_dim1 + 1] * ai[ai_dim1 + 
	    1];
    prodr = br[1] * ar[ar_dim1 + 1] + bi[1] * ai[ai_dim1 + 1];
    prodi = bi[1] * ar[ar_dim1 + 1] - br[1] * ai[ai_dim1 + 1];
    br[1] = prodr / den;
    bi[1] = prodi / den;
    return 0;
/* ----------------------- END OF SUBROUTINE SOLHC ----------------------- */
} /* solhc_ */


/* Subroutine */ int decb_(integer *n, integer *ndim, doublereal *a, integer *
	ml, integer *mu, integer *ip, integer *ier)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, m;
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
    ip[*n] = 1;
    md = *ml + *mu + 1;
    md1 = md + 1;
    ju = 0;
    if (*ml == 0) {
	goto L70;
    }
    if (*n == 1) {
	goto L70;
    }
    if (*n < *mu + 2) {
	goto L7;
    }
    i__1 = *n;
    for (j = *mu + 2; j <= i__1; ++j) {
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L5: */
	    a[i__ + j * a_dim1] = 0.;
	}
    }
L7:
    nm1 = *n - 1;
    i__2 = nm1;
    for (k = 1; k <= i__2; ++k) {
	kp1 = k + 1;
	m = md;
/* Computing MIN */
	i__1 = *ml, i__3 = *n - k;
	mdl = min(i__1,i__3) + md;
	i__1 = mdl;
	for (i__ = md1; i__ <= i__1; ++i__) {
	    if ((d__1 = a[i__ + k * a_dim1], abs(d__1)) > (d__2 = a[m + k * 
		    a_dim1], abs(d__2))) {
		m = i__;
	    }
/* L10: */
	}
	ip[k] = m + k - md;
	t = a[m + k * a_dim1];
	if (m == md) {
	    goto L20;
	}
	ip[*n] = -ip[*n];
	a[m + k * a_dim1] = a[md + k * a_dim1];
	a[md + k * a_dim1] = t;
L20:
	if (t == 0.) {
	    goto L80;
	}
	t = 1. / t;
	i__1 = mdl;
	for (i__ = md1; i__ <= i__1; ++i__) {
/* L30: */
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
	}
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ip[k];
	i__1 = max(i__3,i__4);
	ju = min(i__1,*n);
	mm = md;
	if (ju < kp1) {
	    goto L55;
	}
	i__1 = ju;
	for (j = kp1; j <= i__1; ++j) {
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
	    i__3 = mdl;
	    for (i__ = md1; i__ <= i__3; ++i__) {
		ijk = i__ - jk;
/* L40: */
		a[ijk + j * a_dim1] += a[i__ + k * a_dim1] * t;
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
    k = *n;
    if (a[md + *n * a_dim1] == 0.) {
	goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[*n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECB ------------------------ */
} /* decb_ */



/* Subroutine */ int solb_(integer *n, integer *ndim, doublereal *a, integer *
	ml, integer *mu, doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, m;
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
    nm1 = *n - 1;
    if (*ml == 0) {
	goto L25;
    }
    if (*n == 1) {
	goto L50;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	m = ip[k];
	t = b[m];
	b[m] = b[k];
	b[k] = t;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	mdl = min(i__2,i__3) + md;
	i__2 = mdl;
	for (i__ = md1; i__ <= i__2; ++i__) {
	    imd = i__ + k - md;
/* L10: */
	    b[imd] += a[i__ + k * a_dim1] * t;
	}
/* L20: */
    }
L25:
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[md + k * a_dim1];
	t = -b[k];
	kmd = md - k;
/* Computing MAX */
	i__2 = 1, i__3 = kmd + 1;
	lm = max(i__2,i__3);
	i__2 = mdm;
	for (i__ = lm; i__ <= i__2; ++i__) {
	    imd = i__ - kmd;
/* L30: */
	    b[imd] += a[i__ + k * a_dim1] * t;
	}
/* L40: */
    }
L50:
    b[1] /= a[md + a_dim1];
    return 0;
/* ----------------------- END OF SUBROUTINE SOLB ------------------------ */
} /* solb_ */


/* Subroutine */ int decbc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *ml, integer *mu, integer *ip, integer *ier)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k, m, md, jk, mm;
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
    ip[*n] = 1;
    md = *ml + *mu + 1;
    md1 = md + 1;
    ju = 0;
    if (*ml == 0) {
	goto L70;
    }
    if (*n == 1) {
	goto L70;
    }
    if (*n < *mu + 2) {
	goto L7;
    }
    i__1 = *n;
    for (j = *mu + 2; j <= i__1; ++j) {
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ar[i__ + j * ar_dim1] = 0.;
	    ai[i__ + j * ai_dim1] = 0.;
/* L5: */
	}
    }
L7:
    nm1 = *n - 1;
    i__2 = nm1;
    for (k = 1; k <= i__2; ++k) {
	kp1 = k + 1;
	m = md;
/* Computing MIN */
	i__1 = *ml, i__3 = *n - k;
	mdl = min(i__1,i__3) + md;
	i__1 = mdl;
	for (i__ = md1; i__ <= i__1; ++i__) {
	    if ((d__1 = ar[i__ + k * ar_dim1], abs(d__1)) + (d__2 = ai[i__ + 
		    k * ai_dim1], abs(d__2)) > (d__3 = ar[m + k * ar_dim1], 
		    abs(d__3)) + (d__4 = ai[m + k * ai_dim1], abs(d__4))) {
		m = i__;
	    }
/* L10: */
	}
	ip[k] = m + k - md;
	tr = ar[m + k * ar_dim1];
	ti = ai[m + k * ai_dim1];
	if (m == md) {
	    goto L20;
	}
	ip[*n] = -ip[*n];
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
	i__1 = mdl;
	for (i__ = md1; i__ <= i__1; ++i__) {
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
	    ar[i__ + k * ar_dim1] = -prodr;
	    ai[i__ + k * ai_dim1] = -prodi;
/* L30: */
	}
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ip[k];
	i__1 = max(i__3,i__4);
	ju = min(i__1,*n);
	mm = md;
	if (ju < kp1) {
	    goto L55;
	}
	i__1 = ju;
	for (j = kp1; j <= i__1; ++j) {
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
		i__3 = mdl;
		for (i__ = md1; i__ <= i__3; ++i__) {
		    ijk = i__ - jk;
		    prodr = ar[i__ + k * ar_dim1] * tr;
		    prodi = ai[i__ + k * ai_dim1] * tr;
		    ar[ijk + j * ar_dim1] += prodr;
		    ai[ijk + j * ai_dim1] += prodi;
/* L40: */
		}
		goto L48;
	    }
	    if (tr == 0.) {
		i__3 = mdl;
		for (i__ = md1; i__ <= i__3; ++i__) {
		    ijk = i__ - jk;
		    prodr = -ai[i__ + k * ai_dim1] * ti;
		    prodi = ar[i__ + k * ar_dim1] * ti;
		    ar[ijk + j * ar_dim1] += prodr;
		    ai[ijk + j * ai_dim1] += prodi;
/* L45: */
		}
		goto L48;
	    }
	    i__3 = mdl;
	    for (i__ = md1; i__ <= i__3; ++i__) {
		ijk = i__ - jk;
		prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * 
			ti;
		prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * 
			ti;
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
    k = *n;
    if ((d__1 = ar[md + *n * ar_dim1], abs(d__1)) + (d__2 = ai[md + *n * 
	    ai_dim1], abs(d__2)) == 0.) {
	goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[*n] = 0;
    return 0;
/* ----------------------- END OF SUBROUTINE DECBC ------------------------ */
} /* decbc_ */



/* Subroutine */ int solbc_(integer *n, integer *ndim, doublereal *ar, 
	doublereal *ai, integer *ml, integer *mu, doublereal *br, doublereal *
	bi, integer *ip)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, m, kb, md, lm;
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
    nm1 = *n - 1;
    if (*ml == 0) {
	goto L25;
    }
    if (*n == 1) {
	goto L50;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	m = ip[k];
	tr = br[m];
	ti = bi[m];
	br[m] = br[k];
	bi[m] = bi[k];
	br[k] = tr;
	bi[k] = ti;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	mdl = min(i__2,i__3) + md;
	i__2 = mdl;
	for (i__ = md1; i__ <= i__2; ++i__) {
	    imd = i__ + k - md;
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
	    br[imd] += prodr;
	    bi[imd] += prodi;
/* L10: */
	}
/* L20: */
    }
L25:
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	den = ar[md + k * ar_dim1] * ar[md + k * ar_dim1] + ai[md + k * 
		ai_dim1] * ai[md + k * ai_dim1];
	prodr = br[k] * ar[md + k * ar_dim1] + bi[k] * ai[md + k * ai_dim1];
	prodi = bi[k] * ar[md + k * ar_dim1] - br[k] * ai[md + k * ai_dim1];
	br[k] = prodr / den;
	bi[k] = prodi / den;
	tr = -br[k];
	ti = -bi[k];
	kmd = md - k;
/* Computing MAX */
	i__2 = 1, i__3 = kmd + 1;
	lm = max(i__2,i__3);
	i__2 = mdm;
	for (i__ = lm; i__ <= i__2; ++i__) {
	    imd = i__ - kmd;
	    prodr = ar[i__ + k * ar_dim1] * tr - ai[i__ + k * ai_dim1] * ti;
	    prodi = ai[i__ + k * ai_dim1] * tr + ar[i__ + k * ar_dim1] * ti;
	    br[imd] += prodr;
	    bi[imd] += prodi;
/* L30: */
	}
/* L40: */
    }
    den = ar[md + ar_dim1] * ar[md + ar_dim1] + ai[md + ai_dim1] * ai[md + 
	    ai_dim1];
    prodr = br[1] * ar[md + ar_dim1] + bi[1] * ai[md + ai_dim1];
    prodi = bi[1] * ar[md + ar_dim1] - br[1] * ai[md + ai_dim1];
    br[1] = prodr / den;
    bi[1] = prodi / den;
L50:
    return 0;
/* ----------------------- END OF SUBROUTINE SOLBC ------------------------ */
} /* solbc_ */



/* Subroutine */ int elmhes_(integer *nm, integer *n, integer *low, integer *
	igh, doublereal *a, integer *int__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, m;
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

    i__1 = la;
    for (m = kp1; m <= i__1; ++m) {
	mm1 = m - 1;
	x = 0.;
	i__ = m;

	i__2 = *igh;
	for (j = m; j <= i__2; ++j) {
	    if ((d__1 = a[j + mm1 * a_dim1], abs(d__1)) <= abs(x)) {
		goto L100;
	    }
	    x = a[j + mm1 * a_dim1];
	    i__ = j;
L100:
	    ;
	}

	int__[m] = i__;
	if (i__ == m) {
	    goto L130;
	}
/*    :::::::::: interchange rows and columns of a :::::::::: */
	i__2 = *n;
	for (j = mm1; j <= i__2; ++j) {
	    y = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = a[m + j * a_dim1];
	    a[m + j * a_dim1] = y;
/* L110: */
	}

	i__2 = *igh;
	for (j = 1; j <= i__2; ++j) {
	    y = a[j + i__ * a_dim1];
	    a[j + i__ * a_dim1] = a[j + m * a_dim1];
	    a[j + m * a_dim1] = y;
/* L120: */
	}
/*    :::::::::: end interchange :::::::::: */
L130:
	if (x == 0.) {
	    goto L180;
	}
	mp1 = m + 1;

	i__2 = *igh;
	for (i__ = mp1; i__ <= i__2; ++i__) {
	    y = a[i__ + mm1 * a_dim1];
	    if (y == 0.) {
		goto L160;
	    }
	    y /= x;
	    a[i__ + mm1 * a_dim1] = y;

	    i__3 = *n;
	    for (j = m; j <= i__3; ++j) {
/* L140: */
		a[i__ + j * a_dim1] -= y * a[m + j * a_dim1];
	    }

	    i__3 = *igh;
	    for (j = 1; j <= i__3; ++j) {
/* L150: */
		a[j + m * a_dim1] += y * a[j + i__ * a_dim1];
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

/* Subroutine */ int decomr_(integer *n, doublereal *fjac, integer *ldjac, 
	doublereal *fmas, integer *ldmas, integer *mlmas, integer *mumas, 
	integer *m1, integer *m2, integer *nm1, doublereal *fac1, doublereal *
	e1, integer *lde1, integer *ip1, integer *ier, integer *ijob, logical 
	*calhes, integer *iphes)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, j1, ib, mm, jm1;
    extern /* Subroutine */ int dec_(integer *, integer *, doublereal *, 
	    integer *, integer *);
    static doublereal sum;
    extern /* Subroutine */ int decb_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *), dech_(integer *, 
	    integer *, doublereal *, integer *, integer *, integer *), 
	    elmhes_(integer *, integer *, integer *, integer *, doublereal *, 
	    integer *);


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
	case 8:  goto L55;
	case 9:  goto L55;
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
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e1[i__ + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
	}
	e1[j + j * e1_dim1] += *fac1;
    }
    dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e1[i__ + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1];
	}
	e1[j + j * e1_dim1] += *fac1;
    }
L45:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = 0.;
	    i__3 = mm - 1;
	    for (k = 0; k <= i__3; ++k) {
		sum = (sum + fjac[i__ + (j + k * *m2) * fjac_dim1]) / *fac1;
	    }
	    e1[i__ + j * e1_dim1] -= sum;
	}
    }
    dec_(nm1, lde1, &e1[e1_offset], &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
	}
	e1[linal_1.mdiag + j * e1_dim1] += *fac1;
    }
    decb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1]
		    ;
	}
	e1[linal_1.mdiag + j * e1_dim1] += *fac1;
    }
L46:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = 0.;
	    i__3 = mm - 1;
	    for (k = 0; k <= i__3; ++k) {
		sum = (sum + fjac[i__ + (j + k * *m2) * fjac_dim1]) / *fac1;
	    }
	    e1[i__ + linal_1.mle + j * e1_dim1] -= sum;
	}
    }
    decb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier)
	    ;
    return 0;

/* ----------------------------------------------------------- */

L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e1[i__ + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
	}
/* Computing MAX */
	i__2 = 1, i__3 = j - *mumas;
/* Computing MIN */
	i__5 = *n, i__6 = j + *mlmas;
	i__4 = min(i__5,i__6);
	for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
	    e1[i__ + j * e1_dim1] += *fac1 * fmas[i__ - j + linal_1.mbdiag + 
		    j * fmas_dim1];
	}
    }
    dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__4 = *nm1;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    e1[i__ + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1];
	}
/* Computing MAX */
	i__4 = 1, i__2 = j - *mumas;
/* Computing MIN */
	i__5 = *nm1, i__6 = j + *mlmas;
	i__3 = min(i__5,i__6);
	for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
	    e1[i__ + j * e1_dim1] += *fac1 * fmas[i__ - j + linal_1.mbdiag + 
		    j * fmas_dim1];
	}
    }
    goto L45;

/* ----------------------------------------------------------- */

L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__3 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
	}
	i__3 = linal_1.mbb;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ib = i__ + linal_1.mdiff;
	    e1[ib + j * e1_dim1] += *fac1 * fmas[i__ + j * fmas_dim1];
	}
    }
    decb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L14:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__3 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    e1[i__ + linal_1.mle + j * e1_dim1] = -fjac[i__ + jm1 * fjac_dim1]
		    ;
	}
	i__3 = linal_1.mbb;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    ib = i__ + linal_1.mdiff;
	    e1[ib + j * e1_dim1] += *fac1 * fmas[i__ + j * fmas_dim1];
	}
    }
    goto L46;

/* ----------------------------------------------------------- */

L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__3 = *n;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    e1[i__ + j * e1_dim1] = fmas[i__ + j * fmas_dim1] * *fac1 - fjac[
		    i__ + j * fjac_dim1];
	}
    }
    dec_(n, lde1, &e1[e1_offset], &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__3 = *nm1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    e1[i__ + j * e1_dim1] = fmas[i__ + j * fmas_dim1] * *fac1 - fjac[
		    i__ + jm1 * fjac_dim1];
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
	elmhes_(ldjac, n, &c__1, n, &fjac[fjac_offset], &iphes[1]);
    }
    *calhes = FALSE_;
    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
	j1 = j + 1;
	e1[j1 + j * e1_dim1] = -fjac[j1 + j * fjac_dim1];
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__3 = j;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    e1[i__ + j * e1_dim1] = -fjac[i__ + j * fjac_dim1];
	}
	e1[j + j * e1_dim1] += *fac1;
    }
    dech_(n, lde1, &e1[e1_offset], &c__1, &ip1[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* decomr_ */


/*     END OF SUBROUTINE DECOMR */

/* *********************************************************** */

/* Subroutine */ int decomc_(integer *n, doublereal *fjac, integer *ldjac, 
	doublereal *fmas, integer *ldmas, integer *mlmas, integer *mumas, 
	integer *m1, integer *m2, integer *nm1, doublereal *alphn, doublereal 
	*betan, doublereal *e2r, doublereal *e2i, integer *lde1, integer *ip2,
	 integer *ier, integer *ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e2r_dim1, 
	    e2r_offset, e2i_dim1, e2i_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, j1;
    static doublereal bb;
    static integer ib, mm, jm1;
    static doublereal bet, alp;
    extern /* Subroutine */ int decc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *);
    static doublereal ffma, abno;
    static integer imle;
    static doublereal sumi, sumr, sums;
    extern /* Subroutine */ int decbc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *), dechc_(
	    integer *, integer *, doublereal *, doublereal *, integer *, 
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
	case 8:  goto L55;
	case 9:  goto L55;
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
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
	    e2i[i__ + j * e2i_dim1] = 0.;
	}
	e2r[j + j * e2r_dim1] += *alphn;
	e2i[j + j * e2i_dim1] = *betan;
    }
    decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + jm1 * fjac_dim1];
	    e2i[i__ + j * e2i_dim1] = 0.;
	}
	e2r[j + j * e2r_dim1] += *alphn;
	e2i[j + j * e2i_dim1] = *betan;
    }
L45:
    mm = *m1 / *m2;
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
    alp = *alphn / abno;
    bet = *betan / abno;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sumr = 0.;
	    sumi = 0.;
	    i__3 = mm - 1;
	    for (k = 0; k <= i__3; ++k) {
		sums = sumr + fjac[i__ + (j + k * *m2) * fjac_dim1];
		sumr = sums * alp + sumi * bet;
		sumi = sumi * alp - sums * bet;
	    }
	    e2r[i__ + j * e2r_dim1] -= sumr;
	    e2i[i__ + j * e2i_dim1] -= sumi;
	}
    }
    decc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    imle = i__ + linal_1.mle;
	    e2r[imle + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
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
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e2r[i__ + linal_1.mle + j * e2r_dim1] = -fjac[i__ + jm1 * 
		    fjac_dim1];
	    e2i[i__ + linal_1.mle + j * e2i_dim1] = 0.;
	}
	e2r[linal_1.mdiag + j * e2r_dim1] += *alphn;
	e2i[linal_1.mdiag + j * e2i_dim1] += *betan;
    }
L46:
    mm = *m1 / *m2;
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
    alp = *alphn / abno;
    bet = *betan / abno;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sumr = 0.;
	    sumi = 0.;
	    i__3 = mm - 1;
	    for (k = 0; k <= i__3; ++k) {
		sums = sumr + fjac[i__ + (j + k * *m2) * fjac_dim1];
		sumr = sums * alp + sumi * bet;
		sumi = sumi * alp - sums * bet;
	    }
	    imle = i__ + linal_1.mle;
	    e2r[imle + j * e2r_dim1] -= sumr;
	    e2i[imle + j * e2i_dim1] -= sumi;
	}
    }
    decbc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
	    e2i[i__ + j * e2i_dim1] = 0.;
	}
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = 1, i__3 = j - *mumas;
/* Computing MIN */
	i__5 = *n, i__6 = j + *mlmas;
	i__4 = min(i__5,i__6);
	for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
	    e2r[i__ + j * e2r_dim1] += *alphn * bb;
	    e2i[i__ + j * e2i_dim1] = *betan * bb;
	}
    }
    decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__4 = *nm1;
	for (i__ = 1; i__ <= i__4; ++i__) {
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + jm1 * fjac_dim1];
	    e2i[i__ + j * e2i_dim1] = 0.;
	}
/* Computing MAX */
	i__4 = 1, i__2 = j - *mumas;
/* Computing MIN */
	i__5 = *nm1, i__6 = j + *mlmas;
	i__3 = min(i__5,i__6);
	for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
	    ffma = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
	    e2r[i__ + j * e2r_dim1] += *alphn * ffma;
	    e2i[i__ + j * e2i_dim1] += *betan * ffma;
	}
    }
    goto L45;

/* ----------------------------------------------------------- */

L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__3 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    imle = i__ + linal_1.mle;
	    e2r[imle + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
	    e2i[imle + j * e2i_dim1] = 0.;
	}
/* Computing MAX */
	i__3 = 1, i__4 = *mumas + 2 - j;
/* Computing MIN */
	i__5 = linal_1.mbb, i__6 = *mumas + 1 - j + *n;
	i__2 = min(i__5,i__6);
	for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
	    ib = i__ + linal_1.mdiff;
	    bb = fmas[i__ + j * fmas_dim1];
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
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__2 = linal_1.mbjac;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e2r[i__ + linal_1.mle + j * e2r_dim1] = -fjac[i__ + jm1 * 
		    fjac_dim1];
	    e2i[i__ + linal_1.mle + j * e2i_dim1] = 0.;
	}
	i__2 = linal_1.mbb;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ib = i__ + linal_1.mdiff;
	    ffma = fmas[i__ + j * fmas_dim1];
	    e2r[ib + j * e2r_dim1] += *alphn * ffma;
	    e2i[ib + j * e2i_dim1] += *betan * ffma;
	}
    }
    goto L46;

/* ----------------------------------------------------------- */

L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    bb = fmas[i__ + j * fmas_dim1];
	    e2r[i__ + j * e2r_dim1] = bb * *alphn - fjac[i__ + j * fjac_dim1];
	    e2i[i__ + j * e2i_dim1] = bb * *betan;
	}
    }
    decc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *nm1;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j + *m1;
	i__2 = *nm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e2r[i__ + j * e2r_dim1] = *alphn * fmas[i__ + j * fmas_dim1] - 
		    fjac[i__ + jm1 * fjac_dim1];
	    e2i[i__ + j * e2i_dim1] = *betan * fmas[i__ + j * fmas_dim1];
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
    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
	j1 = j + 1;
	e2r[j1 + j * e2r_dim1] = -fjac[j1 + j * fjac_dim1];
	e2i[j1 + j * e2i_dim1] = 0.;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    e2i[i__ + j * e2i_dim1] = 0.;
	    e2r[i__ + j * e2r_dim1] = -fjac[i__ + j * fjac_dim1];
	}
	e2r[j + j * e2r_dim1] += *alphn;
	e2i[j + j * e2i_dim1] = *betan;
    }
    dechc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &ip2[1], ier);
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* decomc_ */


/*     END OF SUBROUTINE DECOMC */

/* *********************************************************** */

/* Subroutine */ int slvrar_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *fac1, doublereal *e1, integer *lde1, doublereal *z1, 
	doublereal *f1, integer *ip1, integer *iphes, integer *ier, integer *
	ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s1;
    static integer mm, mp, im1, mp1, jkm;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum1;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solh_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal zsafe;


    /* Parameter adjustments */
    --iphes;
    --f1;
    --z1;
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
	case 8:  goto L55;
	case 9:  goto L55;
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z1[i__] -= f1[i__] * *fac1;
    }
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z1[i__] -= f1[i__] * *fac1;
    }
L48:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum1 = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    jkm = j + k * *m2;
	    sum1 = (z1[jkm] + sum1) / *fac1;
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		im1 = i__ + *m1;
		z1[im1] += fjac[i__ + jkm * fjac_dim1] * sum1;
	    }
	}
    }
    sol_(nm1, lde1, &e1[e1_offset], &z1[*m1 + 1], &ip1[1]);
L49:
    for (i__ = *m1; i__ >= 1; --i__) {
	z1[i__] = (z1[i__] + z1[*m2 + i__]) / *fac1;
    }
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z1[i__] -= f1[i__] * *fac1;
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	    );
    return 0;

/* ----------------------------------------------------------- */

L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z1[i__] -= f1[i__] * *fac1;
    }
L45:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum1 = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    jkm = j + k * *m2;
	    sum1 = (z1[jkm] + sum1) / *fac1;
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		im1 = i__ + *m1;
		z1[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * 
			sum1;
	    }
	}
    }
    solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[*m1 + 1],
	     &ip1[1]);
    goto L49;

/* ----------------------------------------------------------- */

L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1 = 0.;
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
	    s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
	}
	z1[i__] += s1 * *fac1;
    }
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
    return 0;

/* ----------------------------------------------------------- */

L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z1[i__] -= f1[i__] * *fac1;
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	im1 = i__ + *m1;
	s1 = 0.;
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
	    s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *m1]
		    ;
	}
	z1[im1] += s1 * *fac1;
    }
    if (*ijob == 14) {
	goto L45;
    }
    goto L48;

/* ----------------------------------------------------------- */

L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1 = 0.;
/* Computing MAX */
	i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__4 = min(i__5,i__6);
	for (j = max(i__2,i__3); j <= i__4; ++j) {
	    s1 -= fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
	}
	z1[i__] += s1 * *fac1;
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	    );
    return 0;

/* ----------------------------------------------------------- */

L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1 = 0.;
	i__4 = *n;
	for (j = 1; j <= i__4; ++j) {
	    s1 -= fmas[i__ + j * fmas_dim1] * f1[j];
	}
	z1[i__] += s1 * *fac1;
    }
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
    return 0;

/* ----------------------------------------------------------- */

L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z1[i__] -= f1[i__] * *fac1;
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	im1 = i__ + *m1;
	s1 = 0.;
	i__4 = *nm1;
	for (j = 1; j <= i__4; ++j) {
	    s1 -= fmas[i__ + j * fmas_dim1] * f1[j + *m1];
	}
	z1[im1] += s1 * *fac1;
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z1[i__] -= f1[i__] * *fac1;
    }
    for (mm = *n - 2; mm >= 1; --mm) {
	mp = *n - mm;
	mp1 = mp - 1;
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L746;
	}
	zsafe = z1[mp];
	z1[mp] = z1[i__];
	z1[i__] = zsafe;
L746:
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
	    z1[i__] -= fjac[i__ + mp1 * fjac_dim1] * z1[mp];
	}
    }
    solh_(n, lde1, &e1[e1_offset], &c__1, &z1[1], &ip1[1]);
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
	mp = *n - mm;
	mp1 = mp - 1;
	i__4 = *n;
	for (i__ = mp + 1; i__ <= i__4; ++i__) {
	    z1[i__] += fjac[i__ + mp1 * fjac_dim1] * z1[mp];
	}
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L750;
	}
	zsafe = z1[mp];
	z1[mp] = z1[i__];
	z1[i__] = zsafe;
L750:
	;
    }
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* slvrar_ */


/*     END OF SUBROUTINE SLVRAR */

/* *********************************************************** */

/* Subroutine */ int slvrai_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *alphn, doublereal *betan, doublereal *e2r, 
	doublereal *e2i, integer *lde1, doublereal *z2, doublereal *z3, 
	doublereal *f2, doublereal *f3, doublereal *cont, integer *ip2, 
	integer *iphes, integer *ier, integer *ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e2r_dim1, 
	    e2r_offset, e2i_dim1, e2i_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s2, s3, bb;
    static integer mm, mp, im1, jm1, mp1;
    static doublereal z2i, z3i;
    static integer jkm, mpi;
    static doublereal sum2, sum3, abno;
    extern /* Subroutine */ int solc_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer iimu;
    static doublereal sumh, e1imp;
    extern /* Subroutine */ int solbc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal zsafe;
    extern /* Subroutine */ int solhc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);


    /* Parameter adjustments */
    --iphes;
    --f3;
    --f2;
    --z3;
    --z2;
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
	case 8:  goto L55;
	case 9:  goto L55;
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
L48:
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum2 = 0.;
	sum3 = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    jkm = j + k * *m2;
	    sumh = (z2[jkm] + sum2) / abno;
	    sum3 = (z3[jkm] + sum3) / abno;
	    sum2 = sumh * *alphn + sum3 * *betan;
	    sum3 = sum3 * *alphn - sumh * *betan;
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		im1 = i__ + *m1;
		z2[im1] += fjac[i__ + jkm * fjac_dim1] * sum2;
		z3[im1] += fjac[i__ + jkm * fjac_dim1] * sum3;
	    }
	}
    }
    solc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[*m1 + 1], &z3[*
	    m1 + 1], &ip2[1]);
L49:
    for (i__ = *m1; i__ >= 1; --i__) {
	mpi = *m2 + i__;
	z2i = z2[i__] + z2[mpi];
	z3i = z3[i__] + z3[mpi];
	z3[i__] = (z3i * *alphn - z2i * *betan) / abno;
	z2[i__] = (z2i * *alphn + z3i * *betan) / abno;
    }
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[1], &z3[1], &ip2[1]);
    return 0;

/* ----------------------------------------------------------- */

L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
L45:
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum2 = 0.;
	sum3 = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    jkm = j + k * *m2;
	    sumh = (z2[jkm] + sum2) / abno;
	    sum3 = (z3[jkm] + sum3) / abno;
	    sum2 = sumh * *alphn + sum3 * *betan;
	    sum3 = sum3 * *alphn - sumh * *betan;
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		im1 = i__ + *m1;
		iimu = i__ + *mujac + 1 - j;
		z2[im1] += fjac[iimu + jkm * fjac_dim1] * sum2;
		z3[im1] += fjac[iimu + jkm * fjac_dim1] * sum3;
	    }
	}
    }
    solbc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
    goto L49;

/* ----------------------------------------------------------- */

L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = 0.;
	s3 = 0.;
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
	    s2 -= bb * f2[j];
	    s3 -= bb * f3[j];
	}
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
    return 0;

/* ----------------------------------------------------------- */

L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	im1 = i__ + *m1;
	s2 = 0.;
	s3 = 0.;
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
	    jm1 = j + *m1;
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
	    s2 -= bb * f2[jm1];
	    s3 -= bb * f3[jm1];
	}
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = 0.;
	s3 = 0.;
/* Computing MAX */
	i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__4 = min(i__5,i__6);
	for (j = max(i__2,i__3); j <= i__4; ++j) {
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
	    s2 -= bb * f2[j];
	    s3 -= bb * f3[j];
	}
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[1], &z3[1], &ip2[1]);
    return 0;

/* ----------------------------------------------------------- */

L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = 0.;
	s3 = 0.;
	i__4 = *n;
	for (j = 1; j <= i__4; ++j) {
	    bb = fmas[i__ + j * fmas_dim1];
	    s2 -= bb * f2[j];
	    s3 -= bb * f3[j];
	}
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
    return 0;

/* ----------------------------------------------------------- */

L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	im1 = i__ + *m1;
	s2 = 0.;
	s3 = 0.;
	i__4 = *nm1;
	for (j = 1; j <= i__4; ++j) {
	    jm1 = j + *m1;
	    bb = fmas[i__ + j * fmas_dim1];
	    s2 -= bb * f2[jm1];
	    s3 -= bb * f3[jm1];
	}
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    for (mm = *n - 2; mm >= 1; --mm) {
	mp = *n - mm;
	mp1 = mp - 1;
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L746;
	}
	zsafe = z2[mp];
	z2[mp] = z2[i__];
	z2[i__] = zsafe;
	zsafe = z3[mp];
	z3[mp] = z3[i__];
	z3[i__] = zsafe;
L746:
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
	    e1imp = fjac[i__ + mp1 * fjac_dim1];
	    z2[i__] -= e1imp * z2[mp];
	    z3[i__] -= e1imp * z3[mp];
	}
    }
    solhc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &z2[1], &z3[1],
	     &ip2[1]);
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
	mp = *n - mm;
	mp1 = mp - 1;
	i__4 = *n;
	for (i__ = mp + 1; i__ <= i__4; ++i__) {
	    e1imp = fjac[i__ + mp1 * fjac_dim1];
	    z2[i__] += e1imp * z2[mp];
	    z3[i__] += e1imp * z3[mp];
	}
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L750;
	}
	zsafe = z2[mp];
	z2[mp] = z2[i__];
	z2[i__] = zsafe;
	zsafe = z3[mp];
	z3[mp] = z3[i__];
	z3[i__] = zsafe;
L750:
	;
    }
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* slvrai_ */


/*     END OF SUBROUTINE SLVRAI */

/* *********************************************************** */

/* Subroutine */ int slvrad_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *fac1, doublereal *alphn, doublereal *betan, 
	doublereal *e1, doublereal *e2r, doublereal *e2i, integer *lde1, 
	doublereal *z1, doublereal *z2, doublereal *z3, doublereal *f1, 
	doublereal *f2, doublereal *f3, doublereal *cont, integer *ip1, 
	integer *ip2, integer *iphes, integer *ier, integer *ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, e2r_dim1, e2r_offset, e2i_dim1, e2i_offset, i__1, i__2,
	     i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s1, s2, s3, bb;
    static integer mm, mp, j1b, j2b, im1, jm1, mp1;
    static doublereal z2i, z3i;
    static integer jkm, mpi;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum1, sum2, sum3, ffja, abno;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solc_(integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *), solh_(integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal sumh, e1imp;
    extern /* Subroutine */ int solbc_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal zsafe;
    extern /* Subroutine */ int solhc_(integer *, integer *, doublereal *, 
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
	case 8:  goto L55;
	case 9:  goto L55;
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z1[i__] -= f1[i__] * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z1[i__] -= f1[i__] * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
L48:
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
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
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		im1 = i__ + *m1;
		z1[im1] += fjac[i__ + jkm * fjac_dim1] * sum1;
		z2[im1] += fjac[i__ + jkm * fjac_dim1] * sum2;
		z3[im1] += fjac[i__ + jkm * fjac_dim1] * sum3;
	    }
	}
    }
    sol_(nm1, lde1, &e1[e1_offset], &z1[*m1 + 1], &ip1[1]);
    solc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[*m1 + 1], &z3[*
	    m1 + 1], &ip2[1]);
L49:
    for (i__ = *m1; i__ >= 1; --i__) {
	mpi = *m2 + i__;
	z1[i__] = (z1[i__] + z1[mpi]) / *fac1;
	z2i = z2[i__] + z2[mpi];
	z3i = z3[i__] + z3[mpi];
	z3[i__] = (z3i * *alphn - z2i * *betan) / abno;
	z2[i__] = (z2i * *alphn + z3i * *betan) / abno;
    }
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z1[i__] -= f1[i__] * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	    );
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[1], &z3[1], &ip2[1]);
    return 0;

/* ----------------------------------------------------------- */

L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z1[i__] -= f1[i__] * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
L45:
/* Computing 2nd power */
    d__1 = *alphn;
/* Computing 2nd power */
    d__2 = *betan;
    abno = d__1 * d__1 + d__2 * d__2;
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
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
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		im1 = i__ + *m1;
		ffja = fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1];
		z1[im1] += ffja * sum1;
		z2[im1] += ffja * sum2;
		z3[im1] += ffja * sum3;
	    }
	}
    }
    solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[*m1 + 1],
	     &ip1[1]);
    solbc_(nm1, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[*m1 + 1], &z3[*m1 + 1], &ip2[1]);
    goto L49;

/* ----------------------------------------------------------- */

L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1 = 0.;
	s2 = 0.;
	s3 = 0.;
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
	    s1 -= bb * f1[j];
	    s2 -= bb * f2[j];
	    s3 -= bb * f3[j];
	}
	z1[i__] += s1 * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
    return 0;

/* ----------------------------------------------------------- */

L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z1[i__] -= f1[i__] * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	im1 = i__ + *m1;
	s1 = 0.;
	s2 = 0.;
	s3 = 0.;
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
	j1b = max(i__3,i__4);
/* Computing MIN */
	i__3 = *nm1, i__4 = i__ + *mumas;
	j2b = min(i__3,i__4);
	i__3 = j2b;
	for (j = j1b; j <= i__3; ++j) {
	    jm1 = j + *m1;
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1 = 0.;
	s2 = 0.;
	s3 = 0.;
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
	    bb = fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1];
	    s1 -= bb * f1[j];
	    s2 -= bb * f2[j];
	    s3 -= bb * f3[j];
	}
	z1[i__] += s1 * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &z1[1], &ip1[1]
	    );
    solbc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &linal_1.mle, &
	    linal_1.mue, &z2[1], &z3[1], &ip2[1]);
    return 0;

/* ----------------------------------------------------------- */

L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s1 = 0.;
	s2 = 0.;
	s3 = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    bb = fmas[i__ + j * fmas_dim1];
	    s1 -= bb * f1[j];
	    s2 -= bb * f2[j];
	    s3 -= bb * f3[j];
	}
	z1[i__] += s1 * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    sol_(n, lde1, &e1[e1_offset], &z1[1], &ip1[1]);
    solc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &z2[1], &z3[1], &ip2[1]
	    );
    return 0;

/* ----------------------------------------------------------- */

L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z1[i__] -= f1[i__] * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	im1 = i__ + *m1;
	s1 = 0.;
	s2 = 0.;
	s3 = 0.;
	i__2 = *nm1;
	for (j = 1; j <= i__2; ++j) {
	    jm1 = j + *m1;
	    bb = fmas[i__ + j * fmas_dim1];
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s2 = -f2[i__];
	s3 = -f3[i__];
	z1[i__] -= f1[i__] * *fac1;
	z2[i__] = z2[i__] + s2 * *alphn - s3 * *betan;
	z3[i__] = z3[i__] + s3 * *alphn + s2 * *betan;
    }
    for (mm = *n - 2; mm >= 1; --mm) {
	mp = *n - mm;
	mp1 = mp - 1;
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L746;
	}
	zsafe = z1[mp];
	z1[mp] = z1[i__];
	z1[i__] = zsafe;
	zsafe = z2[mp];
	z2[mp] = z2[i__];
	z2[i__] = zsafe;
	zsafe = z3[mp];
	z3[mp] = z3[i__];
	z3[i__] = zsafe;
L746:
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
	    e1imp = fjac[i__ + mp1 * fjac_dim1];
	    z1[i__] -= e1imp * z1[mp];
	    z2[i__] -= e1imp * z2[mp];
	    z3[i__] -= e1imp * z3[mp];
	}
    }
    solh_(n, lde1, &e1[e1_offset], &c__1, &z1[1], &ip1[1]);
    solhc_(n, lde1, &e2r[e2r_offset], &e2i[e2i_offset], &c__1, &z2[1], &z3[1],
	     &ip2[1]);
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
	mp = *n - mm;
	mp1 = mp - 1;
	i__2 = *n;
	for (i__ = mp + 1; i__ <= i__2; ++i__) {
	    e1imp = fjac[i__ + mp1 * fjac_dim1];
	    z1[i__] += e1imp * z1[mp];
	    z2[i__] += e1imp * z2[mp];
	    z3[i__] += e1imp * z3[mp];
	}
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L750;
	}
	zsafe = z1[mp];
	z1[mp] = z1[i__];
	z1[i__] = zsafe;
	zsafe = z2[mp];
	z2[mp] = z2[i__];
	z2[i__] = zsafe;
	zsafe = z3[mp];
	z3[mp] = z3[i__];
	z3[i__] = zsafe;
L750:
	;
    }
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* slvrad_ */


/*     END OF SUBROUTINE SLVRAD */

/* *********************************************************** */

/* Subroutine */ int estrad_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, doublereal *h__, doublereal *dd1, 
	doublereal *dd2, doublereal *dd3, FP_CB_f fcn, void* fcn_PY, integer *nfcn, doublereal 
	*y0, doublereal *y, integer *ijob, doublereal *x, integer *m1, 
	integer *m2, integer *nm1, doublereal *e1, integer *lde1, doublereal *
	z1, doublereal *z2, doublereal *z3, doublereal *cont, doublereal *
	werr, doublereal *f1, doublereal *f2, integer *ip1, integer *iphes, 
	doublereal *scal, doublereal *err, logical *first, logical *reject, 
	doublereal *fac1, doublereal *rpar, integer *ipar)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, mm, mp, im1;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum, hee1, hee2, hee3, sum1;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solh_(integer *, 
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
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
    }

L1:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	cont[i__] = f2[i__] + y0[i__];
    }
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L77;

L11:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	cont[i__] = f2[i__] + y0[i__];
    }
L48:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum1 = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    sum1 = (cont[j + k * *m2] + sum1) / *fac1;
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		im1 = i__ + *m1;
		cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
	    }
	}
    }
    sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
	cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
    }
    goto L77;

L2:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	cont[i__] = f2[i__] + y0[i__];
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
	    1]);
    goto L77;

L12:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	cont[i__] = f2[i__] + y0[i__];
    }
L45:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum1 = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		im1 = i__ + *m1;
		cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) * 
			fjac_dim1] * sum1;
	    }
	}
    }
    solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 
	    1], &ip1[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
	cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
    }
    goto L77;

L3:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
	}
	f2[i__] = sum;
	cont[i__] = sum + y0[i__];
    }
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L77;

L13:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	cont[i__] = f2[i__] + y0[i__];
    }
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *
		    m1];
	}
	im1 = i__ + *m1;
	f2[im1] = sum;
	cont[im1] = sum + y0[im1];
    }
    goto L48;

L4:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
/* Computing MAX */
	i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__4 = min(i__5,i__6);
	for (j = max(i__2,i__3); j <= i__4; ++j) {
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j];
	}
	f2[i__] = sum;
	cont[i__] = sum + y0[i__];
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
	    1]);
    goto L77;

L14:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	cont[i__] = f2[i__] + y0[i__];
    }
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * f1[j + *
		    m1];
	}
	im1 = i__ + *m1;
	f2[im1] = sum;
	cont[im1] = sum + y0[im1];
    }
    goto L45;

L5:
/* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *n;
	for (j = 1; j <= i__3; ++j) {
	    sum += fmas[i__ + j * fmas_dim1] * f1[j];
	}
	f2[i__] = sum;
	cont[i__] = sum + y0[i__];
    }
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L77;

L15:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	cont[i__] = f2[i__] + y0[i__];
    }
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
	f1[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *nm1;
	for (j = 1; j <= i__3; ++j) {
	    sum += fmas[i__ + j * fmas_dim1] * f1[j + *m1];
	}
	im1 = i__ + *m1;
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2[i__] = hee1 * z1[i__] + hee2 * z2[i__] + hee3 * z3[i__];
	cont[i__] = f2[i__] + y0[i__];
    }
    for (mm = *n - 2; mm >= 1; --mm) {
	mp = *n - mm;
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L310;
	}
	zsafe = cont[mp];
	cont[mp] = cont[i__];
	cont[i__] = zsafe;
L310:
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
	    cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
	}
    }
    solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
	mp = *n - mm;
	i__3 = *n;
	for (i__ = mp + 1; i__ <= i__3; ++i__) {
	    cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
	}
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L440;
	}
	zsafe = cont[mp];
	cont[mp] = cont[i__];
	cont[i__] = zsafe;
L440:
	;
    }

/* -------------------------------------- */

L77:
    *err = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	werr[i__] = cont[i__] / scal[i__];
/* Computing 2nd power */
	d__1 = werr[i__];
	*err += d__1 * d__1;
    }
/* Computing MAX */
    d__1 = sqrt(*err / *n);
    *err = max(d__1,1e-10);

    if (*err < 1.) {
	return 0;
    }
    if (*first || *reject) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    cont[i__] = y[i__] + cont[i__];
	}
    (*fcn)(n, x, &cont[1], &f1[1], &rpar[1], &ipar[1], fcn_PY);
	++(*nfcn);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    cont[i__] = f1[i__] + f2[i__];
	}
	switch (*ijob) {
	    case 1:  goto L31;
	    case 2:  goto L32;
	    case 3:  goto L31;
	    case 4:  goto L32;
	    case 5:  goto L31;
	    case 6:  goto L32;
	    case 7:  goto L33;
	    case 8:  goto L55;
	    case 9:  goto L55;
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
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
	    sum1 = 0.;
	    for (k = mm - 1; k >= 0; --k) {
		sum1 = (cont[j + k * *m2] + sum1) / *fac1;
		i__3 = *nm1;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    im1 = i__ + *m1;
		    cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
		}
	    }
	}
	sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
	    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
	}
	goto L88;
/* ------ BANDED MATRIX OPTION */
L32:
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &
		ip1[1]);
	goto L88;
/* ------ BANDED MATRIX OPTION, SECOND ORDER */
L42:
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
	    sum1 = 0.;
	    for (k = mm - 1; k >= 0; --k) {
		sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/* Computing MAX */
		i__3 = 1, i__4 = j - *mujac;
/* Computing MIN */
		i__5 = *nm1, i__6 = j + *mljac;
		i__2 = min(i__5,i__6);
		for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
		    im1 = i__ + *m1;
		    cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) * 
			    fjac_dim1] * sum1;
		}
	    }
	}
	solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*
		m1 + 1], &ip1[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
	    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
	}
	goto L88;
/* ------ HESSENBERG MATRIX OPTION */
L33:
	for (mm = *n - 2; mm >= 1; --mm) {
	    mp = *n - mm;
	    i__ = iphes[mp];
	    if (i__ == mp) {
		goto L510;
	    }
	    zsafe = cont[mp];
	    cont[mp] = cont[i__];
	    cont[i__] = zsafe;
L510:
	    i__1 = *n;
	    for (i__ = mp + 1; i__ <= i__1; ++i__) {
		cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
	    }
	}
	solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
	i__1 = *n - 2;
	for (mm = 1; mm <= i__1; ++mm) {
	    mp = *n - mm;
	    i__2 = *n;
	    for (i__ = mp + 1; i__ <= i__2; ++i__) {
		cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
	    }
	    i__ = iphes[mp];
	    if (i__ == mp) {
		goto L640;
	    }
	    zsafe = cont[mp];
	    cont[mp] = cont[i__];
	    cont[i__] = zsafe;
L640:
	    ;
	}
/* ----------------------------------- */
L88:
	*err = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    werr[i__] = cont[i__] / scal[i__];
/* Computing 2nd power */
	    d__1 = werr[i__];
	    *err += d__1 * d__1;
	}
/* Computing MAX */
	d__1 = sqrt(*err / *n);
	*err = max(d__1,1e-10);
    }
    return 0;
/* ----------------------------------------------------------- */
L55:
    return 0;
} /* estrad_ */


/*     END OF SUBROUTINE ESTRAD */

/* *********************************************************** */

/* Subroutine */ int estrav_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, doublereal *h__, doublereal *dd, S_fp 
	fcn, integer *nfcn, doublereal *y0, doublereal *y, integer *ijob, 
	doublereal *x, integer *m1, integer *m2, integer *nm1, integer *ns, 
	integer *nns, doublereal *e1, integer *lde1, doublereal *zz, 
	doublereal *cont, doublereal *ff, integer *ip1, integer *iphes, 
	doublereal *scal, doublereal *err, logical *first, logical *reject, 
	doublereal *fac1, doublereal *rpar, integer *ipar)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e1_dim1, 
	    e1_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, mm, mp, im1;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum, sum1;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solh_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal zsafe;

    /* Parameter adjustments */
    --scal;
    --iphes;
    --cont;
    --y;
    --y0;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip1;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    --dd;
    --ff;
    --zz;
    e1_dim1 = *lde1;
    e1_offset = 1 + e1_dim1;
    e1 -= e1_offset;
    --rpar;
    --ipar;

    /* Function Body */
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L13;
	case 14:  goto L14;
	case 15:  goto L15;
    }

L1:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__ + *n] = sum / *h__;
	cont[i__] = ff[i__ + *n] + y0[i__];
    }
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L77;

L11:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__ + *n] = sum / *h__;
	cont[i__] = ff[i__ + *n] + y0[i__];
    }
L48:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum1 = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    sum1 = (cont[j + k * *m2] + sum1) / *fac1;
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		im1 = i__ + *m1;
		cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
	    }
	}
    }
    sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
	cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
    }
    goto L77;

L2:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__ + *n] = sum / *h__;
	cont[i__] = ff[i__ + *n] + y0[i__];
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
	    1]);
    goto L77;

L12:
/* ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__ + *n] = sum / *h__;
	cont[i__] = ff[i__ + *n] + y0[i__];
    }
L45:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum1 = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		im1 = i__ + *m1;
		cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) * 
			fjac_dim1] * sum1;
	    }
	}
    }
    solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*m1 + 
	    1], &ip1[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
	cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
    }
    goto L77;

L3:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__4 = *ns;
	for (k = 1; k <= i__4; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__] = sum / *h__;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j];
	}
	ff[i__ + *n] = sum;
	cont[i__] = sum + y0[i__];
    }
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L77;

L13:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__ + *n] = sum / *h__;
	cont[i__] = ff[i__ + *n] + y0[i__];
    }
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__] = sum / *h__;
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
/* Computing MAX */
	i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__2 = min(i__5,i__6);
	for (j = max(i__3,i__4); j <= i__2; ++j) {
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j + *
		    m1];
	}
	im1 = i__ + *m1;
	ff[im1 + *n] = sum;
	cont[im1] = sum + y0[im1];
    }
    goto L48;

L4:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *ns;
	for (k = 1; k <= i__2; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__] = sum / *h__;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
/* Computing MAX */
	i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *n, i__6 = i__ + *mumas;
	i__4 = min(i__5,i__6);
	for (j = max(i__2,i__3); j <= i__4; ++j) {
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j];
	}
	ff[i__ + *n] = sum;
	cont[i__] = sum + y0[i__];
    }
    solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &ip1[
	    1]);
    goto L77;

L14:
/* ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__4 = *ns;
	for (k = 1; k <= i__4; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__ + *n] = sum / *h__;
	cont[i__] = ff[i__ + *n] + y0[i__];
    }
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__4 = *ns;
	for (k = 1; k <= i__4; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__] = sum / *h__;
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
/* Computing MAX */
	i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	i__5 = *nm1, i__6 = i__ + *mumas;
	i__3 = min(i__5,i__6);
	for (j = max(i__4,i__2); j <= i__3; ++j) {
	    sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ff[j + *
		    m1];
	}
	im1 = i__ + *m1;
	ff[im1 + *n] = sum;
	cont[im1] = sum + y0[im1];
    }
    goto L45;

L5:
/* ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__] = sum / *h__;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *n;
	for (j = 1; j <= i__3; ++j) {
	    sum += fmas[i__ + j * fmas_dim1] * ff[j];
	}
	ff[i__ + *n] = sum;
	cont[i__] = sum + y0[i__];
    }
    sol_(n, lde1, &e1[e1_offset], &cont[1], &ip1[1]);
    goto L77;

L15:
/* ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    i__1 = *m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__ + *n] = sum / *h__;
	cont[i__] = ff[i__ + *n] + y0[i__];
    }
    i__1 = *n;
    for (i__ = *m1 + 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__] = sum / *h__;
    }
    i__1 = *nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *nm1;
	for (j = 1; j <= i__3; ++j) {
	    sum += fmas[i__ + j * fmas_dim1] * ff[j + *m1];
	}
	im1 = i__ + *m1;
	ff[im1 + *n] = sum;
	cont[im1] = sum + y0[im1];
    }
    goto L48;

L6:
/* ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ------  THIS OPTION IS NOT PROVIDED */
    return 0;

L7:
/* ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__3 = *ns;
	for (k = 1; k <= i__3; ++k) {
	    sum += dd[k] * zz[i__ + (k - 1) * *n];
	}
	ff[i__ + *n] = sum / *h__;
	cont[i__] = ff[i__ + *n] + y0[i__];
    }
    for (mm = *n - 2; mm >= 1; --mm) {
	mp = *n - mm;
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L310;
	}
	zsafe = cont[mp];
	cont[mp] = cont[i__];
	cont[i__] = zsafe;
L310:
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
	    cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
	}
    }
    solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
    i__1 = *n - 2;
    for (mm = 1; mm <= i__1; ++mm) {
	mp = *n - mm;
	i__3 = *n;
	for (i__ = mp + 1; i__ <= i__3; ++i__) {
	    cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
	}
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L440;
	}
	zsafe = cont[mp];
	cont[mp] = cont[i__];
	cont[i__] = zsafe;
L440:
	;
    }

/* -------------------------------------- */

L77:
    *err = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = cont[i__] / scal[i__];
	*err += d__1 * d__1;
    }
/* Computing MAX */
    d__1 = sqrt(*err / *n);
    *err = max(d__1,1e-10);

    if (*err < 1.) {
	return 0;
    }
    if (*first || *reject) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    cont[i__] = y[i__] + cont[i__];
	}
	(*fcn)(n, x, &cont[1], &ff[1], &rpar[1], &ipar[1]);
	++(*nfcn);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    cont[i__] = ff[i__] + ff[i__ + *n];
	}
	switch (*ijob) {
	    case 1:  goto L31;
	    case 2:  goto L32;
	    case 3:  goto L31;
	    case 4:  goto L32;
	    case 5:  goto L31;
	    case 6:  goto L32;
	    case 7:  goto L33;
	    case 8:  goto L55;
	    case 9:  goto L55;
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
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
	    sum1 = 0.;
	    for (k = mm - 1; k >= 0; --k) {
		sum1 = (cont[j + k * *m2] + sum1) / *fac1;
		i__3 = *nm1;
		for (i__ = 1; i__ <= i__3; ++i__) {
		    im1 = i__ + *m1;
		    cont[im1] += fjac[i__ + (j + k * *m2) * fjac_dim1] * sum1;
		}
	    }
	}
	sol_(nm1, lde1, &e1[e1_offset], &cont[*m1 + 1], &ip1[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
	    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
	}
	goto L88;
/* ------ BANDED MATRIX OPTION */
L32:
	solb_(n, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[1], &
		ip1[1]);
	goto L88;
/* ------ BANDED MATRIX OPTION, SECOND ORDER */
L42:
	i__1 = *m2;
	for (j = 1; j <= i__1; ++j) {
	    sum1 = 0.;
	    for (k = mm - 1; k >= 0; --k) {
		sum1 = (cont[j + k * *m2] + sum1) / *fac1;
/* Computing MAX */
		i__3 = 1, i__4 = j - *mujac;
/* Computing MIN */
		i__5 = *nm1, i__6 = j + *mljac;
		i__2 = min(i__5,i__6);
		for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
		    im1 = i__ + *m1;
		    cont[im1] += fjac[i__ + *mujac + 1 - j + (j + k * *m2) * 
			    fjac_dim1] * sum1;
		}
	    }
	}
	solb_(nm1, lde1, &e1[e1_offset], &linal_1.mle, &linal_1.mue, &cont[*
		m1 + 1], &ip1[1]);
	for (i__ = *m1; i__ >= 1; --i__) {
	    cont[i__] = (cont[i__] + cont[*m2 + i__]) / *fac1;
	}
	goto L88;
/* ------ HESSENBERG MATRIX OPTION */
L33:
	for (mm = *n - 2; mm >= 1; --mm) {
	    mp = *n - mm;
	    i__ = iphes[mp];
	    if (i__ == mp) {
		goto L510;
	    }
	    zsafe = cont[mp];
	    cont[mp] = cont[i__];
	    cont[i__] = zsafe;
L510:
	    i__1 = *n;
	    for (i__ = mp + 1; i__ <= i__1; ++i__) {
		cont[i__] -= fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
	    }
	}
	solh_(n, lde1, &e1[e1_offset], &c__1, &cont[1], &ip1[1]);
	i__1 = *n - 2;
	for (mm = 1; mm <= i__1; ++mm) {
	    mp = *n - mm;
	    i__2 = *n;
	    for (i__ = mp + 1; i__ <= i__2; ++i__) {
		cont[i__] += fjac[i__ + (mp - 1) * fjac_dim1] * cont[mp];
	    }
	    i__ = iphes[mp];
	    if (i__ == mp) {
		goto L640;
	    }
	    zsafe = cont[mp];
	    cont[mp] = cont[i__];
	    cont[i__] = zsafe;
L640:
	    ;
	}
/* ----------------------------------- */
L88:
	*err = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = cont[i__] / scal[i__];
	    *err += d__1 * d__1;
	}
/* Computing MAX */
	d__1 = sqrt(*err / *n);
	*err = max(d__1,1e-10);
    }
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* estrav_ */


/*     END OF SUBROUTINE ESTRAV */

/* *********************************************************** */

/* Subroutine */ int slvrod_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *fac1, doublereal *e, integer *lde, integer *ip, 
	doublereal *dy, doublereal *ak, doublereal *fx, doublereal *ynew, 
	doublereal *hd, integer *ijob, logical *stage1)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e_dim1, e_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, mm, im1, jkm;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *);


    /* Parameter adjustments */
    --ynew;
    --fx;
    --ak;
    --dy;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;

    /* Function Body */
    if (*hd == 0.) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ak[i__] = dy[i__];
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ak[i__] = dy[i__] + *hd * fx[i__];
	}
    }

    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
	case 6:  goto L6;
	case 7:  goto L55;
	case 8:  goto L55;
	case 9:  goto L55;
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
    if (*stage1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ak[i__] += ynew[i__];
	}
    }
    sol_(n, lde, &e[e_offset], &ak[1], &ip[1]);
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    if (*stage1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ak[i__] += ynew[i__];
	}
    }
L48:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    jkm = j + k * *m2;
	    sum = (ak[jkm] + sum) / *fac1;
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		im1 = i__ + *m1;
		ak[im1] += fjac[i__ + jkm * fjac_dim1] * sum;
	    }
	}
    }
    sol_(nm1, lde, &e[e_offset], &ak[*m1 + 1], &ip[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
	ak[i__] = (ak[i__] + ak[*m2 + i__]) / *fac1;
    }
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    if (*stage1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ak[i__] += ynew[i__];
	}
    }
    solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]);
    return 0;

/* ----------------------------------------------------------- */

L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    if (*stage1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ak[i__] += ynew[i__];
	}
    }
L45:
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    jkm = j + k * *m2;
	    sum = (ak[jkm] + sum) / *fac1;
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		im1 = i__ + *m1;
		ak[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * sum;
	    }
	}
    }
    solb_(nm1, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[*m1 + 1], &
	    ip[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
	ak[i__] = (ak[i__] + ak[*m2 + i__]) / *fac1;
    }
    return 0;

/* ----------------------------------------------------------- */

L3:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX */
    if (*stage1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = 0.;
/* Computing MAX */
	    i__4 = 1, i__2 = i__ - *mlmas;
/* Computing MIN */
	    i__5 = *n, i__6 = i__ + *mumas;
	    i__3 = min(i__5,i__6);
	    for (j = max(i__4,i__2); j <= i__3; ++j) {
		sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[
			j];
	    }
	    ak[i__] += sum;
	}
    }
    sol_(n, lde, &e[e_offset], &ak[1], &ip[1]);
    return 0;

/* ----------------------------------------------------------- */

L13:
/* ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    if (*stage1) {
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ak[i__] += ynew[i__];
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = 0.;
/* Computing MAX */
	    i__3 = 1, i__4 = i__ - *mlmas;
/* Computing MIN */
	    i__5 = *nm1, i__6 = i__ + *mumas;
	    i__2 = min(i__5,i__6);
	    for (j = max(i__3,i__4); j <= i__2; ++j) {
		sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[
			j + *m1];
	    }
	    im1 = i__ + *m1;
	    ak[im1] += sum;
	}
    }
    if (*ijob == 14) {
	goto L45;
    }
    goto L48;

/* ----------------------------------------------------------- */

L4:
/* ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX */
    if (*stage1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = 0.;
/* Computing MAX */
	    i__2 = 1, i__3 = i__ - *mlmas;
/* Computing MIN */
	    i__5 = *n, i__6 = i__ + *mumas;
	    i__4 = min(i__5,i__6);
	    for (j = max(i__2,i__3); j <= i__4; ++j) {
		sum += fmas[i__ - j + linal_1.mbdiag + j * fmas_dim1] * ynew[
			j];
	    }
	    ak[i__] += sum;
	}
    }
    solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]);
    return 0;

/* ----------------------------------------------------------- */

L5:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX */
    if (*stage1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = 0.;
	    i__4 = *n;
	    for (j = 1; j <= i__4; ++j) {
		sum += fmas[i__ + j * fmas_dim1] * ynew[j];
	    }
	    ak[i__] += sum;
	}
    }
    sol_(n, lde, &e[e_offset], &ak[1], &ip[1]);
    return 0;

/* ----------------------------------------------------------- */

L15:
/* ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER */
    if (*stage1) {
	i__1 = *m1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ak[i__] += ynew[i__];
	}
	i__1 = *nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = 0.;
	    i__4 = *nm1;
	    for (j = 1; j <= i__4; ++j) {
		sum += fmas[i__ + j * fmas_dim1] * ynew[j + *m1];
	    }
	    im1 = i__ + *m1;
	    ak[im1] += sum;
	}
    }
    goto L48;

/* ----------------------------------------------------------- */

L6:
/* ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX */
/* ---  THIS OPTION IS NOT PROVIDED */
    if (*stage1) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sum = 0.;
	    i__4 = *n;
	    for (j = 1; j <= i__4; ++j) {
/* L623: */
		sum += fmas[i__ + j * fmas_dim1] * ynew[j];
	    }
/* L624: */
	    ak[i__] += sum;
	}
	solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &ak[1], &ip[1]
		);
    }
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* slvrod_ */


/*     END OF SUBROUTINE SLVROD */


/* *********************************************************** */

/* Subroutine */ int slvseu_(integer *n, doublereal *fjac, integer *ldjac, 
	integer *mljac, integer *mujac, doublereal *fmas, integer *ldmas, 
	integer *mlmas, integer *mumas, integer *m1, integer *m2, integer *
	nm1, doublereal *fac1, doublereal *e, integer *lde, integer *ip, 
	integer *iphes, doublereal *del, integer *ijob)
{
    /* System generated locals */
    integer fjac_dim1, fjac_offset, fmas_dim1, fmas_offset, e_dim1, e_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static integer i__, j, k, mm, mp, im1, mp1, jkm, mmm;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static doublereal sum;
    extern /* Subroutine */ int solb_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *), solh_(integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    static doublereal zsafe;


    /* Parameter adjustments */
    --del;
    --iphes;
    fjac_dim1 = *ldjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --ip;
    fmas_dim1 = *ldmas;
    fmas_offset = 1 + fmas_dim1;
    fmas -= fmas_offset;
    e_dim1 = *lde;
    e_offset = 1 + e_dim1;
    e -= e_offset;

    /* Function Body */
    switch (*ijob) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L1;
	case 4:  goto L2;
	case 5:  goto L1;
	case 6:  goto L55;
	case 7:  goto L7;
	case 8:  goto L55;
	case 9:  goto L55;
	case 10:  goto L55;
	case 11:  goto L11;
	case 12:  goto L12;
	case 13:  goto L11;
	case 14:  goto L12;
	case 15:  goto L11;
    }

/* ----------------------------------------------------------- */

L1:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX */
    sol_(n, lde, &e[e_offset], &del[1], &ip[1]);
    return 0;

/* ----------------------------------------------------------- */

L11:
/* ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER */
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    jkm = j + k * *m2;
	    sum = (del[jkm] + sum) / *fac1;
	    i__2 = *nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		im1 = i__ + *m1;
		del[im1] += fjac[i__ + jkm * fjac_dim1] * sum;
	    }
	}
    }
    sol_(nm1, lde, &e[e_offset], &del[*m1 + 1], &ip[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
	del[i__] = (del[i__] + del[*m2 + i__]) / *fac1;
    }
    return 0;

/* ----------------------------------------------------------- */

L2:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX */
    solb_(n, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &del[1], &ip[1]);
    return 0;

/* ----------------------------------------------------------- */

L12:
/* ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER */
    mm = *m1 / *m2;
    i__1 = *m2;
    for (j = 1; j <= i__1; ++j) {
	sum = 0.;
	for (k = mm - 1; k >= 0; --k) {
	    jkm = j + k * *m2;
	    sum = (del[jkm] + sum) / *fac1;
/* Computing MAX */
	    i__2 = 1, i__3 = j - *mujac;
/* Computing MIN */
	    i__5 = *nm1, i__6 = j + *mljac;
	    i__4 = min(i__5,i__6);
	    for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		im1 = i__ + *m1;
		del[im1] += fjac[i__ + *mujac + 1 - j + jkm * fjac_dim1] * 
			sum;
	    }
	}
    }
    solb_(nm1, lde, &e[e_offset], &linal_1.mle, &linal_1.mue, &del[*m1 + 1], &
	    ip[1]);
    for (i__ = *m1; i__ >= 1; --i__) {
	del[i__] = (del[i__] + del[*m2 + i__]) / *fac1;
    }
    return 0;

/* ----------------------------------------------------------- */

L7:
/* ---  HESSENBERG OPTION */
    for (mmm = *n - 2; mmm >= 1; --mmm) {
	mp = *n - mmm;
	mp1 = mp - 1;
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L110;
	}
	zsafe = del[mp];
	del[mp] = del[i__];
	del[i__] = zsafe;
L110:
	i__1 = *n;
	for (i__ = mp + 1; i__ <= i__1; ++i__) {
	    del[i__] -= fjac[i__ + mp1 * fjac_dim1] * del[mp];
	}
    }
    solh_(n, lde, &e[e_offset], &c__1, &del[1], &ip[1]);
    i__1 = *n - 2;
    for (mmm = 1; mmm <= i__1; ++mmm) {
	mp = *n - mmm;
	mp1 = mp - 1;
	i__4 = *n;
	for (i__ = mp + 1; i__ <= i__4; ++i__) {
	    del[i__] += fjac[i__ + mp1 * fjac_dim1] * del[mp];
	}
	i__ = iphes[mp];
	if (i__ == mp) {
	    goto L240;
	}
	zsafe = del[mp];
	del[mp] = del[i__];
	del[i__] = zsafe;
L240:
	;
    }
    return 0;

/* ----------------------------------------------------------- */

L55:
    return 0;
} /* slvseu_ */