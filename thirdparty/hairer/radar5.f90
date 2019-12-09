! ---------------------------------------------------------- 
!     NUMERICAL SOLUTION OF A STIFF DIFFERENTIAL  
!     (OR DIFFERENTIAL ALGEBRAIC) SYSTEM OF FIRST 0RDER  
!     DELAY DIFFERENTIAL EQUATIONS   
!                     M*Y'(X)=F(X,Y(X),Y(X-A),...). 
!     THE SYSTEM CAN BE (LINEARLY) IMPLICIT (MASS-MATRIX M .NE. I) 
!     OR EXPLICIT (M=I). 
!     NOTE: THIS FORM ALSO ALLOWS TO SOLVE NEUTRAL DIFFERENTIAL PROBLEMS 
! 
!     NOTE: THIS VERSION ALLOWS ARBITRARILY LARGE STEPSIZES 
!     (POSSIBLY LARGER THAN THE DELAYS) 
! 
!     THE METHOD USED IS AN IMPLICIT RUNGE-KUTTA METHOD (3 STAGE  
!     RADAU IIA) OF ORDER 5 WITH STEP SIZE CONTROL AND CONTINUOUS  
!     EXTENSION OF ORDER 3  (C.F. SECTION IV.8 OF (HW)) 
! 
!     AUTHORS: N. GUGLIELMI(*) AND E. HAIRER($)  
!          (*) UNIVERSITA` DELL'AQUILA, DIP. DI MATEMATICA 
!              VIA VETOIO (COPPITO), 67010 L'AQUILA, ITALY 
!          ($) UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES 
!              CH-1211 GENEVE 24, SWITZERLAND  
!                        ---------------------- 
!              E-MAIL ADRESSES:   
!                                         guglielm@univaq.it 
!                                     Ernst.Hairer@math.unige.ch 
!      
!     THIS PROGRAM EXTENDS THE CODE RADAU5 (BY E. HAIRER AND G. WANNER) 
!     TO THE CASE OF DELAY DIFFERENTIAL EQUATIONS.  
!     DETAILS ABOUT RADAU5 CAN BE FOUND IN THE BOOK: 
!     (HW)  E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL 
!           EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS. 
!           SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14, 
!           SPRINGER-VERLAG 1991, SECOND EDITION 1996. 
!     DETAILS ABOUT RADAR5 CAN BE FOUND IN THE PAPERS:
!     (GH1)
!     (GH2)       
!
!     VERSION 2.1 OF JULY 21, 2005 
! 
! ---------------------------------------------------------- 
      MODULE IP_ARRAY 
!        THIS VECTOR IPOSV HAS THE DIMENSION OF THE MAXIMUM NUMBER OF 
!        ALLOWED DELAYS; IF A LARGER NUMBER OF RETARDED ARGUMENT IS  
!        REQUIRED CHANGE THE DIMENSION TO THE DESIRED VALUE AND RECOMPILE 
!        HERE THE DIMENSION IS SET TO 100 
         INTEGER, dimension(10000) :: IPOSV 
      END MODULE IP_ARRAY 
! 
      SUBROUTINE RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H, &
                        RTOL,ATOL,ITOL, &
                        JAC,IJAC,MLJAC,MUJAC, &
                        JACLAG,NLAGS,NJACL, &
                        IMAS,SOLOUT,IOUT, &
                        WORK,IWORK,RPAR,IPAR,IDID, &
                        GRID,IPAST,MAS,MLMAS,MUMAS, &
                        PAST,LRPAST) 
! ---------------------------------------------------------- 
!     INPUT PARAMETERS   
! --------------------   
!     N           DIMENSION OF THE SYSTEM  
! 
!     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE RIGHT- 
!                 HAND-SIDE OF THE DELAY EQUATION, E.G., 
!                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR,...) 
!                    DOUBLE PRECISION X,Y(N),F(N) 
!                    EXTERNAL PHI 
!                    F(1)=G1(X,Y(*),YLAGR5(*,X-TAU(X,Y(*))),PHI,...)) 
!                    F(2)=G2(X,Y(*),YLAGR5(*,X-TAU(X,Y(*))),PHI,...)) 
!                    ETC. 
!                    (*) MEANS ALL POSSIBLE COMPONENTS 
!                 FOR AN EXPLICATION OF YLAGR5 SEE BELOW. 
!                 DO NOT USE YLAGR5(I,X-0.D0,PHI,RPAR,IPAR,...) ! 
!                 Note: 
!                 THE INITIAL FUNCTION HAS TO BE SUPPLIED BY: 
!                    FUNCTION PHI(I,X,RPAR,IPAR) 
!                    DOUBLE PRECISION PHI,X 
!                 WHERE I IS THE COMPONENT AND X THE ARGUMENT 
!                 RPAR, IPAR (SEE BELOW) 
! 
!     X           INITIAL X-VALUE 
! 
!     Y(N)        INITIAL VALUES FOR Y (MAY BE DIFFERENT FROM PHI (I,X), 
!                 IN THIS CASE IT IS HIGHLY RECOMMENDED TO SET IWORK(13) 
!                 AND GRID(1),..., (SEE BELOW) 
! 
!     XEND        FINAL X-VALUE (XEND-X HAS TO BE POSITIVE) 
! 
!     H           INITIAL STEP SIZE GUESS; 
!                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,  
!                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD. 
!                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS 
!                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6). 
! 
!     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY 
!                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N. 
! 
!     ITOL        SWITCH FOR RTOL AND ATOL: 
!                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS. 
!                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF 
!                     Y(I) OVER RTOL*ABS(Y(I))+ATOL 
!                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS. 
!                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) OVER 
!                     RTOL(I)*ABS(Y(I))+ATOL(I). 
! 
!     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES 
!                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y 
!                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1;  
!                 THE USER HAS TO SUPPLY A DUMMY SUBROUTINE  
!                 IN THE CASE IJAC=0). 
!                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM 
!                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR,...) 
!                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N) 
!                    DFY(1,1)= ... 
!                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS 
!                 FURNISHED BY THE CALLING PROGRAM. 
!                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO 
!                    BE FULL AND THE PARTIAL DERIVATIVES ARE 
!                    STORED IN DFY AS 
!                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J) 
!                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND 
!                    THE PARTIAL DERIVATIVES ARE STORED 
!                    DIAGONAL-WISE AS 
!                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J). 
! 
!     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN: 
!                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE 
!                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED. 
!                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC. 
! 
!     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN: 
!                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR 
!                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. 
!                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN  
!                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW 
!                       THE MAIN DIAGONAL). 
! 
!     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON- 
!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). 
!                 DOES NOT NEED TO BE DEFINED IF MLJAC=N. 
! 
!     JACLAG      NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES 
!                 THE PARTIAL DERIVATIVES OF F(X,Y,YLAG) WITH RESPECT TO  
!                 YLAG(*) (YLAG DENOTE THE DELAYED VARIABLES) 
! 
!     NLAGS       DENOTES THE NUMBER OF DELAY ARGUMENTS.  
!                 THIS PARAMETER IS OF INTEREST FOR THE COMPUTATION OF THE 
!                 JACOBIAN. 
!                 TO BE SET = 0 IF ONE DOES WANT TO COMPUTE THE TRADITIONAL 
!                 JACOBIAN;  
!                 TO BE SET = NUMBER OF DISTINCT DELAY ARGUMENTS 
!                 IF ONE WANTS TO CORRECT THE STANDARD JACOBIAN (THROUGH 
!                 THE SUBROUTINE JACLAG) WHEN ADVANCED ARGUMENTS ARE USED. 
! 
!     NJACL       NUMBER OF TERMS IN THE JACOBIAN W.R.T. 
!                 RETARDED COMPONENTS (WHICH IS THOUGHT AS A SPARSE MATRIX). 
! 
!     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      ----- 
!     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): - 
! 
!     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS- 
!                 MATRIX M. 
!                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY 
!                 MATRIX AND NEEDS NOT TO BE DEFINED; 
!                 THE USER HAS TO SUPPLY A DUMMY SUBROUTINE IN THIS CASE. 
!                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM 
!                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR) 
!                    DOUBLE PRECISION AM(LMAS,N) 
!                    AM(1,1)= .... 
!                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED 
!                    AS FULL MATRIX LIKE 
!                         AM(I,J) = M(I,J) 
!                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED 
!                    DIAGONAL-WISE AS 
!                         AM(I-J+MUMAS+1,J) = M(I,J). 
! 
!     IMAS       GIVES INFORMATION ON THE MASS-MATRIX: 
!                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY 
!                       MATRIX, MAS IS NEVER CALLED. 
!                    IMAS=1: MASS-MATRIX  IS SUPPLIED. 
! 
!     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX: 
!                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR 
!                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. 
!                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE 
!                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW 
!                       THE MAIN DIAGONAL). 
!                 MLMAS IS SUPPOSED TO BE <= MLJAC. 
! 
!     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON- 
!                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). 
!                 DOES NOT NEED TO BE DEFINED IF MLMAS=N. 
!                 MUMAS IS SUPPOSED TO BE <= MUJAC. 
! 
!     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE 
!                 NUMERICAL SOLUTION DURING INTEGRATION.  
!                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP. 
!                 THE USER HAS TO SUPPLY A DUMMY SUBROUTINE IF IOUT=0.  
!                 IT MUST HAVE THE FORM 
!                    SUBROUTINE SOLOUT (NR,XOLD,X,HSOL,Y,CONT,LRC,N, 
!                                       RPAR,IPAR,IRTRN) 
!                    DOUBLE PRECISION X,Y(N),CONT(LRC) 
!                    ....   
!                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH 
!                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS 
!                    THE FIRST GRID-POINT). 
!                 "XOLD" IS THE PRECEEDING GRID-POINT. 
!                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN 
!                    IS SET <0, RADAR5 RETURNS TO THE CALLING PROGRAM. 
!            
!          -----  CONTINUOUS OUTPUT: ----- 
!                 DURING CALLS TO "SOLOUT" AS WELL AS TO "FCN", A 
!                 CONTINUOUS SOLUTION IS AVAILABLE THROUGH HTHE FUNCTION 
!                        >>>   YLAGR5(I,S,PHI,RPAR,IPAR,...)   <<< 
!                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH 
!                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE S 
!                 HAS TO LIE IN AN INTERVAL WHERE THE NUMERICAL SOLUTION 
!                 IS ALREADY COMPUTED. IT DEPENDS ON THE SIZE OF LRPAST 
!                 (SEE BELOW) HOW FAR BACK THE SOLUTION IS AVAILABLE. 
! 
!     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT: 
!                    IOUT=0: SUBROUTINE IS NEVER CALLED 
!                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT. 
! 
!     WORK        ARRAY OF STATE VARIABLES OF REAL TYPE FOR EXECUTION. 
!                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS 
!                 FOR THE CODE. FOR STANDARD USE OF THE CODE 
!                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE 
! 
!     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK". 
!                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS 
!                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),.., 
!                 IWORK(20) TO ZERO BEFORE CALLING. 
! 
!     GRID        CONTAINS PRESCRIBED GRID POINTS, WHICH THE 
!                 INTEGRATION METHOD HAS TO TAKE AS GRID-POINTS 
!                 NORMALLY, IF GRID(1) > X, THEN ONE HAS 
!                 X < GRID(1) < GRID(2) < ... < GRID(NGRID) <= XEND 
!                 IN SOME CASES IF THERE ARE DISCONTINUITIES IN THE  
!                 INITIAL FUNCTIONS THEY ARE ALSO SET IN THE GRID  
!                 VECTOR; THEN X < GRID(J) < GRID(J+1) ... < XEND 
!                 WHERE J ADDRESSES THE FIRST ASCISSA IN GRID > X 
! 
!     LGRID       DECLARED LENGTH OF GRID VECTOR, 
!                 GRID(LGRID), 
!                 WHICH MUST BE DECLARED IN THE CALLING PROGRAM. 
!                 "LGRID" MUST BE AT LEAST 
!                              NGRID + 1 
! 
!     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH   
!                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING 
!                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.  
! 
! ---------------------------------------------------------------------- 
!  
!     SOPHISTICATED SETTING OF PARAMETERS 
!     ----------------------------------- 
!              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK  
!              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),... 
!              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO. 
!              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES: 
! 
!    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN 
!              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY 
!              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN. 
!              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N) 
!              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1). 
! 
!    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. 
!              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000. 
! 
!    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE 
!              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP. 
!              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7. 
! 
!    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION 
!              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD. 
!              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED. 
!              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS 
!              DIFFICULTIES WITH CONVERGENCE (THIS IS SEEN IN THE CASE WHEN 
!              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.). 
!              DEFAULT IS IWORK(4)=0. 
! 
!       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR 
!       DELAY DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1. 
!       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT 
!       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER.  
!       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE 
!       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2. 
! 
!    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR  
!              DDE'S THIS EQUALS THE DIMENSION OF THE SYSTEM. 
!              DEFAULT IWORK(5)=N. 
! 
!    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0. 
! 
!    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0. 
! 
!    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY 
!              IF IWORK(8).EQ.1  MODIFIED PREDICTIVE CONTROLLER  
!              (GUSTAFSSON) 
!              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL 
!              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1. 
!              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS; 
!              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES 
!              OFTEN SLIGHTLY FASTER RUNS 
! 
!       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT 
!            Y(I)' = Y(I+M2)   FOR  I=1,...,M1, 
!       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME 
!       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10). 
!       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE  
!       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2. 
!       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS: 
!       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE 
!              JACOBIAN HAVE TO BE STORED 
!              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL 
!                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J) 
!                FOR I=1,N-M1 AND J=1,N. 
!              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM ) 
!                 DFY(I-J+MUJAC+1,J+K*M2) = PARTIAL F(I+M1) / PARTIAL Y(J+K*M2) 
!                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM. 
!       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL 
!                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM) 
!                     PARTIAL F(I+M1) / PARTIAL Y(J+K*M2),  I,J=1,M2 
!                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH 
!                    OF THESE MM+1 SUBMATRICES 
!       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES; 
!                DOES NOT NEED TO BE DEFINED IF MLJAC=N-M1 
!       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND 
!              DOES NOT NEED TO BE DEFINED.  
!              THE USER HAS TO SUPPLY A DUMMY SUBROUTINE IN THIS CASE. 
!              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF 
!              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX. 
!              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL 
!                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1. 
!              ELSE, THE MASS MATRIX IS BANDED 
!                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1) 
!       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL 
!                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX 
!       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX 
!                DOES NOT NEED TO BE DEFINED IF MLMAS=N-M1 
! 
!    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0. 
! 
!    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1. 
! 
! 
! 
!    IWORK(11) SELECT THE TYPE OF ERROR CONTROL:  
!              -1: FOR A PURE CONTROL OF THE DENSE OUTPUT  
!                  (MAKES USE OF A QUADRATIC AND A LINEAR INTERPOLATING 
!                   POLYNOMIALS);  
!               1: FOR A MIXED CONTROL OF DENSE OUTPUT AND DISCRETE OUTPUT 
!               2: FOR A PURE CONTROL OF THE DISCRETE OUTPUT 
!                  (ERROR CONTROL PROVIDED BY THE SUBROUTINE ESTRAD). 
!               3: FOR A SIMPLER MIXED CONTROL OF DENSE OUTPUT AND   
!                  DISCRETE OUTPUT 
!               DEFAULT VALUE IWORK(11)=2. 
! 
!    IWORK(12) = MXST = (ON INPUT)  
!              DECLARED NUMBER OF STEPS STORED IN THE ``PAST VECTOR'',  
!              PAST(LRPAST), 
!              WHICH MUST BE DECLARED IN THE CALLING PROGRAM. 
!              "MXST" MUST BE SUFFICIENTLY LARGE. IF THE DENSE 
!              OUTPUT OF MXST BACK STEPS HAS TO BE STORED,  
!              THE DIMENSION OF PAST MUST BE  
!                       LRPAST=MXST*(4*NRDENS+2) 
!              WHERE NRDENS=IWORK(15) (SEE BELOW). 
! 
!    IWORK(13) = NGRID = (ON INPUT) 
!              NUMBER OF PRESCRIBED POINTS IN THE 
!              INTEGRATION INTERVAL WHICH HAVE TO BE GRID-POINTS 
!              IN THE INTEGRATION. USUALLY, AT THESE POINTS THE 
!              SOLUTION OR ONE OF ITS DERIVATIVE HAS A DISCONTINUITY. 
!              DEFINE THESE POINTS IN GRID(1),...,GRID(NGRID) 
!              DEFAULT VALUE:  IWORK(13)=0 
! 
!    IWORK(14) = SELECTOR FOR FULL ITERATION (2) OR SIMPLIFIED  
!              ITERATION (1) (TAKING INTO ACCOUNT POSSIBLE  
!              ADVANCED ARGUMENTS BUT PRESERVING TENSOR STRUCTURE  
!              OF THE JACOBIAN. 
!              DEFAULT VALUE:  IWORK(14)=1 
! 
!    IWORK(15) = NRDENS = (ON INPUT)  
!              NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT 
!              IS REQUIRED (EITHER BY "SOLOUT" OR BY "FCN"); 
!              DEFAULT VALUE (FOR IWORK(15)=0) IS IWORK(15)=N; 
!              FOR   0 < NRDENS < N   THE COMPONENTS (FOR WHICH DENSE 
!              OUTPUT IS REQUIRED) HAVE TO BE SPECIFIED IN 
!              IPAST(1),...,IPAST(NRDENS); 
!              FOR  NRDENS=N  THIS IS DONE BY THE CODE. 
!    IWORK(16) = NDIMN = (ON INPUT)
!              OPTION VALID FOR NEUTRAL PROBLEMS
!              NUMBER OF DERIVATIVE COMPONENTS (Z) OF THE NEUTRAL PROBLEM 
!              EXCLUDED TRUE SOLUTION COMPONENTS
! 
! ---------- 
! 
!    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16. 
! 
!    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION, 
!              DEFAULT 0.9D0. 
! 
!    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; 
!              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS 
!              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER  
!              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO 
!              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP.      
!              DEFAULT 0.001D0. 
! 
!    WORK(4)   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. 
!              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER. 
!              DEFAULT MIN(0.03D0,RTOL(1)**0.5D0) 
! 
!    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE 
!              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A 
!              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR 
!              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE 
!              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS 
!              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD. 
!              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 . 
! 
!    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X. 
! 
!    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION 
!              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION 
!                 WORK(8) <= HNEW/HOLD <= WORK(9) 
!              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0 
! 
!    WORK(10)  PARAMETER FOR CONTROLLING THE ERROR CONTROL OF DENSE 
!              OUTPUT (0 <= WORK(10) <= 1). (0: STRONG CONTROL, 1: WEAKER) 
!              SUGGESTED VALUES: 
!              FOR PROBLEMS WITH `ALMOST DISCONTINUOUS' SOLUTIONS  
!              (LIKE SHOCKS):  WORK(10)=0.D0 
!              FOR PROBLEMS WITH FAIRLY SMOOTH SOLUTION:  WORK(10)=1.D0 
!              FOR INTERMEDIATE PROBLEMS:  WORK(10)=1.D-M (M=1,2,3,..,) 
!              DEFAULT VALUE: WORK(10)=0.D0  
!    WORK(11)  PARAMETER FOR CONTROLLING THE SEARCH OF BREAKING POINTS:
!              IF THE ERROR INCREASES OF A FACTOR LARGER THAN WORK(11)
!              FROM A STEP TO THE SUBSEQUENT, THE ROUTINE SEARCHING BREAKING
!              POINTS ACTIVATES.
!              DEFAULT VALUE: WORK(11)=5.D0
!----------------------------------------------------------------------- 
! 
!     OUTPUT PARAMETERS  
!     -----------------  
!     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED 
!                 (AFTER SUCCESSFUL RETURN X=XEND). 
! 
!     Y(N)        NUMERICAL SOLUTION AT X 
!  
!     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP 
! 
!     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN: 
!                   IDID= 1  COMPUTATION SUCCESSFUL, 
!                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT) 
!                   IDID=-1  INPUT IS NOT CONSISTENT, 
!                   IDID=-2  LARGER NMAX IS NEEDED, 
!                   IDID=-3  STEP SIZE BECOMES TOO SMALL, 
!                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR. 
!                   IDID=-5  COMPUTATION INTERRUPTED BY YLAGR5.    
!                   IDID=-6  THE EQUATION USES ADVANCED ARGUMENTS 
! 
!   IWORK(13)  NFULL   NUMBER OF FULL NEWTON ITERATIONS                    
!   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL 
!                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)   
!   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY 
!                      OR NUMERICALLY) 
!   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS 
!   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS 
!   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST), 
!                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED) 
!   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF BOTH MATRICES 
!   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF BOTH 
!                      SYSTEMS; THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS, 
!                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED 
!----------------------------------------------------------------------- 
! *** *** *** *** *** *** *** *** *** *** *** *** *** 
!          DECLARATIONS  
! *** *** *** *** *** *** *** *** *** *** *** *** *** 
      USE IP_ARRAY 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
      REAL(kind=DP), dimension(N), intent(inout) ::  &
                                  Y 
      REAL(kind=DP), dimension(30), intent(inout) ::  &
                                  WORK 
      REAL(kind=DP), dimension(1), intent(inout) ::  &
                                  ATOL,RTOL 
      INTEGER, dimension(30), intent(inout) :: IWORK 
      REAL(kind=DP), dimension(1), intent(inout) :: GRID 
      INTEGER, dimension(1), intent(inout) :: IPAST 
      REAL(kind=DP), dimension(1), intent(in) :: RPAR 
      INTEGER, dimension(1), intent(in) :: IPAR 
 
      REAL(kind=DP), dimension(LRPAST), intent(inout) :: PAST
      INTEGER, intent(in) :: LRPAST
 
      LOGICAL IMPLCT,NEUTRAL,JBAND,ARRET,STARTN,PRED 
      LOGICAL FLAGS,FLAGN 
 
      EXTERNAL FCN,PHI,ARGLAG,JAC,JACLAG,MAS,SOLOUT 
! ----> COMMON BLOCKS <---- 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
 
! *** *** *** *** *** *** *** 
!        SETTING THE PARAMETERS  
! *** *** *** *** *** *** *** 
!      IF (NLAGS.GT.0) THEN 
!       ALLOCATE (IPOSV(NLAGS)) 
!      ELSE 
!       ALLOCATE (IPOSV(5)) 
!      END IF 
       NN=N 
 
       NFCN=0 
       NJAC=0 
       NSTEP=0 
       NACCPT=0 
       NREJCT=0 
       NDEC=0 
       NSOL=0 
       ARRET=.FALSE. 
       FLAGS=.FALSE. 
       FLAGN=.FALSE. 
 
!       IF (IOUT.EQ.1) WRITE (6,*) 'STARTING INTEGRATION...' 
!        
! ------> OPERATIONS RELEVANT TO THE DELAY DEPENDENCE <------ 
! 
! -------- ERROR CONTROL 
      IF (IWORK(11).EQ.0) THEN 
       IEFLAG=2 
      ELSE 
       IEFLAG=IWORK(11) 
      END IF 
      IF (IEFLAG.EQ.2) WORK(10)=1.D0 
! -------- NGRID   NUMBER OF PRESCRIBED GRID-POINTS 
      NGRID=IWORK(13) 
      IF (NGRID.LT.0) NGRID=0 
!      IF (IOUT.EQ.1) WRITE(6,*)  
!     &           'NUMBER OF PRESCRIBED GRID POINTS: ',NGRID 
! ------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS 
      NRDENS=IWORK(15)
! ------- NDIMN   NUMBER OF COMPONENTS OF A NEUTRAL PROBLEM
      IF (IMAS.EQ.2) THEN
      IF (IWORK(16).EQ.0) THEN
       WRITE(6,*) 'NUMBER OF Y COMPONENTS HAS TO BE SPECIFIED'
       ARRET=.TRUE.
      END IF
       NDIMN=IWORK(16)
      ELSE 
       NDIMN=N
      END IF     
! ------- LIPAST   DIMENSION OF VECTOR IPAST 
      LIPAST=NRDENS+1  
      IF(NRDENS.LT.0.OR.NRDENS.GT.N) THEN 
         IF (IOUT.GT.0) WRITE(6,*) &
                ' CURIOUS INPUT IWORK(15)=',IWORK(15) 
         ARRET=.TRUE. 
      ELSE IF (NRDENS.EQ.0) THEN 
            NRDS=N 
      ELSE 
            NRDS=NRDENS 
      END IF 
      IF (NRDS.EQ.N) THEN 
            DO 16 I=1,NRDS 
  16           IPAST(I)=I 
      END IF 
!      IF (IOUT.EQ.1) WRITE(6,*) 'NUMBER OF DELAYED COMPONENTS: ',NRDS 
! ------- LRPAST   DIMENSION OF VECTOR PAST 
      MXST=IWORK(12) 
! ------- CONTROL OF LENGTH OF PAST  ------- 
      IF(MXST.LT.1)THEN 
         IF (IOUT.GT.0) WRITE(6,*) &
      ' INSUFFICIENT STORAGE FOR PAST, MIN. LRPAST=',1 
         ARRET=.TRUE. 
      END IF 
! ------- DIM. of PAST  -------- 
      IDIF=4*NRDS+2 
!      LRPAST=MXST*IDIF 
!      WRITE (6,*) ' AT INITIALIZATION, LRPAST = ',LRPAST 
! -------------------------------------------------  
! ------- CONTROL OF SIMPLE NEWTON ITERATION  ------- 
      ISWJL=IWORK(14) 
      IF (ISWJL.EQ.0) ISWJL=1 
 
! -------- UROUND : SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0   
      IF (WORK(1).EQ.0.0D0) THEN 
         UROUND=1.0D-16 
      ELSE 
         UROUND=WORK(1) 
         IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN 
            WRITE(6,*)' COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
! -------> CHECK AND CHANGE THE TOLERANCES <------ 
      EXPM=2.0D0/3.0D0 
      IF (ITOL.EQ.0) THEN 
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN 
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL' 
              ARRET=.TRUE. 
          ELSE 
              QUOT=ATOL(1)/RTOL(1) 
              RTOL(1)=0.1D0*RTOL(1)**EXPM 
              ATOL(1)=RTOL(1)*QUOT 
          END IF 
      ELSE 
          DO I=1,N 
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN 
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL' 
              ARRET=.TRUE. 
          ELSE 
              QUOT=ATOL(I)/RTOL(I) 
              RTOL(I)=0.1D0*RTOL(I)**EXPM 
              ATOL(I)=RTOL(I)*QUOT 
          END IF 
          END DO 
      END IF 
 
! -------> NMAX : THE MAXIMAL NUMBER OF STEPS <------- 
      IF (IWORK(2).EQ.0) THEN 
         NMAX=100000 
      ELSE 
         NMAX=IWORK(2) 
         IF (NMAX.LE.0) THEN 
            WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
! -------> NIT :  MAXIMAL NUMBER OF NEWTON ITERATIONS <------- 
      IF (IWORK(3).EQ.0) THEN 
         NIT=7 
      ELSE 
         NIT=IWORK(3) 
         IF (NIT.LE.0) THEN 
            WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3) 
            ARRET=.TRUE. 
         END IF 
      END IF 
! -------- STARTN : SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS 
      IF(IWORK(4).EQ.0)THEN 
         STARTN=.FALSE. 
      ELSE 
         STARTN=.TRUE. 
      END IF 
 
! -------> PARAMETERS (IF ANY) FOR DIFFERENTIAL-ALGEBRAIC COMPONENTS <------- 
      NIND1=IWORK(5) 
      NIND2=IWORK(6) 
      NIND3=IWORK(7) 
      IF (NIND1.EQ.0) NIND1=N 
      IF (NIND1+NIND2+NIND3.NE.N) THEN 
       WRITE(6,*)' CURIOUS INPUT FOR IWORK(5,6,7)=',NIND1,NIND2,NIND3 
       ARRET=.TRUE. 
      END IF 
 
! -------> PRED   STEP SIZE CONTROL <------- 
      IF(IWORK(8).LE.1)THEN 
         PRED=.TRUE. 
      ELSE 
         PRED=.FALSE. 
      END IF 
 
! -------> PARAMETER FOR SECOND ORDER EQUATIONS <------- 
      M1=IWORK(9) 
      M2=IWORK(10) 
      NM1=N-M1 
      IF (M1.EQ.0) M2=N 
      IF (M2.EQ.0) M2=M1 
      IF (M1.LT.0.OR.M2.LT.0.OR.M1+M2.GT.N) THEN 
       WRITE(6,*)' CURIOUS INPUT FOR IWORK(9,10)=',M1,M2 
       ARRET=.TRUE. 
      END IF 
 
! -------> SAFE  :  SAFETY FACTOR IN STEP SIZE PREDICTION <------- 
      IF (WORK(2).EQ.0.0D0) THEN 
         SAFE=0.9D0 
      ELSE 
         SAFE=WORK(2) 
         IF (SAFE.LE.0.001D0.OR.SAFE.GE.1.0D0) THEN 
            WRITE(6,*)' CURIOUS INPUT FOR WORK(2)=',WORK(2) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
! ------> THET : DETERMINES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; 
      IF (WORK(3).EQ.0.D0) THEN 
         THET=0.001D0 
      ELSE 
         THET=WORK(3) 
         IF (THET.GE.1.0D0) THEN 
            WRITE(6,*)' CURIOUS INPUT FOR WORK(3)=',WORK(3) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
! ---> FNEWT : STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. <--- 
      TOLST=RTOL(1) 
      IF (WORK(4).EQ.0.D0) THEN 
         FNEWT=MAX(10*UROUND/TOLST,MIN(0.03D0,TOLST**0.5D0)) 
      ELSE 
         FNEWT=WORK(4) 
         IF (FNEWT.LE.UROUND/TOLST) THEN 
            WRITE(6,*)' CURIOUS INPUT FOR WORK(4)=',WORK(4) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
! ---> QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST. <--- 
      IF (WORK(5).EQ.0.D0) THEN 
         QUOT1=1.D0 
      ELSE 
         QUOT1=WORK(5) 
      END IF 
      IF (WORK(6).EQ.0.D0) THEN 
         QUOT2=1.2D0 
      ELSE 
         QUOT2=WORK(6) 
      END IF 
      IF (QUOT1.GT.1.0D0.OR.QUOT2.LT.1.0D0) THEN 
         WRITE(6,*)' CURIOUS INPUT FOR WORK(5,6)=',QUOT1,QUOT2 
         ARRET=.TRUE. 
      END IF 
! -------------------------------------------------------  
 
! ---->    GRID WITH DISCONTINUITIES  <---- 
      XURO=100*UROUND*ABS(XEND) 
      IF (NGRID.GT.0) THEN 
         IF (GRID(NGRID)-XEND.GE.XURO) THEN 
            IF(IOUT.GT.0) WRITE(6,*) &
                ' GRID(NGRID) HAS TO BE <= XEND' 
            ARRET=.TRUE. 
         ENDIF 
         IF (ABS(GRID(NGRID)-XEND).GE.XURO) NGRID=NGRID+1 
      ELSE 
         NGRID=NGRID+1 
      END IF 
      GRID(NGRID)=XEND 
      
!      WRITE(6,*) 'FINAL GRID: '
!      DO I = 1, NGRID
!        WRITE(6,*) GRID(I)
!      END DO
! -------------------------------------------------------  
 
! -------> MAXIMAL STEP SIZE <------- 
      IF (WORK(7).EQ.0.D0) THEN 
         DO I=1,NGRID 
          IF (GRID(I).GT.X) THEN  
           IGRID=I 
           GO TO 2 
          END IF 
         END DO 
 2       CONTINUE 
         HMAX=GRID(IGRID)-X 
      ELSE 
         HMAX=WORK(7) 
      END IF  
 
! ------->  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION <------- 
      IF(WORK(8).EQ.0.D0)THEN 
         FACL=5.D0 
      ELSE 
         FACL=1.D0/WORK(8) 
      END IF 
      IF(WORK(9).EQ.0.D0)THEN 
         FACR=1.D0/8.0D0 
      ELSE 
         FACR=1.D0/WORK(9) 
      END IF 
      IF (FACL.LT.1.0D0.OR.FACR.GT.1.0D0) THEN 
            WRITE(6,*)' CURIOUS INPUT WORK(8,9)=',WORK(8),WORK(9) 
            ARRET=.TRUE. 
      END IF 
! ------->  PARAMETER FOR THE CONTROL OF DENSE OUTPUT <------- 
      ALPHA=WORK(10) 
      IF (ALPHA.LT.0.D0.OR.ALPHA.GT.1.D0) THEN 
            WRITE(6,*)' CURIOUS INPUT WORK(10)=',WORK(10) 
            ARRET=.TRUE. 
      END IF 
! ------->   PARAMETER FOR CONTROLLING THE SEARCH OF BP <-------
      TCKBP=WORK(11)
      IF (TCKBP.LE.0.D0) THEN
           TCKBP=5.D0
      END IF           
   
! *** *** *** *** *** *** *** *** *** *** *** *** *** 
!         COMPUTATION OF ARRAY ENTRIES 
! *** *** *** *** *** *** *** *** *** *** *** *** *** 
! ---- IMPLICIT, BANDED OR NOT ? 
      IMPLCT=IMAS.NE.0 
       NEUTRAL=IMAS.EQ.2
      JBAND=MLJAC.LT.NM1 
! -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- 
! -- JACOBIAN  AND  MATRICES E1, E2 
      IF (JBAND) THEN 
         LDJAC=MLJAC+MUJAC+1 
         LDE1=MLJAC+LDJAC 
      ELSE 
         MLJAC=NM1 
         MUJAC=NM1 
         LDJAC=NM1 
         LDE1=NM1 
      END IF 
! -- MASS MATRIX 
      IF (IMPLCT) THEN 
          IF (MLMAS.NE.NM1) THEN 
              LDMAS=MLMAS+MUMAS+1 
              IF (JBAND) THEN 
                 IJOB=4 
              ELSE 
                 IJOB=3 
              END IF 
          ELSE 
              LDMAS=NM1 
              IJOB=5 
          END IF 
! ------ BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC" 
          IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN 
             WRITE (6,*) 'BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF &
      "JAC"' 
            ARRET=.TRUE. 
          END IF 
      ELSE 
          LDMAS=0 
          IF (JBAND) THEN 
             IJOB=2 
          ELSE 
             IJOB=1 
             IF (N.GT.2.AND.IWORK(1).NE.0) IJOB=7 
          END IF 
      END IF 
      LDMAS2=MAX(1,LDMAS) 
! ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN 
      IF ((IMPLCT.OR.JBAND).AND.IJOB.EQ.7) THEN 
         WRITE(6,*)' HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS WITH  &
     FULL JACOBIAN' 
         ARRET=.TRUE. 
      END IF 
 
! ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 
      IF (ARRET) THEN 
         IDID=-1 
         RETURN 
      END IF 
!      WRITE(6,*), (IWORK(j), j = 1,30)
!      WRITE(6,*), (WORK(j), j = 1,30)
!  136 format (30I3.1)
!     NUMERICAL KERNEL
!      WRITE(6,*) 'INTEGRATION...'
! -------- CALL TO CORE INTEGRATOR ------------ 
      CALL RADCOR(N,X,Y,XEND,H,FCN,PHI,ARGLAG,RTOL,ATOL,ITOL, &
         JAC,IJAC,MLJAC,MUJAC,JACLAG,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID, &
         NMAX,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,IJOB,STARTN, &
         NIND1,NIND2,NIND3,PRED,FACL,FACR,M1,M2,NM1, &
         IMPLCT,NEUTRAL,NDIMN,JBAND,LDJAC,LDE1,LDMAS2, &
         NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,NFULL,RPAR,IPAR, &
         IPAST,GRID,NRDS,NLAGS,NJACL, &
         NGRID,IEFLAG,WORK(7),TCKBP,ALPHA,ISWJL, &
         PAST,LRPAST) 
      IWORK(13)=NFULL 
      IWORK(14)=NFCN 
      IWORK(15)=NJAC 
      IWORK(16)=NSTEP 
      IWORK(17)=NACCPT 
      IWORK(18)=NREJCT 
      IWORK(19)=NDEC 
      IWORK(20)=NSOL 
! -------- RESTORE TOLERANCES 
      EXPM=1.0D0/EXPM 
      IF (ITOL.EQ.0) THEN 
              QUOT=ATOL(1)/RTOL(1) 
              RTOL(1)=(10.0D0*RTOL(1))**EXPM 
              ATOL(1)=RTOL(1)*QUOT 
      ELSE 
          DO I=1,N 
              QUOT=ATOL(I)/RTOL(I) 
              RTOL(I)=(10.0D0*RTOL(I))**EXPM 
              ATOL(I)=RTOL(I)*QUOT 
          END DO 
      END IF 
! ----------- RETURN ----------- 
!     DEALLOCATE (IPOSV) 
      RETURN 
      END 
! 
!     END OF SUBROUTINE RADAR5 
! 
! *********************************************************** 
! 
 
      SUBROUTINE RADCOR(N,X,Y,XEND,H,FCN,PHI,ARGLAG,RTOL,ATOL,ITOL, & 
         JAC,IJAC,MLJAC,MUJAC,JACLAG,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID, &
         NMAX,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,IJOB,STARTN, &
         NIND1,NIND2,NIND3,PRED,FACL,FACR,M1,M2,NM1, &
         IMPLCT,NEUTRAL,NDIMN,BANDED,LDJAC,LDE1,LDMAS, &
         NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,NFULL,RPAR,IPAR, &
         IPAST,GRID,NRDS,NLAGS,NJACL, &
         NGRID,IEFLAG,WORK7,TCKBP,ALPHA,ISWJL, &
         PAST,LRPAST) 
! ---------------------------------------------------------- 
!     CORE INTEGRATOR FOR RADAR5 
!     PARAMETERS SAME AS IN RADAR5 WITH WORKSPACE ADDED  
! ----------------------------------------------------------  
!         DECLARATIONS  
! ----------------------------------------------------------  
!     use definitions 
      USE IP_ARRAY 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
      REAL(kind=DP), dimension(1), intent(inout) :: Y 
      REAL(kind=DP), dimension(:), allocatable :: & 
                                   Z1,Z2,Z3,Y0,SCAL,F1,F2,F3,CONT
      REAL(kind=DP), dimension(LRPAST), intent(inout) :: PAST 
      REAL(kind=DP), dimension(:), allocatable :: BPV,UCONT 
      REAL(kind=DP), dimension(:,:), allocatable ::  &
                                   FJAC,FJACS,FMAS,E1,E2R,E2I 
      REAL(kind=DP), dimension(1), intent(inout) ::  &
                                   ATOL,RTOL 
      REAL(kind=DP), dimension(1), intent(in) ::  &
                                   RPAR 
      INTEGER, dimension(1), intent(in) ::  &
                                   IPAR,IPAST 
      REAL(kind=DP), dimension(1), intent(in) ::  &
                                   GRID 
      REAL(kind=DP), dimension(:,:), allocatable :: & 
                                   FJACL,XLAG 
      REAL(kind=DP), dimension(:), allocatable :: &
                                   FJACLAG,ZL  
      INTEGER, dimension(:,:), allocatable :: ICOUN 
      INTEGER, dimension(:), allocatable ::  &
                                   IVL,IVE,IVC,ILS 
      INTEGER, dimension(:), allocatable ::  &
                                   IP1,IP2,IPHES,IPJ 
!      INTEGER :: LRPAST
       INTEGER, intent(in) :: LRPAST
 
      LOGICAL FLAGS,FLAGN,FLAGUS 
      LOGICAL QUADR 
      LOGICAL BPC,BPD,BPDMEM,LEFT 
      LOGICAL REPEAT
! ----> COMMON BLOCKS <---- 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
      COMMON /BPLOG/FIRST,LAST,REJECT,BPD 
      COMMON /BPCOM/BPP,ILBP,LEFT 
      COMMON /LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG 
 
! ---- 
      LOGICAL REJECT,FIRST,IMPLCT,NEUTRAL,PROJECT,BANDED,CALJAC 
      LOGICAL STARTN,CALHES,CALJACL,CALLAG 
      LOGICAL INDEX1,INDEX2,INDEX3,LAST,PRED 
      EXTERNAL FCN,PHI,ARGLAG 
! *** *** *** *** *** *** *** 
!  INITIALISATIONS 
! *** *** *** *** *** *** *** 
       
      ALLOCATE (Z1(N),Z2(N),Z3(N),Y0(N),SCAL(N),F1(N),F2(N),F3(N)) 
      ALLOCATE(BPV(10000)) 
      ALLOCATE (FJAC(LDJAC,N),ZL(3*N)) 
      IF (IMPLCT) ALLOCATE(FMAS(LDMAS,NM1)) 
      ALLOCATE (IP1(NM1),IP2(NM1),IPHES(NM1)) 
      ALLOCATE (E1(LDE1,NM1),E2R(LDE1,NM1),E2I(LDE1,NM1)) 
!      LRPAST = MXST*IDIF
!      ALLOCATE (PAST(MXST*IDIF)) 
!      WRITE (6,*) ' IN RADCOR, LRPAST = ',LRPAST
      IF (NLAGS.GT.0) THEN 
       ALLOCATE (FJACS(LDJAC,N),FJACLAG(NJACL)) 
       ALLOCATE (IVL(NJACL),IVE(NJACL),IVC(NJACL), &
                ILS(2*NLAGS+NJACL),ICOUN(3,NLAGS)) 
       IF (ISWJL.NE.1) THEN 
        ALLOCATE (IPJ(3*N),FJACL(3*N,3*N)) 
       END IF 
       ALLOCATE (XLAG(3,NLAGS)) 
      END IF 
  
!     AMPLITUDE OF CONT 
      LRC=4*N 
      ALLOCATE (CONT(LRC)) 
      ALLOCATE (UCONT(LRC+2)) 
! --- 
!       OPEN(8,FILE='radar5.log')
!       REWIND 8

! -------------------------------------------------  
      BPC=.FALSE. 
      BPD=.FALSE. 
      QUADR=.FALSE. 
       
! --- INITIAL PREPARATIONS 
      CALLAG=.FALSE. 
      IACT=1 
      IPOS=1 
!      DO I=1,NLAGS 
       DO I=1, SIZE(IPOSV)
       IPOSV(I)=1 
      END DO 
!      WRITE(6,*) 'INIT IPOSV', IPOSV(1),NLAGS
      X0B=X 
      DO I=1,NGRID 
       IF (GRID(I).GT.X0B) THEN  
        IGRID=I 
        GO TO 2 
       END IF 
      END DO 
 2    XEND=GRID(IGRID) 
      IBP=1 
      BPV(1)=X0B 
      BTOL=1.D1*UROUND 
      KMAX=10
      IMANT=0
       
! --- GUSTAFFSON TECHNIQUE AFTER BREAKING POINTS IS NOT APPLIED        
      FLAGUS=.FALSE. 
      ERRACC=1.D0
 
      IRTRN=2 
!      WRITE(6,*) 'INIT FCN6'
      CALL FCN(N,X,Y,Y0,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS, &
       LRPAST) 
!      WRITE(6,*) 'END INIT FNC6'
      IRTRN=1 
 
!     TOLMIN
      IF (ITOL.EQ.0) THEN
       RTOLM=RTOL(1)
      ELSE
       RTOLM=RTOL(1)
       DO I=2,N
        IF (RTOL(I).LT.RTOLM) RTOLM=RTOL(I)
       END DO
      END IF
 
! -------- CHECK THE INDEX OF THE PROBLEM -----  
      INDEX1=NIND1.NE.0 
      INDEX2=NIND2.NE.0 
      INDEX3=NIND3.NE.0 
! ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- 
      IF (IMPLCT) CALL MAS(NM1,FMAS,LDMAS,RPAR,IPAR) 
! ---------> REQUIRED CONSTANTS <--------- 
      SQ6=DSQRT(6.D0) 
      C1=(4.D0-SQ6)/10.D0 
      C2=(4.D0+SQ6)/10.D0 
      C1M1=C1-1.D0 
      C2M1=C2-1.D0 
      C1MC2=C1-C2 
      CQ1=(2.D0+3.D0*SQ6)/6.D0 
      CQ2=(2.D0-3.D0*SQ6)/6.D0 
      CQ3=1.D0/3.D0 
      CL1=10.D0/(6.D0+SQ6)                   
      CL2=0.D0                 
      CL3=(-4.D0+SQ6)/(6.D0+SQ6)  
      CERS=5.D-1 
      CERC=5.D-1 
      CERLQ=1.D-2 
      THRS=100.D0 
      DD1=-(13.D0+7.D0*SQ6)/3.D0 
      DD2=(-13.D0+7.D0*SQ6)/3.D0 
      DD3=-1.D0/3.D0 
      U1=(6.D0+81.D0**(1.D0/3.D0)-9.D0**(1.D0/3.D0))/30.D0 
      ALPH=(12.D0-81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))/60.D0 
      BETA=(81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))*DSQRT(3.D0)/60.D0 
      CNO=ALPH**2+BETA**2 
      U1=1.0D0/U1 
      ALPH=ALPH/CNO 
      BETA=BETA/CNO 
      T11=9.1232394870892942792D-02 
      T12=-0.14125529502095420843D0 
      T13=-3.0029194105147424492D-02 
      T21=0.24171793270710701896D0 
      T22=0.20412935229379993199D0 
      T23=0.38294211275726193779D0 
      T31=0.96604818261509293619D0 
      TI11=4.3255798900631553510D0 
      TI12=0.33919925181580986954D0 
      TI13=0.54177053993587487119D0 
      TI21=-4.1787185915519047273D0 
      TI22=-0.32768282076106238708D0 
      TI23=0.47662355450055045196D0 
      TI31=-0.50287263494578687595D0 
      TI32=2.5719269498556054292D0 
      TI33=-0.59603920482822492497D0 
! 
!     INVERSE OF A 
      IF (NLAGS.GT.0.AND.ISWJL.NE.1) THEN 
       AI11= 3.22474487139158904909864D0 
       AI12= 1.16784008469040549492404D0 
       AI13=-0.25319726474218082618594D0 
       AI21=-3.56784008469040549492404D0 
       AI22= 0.77525512860841095090136D0 
       AI23= 1.05319726474218082618594D0 
       AI31= 5.53197264742180826185942D0 
       AI32=-7.53197264742180826185942D0 
       AI33= 5.00000000000000000000000D0 
      END IF 
! 
      IF (M1.GT.0) IJOB=IJOB+10 
      HMAXN=MIN(HMAX,XEND-X)  
      IF (H.LE.10.D0*UROUND) H=1.0D-6  
      H=MIN(H,HMAXN) 
      HOLD=H 
      REJECT=.FALSE. 
      FIRST=.TRUE. 
      LAST=.FALSE. 
      NITER=0
      IF ((X+H*1.0001D0-XEND).GE.0.D0) THEN 
         H=XEND-X 
         LAST=.TRUE. 
      END IF 
! ---  INITIALIZATION FOR THE ARRAY PAST
        DO I=1,MXST*IDIF
           PAST(I)=0.D0
        ENDDO
       DO 3 I=0,MXST-1 
          PAST(1+IDIF*I)=X  
   3   CONTINUE 
       IPA=(MXST-1)*IDIF+1 
       DO J=1,NRDS 
          K=IPAST(J) 
          PAST(J+IPA)=Y(K) 
          PAST(J+1*NRDS+IPA)=0.D0 
          PAST(J+2*NRDS+IPA)=0.D0 
          PAST(J+3*NRDS+IPA)=0.D0 
       ENDDO 
!       WRITE (6,*) 
!     &        ' THE INITIALIZATION WAS RUN.', X, Y(1), H
          PAST(IPA+IDIF-1)=H  
! ---  END OF THE INITIALIZATION      
      FACCON=1.D0 
      CFAC=SAFE*(1+2*NIT) 
      NSING=0 
      XOLD=X 
      IF (IOUT.NE.0) THEN 
          IRTRN=1 
          NRSOL=1 
          XOSOL=XOLD 
          XSOL=X 
          CONT(1:N)=Y(1:N) 
          NSOLU=N 
          HSOL=HOLD
          CALL SOLOUT(NRSOL,XOSOL,XSOL,HSOL,Y,CONT,LRC,NSOLU, &
                     RPAR,IPAR,IRTRN) 
          IF (IRTRN.LT.0) GOTO 179 
      END IF 
      MLE=MLJAC 
      MUE=MUJAC 
      MBJAC=MLJAC+MUJAC+1 
      MBB=MLMAS+MUMAS+1 
      MDIAG=MLE+MUE+1 
      MDIFF=MLE+MUE-MUMAS 
      MBDIAG=MUMAS+1 
      N2=2*N 
      N3=3*N 
      IF (ITOL.EQ.0) THEN 
          DO I=1,N 
             SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I)) 
          END DO 
      ELSE 
          DO I=1,N 
             SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I)) 
          END DO 
      END IF 
      HHFAC=H 
      NFCN=NFCN+1 
 
      NFULL=0
      
!      WRITE(6,*) 'COMMON BLOCKS'
!      WRITE(6,*) XEND,IGRID
!      PRINT *, GRID 
!      PRINT *, IVL
!      ,IVE,IVC,ILS
!      WRITE(6,*) FLAGS,FLAGN,FLAGUS,QUADR BPC,BPD,BPDMEM,LEFT,REPEAT
!      WRITE(6,*) REJECT,FIRST,IMPLCT,NEUTRAL,PROJECT,BANDED,CALJAC 
!      WRITE(6,*) STARTN,CALHES,CALJACL,CALLAG 
!      WRITE(6,*) INDEX1,INDEX2,INDEX3,LAST,PRED
!      PRINT *, ATOL, RTOL
!      WRITE(6,*) X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
!      WRITE(6,*) C1,C2,C1M1,C2M1,C1MC2 
!      WRITE(6,*) FIRST,LAST,REJECT,BPD 
!      WRITE(6,*) BPP,ILBP,LEFT 
!      WRITE(6,*) MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
!      WRITE(6,137) ( IPOSV(j), j=1,100 ) 

!  137 format (100I3.1)
! -------------------------- 
! --- BASIC INTEGRATION STEP   
! -------------------------- 
  10  CONTINUE 
! *** *** *** *** *** *** *** 
!  COMPUTATION OF THE JACOBIAN 
! *** *** *** *** *** *** *** 
! ----------------------- 
      FLAGS=.FALSE. 
      FLAGN=.FALSE. 
! ----------------------- 
      ALOPT=0.D0 
      NJAC=NJAC+1 
      IF (BPD) THEN
       BPDMEM=.TRUE.
       BPD=.FALSE.
      END IF
      IF (IJAC.EQ.0) THEN 
! --- COMPUTE JACOBIAN MATRIX NUMERICALLY 
         IF (BANDED) THEN 
! --- JACOBIAN IS BANDED 
            MUJACP=MUJAC+1 
            MD=MIN(MBJAC,M2) 
            DO MM=1,M1/M2+1 
               DO K=1,MD 
                  J=K+(MM-1)*M2 
 12               F1(J)=Y(J) 
                  F2(J)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J)))) 
                  Y(J)=Y(J)+F2(J) 
                  J=J+MD 
                  IF (J.LE.MM*M2) GOTO 12
!                  WRITE(6,*) 'BANDED JAC FCN9'
                  CALL FCN(N,X,Y,CONT,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST, &
                          NRDS, &
       LRPAST) 
!                  WRITE(6,*) 'END BANDED JAC FCN9'
                  J=K+(MM-1)*M2 
                  J1=K 
                  LBEG=MAX(1,J1-MUJAC)+M1 
 14               LEND=MIN(M2,J1+MLJAC)+M1 
                  Y(J)=F1(J) 
                  MUJACJ=MUJACP-J1-M1 
                  DO L=LBEG,LEND 
                     FJAC(L+MUJACJ,J)=(CONT(L)-Y0(L))/F2(J)  
                  END DO 
                  J=J+MD 
                  J1=J1+MD 
                  LBEG=LEND+1 
                  IF (J.LE.MM*M2) GOTO 14 
               END DO 
            END DO 
         ELSE 
! --- JACOBIAN IS FULL 
            DO I=1,N 
               YSAFE=Y(I) 
              DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE))) 
                Y(I)=YSAFE+DELT
!                 WRITE(6,*) 'JAC FCN8'
                 CALL FCN(N,X,Y,CONT,ARGLAG,PHI, &
                       RPAR,IPAR,PAST,IPAST,NRDS, &
       LRPAST) 
!                 WRITE(6,*) 'END JAC FCN8'
              DO J=M1+1,N 
                 FJAC(J-M1,I)=(CONT(J)-Y0(J))/DELT 
               END DO 
               Y(I)=YSAFE 
            END DO
         END IF 
      ELSE 
! --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
         CALL JAC(N,X,Y,FJAC,LDJAC,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST, &
        NRDS) 
      END IF
      IF (BPDMEM) THEN
      BPDMEM=.FALSE.
      BPD=.TRUE.
      END IF 
      CALJAC=.TRUE. 
      CALHES=.TRUE. 
      JLFLAG=0 
! --- SAVE FJAC 
      IF (NLAGS.GT.0) FJACS=FJAC 
 
! ------------------------------------------------ 
! --- GLOBAL ITERATION BEGINS HERE
! ------------------------------------------------ 
  20  CONTINUE 
 
! *** *** *** *** *** *** *** 
!  STARTING VALUES FOR NEWTON ITERATION 
! *** *** *** *** *** *** *** 
    
      IF (FIRST.OR.STARTN) THEN 
            DO I=1,N 
             Z1(I)=0.D0 
             Z2(I)=0.D0 
             Z3(I)=0.D0 
             F1(I)=0.D0 
             F2(I)=0.D0 
             F3(I)=0.D0 
! 
             A1=Y(I) 
             ZL(I)=A1 
             ZL(I+N)=A1 
             ZL(I+N2)=A1 
            END DO 
      ELSE 
         C3Q=H/HOLD 
         C1Q=C1*C3Q 
         C2Q=C2*C3Q 
         DO I=1,N 
            A1=Y(I) 
            AK1=CONT(I+N) 
            AK2=CONT(I+N2) 
            AK3=CONT(I+N3) 
            Z1I=C1Q*(AK1+(C1Q-C2M1)*(AK2+(C1Q-C1M1)*AK3)) 
            Z2I=C2Q*(AK1+(C2Q-C2M1)*(AK2+(C2Q-C1M1)*AK3)) 
            Z3I=C3Q*(AK1+(C3Q-C2M1)*(AK2+(C3Q-C1M1)*AK3)) 
            Z1(I)=Z1I 
            Z2(I)=Z2I 
            Z3(I)=Z3I 
! 
            ZL(I)=A1+Z1I 
            F1(I)=TI11*Z1I+TI12*Z2I+TI13*Z3I 
            F2(I)=TI21*Z1I+TI22*Z2I+TI23*Z3I 
            F3(I)=TI31*Z1I+TI32*Z2I+TI33*Z3I 
         END DO 
            IF (NLAGS.GT.0) THEN 
             DO I=1,N 
              A1=Y(I) 
              ZL(I+N)=A1+Z2(I) 
              ZL(I+N2)=A1+Z3(I) 
             END DO 
            END IF 
      END IF 
! --- 
      IF (JLFLAG.EQ.0) THEN 
       IF (NLAGS.EQ.0) THEN 
        IF (.NOT.(REJECT.OR.LAST)) THEN  
         IF (IMANT.EQ.1) GOTO 30 
        END IF 
        GOTO 22  
       END IF 
! --------------------- 
!      INITIALIZATION   
! 
! ---  LAGS CYCLE --- ! 
! --- 
       CALJACL=.FALSE. 
       X1=X+C1*H 
       X2=X+C2*H 
       X3=X+H 
       DO IL=1,NLAGS 
        XLAG(1,IL)=0.D0 
        XLAG(2,IL)=0.D0 
        XLAG(3,IL)=0.D0 
        ICOUN(1,IL)=0 
        ICOUN(2,IL)=0 
        ICOUN(3,IL)=0 
       END DO 
! ---  LOOP ON LAG TERMS 
       DO IL=1,NLAGS 
! ---   DELAYED ARGUMENTS ARE COMPUTED 
        XLAG(1,IL)=ARGLAG(IL,X1,ZL,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST, N) 
           IF (XLAG(1,IL).GT.X) ICOUN(1,IL)=1 
        XLAG(2,IL)=ARGLAG(IL,X2,ZL(N+1),RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST, N)  
           IF (XLAG(2,IL).GT.X) ICOUN(2,IL)=1 
        XLAG(3,IL)=ARGLAG(IL,X3,ZL(N2+1),RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST, N)  
           IF (XLAG(3,IL).GT.X) ICOUN(3,IL)=1 
        IF (ICOUN(1,IL)+ICOUN(2,IL)+ICOUN(3,IL).GE.1) CALJACL=.TRUE.
       END DO 
 
       IF (CALJACL) THEN 
        CALL JACLAG(N,X,Y,FJACLAG,ARGLAG,PHI,IVE,IVC,IVL, &
                       RPAR,IPAR,PAST,IPAST,NRDS) 
        IF (.NOT.CALLAG) THEN 
         CALLAG=.TRUE. 
! --     ORDERING STEP 
         LL=2*NLAGS+1 
         DO L=1,NLAGS 
          NL=0 
          DO I=1,NJACL 
           IF (IVL(I).EQ.L) THEN 
            ILS(LL)=I 
            NL=NL+1 
            LL=LL+1 
           END IF 
          END DO 
          ILS(2*L-1)=NL 
          ILS(2*L)=LL-1 
         END DO 
        END IF 
! 
        DO IL=1,NLAGS 
          DO IS=1,3 
           SELECT CASE (IS) 
            CASE (1) 
             XACT=X1 
            CASE (2) 
             XACT=X2 
            CASE (3) 
             XACT=X3 
           END SELECT  
           IF (XLAG(IS,IL).GT.XACT) THEN 
            IF (IOUT.EQ.2) WRITE (6,*) &
                ' WARNING!: ADVANCED ARGUMENTS ARE USED AT X= ',XACT 
            XLAG(IS,IL)=XACT 
           END IF 
          END DO 
 
! ------ JACOBIAN MAINTAINS THE TENSOR STRUCTURE
! ------ UPDATING CONDITION
          ALOPT=0.D0 
          IF (ICOUN(1,IL).EQ.1) THEN 
             S1=DIM(XLAG(1,IL),X)/H 
             ALOPT=(-1.D0+S1)*S1*(-13.D0-7.D0*Sqrt(6.D0) + &
              5.D0*(2.D0+3.D0*Sqrt(6.D0))*S1) 
          END IF 
          IF (ICOUN(2,IL).EQ.1) THEN 
             S2=DIM(XLAG(2,IL),X)/H 
             ALOPT=ALOPT-(-1+S2)*S2* &
             (13.D0-7.D0*Sqrt(6.D0)+5.D0*(-2.D0+3.D0*Sqrt(6.D0))*S2) 
          END IF 
          IF (ICOUN(3,IL).EQ.1) THEN 
             S3=DIM(XLAG(3,IL),X)/H 
             ALOPT=ALOPT+S3*(1.D0 - 8.D0*S3 + 10.D0*S3**2) 
          END IF 
          ALOPT=ALOPT/9.D0         
          
!         OPTIMAL COEFFICIENT (W.R.T. FROBENIUS NORM)
!         JACLAG ~= ALOPT*I 
! 
!         ACTIVATES IF ALOPT DIFFERENT FROM ZERO
          IF (ABS(ALOPT).GE.1.D-8) THEN 
           NL =ILS(2*IL-1) 
           ILE=ILS(2*IL) 
           DO K=1,NL 
            KK=ILS(ILE-K+1) 
            IK=IVE(KK) 
            JK=IVC(KK) 
            FJAC(IK,JK)=FJAC(IK,JK)+ALOPT*FJACLAG(KK)
!                       ----S         
           END DO 
          ELSE 
           CALJACL=.FALSE. 
          END IF 
        END DO 
! 
       ELSE 
! ---   FACTORIZATION IS CONSERVED  
       END IF 
       IF (.NOT.(REJECT.OR.LAST.OR.CALJACL)) THEN  
        IF (IMANT.EQ.1) GOTO 30 
       END IF 
       GOTO 22  
! --- 
      ELSE IF (JLFLAG.EQ.1) THEN 
! ---  THE FULL ITERATION MAKES USE OF THE SAME STEPSIZE OF THE SIMPLIFIED
       GOTO 23 
! --- 
      END IF 
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! - - THE DIFFERENCE IN THE SOLUTION OF THE SYSTEM STARTS HERE 
! ------------------------------- - - - - - - - - - - - - - - - 
! --- SIMPLE NEWTON ITERATION 
! ------------------------------- 
  22  CONTINUE 
! 
! --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS 
      FAC1=U1/H 
      ALPHN=ALPH/H 
      BETAN=BETA/H 
! --- RK JACOBIAN FACTORIZATION 
      CALL DECOMR(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS, &
                  M1,M2,NM1,FAC1,E1,LDE1,IP1,IER,IJOB,CALHES,IPHES) 
      IF (IER.NE.0) THEN 
          GOTO 78 
      END IF 
      CALL DECOMC(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS, &
                  M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,IP2,IER,IJOB) 
      IF (IER.NE.0) THEN 
          GOTO 78 
      END IF 
      NDEC=NDEC+1 
! --- 
  30  CONTINUE 
! --- UPDATE NSTEP 
      NSTEP=NSTEP+1 
      IF (NSTEP.GT.NMAX) GOTO 178 
      IF (0.1D0*H.LE.ABS(X)*UROUND) GOTO 177 
          IF (INDEX2) THEN 
             DO I=NIND1+1,NIND1+NIND2 
                SCAL(I)=SCAL(I)/HHFAC 
             END DO 
          END IF 
          IF (INDEX3) THEN 
             DO I=NIND1+NIND2+1,NIND1+NIND2+NIND3 
                SCAL(I)=SCAL(I)/(HHFAC*HHFAC) 
             END DO 
          END IF 
      XPH=X+H 
 
! *** *** *** *** *** *** *** 
!  LOOP FOR THE SIMPLE NEWTON ITERATION 
! *** *** *** *** *** *** *** 
!              -----------      
            NEWT=0 
            IMANT=0 
! ----------------------- 
            FLAGS=.FALSE. 
            FLAGN=.FALSE. 
! ----------------------- 
            FACCON=MAX(FACCON,UROUND)**0.8D0 
            THETA=ABS(THET) 
! ------------------------------------------------------- 
!           REFERENCE POINT FOR THE SIMPLE AZIONE SIMPLE 
! ------------------------------------------------------- 
  40        CONTINUE  
! ------------------------- 
            IF (FLAGS) THEN 
             FLAGN=.TRUE. 
! ---        CALLS SUBROUTINE LAGR5 
! ------------------------------------------------------------- 
! ---------- DYNAMIC UPDATE OF INTERPOLANT (in PAST).  
! ------------------------------------------------------------- 
             DO J=1,NRDS 
                I=IPAST(J) 
                PAST(J+IACT)=Y(I)+Z3(I) 
                 Z2I=Z2(I) 
                 Z1I=Z1(I) 
                A1=(Z2I-Z3(I))/C2M1 
                PAST(J+1*NRDS+IACT)=A1 
                 AK=(Z1I-Z2I)/C1MC2 
                 ACONT3=Z1I/C1 
                 ACONT3=(AK-ACONT3)/C2 
                A2=(AK-A1)/C1M1  
                PAST(J+2*NRDS+IACT)=A2 
                PAST(J+3*NRDS+IACT)=A2-ACONT3 
             ENDDO 
! ---        UPDATE DI PAST 
             PAST(IACT)=X 
!            LEFT ENDPOINT 
             PAST(IACT+IDIF-1)=H 
!            STEPSIZE 
            END IF 
! ---------------- 
            IF (NEWT.GE.NIT) THEN 
             INREJ=2 
             GOTO 421 
            END IF 
! ----------------------------------- 
! ---     COMPUTE THE RIGHT-HAND SIDE 
            DO I=1,N 
             A1=Y(I) 
             ZL(I)=A1+Z1(I) 
             ZL(I+N)=A1+Z2(I) 
             ZL(I+N2)=A1+Z3(I) 
            END DO 
!           COMPUTATION OF STAGE VALUES
!           WRITE(6,*) 'STAGE FCN3'
            CALL FCN(N,X+C1*H,ZL,Z1,ARGLAG,PHI,RPAR,IPAR,PAST, &
                     IPAST,NRDS, &
        LRPAST) 
            CALL FCN(N,X+C2*H,ZL(N+1),Z2,ARGLAG,PHI,RPAR,IPAR,PAST, &
                     IPAST,NRDS, &
        LRPAST) 
!            WRITE(6,*) 'END STAGE FCN3'
            IF (BPD) THEN
! ---------------------------------------------------------------------- 
! ---         A BREAKING POINT HAS BEEN DETECTED
! ---------------------------------------------------------------- 
              LEFT=.FALSE. 
!            WRITE(6,*) 'Calling arglag (1)'
              DBP=ARGLAG(ILBP,X+C2*H,ZL(N+1),RPAR,IPAR,PHI, &
                         PAST,IPAST,NRDS, &
        LRPAST,N) 
              IF (DBP.LT.BPP) LEFT=.TRUE. 
            ELSE IF (LAST) THEN 
! ---         DISC. FLAG 
              XX=(X+H)*(1.D0-BTOL) 
              DO I=1,N 
                 ZL(N2+I)=(1.D0-BTOL)*ZL(N2+I)+BTOL*ZL(N+I) 
              END DO
!              WRITE(6,*) 'BREAKING FCN4'
              CALL FCN(N,XX,ZL(N2+1),Z3,ARGLAG,PHI,RPAR,IPAR,PAST, &
                       IPAST,NRDS, &
        LRPAST)
!              WRITE(6,*) 'END BREAKING FCN4'
              GO TO 42 
            END IF
!            WRITE(6,*) 'DETECTED FCN5'
            CALL FCN(N,X+H,ZL(N2+1),Z3,ARGLAG,PHI,RPAR,IPAR,PAST, &
                     IPAST,NRDS, &
        LRPAST) 
!            WRITE(6,*) 'END DETECTED FCN5'
 42         CONTINUE
    
            NFCN=NFCN+3 
! ----------------------------------- 
! ---     SOLVE THE LINEAR SYSTEMS 
           DO I=1,N 
              A1=Z1(I) 
              A2=Z2(I) 
              A3=Z3(I) 
              Z1(I)=TI11*A1+TI12*A2+TI13*A3 
              Z2(I)=TI21*A1+TI22*A2+TI23*A3 
              Z3(I)=TI31*A1+TI32*A2+TI33*A3 
           END DO 
            
        CALL SLVRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS, &
                M1,M2,NM1,FAC1,ALPHN,BETAN,E1,E2R,E2I,LDE1,Z1,Z2,Z3, &
                F1,F2,F3,CONT,IP1,IP2,IPHES,IER,IJOB) 
            NSOL=NSOL+1 
            NEWT=NEWT+1 
! ---       NORM OF DY 
            DYNO=0.D0 
           DO I=1,N
            END DO
            DO I=1,N 
               DENOM=SCAL(I) 
               DYNO=DYNO+(Z1(I)/DENOM)**2+(Z2(I)/DENOM)**2 &
                +(Z3(I)/DENOM)**2 
            END DO 
            DYNO=DSQRT(DYNO/N3) 
! -------------------------------------------------------- 
! ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE 
! --------------------------------------------------------
            IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN 
!               WRITE(6,*) 'DYNOLD:', DYNOLD
                THQ=DYNO/DYNOLD 
                IF (NEWT.EQ.2) THEN 
                   THETA=THQ 
                ELSE 
!                   WRITE(6,*) 'THQOLD:', THQOLD
                   THETA=SQRT(THQ*THQOLD) 
                END IF 
                THQOLD=THQ 
                INREJ=0 
! --- 1 
                IF (THETA.LT.0.99D0) THEN 
                    FACCON=THETA/(1.0D0-THETA) 
                    DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT 
                    IF (DYTH.GE.1.0D0) INREJ=1 
! ----------------  SLOW CONVERGENCE --- 
                ELSE 
                    INREJ=2 
! ----------------  DIVERGENCE --- 
                END IF 
! --- 1 
            ELSE 
! ---           NEWT=1 
                INREJ=0 
            END IF 
! ---------------------------------------- 
! ------------------ THE STEP IS REPEATED   
! ---------------------------------------- 
! ---- 2 
 421        IF (INREJ.GT.0) THEN 
! 
! ----- 3 
             IF (.NOT.(BPC.OR.FIRST)) THEN    
              IF (BPD) THEN
! ---          BP IS WRONG                      
               IBP=IBP-1
               BPD=.FALSE.
              END IF
              HP=H*0.99D0
              WRITE(6,*) 'Calling BPDTCT from row 1661.'
              CALL BPDTCT(N,X,HP,Y,ARGLAG,RPAR,IPAR,UCONT,GRID,NLAGS, &
                          FIRST,LAST,XEND,IGRID,BPV,IBP,ILBP,BPP,BPD, &
                          KMAX,PHI,PAST,IPAST,NRDS) 
            
              BPC=.TRUE. 
! ------ 4 
              IF (BPD) THEN 
               H=HP 
               LAST=.TRUE.
              ELSE IF (INREJ.EQ.1) THEN 
               QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH)) 
               HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT)) 
               H=HHFAC*H 
               LAST=.FALSE. 
!
              ELSE IF (INREJ.EQ.2) THEN 
               H=H*0.55D0 
               HHFAC=0.55D0 
               LAST=.FALSE. 
!
              END IF 
! ------ 4 
             ELSE   
! ----- 3 
              IF (BPD) THEN  
!              BP WRONG
               IBP=IBP-1 
               BPD=.FALSE. 
               LAST=.FALSE. 
              END IF 
! ------------------------------------------- 
              IF (.NOT.CALJACL.OR.ISWJL.EQ.1) THEN  
! ---          THERE ARE NOT SMALL DELAYS 
               IF (REJECT.AND.CALJAC) THEN  
                H=H*0.12D0 
                HHFAC=0.12D0 
               ELSE IF (INREJ.EQ.1) THEN 
                QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH)) 
                HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT)) 
                H=HHFAC*H 
          LAST=.FALSE. 
               ELSE IF (INREJ.EQ.2) THEN 
                H=H*0.55D0 
                HHFAC=0.55D0 
                LAST=.FALSE. 
               END IF 
              ELSE 
! ---          CALJACL IS TRUE: FULL ITEARATION IS DONE          
! ---          THE STEPSIZE DOES NOT CHANGET
               JLFLAG=1 
              END IF 
             END IF 
! ----- 3 
             REJECT=.TRUE. 
             IF (CALJAC) GOTO 20 
             GOTO 10 
            END IF 
! ---- 2 
! -------------------------------------------------------- 
            DYNOLD=MAX(DYNO,UROUND)
            DO I=1,N 
               F1I=F1(I)+Z1(I) 
               F2I=F2(I)+Z2(I) 
               F3I=F3(I)+Z3(I) 
               F1(I)=F1I 
               F2(I)=F2I 
               F3(I)=F3I 
               Z1(I)=T11*F1I+T12*F2I+T13*F3I 
               Z2(I)=T21*F1I+T22*F2I+T23*F3I 
               Z3I=T31*F1I+    F2I 
               Z3(I)=Z3I 
! ---          APPROX DELLA SOLUZIONE 
               ZL(I+N2)=Y(I)+Z3I 
            END DO 
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
            IF (NEWT.EQ.1.OR.FACCON*DYNO.GT.FNEWT) THEN 
! ---        NEWTON PROCEDE
             GOTO 40 
            END IF 

! ----------------------------------------------------------------
! ---       ITERATIVE CORRECTION OF THE BREAKING POINT
! ----------------------------------------------------------------
            IF (BPD) THEN
          HNEWT=H
          NITER = NITER+1
          IF (NITER.GT.NIT) THEN  
!              BP WRONG                 
               IBP=IBP-1 
               BPD=.FALSE. 
               LAST=.FALSE. 
               H=H*0.55D0 
               HHFAC=0.55D0
            REJECT=.TRUE.
            NITER=0 
               IF (CALJAC) GOTO 20 
               GOTO 10   
             END IF  
              
             CALL BPACC(N,X,H,Y,ARGLAG,RPAR,IPAR,Z1,Z2,Z3,FIRST, &
                        BPV,IBP,ILBP,BPP,KMAX,PHI,PAST,IPAST,NRDS) 
             IF (ABS(H-HNEWT)/HNEWT.GE.MAX(BTOL,RTOLM*1.D-2)) THEN
              GOTO 20
!             REF POINT
             ELSE
             H=HNEWT
             NITER=0
             END IF
            END IF  
! ----------------------------------------------------------------


! *** *** *** *** *** *** *** *** *** *** *** 
! END LOOP 
! *** *** *** *** *** *** *** *** *** *** *** 
      GOTO 55 
! *** *** *** *** *** 
! 
 
! --- FULL NEWTON ITERATION 
  23  CONTINUE 
      NFULL =NFULL +1 
! --- 
 
! ---------------------------------------------------- 
! --- ALTERNATIVE FULL NEWTON JACOB. INTEGRATION STEP   
! ---------------------------------------------------- 
! 
      DO I=1,N 
       DO J=1,N 
        FJACL(I+N,J)=0.D0 
        FJACL(I+2*N,J)=0.D0 
        FJACL(I+2*N,J+N)=0.D0 
        FJACL(I,J+N)=0.D0 
        FJACL(I,J+2*N)=0.D0 
        FJACL(I+N,J+2*N)=0.D0 
        FIJ=FJACS(I,J) 
        FJACL(I,J)=FIJ 
        FJACL(I+N,J+N)=FIJ 
        FJACL(I+2*N,J+2*N)=FIJ 
       END DO 
      END DO 
 
      QUADR=FIRST 
      IF (.NOT.QUADR) THEN 
       DL1=C1*(C1-C2)*(C1-1.D0) 
       DL2=C2*(C2-C1)*(C2-1.D0) 
       DL3=(1.D0-C1)*(1.D0-C2) 
      ELSE 
       DL1=(C1-C2)*(C1-1.D0) 
       DL2=(C2-C1)*(C2-1.D0) 
       DL3=(1.D0-C1)*(1.D0-C2) 
      ENDIF 
      DO IL=1,NLAGS 
       DO I=1,3 
            XL=XLAG(I,IL) 
             IF (XL.GT.X) THEN 
! 
!            DERIVATIVES OF THE COLLOCATION POLYNOMIAL WRT STAGES 
!            d u/d Y_k = L_k(xlag_i) 
!                                    
               IF (.NOT.QUADR) THEN 
                DCOLI1=((XL-X)/H)*((XL-X2)/H)*((XL-X3)/H)/DL1 
                DCOLI2=((XL-X)/H)*((XL-X1)/H)*((XL-X3)/H)/DL2 
                DCOLI3=((XL-X)/H)*((XL-X1)/H)*((XL-X2)/H)/DL3 
               ELSE 
                DCOLI1=((XL-X2)/H)*((XL-X3)/H)/DL1 
                DCOLI2=((XL-X1)/H)*((XL-X3)/H)/DL2 
                DCOLI3=((XL-X1)/H)*((XL-X2)/H)/DL3 
               END IF 
 
! -----------> JACOBIAN IS UPDATED 
               NL=ILS(2*IL-1) 
               ILE=ILS(2*IL) 
 
! -----------> FULL JACOBIAN MATRIX FJACL 
               DO K=1,NL 
                 KK=ILS(ILE-K+1) 
                 IK=IVE(KK) 
                 JK=IVC(KK) 
                 FJLK=FJACLAG(KK) 
                 IKI=IK+(I-1)*N 
! 
                  FJACL(IKI,JK)=FJACL(IKI,JK)+FJLK*DCOLI1 
                  FJACL(IKI,JK+N)=FJACL(IKI,JK+N)+FJLK*DCOLI2    
                  FJACL(IKI,JK+2*N)=FJACL(IKI,JK+2*N)+FJLK*DCOLI3 
               END DO 
             END IF 
       END DO            
! 
! --> NLAGS 
      END DO 
! <-- 
 
      AI11H   =-AI11/H 
      AI12H   =-AI12/H 
      AI13H   =-AI13/H 
      AI21H   =-AI21/H 
      AI22H   =-AI22/H 
      AI23H   =-AI23/H 
      AI31H   =-AI31/H 
      AI32H   =-AI32/H 
      AI33H   =-AI33/H 
 
! --- FJACL 
      IF (IMPLCT) THEN 
       DO I1=1,N 
         DO J1=1,N 
           FJACL(I1,J1)=FJACL(I1,J1)+AI11H*FMAS(I1,J1) 
           FJACL(I1,J1+N)=FJACL(I1,J1+N)+AI12H*FMAS(I1,J1) 
           FJACL(I1,J1+2*N)=FJACL(I1,J1+2*N)+AI13H*FMAS(I1,J1) 
 
           FJACL(I1+N,J1)=FJACL(I1+N,J1)+AI21H*FMAS(I1,J1) 
           FJACL(I1+N,J1+N)=FJACL(I1+N,J1+N)+AI22H*FMAS(I1,J1) 
           FJACL(I1+N,J1+2*N)=FJACL(I1+N,J1+2*N)+AI23H*FMAS(I1,J1) 
 
           FJACL(I1+2*N,J1)=FJACL(I1+2*N,J1)+AI31H*FMAS(I1,J1) 
           FJACL(I1+2*N,J1+N)=FJACL(I1+2*N,J1+N)+AI32H*FMAS(I1,J1) 
           FJACL(I1+2*N,J1+2*N)=FJACL(I1+2*N,J1+2*N)+AI33H*FMAS(I1,J1) 
         END DO 
       END DO 
      ELSE 
! ---  EXPLICIT CASE 
       DO I1=1,N 
         FJACL(I1,I1)=FJACL(I1,I1)+AI11H 
         FJACL(I1,I1+N)=FJACL(I1,I1+N)+AI12H 
         FJACL(I1,I1+2*N)=FJACL(I1,I1+2*N)+AI13H 
 
         FJACL(I1+N,I1)=FJACL(I1+N,I1)+AI21H 
         FJACL(I1+N,I1+N)=FJACL(I1+N,I1+N)+AI22H 
         FJACL(I1+N,I1+2*N)=FJACL(I1+N,I1+2*N)+AI23H 
 
         FJACL(I1+2*N,I1)=FJACL(I1+2*N,I1)+AI31H 
         FJACL(I1+2*N,I1+N)=FJACL(I1+2*N,I1+N)+AI32H 
         FJACL(I1+2*N,I1+2*N)=FJACL(I1+2*N,I1+2*N)+AI33H 
       END DO 
      END IF 
       
 
  
! ----- FACTORIZATION OF THE FULL JACOBIAN 
      CALL DEC(3*N,3*LDJAC,FJACL,IPJ,IER)  
      IF (IER.NE.0) THEN 
       GOTO 78 
      END IF 
! --->                                
      NDEC=NDEC+1 
  33  CONTINUE 
! --- 
! --- R: EVERY STEP STARTS WITH A SIMPLE ITERATION
      NSTEP=NSTEP+1 
      IF (NSTEP.GT.NMAX) GOTO 178 
      IF (0.1D0*H.LE.ABS(X)*UROUND) GOTO 177 
! --- 
      XPH=X+H 
 
! *** *** *** *** *** *** *** 
!  LOOP FOR NEWTON ITERATION 
! *** *** *** *** *** *** *** 
!              -----------      
            NEWT=0 
! ----------------------- 
            FLAGS=.FALSE. 
            FLAGN=.FALSE. 
! ----------------------- 
            FACCON=MAX(FACCON,UROUND)**0.8D0 
            THETA=ABS(THET) 
! --- --- --- --- --- --- --- --- --- --- --- --- --- 
!         REFERENCE POINT FOR FULL ITERATION          
! --- --- --- --- --- --- --- --- --- --- --- --- --- 
  43        CONTINUE 
! --- --- --- --- --- --- --- --- --- --- --- --- 
            IF (FLAGS) THEN 
             FLAGN=.TRUE. 
! ***************** 
! ---      DYNAMIC UPDATE OF THE CURRENT INTETPOLANT (in PAST)        
             DO J=1,NRDS 
                I=IPAST(J) 
                PAST(J+IACT)=Y(I)+Z3(I) 
                 Z2I=Z2(I) 
                 Z1I=Z1(I) 
                PAST(J+1*NRDS+IACT)=(Z2I-Z3(I))/C2M1 
                 AK=(Z1I-Z2I)/C1MC2 
                 ACONT3=Z1I/C1 
                 ACONT3=(AK-ACONT3)/C2 
                PAST(J+2*NRDS+IACT)=(AK-PAST(J+1*NRDS+IACT))/C1M1 
                PAST(J+3*NRDS+IACT)=PAST(J+2*NRDS+IACT)-ACONT3 
             ENDDO 
!            UPDATE 
             PAST(IACT)=X 
!            LEFT ENDPOINT OF CURRENT INTERVAL 
             PAST(IACT+IDIF-1)=H 
!            USED STEPSIZE                                                 
            END IF 
! --- --- --- --- --- --- --- --- --- --- --- --- 
            IF (NEWT.GE.NIT) THEN  
             INREJ=2 
             GOTO 431 
! ---------> UNEXPECTED STEP-REJECTION 
            END IF 
! ---     COMPUTE THE RIGHT-HAND SIDE 
            CONT(1:N)=Y(1:N)+Z1(1:N) 
!            WRITE(6,*) 'FCN1'
            CALL FCN(N,X+C1*H,CONT,F1,ARGLAG,PHI,RPAR,IPAR,PAST, &
                     IPAST,NRDS, &
        LRPAST) 
            CONT(1:N)=Y(1:N)+Z2(1:N) 
            CALL FCN(N,X+C2*H,CONT,F2,ARGLAG,PHI,RPAR,IPAR,PAST, &
                     IPAST,NRDS, &
        LRPAST) 
            CONT(1:N)=Y(1:N)+Z3(1:N) 
            CALL FCN(N,XPH,CONT,F3,ARGLAG,PHI,RPAR,IPAR,PAST, &
                     IPAST,NRDS, &
        LRPAST)
!            WRITE(6,*) 'FCN2'
            NFCN=NFCN+3 
! 
! --->      RHS COMPUTATION 
            IF (IMPLCT) THEN 
             ZL=0.D0 
             DO I1=1,N 
               DO J1=1,N 
                 ZL(I1)=ZL(I1)+AI11H*FMAS(I1,J1)*Z1(J1) 
                 ZL(I1)=ZL(I1)+AI12H*FMAS(I1,J1)*Z2(J1) 
                 ZL(I1)=ZL(I1)+AI13H*FMAS(I1,J1)*Z3(J1) 
 
                 ZL(I1+N)=ZL(I1+N)+AI21H*FMAS(I1,J1)*Z1(J1) 
                 ZL(I1+N)=ZL(I1+N)+AI22H*FMAS(I1,J1)*Z2(J1) 
                 ZL(I1+N)=ZL(I1+N)+AI23H*FMAS(I1,J1)*Z3(J1) 
 
                 ZL(I1+2*N)=ZL(I1+2*N)+AI31H*FMAS(I1,J1)*Z1(J1) 
                 ZL(I1+2*N)=ZL(I1+2*N)+AI32H*FMAS(I1,J1)*Z2(J1) 
                 ZL(I1+2*N)=ZL(I1+2*N)+AI33H*FMAS(I1,J1)*Z3(J1) 
               END DO 
             END DO 
            ELSE 
             ZL(1:N)=AI11H*Z1(1:N)+AI12H*Z2(1:N)+AI13H*Z3(1:N) 
             ZL(N+1:2*N)=AI21H*Z1(1:N)+AI22H*Z2(1:N)+AI23H*Z3(1:N) 
             ZL(2*N+1:3*N)=AI31H*Z1(1:N)+AI32H*Z2(1:N)+AI33H*Z3(1:N) 
            END IF 
! 
            ZL(1:N)=-ZL(1:N)-F1(1:N) 
            ZL(N+1:2*N)=-ZL(N+1:2*N)-F2(1:N) 
            ZL(2*N+1:3*N)=-ZL(2*N+1:3*N)-F3(1:N) 
               
! --------> SOLVE THE LINEAR SYSTEMS 
            CALL SOL(3*N,3*LDJAC,FJACL,ZL,IPJ) 
            NSOL=NSOL+1 
            NEWT=NEWT+1 
            DYNO=0.D0 
            DO I=1,N 
               DENOM=SCAL(I) 
               DYNO=DYNO+(ZL(I)/DENOM)**2+(ZL(I+N)/DENOM)**2 &
                +(ZL(I+2*N)/DENOM)**2 
            END DO 
            DYNO=DSQRT(DYNO/N3) 
! -------------------------------------------------------- 
! ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE 
! -------------------------------------------------------- 
              IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN 
                THQ=DYNO/DYNOLD 
                IF (NEWT.EQ.2) THEN 
                   THETA=THQ 
                ELSE 
                   THETA=SQRT(THQ*THQOLD) 
                END IF 
                THQOLD=THQ 
                INREJ=0 
! --- 1 
                IF (THETA.LT.0.99D0) THEN 
                    FACCON=THETA/(1.0D0-THETA) 
                    DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT 
                    IF (DYTH.GE.1.0D0) INREJ=1 
! ----------------  CONVERGENZA LENTA --- 
                ELSE 
                    INREJ=2 
! ----------------  DIVERGENZA --- 
                END IF 
! --- 1 
            ELSE 
                INREJ=0 
            END IF 
! ---------------------------------------- 
! ------------------ THE STEP IS REPEATED  
! ---------------------------------------- 
 431        IF (INREJ.GT.0) THEN 
             REJECT=.TRUE. 
             LAST=.FALSE.
             JLFLAG=0 
! -------------------------------------------- 
! ---        TURN TO A SIMPLE ITERATION         
             IF (INREJ.EQ.1) THEN 
                QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH)) 
                HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT)) 
                H=HHFAC*H 
             ELSE IF (INREJ.EQ.2) THEN 
                H=H*0.55D0 
                HHFAC=0.55D0 
             END IF 
 
             IF (CALJAC) GOTO 20 
             GOTO 10 
            END IF 
! -------------------------------------------------------- 
            DYNOLD=MAX(DYNO,UROUND) 
! --        UPDATE OF Z VALUES  
            Z1(1:N)=Z1(1:N)+ZL(1:N) 
            Z2(1:N)=Z2(1:N)+ZL(N+1:2*N) 
            Z3(1:N)=Z3(1:N)+ZL(2*N+1:3*N) 
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
            IF (FACCON*DYNO.GT.FNEWT) THEN 
             GOTO 43 
            END IF 
! - - - - - END FULL NEWTON ITERATION
! ----------------------------------------------------------------- 
! -----------------------------------------------------------------
!          RK EQUATIONS SUCCESFULLY SOLVED 
! -----------------------------------------------------------------
! ----------------------------------------------------------------- 
! *** *** *** *** *** *** *** *** *** *** *** 
! END LOOP 
! *** *** *** *** *** *** *** *** *** *** *** 
 
! -----------------------------------------------------------------
  55  CONTINUE 
! ******************** 
! --- ERROR ESTIMATION   
! ******************** 
! 
!      ERROR ESTIMATES                                        
       CALL ESTRAD (N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS, &
                    H,U1,DD1,DD2,DD3,CL1,CL3,CQ1,CQ2,CQ3,CERLQ, &
                    FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,E1,LDE1,ALPHA, & 
                    Z1,Z2,Z3,CONT,F1,F2,F3,IP1,IPHES,SCAL,SERR,CERR, &
                    FIRST,REJECT,FAC1,ARGLAG,PHI,RPAR,IPAR, &
                    IOUT,PAST,IPAST,NRDS,JLFLAG,IEFLAG,LRPAST) 
 
      FAC=MIN(SAFE,CFAC/(NEWT+2*NIT)) 
 
      IF (FIRST) THEN 
! ------------------------------------------------------------ 
! ---  AFTER A GRID OR BREAKING POINT            
! --------------------------------------------- 
       ERR=SERR
! ---  WE REQUIRE .2<=HNEW/H<=8. 
       QUOT=MAX(FACR,MIN(FACL,ERR**0.25D0/FAC)) 
       HNEW=H/QUOT 
      ELSE 
!   ---------------------------------------------------------- 
       IF (IEFLAG.EQ.-1) THEN 
        CERR2=H*CERR 
! ---    
        ERR=CERR2 
       ELSE IF (IEFLAG.EQ.1) THEN 
! -----    
        CERR2=H*CERR 
        ERR=CERS*SERR+CERC*CERR2 
       ELSE IF (IEFLAG.EQ.2) THEN 
! ----- STANDARD DISCRETE ERROR 
        ERR=SERR 
       ELSE IF (IEFLAG.EQ.3) THEN 
! ----- NOT AVAILABLE AT THE MOMENT 
        CERR2=H*CERR 
        ERR=CERS*SERR+CERC*CERR2 
       END IF 
! ------------------------------------ 
! ---  COMPUTATION OF HNEW 
! ---  AS PREVIOUSLY COMPUTED: 
! ------------------------------------ 
! ---  WE REQUIRE .2<=HNEW/H<=8. 
! --- 
! ---  LINEAR COMBINATION OF BOTH ERROR COMPONENTS
       QUOT=MAX(FACR,MIN(FACL,ERR**.25D0/FAC)) 
 
       HNEW=H/QUOT 
      END IF 
! --- 
! *** *** *** *** *** *** *** 
!  DOES THE ERROR INCREASE MUCH? 
! *** *** *** *** *** *** *** 
      REPEAT=.FALSE.
      IF (ERR.GE.1.D0.OR.((ERR/ERRACC.GE.TCKBP).AND.(.NOT.BPD))) THEN 
! ---   KIND OF STEPSIZE
        IF (BPD) THEN
! ---     BP IS WRONG! REPEAT
          IBP=IBP-1
          BPD=.FALSE.
          REPEAT=.TRUE.
          IF (ERR.LT.1.D0) HNEW=H*0.55D0
        ELSE
! ---     LOOK FOR A BP        
          HP=H*0.99D0
!          WRITE(6,*) 'Calling BPDTCT from row 2154.'
          CALL BPDTCT(N,X,HP,Y,ARGLAG,RPAR,IPAR,UCONT,GRID,NLAGS, &
                      FIRST,LAST,XEND,IGRID,BPV,IBP,ILBP,BPP,BPD, &
                      KMAX,PHI,PAST,IPAST,NRDS) 
          IF (BPD) REPEAT=.TRUE.
        END IF
      END IF 
      BPC=.FALSE. 
! ---
! *** *** *** *** *** *** *** 
!  IS THE ERROR SMALL ENOUGH ? 
! *** *** *** *** *** *** *** 
!     IF (ERR.LT.1.D0) THEN 
 56   IF (.NOT.REPEAT.AND.ERR.LT.1.D0) THEN
! -------------------- 
! --- STEP IS ACCEPTED   
! -------------------- 
         FLAGN=.FALSE. 
         CALJACL=.FALSE. 
! 
! ---    UPDATE NUMBER OF ACCEPTED STEPS              
         NACCPT=NACCPT+1 
         IF (PRED) THEN 
! --------> PREDICTIVE CONTROLLER OF GUSTAFSSON 
            IF (NACCPT.GT.1) THEN 
             IF (FLAGUS) THEN 
!               WRITE(6,*) 'HACC:', HACC
               FACGUS=(HACC/H)*(ERR**2/ERRACC)**0.25D0/SAFE 
               FACGUS=MAX(FACR,MIN(FACL,FACGUS)) 
               QUOT=MAX(QUOT,FACGUS) 
               HNEW=H/QUOT 
             ELSE 
               FLAGUS=.TRUE. 
             END IF 
            END IF 
            HACC=H 
         END IF 
         ERRACC=MAX(1.0D-2,ERR)
!        ERRACC=ERR 
         XOLD=X 
         HOLD=H 
         X=XPH  

! ---    AGGIORNAMENTO DELLA SOLUZIONE 
         DO I=1,N 
            Z3I=Z3(I) 
            YI=Y(I)+Z3I   
            Y(I)=YI 
! ------------------------------------------------------------------------- 
! --------- UPDATE THE STAGES AND SOLUTION 
! ------------------------------------------------------------------------- 
! --- 
            CONT(I)=YI 
! ---    
              Z2I=Z2(I) 
              Z1I=Z1(I) 
! ------------------------------------------------------------------- 
! ---         INVERSE DEVIDED DIFFERENCES
! ------------------------------------------------------------------- 
              A1=(Z2I-Z3I)/C2M1 
              CONT(I+N)=A1 
              AK=(Z1I-Z2I)/C1MC2 
              ACONT3=Z1I/C1 
              ACONT3=(AK-ACONT3)/C2 
              A2=(AK-CONT(I+N))/C1M1 
              CONT(I+N2)=A2 
            IF (.NOT.FIRST) THEN 
              CONT(I+N3)=A2-ACONT3 
            ELSE 
! ---         QUADRATIC APPROXIMATION
              CONT(I+N3)=0.D0           
! ---         INVECE DI:  
! ---         CONT(I+N3)=CONT(I+N2)-ACONT3 
            END IF 
         END DO 
! ---------------------------------------- 
! ---    SAVE LAST ACCEPTED STEP INFORMATION   
         DO I=1,LRC 
          UCONT(I)=CONT(I) 
         END DO 
         UCONT(LRC+1)=X 
         UCONT(LRC+2)=H 
! ---    FOR POSSIBLE SEARCH OF BREAKING POINTS     
! ------------------------------------------------- 
 
! ------------------------------------------------------------------ 
! ---    STEP IS ACCEPTED> DENSE OUTPUT IS STORED IN PAST        
! ---    
         DO J=1,NRDS 
            I=IPAST(J) 
            PAST(J+IACT)=CONT(I) 
            PAST(J+1*NRDS+IACT)=CONT(I+N) 
            PAST(J+2*NRDS+IACT)=CONT(I+N2) 
            IF (.NOT.FIRST) THEN 
              PAST(J+3*NRDS+IACT)=CONT(I+N3) 
            ELSE 
! ---       QUADRATIC APPROXIMATION                        
              PAST(J+3*NRDS+IACT)=0.D0 
            END IF 
         ENDDO 
! ---------> AGGIORNAMENTO DI PAST <--------- 
         PAST(IACT)=XOLD 
! ---    
         IACT=IACT+IDIF 
! ---    POINTER TO NEXT STEP 
         PAST(IACT-1)=H 
! ---    
         IF (IACT+IDIF-1.GT.MXST*IDIF) IACT=1  
! ---    CONTROL ON THE MEMORY DIMENSION
! -------------------------------------------------------------- 

! ----------------------------------------------------------------- 
! ---    COMPUTATION AT BP FOR NEXT STEP
         IF (LAST) THEN 
! ---     LAST HAS TO BE RE/DEFINED
          IF (BPD) THEN 
           WRITE(6,*) 'Found a BP at ', X, ', decrementing IGRID.'
           IGRID=IGRID-1
       END IF 
! --- 
          IF (.NOT.IMPLCT.OR.NEUTRAL) THEN  ! EXPLICIT PROBLEM
            HE=DMAX1(H/1.D4,10.D0*UROUND) 
! -------------------- 
            LEFT=.TRUE. 
! -------------------- 
! --- 
! ---       EULER STEP
!            WRITE(6,*) 'EULER STEP FCN6'
            CALL FCN(N,X,Y,F2,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS, &
       LRPAST)
!            WRITE(6,*) 'END EULER STEP FCN6'
            IF (NEUTRAL) THEN
             DO I=1,N-NDIMN 
               Z2(I)=Y(I)+HE*F2(I) 
              END DO
              DO I=1,NDIMN
              Z2(N-NDIMN+I)=F2(IPAST(NRDS+I))+HE*F2(N-NDIMN+I)
             END DO
            ELSE
              DO I=1,N 
               Z2(I)=Y(I)+HE*F2(I) 
              END DO
            END IF  
! -------------------- 
            LEFT=.FALSE. 
! -------------------- 
!            WRITE(6,*) 'SECOND EULER STEP FCN6'
            CALL FCN(N,X,Y,F3,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS, &
        LRPAST)
!            WRITE(6,*) 'END SECOND EULER STEP FCN6'
            IF (NEUTRAL) THEN
             DO I=1,N-NDIMN 
               Z3(I)=Y(I)+HE*F3(I) 
              END DO
              DO I=1,NDIMN
              Z3(N-NDIMN+I)=F3(IPAST(NRDS+I))+HE*F3(N-NDIMN+I)
             END DO
            ELSE
              DO I=1,N 
               Z3(I)=Y(I)+HE*F3(I) 
              END DO 
            END IF
!            WRITE(6,*) 'Calling arglag (2), ILBP = ', ILBP
            XL =ARGLAG(ILBP,X+HE,Z2,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST,N)  
!            WRITE(6,*) 'Calling arglag (3)'
            XLR=ARGLAG(ILBP,X+HE,Z3,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST,N)  
            IF (XL.GE.BPP.AND.XLR.GE.BPP) THEN 
             LEFT=.FALSE. 
            ELSE IF (XL.LT.BPP.AND.XLR.LT.BPP) THEN 
             LEFT=.TRUE. 
            ELSE 
             IF (IOUT.EQ.1) THEN 
              IF (XL.GT.BPP) THEN
                WRITE (6,*) &
              ' WARNING!: SOLUTION DOES NOT EXIST AT X= ',X
              ELSE
              WRITE (6,*) &
              ' WARNING!: SOLUTION IS  NOT UNIQUE AT X= ',X
             END IF
              GO TO 980
!             RETURN 
             END IF
            END IF
! ---       PROJECTION FOR DERIVATIVE COMPONENTS OF NEUTRAL EXPLICIT PROBLEMS
            PROJECT=.TRUE.
         IF (NEUTRAL.AND.PROJECT) THEN
            IF (LEFT) THEN
              DO J=1,NDIMN
              Y(N-NDIMN+J)=F2(IPAST(NRDS+J))
             END DO
            ELSE
             DO J=1,NDIMN
              Y(N-NDIMN+J)=F3(IPAST(NRDS+J))
             END DO
            END IF
            END IF
       ELSE ! GENERAL IMPLICIT
           LEFT=.TRUE.
       END IF
! --- 
          BPD=.FALSE. 
         END IF  
! ----------------------------------------------------------------- 
! 
! --- 
         IF (ITOL.EQ.0) THEN 
             DO I=1,N 
                SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I)) 
             END DO 
         ELSE 
             DO I=1,N 
                SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I)) 
             END DO 
         END IF 
         IF (IOUT.NE.0) THEN 
! --- 
             NRSOL=NACCPT+1 
             XSOL=X 
             XOSOL=XOLD 
             NSOLU=N 
             HSOL=HOLD 
! --- 
             CALL SOLOUT(NRSOL,XOSOL,XSOL,HSOL,Y,CONT,LRC,NSOLU, &
                         RPAR,IPAR,IRTRN)
             IF (IRTRN.LT.0) GOTO 179 
         END IF 
         CALJAC=.FALSE. 
! ----------------------------- 
! ---    FIRST IS RESTORED 
         FIRST=.FALSE. 
         IF (LAST) FIRST=.TRUE. 
! ---    COMPUTATION OF Y0 
!        WRITE(6,*) 'ERROR ELLER?' 
         CALL FCN(N,X,Y,Y0,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS, &
        LRPAST) 
!       
         NFCN=NFCN+1 
         FIRST=.FALSE.  
! ----------------------------- 
      
! ------ FINAL POINT                               
!         WRITE(6, *) 'X:', X
!         WRITE(6, *) 'IGRID:', IGRID
         IF (LAST) THEN 
 45         CONTINUE
            IF (IGRID.EQ.NGRID) THEN 
!               WRITE(6,*) 'HOPT:', HOPT
               H=HOPT 
               IDID=1 
! ---          END OF COMPUTATION 
               GOTO 980 
            ELSE 
               IGRID=IGRID+1 
               LAST=.FALSE. 
!              
!              LEFT=.FALSE. 
               FIRST=.TRUE. 
               XEND=GRID(IGRID)
            IF (ABS(XEND-X).LE.(H*1.D-2)) THEN
            IGRID=IGRID+1
            GO TO 45 
            END IF 
               FLAGUS=.FALSE. 
               IF (WORK7.EQ.0.D0) THEN 
                  HMAXN=XEND-X  
                  HMAX=HMAXN 
               END IF 
               HNEW=MIN(HNEW,H)
            END IF 
         END IF 
!         WRITE(6,*) 'Before updating hnew and hopt'
!         WRITE(6, *) 'IGRID:', IGRID
         HNEW=MIN(HNEW,HMAXN) 
         HOPT=MIN(H,HNEW) 
         IF (REJECT) HNEW=MIN(HNEW,H)  
         REJECT=.FALSE. 
         IF ((X+HNEW/QUOT1-XEND).GE.0.D0) THEN 
            H=XEND-X 
            IF (H.LT.0.D0) THEN 
             WRITE(6,*) 'ERROR!: NEGATIVE STEPSIZE! AT ' 
             WRITE(6,*) 'X > XEND = ',X,XEND 
             STOP 
            END IF 
            LAST=.TRUE. 
         ELSE 
! ----------------------------------------------------- 
! --------  IN ORDER TO AVOID VERY SMALL END-STEPSIZES:    
            IF ((X+1.8D0*HNEW-XEND).GT.0.D0) THEN  
              H=(XEND-X)*0.55D0 
            ELSE 
              QT=HNEW/H  
              HHFAC=H 
! ------------------ 
              IMANT=0 
! ---         STEP IS MAINTAINED     
              IF (THETA.LE.THET.AND.QT.GE.QUOT1.AND.QT.LE.QUOT2) THEN 
               IF ((JLFLAG.EQ.0).AND.(.NOT.FIRST)) IMANT=1 
              ELSE  
               H=HNEW  
              END IF 
            END IF 
         END IF 
         HHFAC=H 
! ------------------------------------------- 
! ---    SIMPLE ITERATION FIRST                 
         JLFLAG=0 
! ---    
         IF ((.NOT.FIRST).AND.(THETA.LE.THET)) GOTO 20 
         GOTO 10 
! -------------------- 
! --- STEP IS ACCEPTED   
! -------------------- 
      ELSE 
! -------------------- 
! --- STEP IS REJECTED   
! -------------------- 
! --- 
         IF (BPD) THEN  
          IF (.NOT.FIRST) THEN 
           HNEW=HP 
           LAST=.TRUE. 
          END IF 
         ELSE 
          LAST=.FALSE. 
         END IF 
  
         FLAGN=.FALSE. 
! ------------------------------- 
         IF (IRTRN.LT.0) GOTO 179 
 
         REJECT=.TRUE. 
! --- 
         IF (FIRST) THEN 
             H=H*0.12D0 
             HHFAC=0.12D0 
         ELSE  
             HHFAC=HNEW/H 
             H=HNEW  
         END IF 
         IF (NACCPT.GE.1) NREJCT=NREJCT+1 
! ---> 
         JLFLAG=0  
         IF (CALJAC) GOTO 20 
         GOTO 10 
! -------------------- 
! --- STEP IS REJECTED   
! -------------------- 
      END IF 
! ---------------------------------- 
! --- END OF ERROR CONTROL PROCEDURE 
! ---------------------------------- 
! ----------------------------- 
! --- UNEXPECTED STEP-REJECTION 
  78  CONTINUE 
 
      REJECT=.TRUE. 
      FLAGN=.FALSE. 
! --- 
      LAST=.FALSE. 
! --- 
      IF (IER.NE.0) THEN 
          NSING=NSING+1 
          IF (NSING.GE.5) GOTO 176 
      END IF 
! --- 
      H=H*0.55D0  
      HHFAC=0.55D0 
! --- 
      IF (CALJAC) GOTO 20 
      IF (IRTRN.LT.0) GOTO 175 
      GOTO 10 
! --- FAIL EXIT 
 175  CONTINUE 
      IDID=-5 
      GOTO 980 
 176  CONTINUE 
      WRITE(6,979)X    
      WRITE(6,*) ' MATRIX IS REPEATEDLY SINGULAR, IER=',IER 
      IDID=-4 
      GOTO 980 
 177  CONTINUE 
      WRITE(6,979)X    
      WRITE(6,*) ' STEP SIZE TOO SMALL, H=',H 
      IDID=-3 
      GOTO 980 
 178  CONTINUE 
      WRITE(6,979)X    
      WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED'  
      IDID=-2 
      GOTO 980 
! --- EXIT CAUSED BY SOLOUT 
 179  CONTINUE 
      WRITE(6,979)X 
 979  FORMAT(' EXIT OF RADAR5 AT X=',E18.4)  
      IDID=2 
 
! --- RETURN LABEL 
 980  CONTINUE 
!      WRITE(6,*) IBP-1,' COMPUTED BREAKING POINTS: '
!      WRITE(8,*) 'BREAKING POINTS: '
!      DO I=1,IBP
!      WRITE(8,*) BPV(I)
!      END DO
!      WRITE(8,*) ' -------------- '
!      CLOSE(8)

! --- DEALLOCATION OF THE MEMORY 
      DEALLOCATE (Z1,Z2,Z3,Y0,SCAL,F1,F2,F3) 
      DEALLOCATE(BPV) 
      DEALLOCATE(FJAC,ZL) 
      IF (IMPLCT) DEALLOCATE(FMAS) 
      DEALLOCATE (IP1,IP2,IPHES) 
      DEALLOCATE (E1,E2R,E2I) 
!      DEALLOCATE (PAST) 
      IF (NLAGS.GT.0) THEN 
       DEALLOCATE (FJACS,FJACLAG) 
       DEALLOCATE (IVL,IVE,IVC,ILS,ICOUN) 
       IF (ISWJL.NE.1) DEALLOCATE (IPJ,FJACL,XLAG) 
      END IF 
      DEALLOCATE (CONT,UCONT) 
      RETURN 
      END 
! 
!     END OF SUBROUTINE RADCOR 
! 
 
 
 
! *********************************************************** 
! 
      SUBROUTINE LAGR5(IL,X,Y,ARGLAG,PAST,THETA,IPOS,RPAR,IPAR, &
                       PHI,IPAST,NRDS, &
                       LRPAST,N)
! ---------------------------------------------------------- 
!     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION 
!     WITH THE OUTPUT-SUBROUTINE FOR RADAR5. IT PROVIDES THE 
!     POSITION OF THE DENSE OUTPUT AT THE IL-TH DELAY. 
! ---------------------------------------------------------- 
      USE IP_ARRAY 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
      REAL(kind=DP), dimension(N), intent(in) :: Y
      REAL(kind=DP), dimension(LRPAST), intent(in) :: PAST 
      REAL(kind=DP), dimension(1), intent(in) :: RPAR 
      INTEGER, dimension(1), intent(in) :: IPAR 
      INTEGER, dimension(NRDS), intent(in) :: IPAST 
      EXTERNAL PHI
       INTEGER, intent(in) :: LRPAST 
       INTEGER :: N
       
      LOGICAL FLAGS,FLAGN,FIRST,LAST,REJECT,BPD,LEFT 
! --- COMMON BLOCKS 
      COMMON /BPLOG/FIRST,LAST,REJECT,BPD 
      COMMON /BPCOM/BPP,ILBP,LEFT 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
! 
!      WRITE(6,*) 'iposv(il) info. IL, IPOSV(IL), IPOS', 
!     & IL, IPOSV(IL), IPOS
! --- COMPUTE DEVIATED ARGUMENT FOR IL-TH DELAY 
!      WRITE(6,*) 'Calling arglag (4)'
      XLAG=ARGLAG(IL,X,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST,N)
!      IF (IL .EQ. 10) THEN
!        WRITE(6, *) 'XLAG = ', XLAG
!      END IF
!      WRITE (6,*) ' LRPAST,NRDS,N = ',X,XLAG,X0B, IPOSV(1)
! --- INITIAL PHASE 
      THETA=XLAG 
      IPOS=-1 

!     EPSILON GIVES THE OVERLAPPING
!     MIN VALUE FOR THE SUPER-POSITION NEIGHBOURHOOD OF A BP

      COMPAR=UROUND*MAX(ABS(XLAG),ABS(X0B))
      EPSACT=10.D0*COMPAR
      IF (IACT.GT.1) EPSACT=DMAX1(PAST(IACT-1)*1.D-2,EPSACT)

      IF (XLAG.LE.X0B) THEN 
! ---     DEVIATING ARGUMENT ON THE INITIAL SEGMENT         
        IF (.NOT.((IL.EQ.ILBP).AND.(BPD.OR.FIRST))) THEN 
          IF (XLAG-X0B.LT.0.D0) THEN
            RETURN
          ELSE
            IPOS=1
            THETA=-1.D0 
          END IF
        ELSE
          IF (ABS(XLAG-X0B).LE.EPSACT) THEN
            IF (LEFT) THEN
              IPOS=-1 
              THETA=XLAG
            ELSE
              IPOS=1
              THETA=(XLAG-(PAST(IPOS)+PAST(IPOS+IDIF-1)))/ &
              PAST(IPOS+IDIF-1) 
            END IF
          ELSE IF (ABS(XLAG-BPP).LE.EPSACT) THEN 
            IPOS=-1 
            IF (LEFT) THEN 
              IF (XLAG.GT.BPP) THEN 
                IF (BPP.GT.0.D0) THEN 
                  THETA=BPP*(1.D0-100*UROUND) 
                ELSE 
                  THETA=BPP*(1.D0+100*UROUND) 
                END IF 
              END IF 
            ELSE 
              IF (XLAG.LT.BPP) THEN 
                IF (BPP.GT.0) THEN 
                  THETA=BPP*(1.D0+100*UROUND) 
                ELSE 
                  THETA=BPP*(1.D0-100*UROUND) 
                END IF 
              END IF 
            END IF 
          END IF
        END IF
        
        RETURN
      END IF 

! --- COMPUTE THE POSITION OF XLAG 
      IPA = IACT+IDIF 
      IF (IPA.GT.(MXST-1)*IDIF+1) IPA=1 
      IF (XLAG-PAST(IPA).LT.0.D0) THEN
         WRITE (6,*) ' MEMORY FULL, MXST = ',MXST 
         IRTRN=-1
        STOP 
!        RETURN 
      END IF 

      INEXT=IACT-IDIF 
      IF (INEXT.LT.1) INEXT=(MXST-1)*IDIF+1 
      XRIGHT=PAST(INEXT)+PAST(INEXT+IDIF-1) 

! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
! --- INITIALIZE POSITION INSIDE THE MEMORY 
      IPOS=IPOSV(IL)
!      IPOS = 1
!      WRITE(6,*) 'IPOS=', IPOS
! --- STEP AND DELAYS 
      IF (XLAG-XRIGHT.GT.0.D0) THEN 
        IF (.NOT.FLAGN) THEN 
           IPOS=IACT-IDIF 
           IF (IPOS.LT.1) IPOS=(MXST-1)*IDIF+1 
        ELSE 
           IPOS=IACT 
        END IF 
        FLAGS=.TRUE. 
! ------------------------ 
! ----- COMPUTE THETA (>0)     
! ------------------------ 
        H=PAST(IPOS+IDIF-1) 
        THETA=(XLAG-(PAST(IPOS)+H))/H 
! -----    EXTRAPOLATION USE OF THE COLLOCATION POLYNOMIAL            
! ----- 
! -----     
        RETURN                       
      ELSE 
   1    CONTINUE 
        IF (XLAG-PAST(IPOS).LE.0.D0) THEN 
           IPOS=IPOS-IDIF 
           IF (IPOS.LT.1) IPOS=(MXST-1)*IDIF+1 
           GOTO 1 
        END IF 
   2    CONTINUE 
        INEXT=IPOS+IDIF 
        IF (INEXT.GT.(MXST-1)*IDIF+1) INEXT=1 
        IF (XLAG.GT.PAST(INEXT).AND.INEXT.NE.IACT) THEN 
           IPOS=INEXT 
           GOTO 2 
        END IF
       IF (IPOS.EQ.1) THEN 
            IPREV=(MXST-1)*IDIF+1 
        ELSE 
            IPREV=IPOS-IDIF 
        END IF
! --- 
! ---   IN CORRESPONDENCE OF BREAKING POINTS                
! --- 
 3      CONTINUE 
        IF (.NOT.((IL.EQ.ILBP).AND.(BPD.OR.FIRST))) GOTO 10
 
        IF (BPP.EQ.X0B) THEN
         IF (LEFT) THEN
           IPOS=-1
           THETA=XLAG
           RETURN
          ELSE
           IF (IPOS.EQ.-1) IPOS=1
            GO TO 10
           END IF
        END IF

       IPOSB=0  
       IF (ABS(BPP-PAST(IPOS)).LE.10.D0*UROUND) THEN
           IPOSB=IPOS
       ELSE IF (ABS(BPP-PAST(INEXT)).LE.10.D0*UROUND) THEN
          IPOSB=INEXT
!      ELSE IF (ABS(BPP-PAST(IPREV)).LE.10.D0*UROUND) THEN
!         IPOSB=IPREV
       END IF

        IF (IPOSB.EQ.0) THEN
        GO TO 10
       END IF

        IF (IPOSB.EQ.1) THEN
         EPSILON=(PAST(IPOSB+IDIF)-PAST(IPOSB))
       ELSE IF (IPOSB.EQ.(MXST-1)*IDIF+1) THEN
         EPSILON=(PAST(IPOSB)-PAST(IPOSB-IDIF))
       ELSE
         EPSILON=DMIN1(PAST(IPOSB+IDIF)-PAST(IPOSB), &
                       PAST(IPOSB)-PAST(IPOSB-IDIF))
       END IF
       EPSILON=DMAX1(EPSILON*1.D-2,EPSACT)
        
        IF (ABS(XLAG-BPP).GT.EPSILON) GOTO 10 

        IF (IPOSB.EQ.1) THEN
         IF (LEFT) THEN 
                IPOS=-1 
!              IF (BPP.GT.0.D0) THEN 
!                  THETA=BPP*(1.D0-100*UROUND) 
!               ELSE 
!                  THETA=BPP*(1.D0+100*UROUND) 
!               END IF
!               SE PERO' IL DATO INIZIALE E' ESTENDIBILE
                THETA=XLAG
                RETURN
         ELSE 
                IPOS=1 
             GO TO 10
         END IF
        END IF
      
        IF (LEFT) THEN           
! ---     PREVIOUS INTERVAL HAS TO BE SELECTED
          IPOS=IPOSB-IDIF
        ELSE 
! ---     NEXT INTERVAL HAS TO BE SELECTED
          IPOS=IPOSB
        END IF 
      
! --------------------------------------------------------- 

! ----- COMPUTE THETA (<0): SITUAZIONE PIU' TIPICA     
  10    THETA=(XLAG-(PAST(IPOS)+PAST(IPOS+IDIF-1)))/PAST(IPOS+IDIF-1)
!        WRITE(6,*) 'SOME STUFF', XLAG, PAST(IPOS), PAST(IPOS+IDIF-1),
!     &  PAST(IPOS+IDIF-1)
!        WRITE(6,*) IPOS, IDIF
! ----- REM: THETA IS NEGATIVE                                  
      END IF 
! --- UPDATE POSITION INSIDE THE MEMORY 
!      WRITE(6,*) 'Setting iposv(il) to IPOS. IL, IPOSV(IL), IPOS', 
!     & IL, IPOSV(IL), IPOS
      IPOSV(IL)=IPOS

!      WRITE(6,*) 'IN LAGR5, THETA, IPOS', THETA, IPOS
      RETURN 
      END 
! 
!     END OF FUNCTION LAGR5 
! 
! -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

! *********************************************************** 
! 
      FUNCTION YLAGR5(IC,THETA,IPOS,PHI,RPAR,IPAR,PAST,IPAST,NRDS, &
                      LRPAST) 
! ---------------------------------------------------------- 
!     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION 
!     WITH THE SUBROUTINE LAGR5. IT PROVIDES AN APPROXIMATION  
!     TO THE IC-TH COMPONENT OF THE SOLUTION AT XLAG. 
! ---------------------------------------------------------- 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
! --- REAL(kind=DP) PHI,YLAGR5
      REAL(kind=DP), dimension(1), intent(in) :: RPAR 
      INTEGER, dimension(1), intent(in) :: IPAR 
! --- 
      REAL(kind=DP), dimension(LRPAST), intent(in) :: PAST 
      INTEGER, dimension(NRDS), intent(in) :: IPAST 
      
      INTEGER :: LRPAST
      
      LOGICAL FLAGS,FLAGN 
!---- COMMON BLOCKS 
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
 
! 
! --- INITIAL PHASE 
      IF (IPOS.EQ.-1) THEN 
            YLAGR5=PHI(IC,THETA,RPAR,IPAR)   
! ---       CALL PHI(IC,THETA,YLAGR5,RPAR,IPAR)
            RETURN 
      END IF 
! --- 
! --- COMPUTE PLACE OF IC-TH COMPONENT  
      I=0  
      DO 5 J=1,NRDS  
      IF (IPAST(J).EQ.IC) I=J 
   5  CONTINUE 
      IF (I.EQ.0) THEN 
         WRITE (6,*) ' NO DENSE OUTPUT AVAILABLE FOR COMP.',IC 
         RETURN 
      END IF   
! ----- COMPUTE DESIRED APPROXIMATION 
      I=I+IPOS 
      YLAGR5=PAST(I)+THETA*(PAST(NRDS+I)+(THETA-C2M1)*(PAST(2*NRDS+I) &
                  +(THETA-C1M1)*(PAST(3*NRDS+I)))) 
      RETURN 
      END 
! 
!     END OF FUNCTION YLAGR5 
! 
 
! *********************************************************** 
! 
      FUNCTION DLAGR5(IC,THETA,IPOS,PHI,RPAR,IPAR,PAST,IPAST,NRDS, &
                      LRPAST) 
! ---------------------------------------------------------- 
!     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION 
!     WITH THE SUBROUTINE LAGR5. IT PROVIDES AN APPROXIMATION  
!     TO THE IC-TH COMPONENT OF THE SOLUTION DERIVATIVE AT XLAG. 
! ---------------------------------------------------------- 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
! --- REAL(kind=DP) PHI,DLAGR5
      REAL(kind=DP), dimension(1), intent(in) :: RPAR 
      INTEGER, dimension(1), intent(in) :: IPAR 
! 
      REAL(kind=DP), dimension(LRPAST), intent(in) :: PAST 
      INTEGER, dimension(NRDS), intent(in) :: IPAST 
      
       INTEGER :: LRPAST
      
      LOGICAL FLAGS,FLAGN 
!---- COMMON BLOCKS 
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
 
! 
! --- INITIAL PHASE 
      IF (IPOS.EQ.-1) THEN 
!           DLAGR5=DPHI(IC,THETA,RPAR,IPAR)   
            DLAGR5=0.D0 
            RETURN 
      END IF 
! --- 
! --- COMPUTE PLACE OF IC-TH COMPONENT  
      I=0  
      DO 5 J=1,NRDS  
      IF (IPAST(J).EQ.IC) I=J 
   5  CONTINUE 
      IF (I.EQ.0) THEN 
         WRITE (6,*) ' NO DENSE OUTPUT AVAILABLE FOR COMP.',IC 
         RETURN 
      END IF   
! ----- COMPUTE DESIRED APPROXIMATION 
      H=PAST(IPOS+IDIF-1) 
      I=I+IPOS 
      DLAGR5=PAST(NRDS+I)+(THETA-C2M1)*(PAST(2*NRDS+I) &
           +(THETA-C1M1)*PAST(3*NRDS+I)) &
           +THETA*(PAST(2*NRDS+I)+(2.D0*THETA-C2M1-C1M1) &
           *PAST(3*NRDS+I)) 
      DLAGR5=DLAGR5/H 
      RETURN 
      END 
! 
!     END OF FUNCTION DLAGR5 
! 
 
! 
      SUBROUTINE BPDTCT(N,X,H,Y,ARGLAG,RPAR,IPAR,UCONT,GRID,NLAGS, &
                        FIRST,LAST,XEND,IGRID,BPV,IBP,ILBP,BPP,BPD, &
                        KMAX,PHI,PAST,IPAST,NRDS) 
                       
! ---------------------------------------------------------- 
!     THIS SUBROUTINE CAN BE USED FOR DETECTING BREAKING POINTS  
!     WITH THE OUTPUT-SUBROUTINE FOR RADAR5. IT PROVIDES AN 
!     APPROXIMATION TO THE IC-TH COMPONENT OF THE SOLUTION AT X. 
! ---------------------------------------------------------- 
        IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
        INTEGER, PARAMETER :: DP=kind(1D0) 
        REAL(kind=DP), dimension(N) :: Y 
        REAL(kind=DP), dimension(1), intent(in) :: PAST 
        INTEGER, dimension(1), intent(in) :: IPAST 
        REAL(kind=DP), dimension(:), allocatable :: YADV 
        REAL(kind=DP), dimension(1) :: UCONT 
        REAL(kind=DP), dimension(1), intent(inout) :: GRID 
        REAL(kind=DP), dimension(1) :: BPV 
        INTEGER, dimension(NRDS) :: IPAR 
        REAL(kind=DP), dimension(1) :: RPAR 
        EXTERNAL PHI
        LOGICAL FIRST,LAST,BPD,FLAGS,FLAGN 
!----   COMMON BLOCKS 
        COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
        COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
        
        INTEGER LRPAST 
! 
! ---   
!        WRITE(6,*) 'In BPDTCT, ILBP = ', ILBP
        LRPAST = MXST*IDIF
        IF (FIRST) RETURN
        BPD=.FALSE.  
        LRC=4*N 
        EPSILON=1.D-10
        ALLOCATE(YADV(N)) 
        COMPAR=UROUND*MAX(ABS(X),ABS(X+H)) 
!        WRITE(6,*) 'COMPAR:', COMPAR, X, X+H, UROUND
        
      XLAST=UCONT(LRC+1) 
        HLAST=UCONT(LRC+2) 
        DO IL=1,NLAGS 
         ALS = ARGLAG(IL,X,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST,N)  
! -----  DEVIATING ARGUMENT AT X 
! -----  EXTRAPOLATION OF THE COLLOCATION POLYNOMIAL       
         DO IC=1,N 
          YADV(IC)=CONTR5(IC,N,X+H,UCONT,XLAST,HLAST) 
         END DO 
         ALD = ARGLAG(IL,X+H,YADV,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST,N) 
! -----  DEVIATING ARGUMENT AT X+H        
 
         IF (ABS(ALS-ALD).LE.COMPAR) GO TO 33
         DO L=1,IGRID-1 
          BPP=GRID(L) 
          IF ((ALS-BPP)*(ALD-BPP).LT.COMPAR) THEN 
           BPD=.TRUE. 
!          BREAKING POINT! 
           GO TO 33 
          END IF 
         END DO    
         DO L=IBP,1,-1 
          BPP=BPV(L) 
          IF ((ALS-BPP)*(ALD-BPP).LT.COMPAR) THEN 
!          BREAKING POINT! 
           BPD=.TRUE. 
           GO TO 33 
          END IF 
         END DO    
 33      CONTINUE 
         IF (BPD) THEN
! --- 
           THLIM=1.D0 
           THRIGH=THLIM 
! -------------------------------------------------- 
           THLEFT=0.D0
! ---      
           DO K=1,KMAX   
            THNEW = THLEFT - (ALS-BPP)*(THRIGH-THLEFT)/(ALD-ALS)  
! ---       TEST DI CONVERGENZA 
            IF (ABS(THRIGH-THNEW).LE.EPSILON.OR. &
                ABS(THLEFT-THNEW).LE.EPSILON) GOTO 36 
            XA=X+THNEW*H 
            DO IC=1,N 
             YADV(IC)=CONTR5(IC,N,XA,UCONT,XLAST,HLAST) 
            END DO
!            WRITE(6,*) 'Calling arglag (5)' 
            ALN = ARGLAG(IL,XA,YADV,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST,N)  
            IF ((ALS-BPP)*(ALN-BPP).LE.0.D0) THEN 
             ALD=ALN 
             THRIGH=THNEW 
            ELSE 
             ALS=ALN 
             THLEFT=THNEW 
            END IF 
! --- 
           END DO 
 36        CONTINUE 
! ---      BP FOUND!      
           IF ((THNEW.GT.COMPAR).AND.(THNEW.LT.THLIM)) THEN 
            HP=THNEW*H 
! ---
            IF (HP.LE.1.D2*COMPAR) THEN 
             BPD=.FALSE. 
             GO TO 37 
            END IF 
! ---
            IBP=IBP+1 
            BPV(IBP)=X+HP 
            H=HP 
            ILBP=IL 
            GOTO 37 
           ELSE 
! ---       BP ALREADY PRESENT
            BPD=.FALSE. 

           END IF  
         END IF 
        END DO 
 
 37     CONTINUE 
 
        DEALLOCATE(YADV) 
!        WRITE(6,*) 'End of BPDTCT, ILBP = ', ILBP
        RETURN 
        END 
 

      SUBROUTINE BPACC(N,X,H,Y,ARGLAG,RPAR,IPAR,Z1,Z2,Z3,FIRST, &
                       BPV,IBP,ILBP,BPP,KMAX,PHI,PAST,IPAST,NRDS) 
                       
! ---------------------------------------------------------- 
!     THIS SUBROUTINE CAN BE USED FOR APPROXIMATING BREAKING POINTS  
!     IN TANDEM WITH THE SIMPLIFIED NEWTON ITERATION.. 
! ---------------------------------------------------------- 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
      REAL(kind=DP), dimension(N) :: Y,Z1,Z2,Z3 
      REAL(kind=DP), dimension(1), intent(in) :: PAST
      INTEGER, dimension(NRDS), intent(in) :: IPAST
      REAL(kind=DP), dimension(:), allocatable :: YCONT,YAPP
      REAL(kind=DP), dimension(1) :: BPV 
      INTEGER, dimension(1) :: IPAR 
      REAL(kind=DP), dimension(1) :: RPAR 
      EXTERNAL PHI
      LOGICAL FIRST
      LOGICAL FLAGS,FLAGN
!---- COMMON BLOCKS 
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
        
      INTEGER LRPAST
! 
! --- 
        LRPAST = MXST*IDIF
        ALLOCATE(YCONT(4*N),YAPP(N))
        EPSILON=UROUND*1.D3
! ---   DYNAMIC UPDATE
        DO I=1,N 
            Z3I=Z3(I) 
            YI=Y(I)+Z3I   
            YCONT(I)=YI 
            Z2I=Z2(I) 
            Z1I=Z1(I) 
            A1=(Z2I-Z3I)/C2M1 
            YCONT(I+N)=A1 
            AK=(Z1I-Z2I)/C1MC2 
            ACONT3=Z1I/C1 
            ACONT3=(AK-ACONT3)/C2 
            A2=(AK-YCONT(I+N))/C1M1 
            YCONT(I+2*N)=A2 
            IF (.NOT.FIRST) THEN 
              YCONT(I+3*N)=A2-ACONT3 
            ELSE 
              YCONT(I+3*N)=0.D0 
            END IF 
        END DO 

! ---   INITIAL VALUES FOR THE COMPUTATION
        XSOL=X+H
        HSOL=H
        THLEFT=0.9D0
        THRIGH=1.0D0
        XL=X+THLEFT*H
        DO I=1,N 
            YAPP(I)=CONTR5(I,N,XL,YCONT,XSOL,HSOL) 
        END DO
!        WRITE(6,*) 'Calling arglag (6)' 
        ALS = ARGLAG(ILBP,XL,YAPP,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST,N) 
! ---
        XR=X+THRIGH*H
!        WRITE(6,*) 'Calling arglag (7)'
        ALD = ARGLAG(ILBP,XR,YCONT,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST,N) 
        DO K=1,KMAX 
            THNEW = THRIGH - (ALD-BPP)*(THRIGH-THLEFT)/(ALD-ALS)  
            THLEFT= THRIGH
! ---       TEST DI CONVERGENZA 
            IF (ABS(THNEW-THRIGH).LE.EPSILON) GOTO 36 
            IF ((THNEW.LE.0.5D0).OR.(THNEW.GE.1.5D0)) THEN 
              DEALLOCATE(YCONT,YAPP)
              RETURN
            END IF
            THRIGH=THNEW
            XAP=X+THRIGH*H 
            ALS=ALD
            DO I=1,N 
             YAPP(I)=CONTR5(I,N,XAP,YCONT,XSOL,HSOL) 
            END DO
!            WRITE(6,*) 'Calling arglag (8)' 
            ALD = ARGLAG(ILBP,XAP,YAPP,RPAR,IPAR,PHI,PAST,IPAST,NRDS, &
        LRPAST,N)  
            IF (ABS(ALD-ALS).LE.EPSILON) GOTO 36 
        END DO 

 36     CONTINUE
! ---   BP FOUND 
        H=MIN(THRIGH,THLEFT)*H
        BPV(IBP)=X+H 

        DEALLOCATE(YCONT,YAPP)
 
        RETURN 
        END 

