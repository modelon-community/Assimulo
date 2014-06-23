      SUBROUTINE METAN1 (N,FCN,T,Y,TEND,TOL,HMAX,H,KFLAG)
C
C* Begin Prologue METAN1
C
C  ---------------------------------------------------------------------
C
C* Title
C
C    Integrator for stiff systems of autonomous ordinary differential
C    equations.
C
C* Written by        P. Deuflhard, U. Nowak, U. Poehle
C* Purpose           Solution of systems of initial value problems
C* Method            Semi-implicit mid-point rule with
C                    h**2-extrapolation
C* Category          i1a2a. - System of stiff first order differential
C                             equations
C* Keywords          extrapolation, ODE, mid-point rule, stiff
C* Version           1.1 , July 1989
C* Latest Change     February 1991
C* Library           CodeLib
C* Code              Fortran 77
C                    Double Precision
C* Environment       Standard version for FORTRAN77 environments on
C                    PCs, workstations, and hosts
C* Copyright     (c) Konrad-Zuse-Zentrum fuer Informationstechnik
C                    Berlin (ZIB)
C                    Takustrasse 7, D-14195 Berlin-Dahlem
C                    phone : + 49/30/84185-0
C                    fax   : + 49/30/84185-125
C* Contact           Uwe Poehle
C                    ZIB, Scientific Software Group
C                    phone : + 49/30/84185-241
C                    fax   : + 49/30/84185-107
C                    e-mail: poehle@zib.de
C
C  ---------------------------------------------------------------------
C
C* Licence
C  -------
C
C  You may use or modify this code for your own non-commercial
C  purposes for an unlimited time. 
C  In any case you should not deliver this code without a special 
C  permission of ZIB.
C  In case you intend to use the code commercially, we oblige you
C  to sign an according licence agreement with ZIB.
C
C
C* Warranty
C  --------
C 
C  This code has been tested up to a certain level. Defects and
C  weaknesses, which may be included in the code, do not establish
C  any warranties by ZIB. ZIB does not take over any liabilities
C  which may follow from aquisition or application of this code.
C
C
C* Software status 
C  ---------------
C
C  This code is under care of ZIB and belongs to ZIB software
C  class I.
C
C
C  ---------------------------------------------------------------------
C
C* Short doc:
C
C  Jacobian approximation by numerical differences
C  or user supplied subroutine M1JAC (see constant QJACUS)
C
C  The numerical solution of the arising linear equations is done by
C  means of the subroutines DGEFA and DGESL.  For special purposes these
C  routines may be substituted.
C
C
C* References:
C
C /1/ P. Deuflhard:
C     A Semi-Implicit Midpoint Rule for Stiff Systems of Ordinary
C     Differential Equations
C     Num. Math. 41, 373 - 398 (1983) .
C
C /2/ P. Deuflhard:
C     Order and Stepsize Control in Extrapolation Methods
C     Numer. Math. 41, 399-422 (1983)
C
C /3/ P. Deuflhard:
C     Uniqueness Theorems for Stiff ODE Initial Value Problems
C     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin,
C     Preprint SC-87-3 (1987)
C
C
C* External subroutine: (to be supplied by the user)
C
C    FCN           EXT  Subroutine FCN(N,T,Y,DY,IFAIL)
C                       Right-hand side of first-order
C                       differential equations
C                       N      Number of first-order ODE's
C                       T      Actual position
C                       Y(N)   Values at T
C                       DY(N)  Derivatives at T
C                       IFAIL  Error return code
C
C
C* Parameters: (* marks transient parameters)
C
C    N         I   IN   Number of ODE'S
C  * T         D   IN   Starting point of integration
C                       (T .LE. TEND)
C                  OUT  Achieved final point of integration
C  * Y         D   IN   Array of initial values Y(1),...,Y(N)
C                  OUT  Array of final values
C    TEND      D   IN   Prescribed final point of integration
C    TOL       D   IN   Prescribed relative precision (.GT.0)
C    HMAX      D   IN   Maximum permitted stepsize
C  * H         D   IN   Initial stepsize guess
C                  OUT  Stepsize proposal for next integration step
C                       (H .EQ. 0. ,if METAN1 fails to proceed)
C  * KFLAG     I   IN   Print parameter
C                        0   No output
C                        1   Integration monitor
C                        2   Intermediate solution points  T, Y(I),I=1,N
C                        3   Integration monitor and intermediate points
C                  OUT  Error flag
C                       .GE. 0  Successful integration
C                               (KFLAG not altered internally)
C                       -1   TEND .LT. T
C                       -2   More than NSTMAX basic integration steps
C                            per interval have been performed
C                       -3   More than JRMAX stepsize reductions
C                            occurred per basic integration step
C                       -4   Stepsize proposal for next basic
C                            integration too small
C
C
C* End Prologue
C  ------------
C
C
C  COMMON /STAT/ NFCN, NSTEP, NACCPT, NREJCT, NDEC, NSOL
C                       Internally initialized, for statistical
C                       purposes
C    NFCN               Number of FCN-evaluatios
C    NSTEP              Number of integration steps
C    NACCPT             Number of steps accepted (not used)
C    NREJCT             Number of steps rejected (not used)
C    NDEC               Number of LU-decompositions
C    NSOL               Number of forward-(backward-) substitutions
C
C* Type declaration
C
      INTEGER I, IFAIL, IPIVOT, IS, J, JK, JL, JM, JMACT, JOB, JOPT,
     2JRED, JRMAX, J1, K, KFIN, KFLAG, KM, KMACT, KOPT, K1, L, LOUT, M,
     3MAXODE, MDT, M1, N, NACCPT, NDEC, NFCN, NJ, NREJCT, NSOL, NSTEP,
     4NSTMAX
C
      DOUBLE PRECISION A, ALPHA, ANORM1, AWK, B1, C, CA, CB, CD0, CKAP,
     2CKQ, CMY, CMYH, CMYMAX, COSTF, COSTJ, COSTLR, COSTS, CMYRED, D,
     3DABS, DBLE, DEL, DELT0, DELT0A, DMAX1, DMIN1, DMY, DMT, DSQRT,
     4DT, DY, DZ, EPDIFF, EPKAP, EPMACH, EPSAFE, ERR, ETA, ETAD, ETADIF,
     5ETAMIN, ETAMAX, FC, FCK, FCM, FCO, FMIN, FNJ, FNJ1, FN1, FN, FOUR,
     6H, HALF, HJ, HJA, HL, HMAX, HMAXMY, HMAXT, HMAXU, HN, HREST,
     7HRTRN, M1ERRN, OMJ, OMJO, ONE, PCT1, QQ, RED, RMAX, RO, ROW, SAFE,
     8SAFEIN, SD, SMALL, SUMD, T, TAUD, TAUQ, TEND, THETA, THIRD, THMAX,
     9THMIN, THOPT, THQ, TN, TOL, TOLH, TOLMIN, TOLN, TWO, U, V, W, WY,
     AWZ, XQ1, XQ2, Y, YA, YM, YQ1, YQ2, YWGT, YWGTN, ZERO, ZQ1, ZQ2
C
      LOGICAL QDM, QDMA, QINCR, QJACNW, QJACUS, QKONV, QMY, QOPEN,
     2QPRMON, QPRSOL, QRED, QTEST
C
      CHARACTER CHGDAT*20, PRODCT*8
C
      EXTERNAL FCN
C
C* Constants problem oriented: (to be supplied by the user)
C
C    MAXODE    I   K    Maximal number of first-order ODE's
C
      PARAMETER ( MAXODE = 30            )
C
C* Other constants:
C
C    FOUR      D   K    4
C    HALF      D   K    1/2
C    ONE       D   K    1
C    PCT1      D   K    1 percent
C    THIRD     D   K    1/3
C    TWO       D   K    2
C    ZERO      D   K    0
C
      PARAMETER ( FOUR   = 4.0  D0       ,
     2            HALF   = 0.5  D0       ,
     3            ONE    = 1.0  D0       ,
     4            PCT1   = 1.01 D0       ,
     5            THIRD  = 1.0D0 / 3.0D0 ,
     6            TWO    = 2.0  D0       ,
     7            ZERO   = 0.0  D0       )
C
C* Control parameters: (to be supplied by the user)
C  Standard values fixed below
C
C    QJACUS    L   K    .TRUE. for user supplied Jacobian
C    NSTMAX    I   K    Maximum permitted number of integration steps
C                       per interval  =  10000
C    JRMAX     I   K    Maximum permitted number of stepsize reductions
C    KM        I   K    Prescribed maximum column number
C    JM        I   K    Associated maximum row number
C                       (JM = KM + 1)
C    MDT       I   K    Associated dimension of DT
C
      PARAMETER ( QJACUS = .FALSE.       ,
     2            NSTMAX = 10000         ,
     3            JRMAX  = 20            ,
     4            KM     = 6             ,
     5            JM     = KM + 1        ,
     6            MDT    = MAXODE*JM     )
C
C* Internal parameters: (modification not recommended)
C
      PARAMETER ( EPKAP  = 0.03 D0       ,
     2            ETADIF = 1.0  D-6      ,
     3            FMIN   = 1.0  D-3      ,
     4            RMAX   = 0.75 D0       ,
     5            RO     = 0.25 D0       ,
     6            SAFE   = 0.7  D0       ,
     7            SAFEIN = 2.0  D-2      ,
     8            THOPT  = 0.125 D0      )
C
C* Subroutines called
C
C    M1SEQ         EXT  Subroutine M1SEQ(JM,NJ)
C                       Generate stepsize sequence with respect to /1/
C                       JM     Maximum row number
C                       NJ     Array(JM) of stepsize sequence
C    M1JAC         EXT  Subroutine M1JAC( M, N, T, Y, A)
C                       (may be dummy)
C                       Provide Jacobian matrix A = DF/DY if desired
C                       METAN1 uses internally the matrix
C                              A = -DF/DY (scaled)
C                       M      Maximum number of differential equations
C                       N      Actual number of differential equations
C                       T      Actual position
C                       Y      Array(M) of values at T
C                       A      Array(M,M) DF/DY at T
C    M1SCAL        EXT  Subroutine M1SCAL (IFLAG, Y, N, YWGT, THRESH)
C                       Scaling for METAN1
C                       IFLAG   =0    Initial scaling
C                               else  Rescaling
C                       Y       Array of values Y(1),...,Y(N)
C                       N       Length of vectors Y and YWGT
C                       YWGT    Array of scaled values old
C                               Array of scaled values new
C                       THRESH  Threshold value
C    M1ERRN        EXT  Double precision function M1ERRN(Y, N, YWGT)
C                       Scaled root mean square error
C                       Y      Array of values Y(1),...,Y(N)
C                       N      Length of vectors Y and YWGT
C                       YWGT   Array of scaled values
C    DGEFA         EXT  Subroutine DGEFA (A, LDA, N, IPVT, INFO)
C                       LU-decomposition of a system of linear
C                       equations by Gauss-algorithm in the dense
C                       matrix case.
C                       A      IN  Matrix of coefficients
C                              OUT Matrix of LU-decomposition
C                       LDA    Leading dimension of A in calling program
C                       N      Dimension of the linear system
C                       IPVT   Vector of Pivot rows
C                       INFO   Information flag
C                              = 0         Solution successful
C                              = nonzero   A might be singular
C                       Calls subroutines and functions from BLAS:
C                       DAXPY, DSCAL, IDAMAX
C    DGESL         EXT  Subroutine DGESL (A, LDA, N, IPVT, B, JOB)
C                       Solution of a system of linear equations by
C                       Gauss-algorithm (Forward/backward substitution
C                       after LU-decomposition by DGEFA) in the dense
C                       matrix case.
C                       A      IN  Matrix of LU-decomposition
C                       LDA    Leading dimension of A in calling program
C                       N      Dimension of the linear system
C                       IPVT   IN  Vector of Pivot rows
C                       B      IN  Vector of right-hand side
C                              OUT Vector of solution
C                       JOB    =0          to solve A*X = B
C                              =nonzero    to solve trans(A)*X = B
C                       Calls subroutines and functions from BLAS:
C                       DAXPY, DDOT
C
C
C* Local variables: (workspace)
C
C     QDM      (DMY .GE. SMALL), i.e. CMY and CMYH estimated
C     QDMA     (DMY .GE. SMALL), i.e. CMY and CMYH estimated before
C     QJACNW   .TRUE.  Generate new Jacobian at next step
C     QMY      MY estimated
C     QPRMON   Print integration monitor
C     QPRSOL   Print intermediate solution points
C
C* Dimensions:
C
      DIMENSION A(MAXODE,MAXODE), ALPHA(JM,JM), AWK(JM), D(JM,JM),
     2 DEL(MAXODE), DT(MAXODE,JM), DY(MAXODE), DZ(MAXODE), ETA(MAXODE),
     3 FCK(KM), IPIVOT(MAXODE), NJ(JM), QQ(MAXODE,MAXODE), THQ(JM),
     4 Y(MAXODE), YM(MAXODE), YWGT(MAXODE), YWGTN(MAXODE)
C
      COMMON /STATP/ NFCN, NSTEP, NACCPT, NREJCT, NDEC, NSOL
C
C******  Revision 1.1 ******  Latest change:
      DATA      CHGDAT      /'July 17, 1989       '/
      DATA      PRODCT      /'METAN1'/
C***************************
C
C
      DATA  DT/MDT*0.D0/
C
C---1. Initial preparations
      CALL ZIBCONST(EPMACH,SMALL)
      LOUT   = 6
      EPSAFE = EPMACH*10.0D0
      QPRMON = (KFLAG .EQ. 1 .OR. KFLAG .EQ. 3)
      QPRSOL = (KFLAG .GE. 2)
      IF (TEND .LT. T) THEN
C        Error 1
         IF (QPRMON) WRITE (LOUT, 10001) PRODCT, T, TEND
         KFLAG = -1
         GOTO 9
C        EXIT to return
      ENDIF
      HREST = TEND - T
      H = DMIN1(H, HREST)
      HMAXU = HMAX
      IF (HMAX .GT. EPSAFE) THEN
         FCM = DMAX1(H/HMAX, FMIN)
      ELSE
         FCM = FMIN
      ENDIF
      KMACT = KM
      JMACT = JM
      CALL M1SEQ (JM, NJ)
      FN = DBLE(N)
      FN1 = DBLE(NJ(1))
      CMYRED = FN1*THIRD
      CMYMAX = FN1*0.9 D0
      TOLH = RO*TOL
      TOLN = TOL*TOL*FN
      THMIN = 1.0D-6
      CMYH = ZERO
      EPDIFF = DSQRT(EPSAFE)
      ETAMAX = DSQRT(EPDIFF)
      ETAMIN = EPDIFF*ETAMAX
      TOLMIN = EPSAFE*FN
      IF (TOL .LT. TOLMIN) THEN
         WRITE (LOUT, 10002) PRODCT, TOL, TOLMIN
         TOL = TOLMIN
      ENDIF
C
C---  Compute amount of work per row of extrapolation tableau
      AWK(1) = FN1 + ONE
      DO 101 J=2,JM
         J1 = J - 1
         FNJ = DBLE(NJ(J))
         V = AWK(J1) + FNJ
         AWK(J) = V
         DO 1011 K=1,J1
 1011       D(J,K) = (FNJ / DBLE(NJ(K)))*(FNJ / DBLE(NJ(K)))
C        ENDDO
         IF (J .NE. 2) THEN
            W = V - FN1
            DO 1012 K1=2,J1
               K = K1 - 1
               U = (AWK(K1) - V) / (W*DBLE(K + K1))
               U = TOLH**U
 1012          ALPHA(J1,K) = U
C           ENDDO
         ENDIF
 101     CONTINUE
C     ENDDO
C
C---1.1 Evaluation of cost coefficients
      COSTF = ONE
      COSTJ = FN
      COSTS = ZERO
      COSTLR = ZERO
      IF ((COSTS + COSTLR + COSTJ) .NE. ZERO) THEN
         AWK(1) = COSTJ + COSTLR + (COSTF + COSTS)*(FN1 + ONE)
         DO 11 J=2,JMACT
            J1 = J - 1
 11         AWK(J) = AWK(J1) +
     2         (COSTF + COSTS) * DBLE(NJ(J)) +
     3         COSTS + COSTLR
C        ENDDO
      ENDIF
C
C---1.2 Determination of maximum column number in extrapolation
C---    tableau (information theoretic concept, Ref./2/)
      KOPT = 1
      JOPT = 2
 121  CONTINUE
C     DO WHILE (JOPT .LT. JMACT .AND.
C               AWK(JOPT+1)*PCT1 .LE. AWK(JOPT)*ALPHA(JOPT,KOPT))
         IF (JOPT .GE. JMACT .OR.
     2      AWK(JOPT+1)*PCT1 .GT. AWK(JOPT)*ALPHA(JOPT,KOPT)) GOTO 122
C                                                         EXIT 121
         KOPT = JOPT
         JOPT = JOPT + 1
         GOTO  121
C     ENDDO
 122  KMACT = KOPT
      JMACT = JOPT
      IF (QPRMON) WRITE(LOUT, 11221)
     2   PRODCT, CHGDAT,TOL,KMACT,NJ
C
      IF (QPRSOL) WRITE(LOUT, 11222)
      NSTEP = 0
      QOPEN = .TRUE.
      NFCN = 0
      NDEC = 0
      NSOL = 0
      QJACNW = .TRUE.
      KFIN = 0
      OMJO = ZERO
      DO 123 I=1,N
 123     ETA(I) = ETADIF
C     ENDDO
C
C---  Initial scaling
      CALL M1SCAL (0, Y, N, YWGT, EPSAFE)
C
C---2. Basic integration step
 2    CONTINUE
C     DO WHILE (T .NE. TEND)
         IF (QPRMON) WRITE(LOUT, 12001) NSTEP,NFCN,T,KFIN,KOPT
         IF (QPRSOL) WRITE(LOUT, 12002) NSTEP,NFCN,T,H,(Y(I),I=1,N)
         JRED = 0
         CALL  FCN (N, T, Y, DZ, IFAIL)
         NFCN = NFCN + 1
C
C---2.2  Generate Jacobian
         IF (QJACNW) THEN
            IF (QJACUS) THEN
C
C---2.2.1      Analytic expression of Jacobian
               CALL M1JAC (MAXODE, N, T, Y, A)
               DO 221 I=1,N
                  DO 221 K=1,N
 221                 A(I,K) = -YWGT(K)*A(I,K) / YWGT(I)
C                 ENDDO
C              ENDDO
            ELSE
C
C---2.2.2      Numerical difference approximation of Jacobian
C              A = -DF/DY (scaled)
C              (feed-back control of discretization and rounding errors)
               DO 222 K=1,N
                  IS = 0
 2221             CONTINUE
C                 DO UNTIL (SUMD .GE. ETAMIN)
                     W = Y(K)
                     U = YWGT(K)*ETA(K)
                     IF (DZ(K) .GT. ZERO) U = -U
                     Y(K) = W + U
                     CALL FCN (N, T, Y, DY, IFAIL)
                     NFCN = NFCN + 1
                     Y(K) = W
                     U = YWGT(K) / U
                     SUMD = ZERO
                     DO 2222 I=1,N
                        WZ = DZ(I)
                        WY = DY(I)
                        W = DABS(WY)
                        WY = WZ - WY
                        WZ = DABS(WZ)
                        IF (W .LT. WZ) W = WZ
                        IF (W .NE. ZERO) THEN
                           W = WY / W
                           SUMD = SUMD + W*W
                        ENDIF
                        A(I,K) = WY*U / YWGT(I)
 2222                   CONTINUE
C                    ENDDO
                     SUMD = DSQRT(SUMD/FN)
                     IF (SUMD .EQ. ZERO .OR. IS .GT. 0) GOTO 222
C                                                       EXIT 2221
                     ETAD = DMIN1 (DSQRT(EPDIFF/SUMD)*ETA(K), ETAMAX)
                     ETAD = DMAX1 (ETAD, ETAMIN)
                     ETA(K) = ETAD
                     IS = 1
                     IF (SUMD .LT. ETAMIN) GOTO 2221
C                 ENDDO
 222              CONTINUE
C              ENDDO
            ENDIF
         ELSE
            IF (QPRMON) WRITE(LOUT, 12201)
         ENDIF
         DO 223 K=1,N
 223        DZ(K) = DZ(K) / YWGT(K)
C        ENDDO
C
C---2.4  Initial stepsize guess, if H = zero
C---     Estimation of lipschitz constant CL1
         ANORM1 = ZERO
         DO 24 I=1,N
            ROW = ZERO
            DO 241 K=1,N
 241           ROW = DABS(A(I,K)) + ROW
C           ENDDO
 24         ANORM1 = DMAX1(ANORM1,ROW)
C        ENDDO
         IF (ANORM1 .EQ. ZERO) ANORM1 = ONE
         IF (H .EQ. ZERO) THEN
            H = DMIN1 (SAFEIN/ANORM1, HREST)
         ENDIF
C
C---3.   Basic discretization step
 3       CONTINUE
C        DO WHILE (JRED .LE. JRMAX .AND. .NOT. QKONV)
            IF (H .EQ. HREST) THEN
               TN = TEND
            ELSE
               TN = T + H
            ENDIF
            IF (TN .EQ. T) THEN
C              Error 4
               IF (QPRMON) WRITE(LOUT, 13001) PRODCT
               KFLAG = -4
               GOTO  9
C              EXIT to return
            ENDIF
            HL = H*ANORM1
            HMAXT = H / FMIN
            QTEST = .FALSE.
            QJACNW = .FALSE.
            QINCR = .TRUE.
            QDM = .FALSE.
            QMY = .FALSE.
C
C---3.1     Internal discretization
            DO 31 J=1,JMACT
               M = NJ(J)
               M1 = M - 1
               KFIN = J - 1
               IF (J .GT. 1) THEN
                  FNJ1 = FNJ
                  HJA = HJ
               ENDIF
               FNJ = DBLE(M)
               HJ = H / FNJ
C
C---3.1.1      Semi-implicit Euler starting step
               DO 311 I=1,N
                  DO 3111 K=1,N
 3111                QQ(I,K) = HJ*A(I,K)
C                 ENDDO
                  QQ(I,I) = QQ(I,I) + ONE
                  YM(I) = Y(I)
 311              DEL(I) = HJ*DZ(I)
C              ENDDO
               IFAIL = 0
               CALL DGEFA (QQ, MAXODE, N, IPIVOT, IFAIL)
               NDEC = NDEC + 1
C
C---           Occurrence of zero pivot
               IF (IFAIL .NE. 0) THEN
                  RED = HJ*FN1*HALF / H
                  IF (QPRMON) WRITE(LOUT, 13111) PRODCT
                  GOTO  32
C                 EXIT 3.1 to stepsize reduction
               ENDIF
C
               JOB = 0
               CALL DGESL (QQ, MAXODE, N, IPIVOT, DEL, JOB)
               NSOL = NSOL + 1
               DO 3112 I=1,N
 3112             YM(I) = YM(I) + DEL(I)*YWGT(I)
C              ENDDO
C
C              *************************************************
C
C---3.1.2      Computational estimation of one-sided Lipschitz constant
               DMY = ZERO
               DO 3121 I=1,N
 3121             DMY = DMY + DEL(I)*DEL(I)
C              ENDDO
               DELT0A = DELT0
               DELT0 = DSQRT(DMY)
               QDMA = QDM
               QDM = (DMY .GE. TOLN)
               QMY = .FALSE.
               IF (J .GT. 1) THEN
                  IF (QDMA) THEN
                     CKAP = DELT0 / DELT0A
                     CKQ = FNJ1 / FNJ
                     C = ONE - CKAP
                     IF (C .GT. EPKAP) THEN
                        CMYH = FNJ*(CKQ - CKAP) / C
                        CMY = CMYH / H
                        QMY = .TRUE.
                     ENDIF
                     IF (QMY .AND. (CMYH .GE. CMYMAX)) THEN
                        RED = CMYRED / CMYH
                        IF (QPRMON) WRITE (LOUT, 13121) CMYH
                        GOTO  32
C                       EXIT 3.1 to stepsize reduction
                     ENDIF
                  ENDIF
                  IF (QMY) THEN
                     QMY =(CMY*HJA .GT. -ONE) .AND.
     2                        (CMY*HJA .LT. THIRD)
                     ZQ2 = HJ / (ONE - CMY*HJ)
                     ZQ1 = HJA / (ONE - CMY*HJA)
                  ENDIF
               ENDIF
C              *************************************************
C
C---3.1.3      Semi-implicit midpoint steps
               DO 313 K=1,M1
                  CALL  FCN (N, T + HJ*DBLE(K), YM, DY, IFAIL)
                  NFCN = NFCN + 1
                  DO 3132 I=1,N
 3132                DY(I) = HJ*DY(I) / YWGT(I) - DEL(I)
C                 ENDDO
                  JOB = 0
                  CALL DGESL (QQ, MAXODE, N, IPIVOT, DY, JOB)
                  NSOL = NSOL + 1
C
C                 *************************************************
C---              Newton restriction
C---              THETA .LE. THOPT*TWO
                  IF (K .EQ. 1 .AND. QDM) THEN
                     DMT = ZERO
                     DO 3134 I = 1,N
 3134                  DMT = DMT + DY(I)*DEL(I)
C                    ENDDO
                     THQ(J) = DMT / DMY
                     IF (J .EQ. 1) THEN
                        THETA = DABS(THQ(J))
                        IF (THETA .GE. THOPT*TWO) THEN
                           RED = DMAX1(THOPT, THOPT / THETA)
                           IF (QPRMON) WRITE(LOUT,13134) THETA
                           GOTO 32
C                          EXIT to stepsize reduction
                        ELSEIF (THETA .GT. THOPT) THEN
                           QINCR = .FALSE.
                        ENDIF
C
C---                    Monitor for Jacobian regeneration
                        QJACNW = THETA .GT. THMIN
                     ELSEIF (J .GT. 1) THEN
                        IF (DABS(THQ(J-1)) .GT. THMIN) THEN
                           C = THQ(J) / THQ(J-1)
                           IF (C .LT. ZERO .OR. C .GT. (ONE - EPKAP))
     2                        QINCR = .FALSE.
                        ENDIF
                     ENDIF
                     IF (QMY) THEN
                        IF (J .GT. 1) THEN
                           QTEST = .TRUE.
                           XQ2 = THQ(J) / ZQ2
                           YQ2 = XQ2 / ZQ2
                           XQ1 = THQ(J-1) / ZQ1
                           YQ1 = XQ1 / ZQ1
                           CA = (YQ1 - YQ2)*ZQ1*ZQ2 / (ZQ2 - ZQ1)
                           CB = DMAX1(DABS((XQ1-XQ2)/(ZQ1-ZQ2)),SMALL)
                           CD0 = DABS(CA)
                           TAUQ = HALF / DSQRT(CB)
                           CA = TWO*CD0*TAUQ
                           THMAX = 0.125D0
                           TAUD = TAUQ*FOUR*THMAX /
     2                            (CA + DSQRT(FOUR*THMAX + CA*CA))
                           SD = CMY*TAUD
C
C---                       Check for Newton Restriction
C---                       Implicit Euler Restriction
                           HMAXT = TAUD / EPSAFE
                           IF (SD .GE. (-ONE + EPSAFE)) THEN
                              HMAXT = TAUD / (ONE + SD)
                           ENDIF
                           HMAXT = HMAXT*FNJ
                           RED = HMAXT / H
                           IF (QPRMON .AND. QMY) THEN
                              WRITE (LOUT, 13131)
     2                              J, CKAP, CMYH, HL, RED
                           ENDIF
C
C---                       Possible stepsize reduction
                           IF (RED .LE. RMAX) THEN
                              IF (QPRMON) WRITE(LOUT, 13132)
                              RED = RED*SAFE
                              GOTO  32
C                             EXIT 3.1 TO STEPSIZE REDUCTION
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
C
C                 *************************************************
                  DO 3135 I=1,N
                     DEL(I) = DEL(I) + TWO*DY(I)
 3135                YM(I) = YM(I) + DEL(I)*YWGT(I)
C                 ENDDO
 313              CONTINUE
C              ENDDO
C
C---3.1.4      Smoothing Final Step
               CALL FCN (N, TN, YM, DY, IFAIL)
               NFCN = NFCN + 1
               DO 3142 I=1,N
 3142             DY(I) = HJ*DY(I) / YWGT(I) - DEL(I)
C              ENDDO
               JOB = 0
               CALL DGESL (QQ, MAXODE, N, IPIVOT, DY, JOB)
               NSOL = NSOL + 1
               DO 3143 I=1,N
 3143             YM(I) = YM(I) + DY(I)*YWGT(I)
C              ENDDO
C
C---3.1.5      Extrapolation
               DO 315 I=1,N
                  C = YM(I)
                  V = DT(I,1)
                  DT(I,1) = C
                  IF (J .NE. 1) THEN
                     YA = C
                     DO 3151 K=2,J
                        JK = J - K + 1
                        B1 = D(J,JK)
                        W = C - V
                        U = W / (B1 - ONE)
                        C = B1*U
                        V = DT(I,K)
                        DT(I,K) = U
 3151                   YA = U + YA
C                    ENDDO
                     YM(I) = YA
                     DY(I) = U
                  ENDIF
 315              CONTINUE
C              ENDDO
               IF (J .NE. 1) THEN
C
C---3.1.6         Convergence monitor
                  DO 3161 I=1,N
 3161                YWGTN(I) = YWGT(I)
C                 ENDDO
                  CALL M1SCAL (1, YM, N, YWGTN, EPSAFE)
                  ERR = M1ERRN (DY, N, YWGTN)
                  QKONV = ERR .LE. TOL
                  ERR = ERR / TOLH
C
C---              Order control
                  K = J - 1
                  FC = ERR**(ONE / DBLE(K + J))
                  FCK(K) = FC
C
C---              Order window
                  IF (J .GE. KOPT .OR. QOPEN) THEN
                     IF (QKONV) GOTO 25
C                                EXIT 3 for next basic integration step
C
C---                 Check for possible stepsize reduction
                     RED = ONE / FC
                     QRED = .FALSE.
                     IF (K .EQ. KMACT .OR. K .EQ. JOPT) THEN
                        RED = RED*SAFE
                        QRED = .TRUE.
                     ELSE
                        IF (K .EQ. KOPT) THEN
                           RED = RED*ALPHA(JOPT,KOPT)
                           IF (RED .LT. ONE) THEN
                              RED = ONE / FC
                              QRED = .TRUE.
                           ENDIF
                        ELSE
                           IF (KOPT .EQ. KMACT) THEN
                              RED = RED*ALPHA(KMACT,K)
                              IF (RED .LT. ONE) THEN
                                 RED = RED * SAFE
                                 QRED = .TRUE.
                              ENDIF
                           ELSE
                              RED = RED*ALPHA(JOPT,K)
                              IF (RED .LT. ONE) THEN
                                 RED = ALPHA(KOPT,K) / FC
                                 QRED = .TRUE.
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                     IF (QRED) GOTO 32
C                              EXIT 3.1 to stepsize reduction
                  ENDIF
               ENDIF
 31            CONTINUE
C           ENDDO
C
C---3.2     Prepare stepsize reduction
 32         CONTINUE
            IF (NSTEP .EQ. 0) THEN
               RED = DMIN1(RED, SAFE / HL)
            ENDIF
C
C---3.5     Stepsize reduction
            RED = DMIN1(RED, RMAX)
            H = H*RED
            IF (NSTEP .GT. 0) QOPEN = .FALSE.
            JRED = JRED + 1
            IF (QPRMON) WRITE(LOUT, 13501) JRED,RED,
     2         KFIN,KOPT,KMACT,THETA
            IF (JRED .GT. JRMAX) THEN
C              Error 3
               IF (QPRMON) WRITE(LOUT, 13502) PRODCT, JRMAX
               KFLAG = -3
               GOTO  9
C              EXIT to return
            ENDIF
            GOTO  3
C        ENDDO
C
C        ************************************************
C---2.5  Preparations for next basic integration step
 25      NSTEP = NSTEP + 1
         IF (NSTEP .GT. NSTMAX) THEN
C           Error 2
C           Emergency exit, if too many steps taken
            IF (QPRMON) WRITE(LOUT, 12501) PRODCT, NSTMAX
            KFLAG = -2
            GOTO  9
C           EXIT to return
         ENDIF
C
C---     Restoring
         DO 251 I=1, N
 251        Y(I) = YM(I)
C        ENDDO
         T = TN
         IF (T .EQ. TEND) GOTO 9
C                         EXIT to return
C
C---2.6  Rescaling
C
C---     Possible descaling of Jacobian
         IF (.NOT. QJACNW) THEN
            DO 261 I=1,N
               DO 261 K=1,N
 261              A(I,K) = A(I,K)*YWGT(I) / YWGT(K)
C              ENDDO
C           ENDDO
         ENDIF
C
C---     New scaling
         DO 262 I=1,N
 262        YWGT(I) = YWGTN(I)
C        ENDDO
C
C---     Possible rescaling of Jacobian
         IF (.NOT. QJACNW) THEN
            DO 263 I=1,N
               DO 263 K=1,N
 263              A(I,K) = A(I,K)*YWGT(K) / YWGT(I)
C              ENDDO
C           ENDDO
         ENDIF
C
C---2.7  Order and stepsize selection
C
C---2.7.1 Stepsize restrictions
         HMAXMY = H / FMIN
         IF (QMY) THEN
            IF (CMY .GE. SMALL) THEN
               HMAXMY = CMYRED / CMY
            ENDIF
         ENDIF
         HMAX = DMIN1(HMAXU,HMAXMY,HMAXT,H/FMIN)
         FCM = H / HMAX
C
C---2.7.2 Optimal order determination
         KOPT = 1
         JOPT = 2
         FCO = DMAX1(FCK(1), FCM)
         OMJO = FCO*AWK(2)
         IF (KFIN .GE. 2) THEN
            DO 272 L=2,KFIN
               JL = L + 1
               FC = DMAX1 (FCK(L), FCM)
               OMJ = FC*AWK(JL)
               IF (OMJ*PCT1 .LE. OMJO) THEN
                  KOPT = L
                  JOPT = JL
                  OMJO = OMJ
                  FCO = FC
               ENDIF
 272           CONTINUE
C           ENDDO
         ENDIF
         HREST = TEND - T
         HN = H / FCO
C
C---2.7.3 Possible increase of order
         IF (HN .LT. HREST .AND. QINCR .AND. CMYH .LE. CMYRED) THEN
            IF ((JRED .EQ. 0 .OR. NSTEP .EQ. 0) .AND.
     2           KOPT .GE. KFIN .AND. KOPT .NE. KMACT) THEN
               FC = DMAX1(FCO/ALPHA(JOPT,KOPT), FCM)
               JL = JOPT + 1
               IF (AWK(JL)*FC*PCT1 .LE. OMJO) THEN
                  FCO = FC
                  HN = H / FCO
                  KOPT = JOPT
                  JOPT = JOPT + 1
               ENDIF
            ENDIF
         ENDIF
C
C---2.7.4 Stepsize selection
         QOPEN = .FALSE.
         H = HN
         HRTRN = H
         IF (H .GT. HREST) THEN
            H = HREST
            QOPEN = .TRUE.
         ENDIF
         GO TO  2
C     ENDDO
C
C---9. EXIT
 9    HMAX = HMAXU
      IF ( KFLAG .LT. 0) THEN
C        Fail Exit
         H = ZERO
      ELSE
C        Solution Exit
         H = HRTRN
         IF (QPRMON) WRITE(LOUT, 12001) NSTEP,NFCN,T,KFIN,KOPT
         IF (QPRSOL) WRITE(LOUT, 12002) NSTEP,NFCN,T,H,(Y(I),I=1,N)
      ENDIF
      RETURN
C
C
10001 FORMAT(//,' ',A8,'  - ERROR -  '
     2      ,   ' Direction of integration is reverse to convention.')
10002 FORMAT(//,' ',A8,'  - WARNING -'
     2      ,   ' Desired tolerance ', D10.3, ' too small.', /,
     3      22X,' Tolerance set to  ', D10.3, '.')
11221 FORMAT(1H0,A8,' - ',A20,/,
     2       1H0,' rel.prec. TOL ',D10.3,' max.col.',I3,
     3       ' Sequence ',(1H ,13I4))
11222 FORMAT(//,5X,4HStep,3X,7HF-Calls,8X,1HT,25X,1HH,5X,7HY1(T)..,//)
12001 FORMAT(1H ,2I9,D20.11,I9,I6)
12002 FORMAT(1H ,2I9,D20.11,D12.5,4D20.11,/,(1H ,50X,4D20.11))
12201 FORMAT(' ',18H Old Jacobian kept)
12501 FORMAT(//,' ',A8,'  - ERROR -  '
     2      ,18H More than NSTMAX=,I3,18H integration steps,//)
13001 FORMAT(//,' ',A8,'  - ERROR -  '
     2      ,40H Stepsize reduction failed to succeed  ,//)
13111 FORMAT(//,' ',A8,'  - ERROR -  '
     2      ,36H Zero Pivot in linear system solver,/)
13121 FORMAT(1H0,'Stepsize reduction: CMYH ',D10.3,' too large',/)
13131 FORMAT (1H ,I2,3H CK,D10.3,3H MH,D10.3,
     2        3H HL,D10.3,3H RD,D10.3)
13132 FORMAT(1H0,'Stepsize reduction: Newton restriction activated',/)
13134 FORMAT(1H0,'Stepsize reduction: THETA ',D10.3,' too large',/)
13501 FORMAT(1H ,I3,27H Stepsize reduction factor ,D10.3,
     2      ' KFIN',I3,' KOPT',I3,' KMAX',I3,' THETA',D10.3)
13502 FORMAT(//,' ',A8,'  - ERROR -  '
     2      ,17H More then JRMAX=,I3,29H stepsize reductions per step,/)
C
C
C End METAN1
C
      END
      DOUBLE PRECISION FUNCTION M1ERRN(Y, N, YWGT)
C* Title:
C
C  Scaled root mean square error
C
C
C* Parameters:
C
C    Y         D   IN   Array of values Y(1),...,Y(N)
C    N         I   IN   Length of vectors Y and YWGT
C    YWGT      D   IN   Array of scaled values
C
C* Type declaration
C
      INTEGER N, I
      DOUBLE PRECISION Y, YWGT,
     2 DBLE, DSQRT, SUM, ZERO
C
C* Constants:
C
C    ZERO      D   K    0
C
      PARAMETER ( ZERO   = 0.0  D0       )
C
      DIMENSION Y(N), YWGT(N)
C
      SUM = ZERO
      DO 100 I=1,N
 100     SUM = SUM + (Y(I) / YWGT(I)) * (Y(I) / YWGT(I))
C     ENDDO
      M1ERRN = DSQRT(SUM / DBLE(N))
      RETURN
      END
      SUBROUTINE M1JAC (M, N, T, Y, A)
C
C     Provide Jacobian matrix A = DF/DY if desired
C     METAN1 uses internally the matrix
C            A = -DF/DY (scaled)
C
C
C* Parameters:
C
C    M         I   IN   Maximum number of differential equations
C    N         I   IN   Actual number of differential equations
C    T         D   IN   Actual position
C    Y         D   IN   Array(M) of values Y(1),...,Y(N)
C    A         D   OUT  Array(M,M) DF/DY at T
C
C* Type declaration
C
      INTEGER M, N
      DOUBLE PRECISION T, Y, A
      DIMENSION Y(M), A(M,M)
C
      WRITE (6, '(A)')
     2   ' +++ ERROR +++ M1JAC should not have been called.',
     3   '               Provide routine which computes Jacobian'
C
      RETURN
      END
      SUBROUTINE M1SCAL (IFLAG, Y, N, YWGT, THRESH)
C
C  Standard scaling
C  (For real life applications to be altered
C   by the skillful user)
C
C
C* Parameters:
C
C    IFLAG     I   IN   =0    Initial scaling
C                       else  Rescaling
C    Y         D   IN   Array of values Y(1),...,Y(N)
C    N         I   IN   Length of vectors Y and YWGT
C    YWGT      D   IN   Array of scaled values old
C              D   OUT  Array of scaled values new
C    THRESH    D   IN   Threshold value
C
C* Type declaration
C
      INTEGER IFLAG, N, I
      DOUBLE PRECISION Y, YWGT, THRESH,
     2 DABS, DMAX1, ONE, U, ZERO
C
C* Constants:
C
C    ONE       D   K    1
C    ZERO      D   K    0
C
      PARAMETER ( ONE    = 1.0  D0       ,
     2            ZERO   = 0.0  D0       )
C
      DIMENSION Y(N), YWGT(N)
      IF (IFLAG .EQ. 0) THEN
         DO 100 I=1,N
            U = DABS(Y(I))
            IF (U .EQ. ZERO) U = ONE
 100        YWGT(I) = DMAX1(U, THRESH)
C        ENDDO
      ELSE
         DO 200 I=1,N
 200        YWGT(I) = DMAX1(YWGT(I), DABS(Y(I)))
C        ENDDO
      ENDIF
      RETURN
      END
      SUBROUTINE M1SEQ(M,NJ)
      INTEGER M, NJ, I
      REAL ALPHA
      DIMENSION NJ(M)
C
C  Set stepsize sequence for semi-implicit method
C  with respect to Toeplitz condition
C
      PARAMETER (ALPHA = 5.0/7.0)
C
      NJ(1) = 2
      DO 10 I=2,M
        NJ(I) = NJ(I-1) + 4
C       Do While
 20       IF (FLOAT(NJ(I-1))/FLOAT(NJ(I)) .GT. ALPHA) THEN
            NJ(I) = NJ(I) + 4
            GOTO 20
          ENDIF
C       ENDDO
 10     CONTINUE
C     ENDDO
      RETURN
      END
