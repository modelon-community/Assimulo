      SUBROUTINE EULEX (N,FCN,T,Y,TEND,TOL,HMAX,H,KFLAG)
C
C* Begin Prologue EULEX
C
C  ---------------------------------------------------------------------
C
C* Title
C
C    Explicit extrapolation integrator for non-stiff systems of
C    ordinary first-order differential equations.
C
C* Written by        P. Deuflhard, U. Nowak, U. Poehle
C* Purpose           Solution of systems of initial value problems
C* Method            Explicit Euler discretization with
C                    h-extrapolation
C* Category          i1a1c1. - System of nonstiff first order
C                              differential equations
C* Keywords          extrapolation, ODE, explicit Euler, nonstiff
C* Version           1.0 , February 1988
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
C  class II.
C
C
C  ---------------------------------------------------------------------
C
C* References:
C
C /1/ F.L.Bauer,H.Rutishauser,E.Stiefel:
C     New Aspects in Numerical Quadrature
C     Proc. Symp. Appl. Math. AMS 15, 199-218 (1963)
C
C /2/ R. Bulirsch, J. Stoer:
C     Fehlerabschaetzungen und Extrapolation ...
C     Numer. Math. 6, 413-427 (1964)
C
C /3/ P. Deuflhard:
C     Order and Stepsize Control in Extrapolation Methods
C     Numer. Math. 41, 399-422 (1983)
C
C
C* External Subroutine: (to be Supplied by the User)
C
C    FCN           EXT  Subroutine FCN(N,T,Y,DY)
C                       Right-Hand Side of First-Order
C                       Differential Equations
C                       N      Number of First-Order ODE'S
C                       T      Actual Position
C                       Y(N)   Values at T
C                       DY(N)  Derivatives at T
C
C
C* Parameters: (* Marks Transient Parameters)
C
C    N         I   IN   Number of ODE'S
C  * T         D   IN   Starting Point of Integration
C                       (T .LE. TEND)
C                  OUT  Achieved Final Point of Integration
C  * Y         D   IN   Array of Initial Values Y(1),...,Y(N)
C                  OUT  Array of Final Values
C    TEND      D   IN   Prescribed Final Point of Integration
C    TOL       D   IN   Prescribed Relative Precision (.GT.0)
C    HMAX      D   IN   Maximum Permitted Stepsize
C  * H         D   IN   Initial Stepsize Guess
C                  OUT  Stepsize Proposal for Next Integration Step
C                       (H .EQ. 0. ,if EULEX Fails to Proceed)
C  * KFLAG     I   IN   Print Parameter
C                        0   no Output
C                        1   Integration Monitor
C                        2   Intermediate Solution Points  T,Y(I),I=1,N
C                        3   Integration  Monitor and Solution Points
C                  OUT  Error Flag
C                       .GE. 0  Successful Integration
C                               (KFLAG not Altered Internally)
C                       -1   TEND .Lt. T
C                       -2   More Than NSTMAX Basic Integration Steps
C                            per Interval Have Been Performed
C                       -3   More Than JRMAX Stepsize Reductions
C                            Occurred per Basic Integration Step
C                       -4   Stepsize Proposal for Next Basic
C                            Integration too Small
C
C
C* End Prologue
C  ------------
C
C
C    COMMON /STAT/ NFCN, NSTEP, NACCPT, NREJCT, NDEC, NSOL
C                       Internally Initialized, for Statistical
C                       Purposes
C    NFCN               Number of FCN-Evaluatios
C    NSTEP              Number of Integration Steps
C    NACCPT             Number of Steps Accepted
C    NREJCT             Number of Steps Rejected
C    NDEC               Number of Decompositions
C    NSOL               Number of Substitutions
C
C* Type Declaration
C
C
      INTEGER I, J, JK, JL, JM, JMACT, JOPT, JRED, JRMAX, J1, K, KFIN,
     2KFLAG, KM, KMACT, KOPT, K1, L, LOUT, M, MAXODE, MDT, M1, N,
     3NACCPT, NDEC, NFCN, NJ, NREJCT, NSOL, NSTEP, NSTMAX
C
      DOUBLE PRECISION ALPHA, AWK, B1, C, D, DEL, DFLOAT, DMAX1, DMIN1,
     2DT, DUMMY, DY, DZ, EPMACH, EPSAFE, ERR, EXERRN, FC, FCK, FCM, FCO,
     3FMIN, FNJ, FN1, FN, H, HJ, HMAX, HMAXU, HN, HREST, HRTRN, OMJ,
     4OMJO, ONE, PCT101, PCT90, RED, RMAX, RO, SAFE, T, TEND, TN, TOL,
     5TOLH, TOLMIN, U, V, W, Y, YA, YM, YMAX, YWGT, ZERO, SMALL
C
      LOGICAL QFIRST, QINCR, QKONV, QLAST, QPRMON, QPRSOL, QRED
C
      CHARACTER CHGDAT*20, PRODCT*8
C
      EXTERNAL FCN
C
C* Constants Problem Oriented: (to be Supplied by the User)
C
C    MAXODE    I   K    Maximal Number of First-Order ODE'S
C
      PARAMETER ( MAXODE = 51          )
C
C* Other Constants:
C
C    ONE       D   K    1
C    PCT101    D   K    101 Percent
C    PCT90     D   K    90 Percent
C    ZERO      D   K    0
C
      PARAMETER ( ONE    = 1.0  D0       ,
     2            PCT101 = 1.01 D0       ,
     3            PCT90  = 0.9  D0       ,
     4            ZERO   = 0.0  D0       )
C
C* Control Parameters: (to be Supplied by the User)
C  Standard Values Fixed Below
C
C    NSTMAX    I   K    Maximum Permitted Number of Integration Steps
C                       per Interval
C    JRMAX     I   K    Maximum Permitted Number of Stepsize Reductions
C    KM        I   K    Prescribed Maximum Column Number
C    JM        I   K    Associated Maximum Row Number
C                       (JM = KM + 1)
C    MDT       I   K    Associated Dimension of DT
C    EXSEQ         EXT  Subroutine EXSEQ(JM,NJ)
C                       Generate Stepsize Sequence with Respect to /1/
C                       JM     Maximum Row Number
C                       NJ     Array(JM) of Stepsize Sequence
C    EXSCAL        EXT  Subroutine EXSCAL (MODE, Y, N, YOLD, YWGT,
C                                          YMAX, THREL, THABS)
C                       Scaling for DIFEX1
C                       MODE   ='INITIAL '    Initial Scaling
C                              ='INTERNAL'    Scaling during Discret.
C                              ='ACCEPTED'    Rescaling if Step Accepted
C                              Else           Error
C                       Y      Array of Values Y(1),...,Y(N)
C                       N      Length of Vectors Y, YOLD, YWGT, and YMAX
C                       YOLD   Array of Old Values
C                       YWGT   Array of Scaled Values
C                       YMAX   Array of Maximum Values
C                       THREL  Relative Threshold Value
C                       THABS  Absolute Threshold Value
C    EXERRN        EXT  Double Precision Function EXERRN(Y, N, YWGT)
C                       Scaled Root Mean Square Error
C                       Y      Array of Values Y(1),...,Y(N)
C                       N      Length of Vectors Y and YWGT
C                       YWGT   Array of Scaled Values
C
      PARAMETER ( NSTMAX = 10000         ,
     2            JRMAX  = 20            ,
     3            KM     = 12            ,
     4            JM     = KM + 1        ,
     5            MDT    = MAXODE*JM     )
C
C* Internal Parameters: (Modification not Recommended)
C
C
      PARAMETER ( FMIN   = 1.0  D-3      ,
     2            RMAX   = 0.75 D0       ,
     3            RO     = 0.25 D0       ,
     4            SAFE   = 0.7  D0       )
C
C
C* Local Variables: (Workspace)
C
C
C
C     QFIRST   First Integration Step
C     QKONV    Convergence Detected
C     QLAST    Last Integration Step
C     QPRMON   Print Integration Monitor
C     QPRSOL   Print Intermediate Solution Points
C
C* Dimensions:
C
      DIMENSION ALPHA(JM,JM), AWK(JM), D(JM,JM), DEL(MAXODE),
     2 DT(MAXODE,JM), DY(MAXODE), DZ(MAXODE), FCK(KM), NJ(JM),
     3 Y(MAXODE), YM(MAXODE), YMAX(MAXODE), YWGT(MAXODE)
C
      COMMON /STATP/ NFCN, NSTEP, NACCPT, NREJCT, NDEC, NSOL
C
C*******  Revision 1 *******  Latest Change:
      DATA      CHGDAT      /'February 19, 1988   '/
      DATA      PRODCT      /'EULEX'/
C***************************************************
C
C
      DATA  DT/MDT*0.D0/
C
C---1. Initial Preparations
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
C        EXIT TO RETURN
      ENDIF
C
      DO 1001 I = 1, N
         YMAX(I) = ZERO
 1001 CONTINUE
C
      HREST = TEND - T
      H = DMIN1(H, HREST)
      HMAXU = HMAX
      FCM = DMAX1(H/HMAX, FMIN)
      KMACT = KM
      JMACT = JM
      CALL EXSEQ (JM, NJ)
      FN = DFLOAT(N)
      FN1 = DFLOAT(NJ(1))
      TOLH = RO*TOL
      TOLMIN = EPSAFE*FN
      IF (TOL .LT. TOLMIN) THEN
         WRITE (LOUT, 10002) PRODCT, TOL, TOLMIN
         TOL = TOLMIN
      ENDIF
C
C---  Compute Amount of Work per Row of Extrapolation Tableau
      AWK(1) = FN1
      DO 101 J=2,JM
         J1 = J - 1
         FNJ = DFLOAT(NJ(J))
         V = AWK(J1) + FNJ - ONE
         AWK(J) = V
         DO 1011 K=1,J1
 1011       D(J,K) = FNJ / DFLOAT(NJ(K))
C        ENDDO
         IF (J .NE. 2) THEN
            W = V - AWK(1) + ONE
            DO 1012 K1=2,J1
               K = K1 - 1
               U = (AWK(K1) - V) / (W*DFLOAT(K1))
               U = TOLH**U
 1012          ALPHA(J1,K) = U
C           ENDDO
         ENDIF
 101     CONTINUE
C     ENDDO
C
C---1.2 Determination of Maximum Column Number in Extrapolation
C---    Tableau (Information Theoretic Concept, Ref./2/)
      KOPT = 1
      JOPT = 2
 121  CONTINUE
C     DO WHILE (JOPT .LT. JMACT .AND.
C               AWK(JOPT+1)*PCT101 .LE. AWK(JOPT)*ALPHA(JOPT,KOPT))
         IF (JOPT .GE. JMACT .OR.
     @      AWK(JOPT+1)*PCT101 .GT. AWK(JOPT)*ALPHA(JOPT,KOPT)) GOTO 122
C                                                         EXIT 121
         KOPT = JOPT
         JOPT = JOPT + 1
         GOTO  121
C     ENDDO
 122  KMACT = KOPT
      JMACT = JOPT
      IF (QPRMON) WRITE(LOUT, 11221)
     +   PRODCT, CHGDAT, TOL, KMACT, NJ
C
      IF (QPRSOL) WRITE(LOUT, 11222)
      NSTEP = 0
      QFIRST = .TRUE.
      QLAST = .FALSE.
      NFCN = 0
      KFIN = 0
      OMJO = ZERO
C
C---  Initial Scaling
      CALL EXSCAL ('INITIAL ', Y, N, DUMMY, YWGT, YMAX, TOL, ONE)
C
C---2. Basic Integration Step
 2    CONTINUE
C     DO WHILE (T .NE. TEND)
         IF (QPRMON) WRITE(LOUT, 12001) NSTEP,NFCN,T,H,KFIN,KOPT
         IF (QPRSOL) WRITE(LOUT, 12002) NSTEP,NFCN,T,H,(Y(I),I=1,N)
         JRED = 0
         CALL FCN( N, T, Y, DZ)
         NFCN = NFCN + 1
C
C---3.   Basic Discretization Step
 3       CONTINUE
C        DO WHILE (JRED .LE. JRMAX .AND. .NOT. QKONV)
            IF (QLAST) THEN
               TN = TEND
            ELSE
               TN = T + H
            ENDIF
            IF (TN .EQ. T) THEN
C              Error 4
               IF (QPRMON) WRITE(LOUT, 13001) PRODCT
               KFLAG = -4
               GOTO  9
C              EXIT TO RETURN
            ENDIF
            QINCR = .TRUE.
C
C---3.1     Internal Discretization
            DO 31 J=1,JMACT
               M = NJ(J)
               M1 = M - 1
               KFIN = J - 1
               IF (J .GT. 1) THEN
               ENDIF
               FNJ = DFLOAT(M)
               HJ = H / FNJ
               DO 311 I=1,N
 311              YM(I) = Y(I) + HJ*DZ(I)
C              ENDDO
C
C---3.1.3      Explicit Euler Steps
               DO 313 K=1,M1
                  CALL  FCN (N, T + HJ*DFLOAT(K), YM, DEL)
                  NFCN = NFCN + 1
                  DO 3133 I=1,N
 3133                YM(I) = YM(I) + HJ*DEL(I)
C                 ENDDO
 313              CONTINUE
C              ENDDO
C
C---3.1.5      Extrapolation
               ERR = ZERO
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
C---3.1.6         Convergence Monitor
                  CALL EXSCAL ('INTERNAL',YM,N,Y,YWGT,YMAX,TOL,ONE)
                  ERR = EXERRN (DY, N, YWGT)
                  QKONV = ERR .LE. TOL
                  ERR = ERR / TOLH
C
C---              Order Control
                  K = J - 1
                  FC = ERR**(ONE / DFLOAT(J))
                  FCK(K) = FC
C
C---              Order Window
                  IF (J .GE. KOPT .OR. QFIRST .OR. QLAST) THEN
                     IF (QKONV) GOTO 25
C                                EXIT 3 FOR NEXT BASIC INTEGRATION STEP
C
C---                 Check for Possible Stepsize Reduction
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
C                              EXIT 3.1 TO STEPSIZE REDUCTION
                  ENDIF
               ENDIF
 31            CONTINUE
C           ENDDO
C
C---3.2     Prepare Stepsize Reduction
 32         CONTINUE
C
C---3.5     Stepsize Reduction
            RED = DMIN1(RED, RMAX)
            H = H*RED
            IF (NSTEP .GT. 0) QLAST = .FALSE.
            JRED = JRED + 1
            IF (QPRMON) WRITE(LOUT, 13501) JRED,RED,
     2         KFIN,KOPT,KMACT
            IF (JRED .GT. JRMAX) THEN
C              Error 3
               IF (QPRMON) WRITE(LOUT, 13502) JRMAX
               KFLAG = -3
               GOTO  9
C              EXIT TO RETURN
            ENDIF
            GOTO  3
C        ENDDO
C
C        ************************************************
C---2.5  Preparations for Next Basic Integration Step
 25      NSTEP = NSTEP + 1
         QFIRST = .FALSE.
         IF (NSTEP .GT. NSTMAX) THEN
C           Error 2
C           EMERGENCY EXIT, IF TOO MANY STEPS TAKEN
            IF (QPRMON) WRITE(LOUT, 12501) PRODCT, NSTMAX
            KFLAG = -2
            GOTO  9
C           EXIT TO RETURN
         ENDIF
C
C---     Restoring
         DO 251 I=1, N
 251        Y(I) = YM(I)
C        ENDDO
         T = TN
         IF (T .EQ. TEND) GOTO 9
C                         Exit to Return
         CALL EXSCAL ('ACCEPTED', Y, N, DUMMY, YWGT, YMAX, TOL, ONE)
C
C---2.7  Order and Stepsize Selection
C
C---2.7.1 Stepsize Restrictions
         HMAX = DMIN1(HMAXU,H/FMIN)
         FCM = H / HMAX
C
C---2.7.2 Optimal Order Determination
         KOPT = 1
         JOPT = 2
         FCO = DMAX1(FCK(1), FCM)
         OMJO = FCO*AWK(2)
         IF (KFIN .GE. 2) THEN
            DO 272 L=2,KFIN
               JL = L + 1
               FC = DMAX1 (FCK(L), FCM)
               OMJ = FC*AWK(JL)
               IF (OMJ*PCT101 .LE. OMJO) THEN
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
C---2.7.3 Possible Increase of Order
         IF (HN .LT. HREST .AND. QINCR) THEN
            IF ((JRED .EQ. 0 .OR. NSTEP .EQ. 0) .AND.
     @           KOPT .GE. KFIN .AND. KOPT .NE. KMACT) THEN
               FC = DMAX1(FCO/ALPHA(JOPT,KOPT), FCM)
               JL = JOPT + 1
               IF (AWK(JL)*FC*PCT101 .LE. OMJO) THEN
                  FCO = FC
                  HN = H / FCO
                  KOPT = JOPT
                  JOPT = JOPT + 1
               ENDIF
            ENDIF
         ENDIF
C
C---2.7.4 Stepsize Selection
         H = HN
         HRTRN = H
         IF (H .GT. HREST*PCT90) THEN
            H = HREST
            QLAST = .TRUE.
         ENDIF
         GO TO  2
C     ENDDO
C
C---9. EXIT
 9    HMAX = HMAXU
      IF ( KFLAG .LT. 0) THEN
C        FAIL EXIT
         H = ZERO
      ELSE
C        SOLUTION EXIT
         H = HRTRN
         IF (QPRMON) WRITE(LOUT, 12001) NSTEP,NFCN,T,H,KFIN,KOPT
         IF (QPRSOL) WRITE(LOUT, 12002) NSTEP,NFCN,T,H,(Y(I),I=1,N)
      ENDIF
      RETURN
C
C
10001 FORMAT(//,' ',A8,'  - ERROR -  '
     2      ,   ' DIRECTION IF INTEGRATION IS REVERSE TO CONVENTION.')
10002 FORMAT(//,' ',A8,'  - WARNING -'
     2      ,   ' DESIRED TOLERANCE ', D10.3, ' TOO SMALL.', /,
     3      22X,' TOLERANCE SET TO  ', D10.3, '.')
11221 FORMAT(1H0,A8,' - ',A20,/,
     2       1H0,' REL.PREC. TOL ',D10.3,' MAX.COL. ',I3,
     3       ' SEQUENCE ',(1H ,13I4))
11222 FORMAT(//,5X,4HSTEP,3X,7HF-CALLS,8X,1HT,25X,1HH,5X,7HY1(T)..,//)
12001 FORMAT(1H ,2I9,D20.11,D12.5,I9,I6)
12002 FORMAT(1H ,2I9,D20.11,D12.5,4D20.11,/,(1H ,50X,4D20.11))
12501 FORMAT(//,' ',A8,'  - ERROR -  '
     2      ,18H MORE THAN NSTMAX=,I3,18H INTEGRATION STEPS,//)
13001 FORMAT(//,' ',A8,'  - ERROR -  '
     2      ,40H STEPSIZE REDUCTION FAILED TO SUCCEED  ,//)
13501 FORMAT(1H ,I3,27H STEPSIZE REDUCTION FACTOR ,D10.3,
     +      ' KFIN',I3,' KOPT',I3,' KMAX',I3)
13502 FORMAT(//,' ',A8,'  - ERROR -  '
     2      ,17H MORE THAN JRMAX=,I3,29H STEPSIZE REDUCTIONS PER STEP,/)
C
C
C End EULEX
C
      END
      SUBROUTINE EXSEQ(M,NJ)
      INTEGER I, M, NJ
      DIMENSION NJ(M)
C
C  Set Stepsize Sequence for Explicit Euler Method
C
      NJ(1) = 1
      DO 10 I=2,M
        NJ(I) = NJ(I-1) + 1
 10     CONTINUE
C     ENDDO
      RETURN
      END
      SUBROUTINE EXSCAL (MODE, Y, N, YOLD, YWGT, YMAX, THREL, THABS)
C
C     Scaling for EULEX
C
C       (for Real Life Applications to be Altered
C        by the Skillful User)
C
C
C* Parameters:
C
C    MODE      C*8 IN   ='INITIAL '    Initial Scaling
C                       ='INTERNAL'    Scaling during Discretization
C                       ='ACCEPTED'    Rescaling if Step Accepted
C                       Else           Error
C    Y         D   IN   Array of Values Y(1),...,Y(N)
C    N         I   IN   Length of Vectors Y, YOLD, YWGT, and YMAX
C    YOLD      D   IN   Array of Old Values
C    YWGT      D   OUT  Array of Scaled Values New
C    YMAX      D   IN   Array of Maximum Values Old
C                  OUT  Array of Maximum Values New
C    THREL     D   IN   Relative Threshold Value
C    THABS     D   IN   Absolute Threshold Value
C
C* Local Variables:
C
C    YUSER     D   V    User Defined Array of Maximum Values
C
C* Type Declaration
C
      INTEGER I, LOUT, MAXODE, N
C
      DOUBLE PRECISION DABS, DMAX1, EPMACH, ONE, THABS, THREL, U, Y,
     2YMAX, YOLD, YUSER, YWGT, ZERO, SMALL
C
      CHARACTER MODE*8
C
C* Constants:
C
C    EPMACH    D   K    Relative Machine Precision
C    LOUT      I   K    Output is Written on Logical Unit LOUT
C    MAXODE    I   K    Maximal Number of First-Order ODE's
C    ONE       D   K    1.0
C    ZERO      D   K    0.0
C
      PARAMETER (
     2            LOUT   = 6             ,
     3            MAXODE = 51            ,
     4            ONE    = 1.0  D0       ,
     5            ZERO   = 0.0  D0       )
C
      DIMENSION Y(N), YOLD(N), YWGT(N), YMAX(N), YUSER(MAXODE)
      SAVE YUSER
      CALL ZIBCONST(EPMACH,SMALL)
      IF (MODE .EQ.          'INITIAL '         ) THEN
C                             --------
         DO 100 I=1,N
            YUSER(I) = DABS (YMAX(I))
            U = DABS (Y(I))
            IF (U .LT. EPMACH) U = ONE
            YMAX(I) = DMAX1 (U, YUSER(I), THABS)
 100        YWGT(I) = YMAX(I)
C        ENDDO
      ELSE IF (MODE .EQ.     'INTERNAL'         ) THEN
C                             --------
         DO 200 I=1,N
 200        YWGT(I) = DMAX1 (YMAX(I)*THREL, DABS(Y(I)),
     +                       DABS(YOLD(I)), YUSER(I), THABS)
C        ENDDO
      ELSE IF (MODE .EQ.     'ACCEPTED'         ) THEN
C                             --------
         DO 300 I=1,N
 300        YMAX(I) = DMAX1 (YMAX(I), DABS(Y(I)))
C        ENDDO
      ELSE
         WRITE (LOUT, '(//,A,/)')
     +      ' EXSCAL    - ERROR -   Illegal Mode'
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION EXERRN(Y, N, YWGT)
C* Title:
C
C  Scaled Root Mean Square Error
C
C
C* Parameters:
C
C    Y         D   IN   Array of Values Y(1),...,Y(N)
C    N         I   IN   Length of Vectors Y and YWGT
C    YWGT      D   IN   Array of Scaled Values
C
C* Type Declaration
C
      INTEGER I, N
C
      DOUBLE PRECISION DFLOAT, DSQRT, SUM, Y, YWGT, ZERO
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
      EXERRN = DSQRT(SUM / DFLOAT(N))
      RETURN
      END
