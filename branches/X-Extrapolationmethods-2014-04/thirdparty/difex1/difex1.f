      SUBROUTINE DIFEX1 (N,FCN,T,Y,TEND,TOL,HMAX,H,KFLAG)
C
C* Begin prologue DIFEX1
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
C* Method            Explicit mid-point rule discretization with
C                    h**2-extrapolation
C* Category          i1a1c1. - System of nonstiff first order
C                              differential equations
C* Keywords          extrapolation, ODE, explicit mid-point rule,
C                    nonstiff
C* Version           1.2 , August 1991
C* Latest Change     August 1991
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
C  which may follow from acquisition or application of this code.
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
C /1/ W. B. Gragg:
C     On Extrapolation Algorithms for Ordinary
C     Initial Value Problems
C     SIAM J. Numer. Anal. 2, 384-404 (1965)
C
C /2/ R. Bulirsch, J. Stoer:
C     Numerical Treatment of Ordinary Differential Equations
C     by Extrapolation Methods
C     Num. Math. 8, 1-13 (1966)
C
C /3/ P. Deuflhard:
C     Order and Stepsize Control in Extrapolation Methods
C     Numer. Math. 41, 399-422 (1983)
C
C
C* External subroutine: (to be supplied by the user)
C
C    FCN           EXT  Subroutine FCN (N,T,Y,DY)
C                       Right-hand side of first-order
C                       differential equations
C                       N      Number of first-order ODEs
C                       T      Actual position
C                       Y(N)   Values at T
C                       DY(N)  Derivatives at T
C
C
C* Parameters: (* marks input/output parameters)
C
C    N         I   IN   Number of first-order ODEs
C  * T         D   IN   Starting point of integration
C                  OUT  Achieved final point of integration
C  * Y         D   IN   Array of initial values Y(1),...,Y(N)
C                  OUT  Array of final values
C    TEND      D   IN   Prescribed final point of integration
C    TOL       D   IN   Prescribed relative precision (.GT.0)
C    HMAX      D   IN   Maximum permitted stepsize
C  * H         D   IN   Initial stepsize guess
C                  OUT  Stepsize proposal for next integration step
C                       (H .EQ. 0. ,if DIFEX1 fails to proceed)
C  * KFLAG     I   IN   Print parameter
C                         0   no output
C                         1   Integration monitor
C                         2   Intermediate solution points  T,Y(I),I=1,N
C                         3   Integration monitor and solution points
C                       +10   if KFLAG is augmented by 10, a time
C                             monitor is printed additionally (this may
C                             be very expensive in terms of cpu-time)
C                  OUT  Error flag
C                       .GE. 0  Successful integration
C                               (KFLAG not altered internally)
C                       -2   More than NSTMAX basic integration steps
C                            per interval have been performed
C                       -3   More than JRMAX stepsize reductions
C                            occurred per basic integration step
C                       -4   Stepsize proposal for next basic
C                            integration step too small
C
C
C* End prologue
C  ------------
C
C
C    COMMON /STAT/ NFCN, NSTEP, NACCPT, NREJCT, NDEC, NSOL
C                       Internally initialized, for statistical
C                       purposes
C    NFCN               Number of FCN-evaluations
C    NSTEP              Number of integration steps
C    NACCPT             Number of steps accepted
C    NREJCT             Number of steps rejected
C    NDEC               Number of decompositions
C    NSOL               Number of substitutions
C
C* Type declaration
C
      INTEGER I, J, JK, JL, JM, JMACT, JOPT, JRED, JRMAX, J1, K, KFIN,
     $     KFLAG, KM, KMACT, KOPT, K1, L, LOUT, M, MAXODE, MDT, M1, N,
     $     NACCPT, NDEC, NFCN, NJ, NREJCT, NSOL, NSTEP, NSTMAX
C
      DOUBLE PRECISION ALPHA, AWK, D, DMAX1, DMIN1, DT, DUMMY, DY, DYM,
     $     DZ, D1ERRN, EPMACH, EPSAFE, ERR, FC, FCK, FCM, FCO, FMIN,
     $     FNJ, FN1, FN, H, HALF, HJ,  HJ2, HMAX, HMAXU, HN, HREST,
     $     HRTRN, OMJ, OMJO, ONE, PCT101, PCT90, QUART, RED, RMAX, RO,
     $     SAFE, SMALL, T, TEND, TN, TOL, TOLH, TOLMIN, U, Y, YK, YM,
     $     YMAX, YWGT, ZERO
C
      LOGICAL QFIRST, QKONV, QLAST, QPRMON, QPRSOL, QPRTIM, QRED
C
      CHARACTER CHGDAT*20, PRODCT*8
C
      EXTERNAL FCN
C
C* Constants problem oriented: (to be supplied by the user)
C
C    MAXODE    I   K    Maximal number of first-order ODEs
C
      PARAMETER ( MAXODE = 1024          )
C
C* Other Constants:
C
C    HALF      D   K    1/2
C    ONE       D   K    1
C    PCT101    D   K    101 Percent
C    PCT90     D   K    90 Percent
C    QUART     D   K    1/4
C    ZERO      D   K    0
C
      PARAMETER ( HALF   = 0.5  D0       ,
     $            ONE    = 1.0  D0       ,
     $            PCT101 = 1.01 D0       ,
     $            PCT90  = 0.9  D0       ,
     $            QUART  = 0.25 D0       ,
     $            ZERO   = 0.0  D0       )
C
C* Control parameters: (to be supplied by the user)
C  standard values fixed below
C
C    NSTMAX    I   K    Maximum permitted number of integration steps
C                       per interval  =  10000
C    JRMAX     I   K    Maximum permitted number of stepsize reductions
C    KM        I   K    Prescribed maximum column number
C    JM        I   K    Associated maximum row number
C                       (JM = KM + 1)
C    MDT       I   K    Associated dimension of DT
C    D1SEQ         EXT  Subroutine D1SEQ(JM,NJ)
C                       Generate stepsize sequence with respect to /1/
C                       JM     Maximum row number
C                       NJ     Array(JM) of stepsize sequence
C    D1SCAL        EXT  Subroutine D1SCAL (MODE, Y, N, YOLD, YWGT,
C                                          YMAX, THREL, THABS)
C                       Scaling for DIFEX1
C                       MODE   ='INITIAL '    Initial scaling
C                              ='INTERNAL'    Scaling during 
C                                             discretization
C                              ='ACCEPTED'    Rescaling if step accepted
C                              Else           Error
C                       Y      Array of values Y(1),...,Y(N)
C                       N      Length of vectors Y, YOLD, YWGT, and YMAX
C                       YOLD   Array of old values
C                       YWGT   Array of scaled values
C                       YMAX   Array of maximum values
C                       THREL  Relative threshold value
C                       THABS  Absolute threshold value
C    D1ERRN        EXT  Double Precision function D1ERRN(Y, N, YWGT)
C                       Scaled root mean square error
C                       Y      Array of values Y(1),...,Y(N)
C                       N      Length of vectors Y and YWGT
C                       YWGT   Array of scaled values
C
      PARAMETER ( NSTMAX = 10000         ,
     $            JRMAX  = 10            ,
     $            KM     = 8             ,
     $            JM     = KM + 1        ,
     $            MDT    = MAXODE*JM     )
C
C* Internal parameters: (modification not recommended)
C
C
      PARAMETER ( FMIN   = 1.0  D-3      ,
     $            RMAX   = 0.75 D0       ,
     $            RO     = QUART         ,
     $            SAFE   = 0.7  D0       )
C
C
C* Local variables: (workspace)
C
C
C
C    QFIRST    L   V    First integration step
C    QKONV     L   V    Convergence detected
C    QLAST     L   V    Last integration step
C    QPRMON    L   V    Print integration monitor
C    QPRSOL    L   V    Print intermediate solution points
C    QPRTIM    L   V    Print time monitor (this option may be
C                       expensive in terms of cpu-time)
C
C* Dimensions:
C
      DIMENSION ALPHA(JM,JM), AWK(JM), D(JM,JM), DT(MAXODE,JM),
     $     DY(MAXODE), DYM(MAXODE), DZ(MAXODE), FCK(KM), NJ(JM),
     $     Y(MAXODE), YK(MAXODE), YM(MAXODE), YMAX(MAXODE),
     $     YWGT(MAXODE)
C
      COMMON /STATP/ NFCN, NSTEP, NACCPT, NREJCT, NDEC, NSOL
C
C*******  Revision 1 *******  Latest change:
      DATA      CHGDAT      /'August 27, 1991'/
      DATA      PRODCT      /'DIFEX1'/
C***************************
C
C
C* Modification history
C  --------------------
C
C
C  1.0       Feb  9, 1988    First release at ZIB
C  1.1       Mar 27, 1991    Vectorize extrapolation loop,
C                            Time monitor
C  1.2       Aug 27, 1991    Allow reverse integration direction
C
C
      DATA  DT/MDT*0.D0/
C
C---1. Initial preparations
      CALL ZIBCONST(EPMACH,SMALL)
      LOUT   = 6
      EPSAFE = EPMACH*10.0D0
      HRTRN = H
      QPRTIM = (KFLAG .GE. 10)
      IF (QPRTIM) THEN
         KFLAG = KFLAG - 10
         CALL MONINI ('DIFEX1 Version 1.2', LOUT)
         CALL MONDEF (0, 'Integration')
         CALL MONDEF (1, 'FCN evaluation')
         CALL MONDEF (2, 'Extrapolation')
         CALL MONSTR (IFTIM)
         IF (IFTIM .NE. 0) THEN
            QPRTIM = .FALSE.
         ENDIF
      ENDIF
      QPRMON = (KFLAG .EQ. 1 .OR. KFLAG .EQ. 3)
      QPRSOL = (KFLAG .GE. 2)
C
      DO 1001 I = 1, N
         YMAX(I) = ZERO
 1001 CONTINUE
C
      HREST = TEND - T
      H = SIGN (DMIN1 (DABS(H), DABS(HREST)), HREST)
      HMAX = DABS(HMAX)
      HMAXU = HMAX
      FCM = DMAX1 (DABS(H)/HMAX, FMIN)
      KMACT = KM
      JMACT = JM
      CALL D1SEQ (JM, NJ)
      FN = DBLE (N)
      FN1 = DBLE (NJ(1))
      TOLH = RO*TOL
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
         FNJ = DBLE (NJ(J))
         AWK(J) = AWK(J1) + FNJ
         DO 1011 K=1,J1
            D(J,K) = (FNJ / DBLE (NJ(K)))*(FNJ / DBLE (NJ(K)))
 1011    CONTINUE
C
         IF (J .NE. 2) THEN
            DO 1012 K1=2,J1
               K = K1 - 1
               ALPHA(J1,K) = TOLH**((AWK(K1) - AWK(J)) /
     $              ((AWK(J) - AWK(1) + ONE)*DBLE(K + K1)))
 1012       CONTINUE
C
         ENDIF
 101  CONTINUE
C
C---1.2 Determination of maximum column number in extrapolation
C---    tableau (information theoretic concept, ref./3/)
      KOPT = 1
      JOPT = 2
 121  CONTINUE
C     DO WHILE (JOPT .LT. KMACT .AND.
C               AWK(JOPT+1)*PCT101 .LE. AWK(JOPT)*ALPHA(JOPT,KOPT))
         IF (JOPT .GE. KMACT .OR.
     $     AWK(JOPT+1)*PCT101 .GT. AWK(JOPT)*ALPHA(JOPT,KOPT)) GOTO 122
C                                                              Exit 121
         KOPT = JOPT
         JOPT = JOPT + 1
         GOTO  121
C     ENDDO
 122  KMACT = KOPT + 1
      JMACT = JOPT
      IF (QPRMON) WRITE (LOUT, 11221)
     $     PRODCT, CHGDAT, TOL, KMACT, NJ
C
      IF (QPRSOL) WRITE (LOUT, 11222)
      NSTEP = 0
      QFIRST = .TRUE.
      QLAST = .FALSE.
      NFCN = 0
      KFIN = 0
      OMJO = ZERO
      CALL D1SCAL ('INITIAL ', Y, N, DUMMY, YWGT, YMAX, TOL, ONE)
C
C---2. Basic integration step
 2    CONTINUE
C     DO WHILE (T .NE. TEND)
         IF (QPRMON) WRITE (LOUT, 12001) NSTEP,NFCN,T,H,KFIN,KOPT
         IF (QPRSOL) WRITE (LOUT, 12002) NSTEP,NFCN,T,H,(Y(I),I=1,N)
         JRED = 0
C
C---     Explicit euler starting step
         IF (QPRTIM) CALL MONON(1)
         CALL FCN (N, T, Y, DZ)
         NFCN = NFCN + 1
         IF (QPRTIM) CALL MONOFF(1)
C
C---3.   Basic discretization step
 3       CONTINUE
C        DO WHILE (JRED .LE. JRMAX .AND. .NOT. QKONV)
            IF (QLAST) THEN
               TN = TEND
            ELSE
               TN = T + H
            ENDIF
            IF (TN .EQ. T) THEN
C              Error 4
               IF (QPRMON) WRITE (LOUT, 13001) PRODCT
               KFLAG = -4
               GOTO  9
C              Exit to Return
            ENDIF
C
C---3.1     Internal discretization
            DO 31 J=1,JMACT
               M = NJ(J)
               M1 = M - 1
               KFIN = J - 1
               FNJ = DBLE (M)
               HJ = H / FNJ
               HJ2 = HJ + HJ
               DO 3101 I=1,N
                  YK(I) = Y(I)
                  YM(I) = Y(I) + HJ*DZ(I)
 3101          CONTINUE
C
C---3.1.3      Explicit mid-point rule
               DO 313 K=1,M1
                  IF (QPRTIM) CALL MONON(1)
                  CALL FCN (N, T + HJ*DBLE (K), YM, DY)
                  NFCN = NFCN + 1
                  IF (QPRTIM) CALL MONOFF(1)
                  DO 3135 I=1,N
                     U = YK(I) + HJ2*DY(I)
                     YK(I) = YM(I)
                     YM(I) = U
 3135             CONTINUE
 313           CONTINUE
C
C---3.1.4      Smoothing final step
               IF (QPRTIM) CALL MONON(1)
               CALL FCN (N, TN, YM, DY)
               NFCN = NFCN + 1
               IF (QPRTIM) CALL MONOFF(1)
               DO 3141 I = 1,N
                  YM(I) = (YM(I) + YK(I) + HJ*DY(I))*HALF
 3141          CONTINUE
C
C
C---3.1.5      Extrapolation
               IF (QPRTIM) CALL MONON(2)
               DO 3153 I=1,N
                  DY(I) = YM(I)
                  YK(I) = DT(I,1)
                  DT(I,1) = DY(I)
 3153          CONTINUE
C
               DO 3158 K=2,J
                  JK = J - K + 1
C
                  DO 3155 I=1,N
                     DYM(I) = (DY(I) - YK(I)) / (D(J,JK) - ONE)
                     DY(I) = D(J,JK)*DYM(I)
                     YK(I) = DT(I,K)
                     DT(I,K) = DYM(I)
                     YM(I) = DYM(I) + YM(I)
 3155             CONTINUE
C
 3158          CONTINUE
C
               IF (QPRTIM) CALL MONOFF(2)
C
               IF (J .NE. 1) THEN
C
C---3.1.6         Convergence monitor
                  CALL D1SCAL ('INTERNAL',YM,N,Y,YWGT,YMAX,TOL,ONE)
                  ERR = D1ERRN (DYM, N, YWGT)
                  QKONV = ERR .LE. TOL
                  ERR = ERR / TOLH
C
C---              Order control
                  K = J - 1
                  FC = ERR**(ONE / DBLE(K + J))
                  FCK(K) = FC
C
C---              Order window
                  IF (J .GE. KOPT .OR. QFIRST .OR. QLAST) THEN
                     IF (QKONV) GOTO 25
C                                Exit 3 for next basic integration step
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
C                              Exit 3.1 to stepsize reduction
                  ENDIF
               ENDIF
 31         CONTINUE
C
C---3.2     Prepare stepsize reduction
 32         CONTINUE
C
C---3.5     Stepsize reduction
            RED = DMIN1 (RED, RMAX)
            H = H*RED
            IF (NSTEP .GT. 0) QLAST = .FALSE.
            JRED = JRED + 1
            IF (QPRMON) WRITE (LOUT, 13501) JRED,RED,
     $           KFIN,KOPT,KMACT
            IF (JRED .GT. JRMAX) THEN
C              Error 3
               IF (QPRMON) WRITE (LOUT, 13502) JRMAX
               KFLAG = -3
               GOTO  9
C              Exit to Return
            ENDIF
            GOTO  3
C        ENDDO
C
C        ************************************************
C---2.5  Preparations for next basic integration step
 25      NSTEP = NSTEP + 1
         QFIRST = .FALSE.
         IF (NSTEP .GT. NSTMAX) THEN
C           Error 2
C           Emergency exit, if too many steps taken
            IF (QPRMON) WRITE (LOUT, 12501) PRODCT, NSTMAX
            KFLAG = -2
            GOTO  9
C           Exit to return
         ENDIF
C
C---     Restoring
         DO 251 I=1, N
            Y(I) = YM(I)
 251     CONTINUE
C
         T = TN
         IF (T .EQ. TEND) GOTO 9
C                         Exit to return
         CALL D1SCAL ('ACCEPTED', Y, N, DUMMY, YWGT, YMAX, TOL, ONE)
C
C---2.7  Order and stepsize selection
C
C---2.7.1 Stepsize restrictions
         HMAX = DMIN1(HMAXU, DABS(H)/FMIN)
         FCM = DABS(H) / HMAX
C
C---2.7.2 Optimal order determination
         KOPT = 1
         JOPT = 2
         FCO = DMAX1 (FCK(1), FCM)
         OMJO = FCO*AWK(2)
         IF (KFIN .GE. 2) THEN
            DO 272 L=2,KFIN
               JL = L + 1
               FC = DMAX1 (FCK(L), FCM)
               OMJ = FC*AWK(JL)
               IF (OMJ*PCT101 .LE. OMJO .AND. L .LT. KMACT) THEN
                  KOPT = L
                  JOPT = JL
                  OMJO = OMJ
                  FCO = FC
               ENDIF
 272        CONTINUE
         ENDIF
         HREST = TEND - T
         HN = H / FCO
C
C---2.7.3 Possible increase of order
         IF (DABS(HN) .LT. DABS(HREST)) THEN
            IF ((JRED .EQ. 0 .OR. NSTEP .EQ. 0) .AND.
     $           KOPT .GE. KFIN .AND. KOPT .NE. KMACT) THEN
               FC = DMAX1 (FCO/ALPHA(JOPT,KOPT), FCM)
               JL = JOPT + 1
               IF (AWK(JL)*FC*PCT101 .LE. OMJO .AND.
     $              JOPT .LT. KMACT) THEN
                  FCO = FC
                  HN = H / FCO
                  KOPT = JOPT
                  JOPT = JOPT + 1
               ENDIF
            ENDIF
         ENDIF
C
C---2.7.4 Stepsize selection
         H = HN
         HRTRN = H
         IF (DABS(H) .GT. DABS(HREST)*PCT90) THEN
            H = HREST
            QLAST = .TRUE.
         ENDIF
         GOTO  2
C     ENDDO
C
C---9. Exit
 9    HMAX = HMAXU
      IF (KFLAG .LT. 0) THEN
C        Fail exit
         H = ZERO
      ELSE
C        Solution exit
         H = HRTRN
         IF (QPRMON) WRITE (LOUT, 12001) NSTEP,NFCN,T,H,KFIN,KOPT
         IF (QPRSOL) WRITE (LOUT, 12002) NSTEP,NFCN,T,H,(Y(I),I=1,N)
         IF (QPRTIM) THEN
            CALL MONHLT
            CALL MONPRT
         ENDIF
      ENDIF
      RETURN
C
C
10001 FORMAT(//,' ',A8,'  - Error -  '
     $      ,   ' Direction if integration is reverse to convention.')
10002 FORMAT(//,' ',A8,'  - Warning -'
     $      ,   ' Desired tolerance ', D10.3, ' too small.', /,
     $      22X,' tolerance set to  ', D10.3, '.')
C
11221 FORMAT(1H0,A8,' - ',A20,/,
     $       1H0,' Rel.prec. TOL ',D10.3,' max.col. ',I3, /,
     $       ' sequence ',(1H ,13I4))
11222 FORMAT(//,5X,4HStep,3X,7HF-calls,8X,1HT,25X,1HH,5X,7HY1(T)..,//)
12001 FORMAT(1H ,2I9,D20.11,D12.4,I9,I6)
12002 FORMAT(1H ,2I9,D20.11,D12.4,4D20.11,/,(1H ,50X,4D20.11))
12501 FORMAT(//,' ',A8,'  - Error -  '
     $      ,18H more than NSTMAX=,I3,18H integration steps,//)
13001 FORMAT(//,' ',A8,'  - Error -  '
     $      ,40H stepsize reduction failed to succeed  ,//)
13501 FORMAT(1H ,I3,27H Stepsize reduction factor ,D10.3,
     $      ' KFIN',I3,' KOPT',I3,' KMAX',I3)
13502 FORMAT(//,' ',A8,'  - Error -  '
     $      ,17H more than JRMAX=,I3,29H stepsize reductions per step,/)
C
C
C End DIFEX1
C
      END
      SUBROUTINE D1SEQ(M,NJ)
      INTEGER I, M, NJ
      DIMENSION NJ(M)
C
C  Set stepsize sequence for DIFEX1
C
      NJ(1) = 2
      DO 10 I=2,M
         NJ(I) = NJ(I-1) + 2
 10   CONTINUE
C
      RETURN
      END
      SUBROUTINE D1SCAL (MODE, Y, N, YOLD, YWGT, YMAX, THREL, THABS)
C
C     Scaling for DIFEX1
C
C       (May be altered for real life applications
C        by the skillful user)
C
C
C* Parameters:
C
C    MODE      C*8 IN   ='INITIAL '    Initial scaling
C                       ='INTERNAL'    Scaling during discretization
C                       ='ACCEPTED'    Rescaling if step accepted
C                       Else           Error
C    Y         D   IN   Array of values Y(1),...,Y(N)
C    N         I   IN   Length of vectors Y, YOLD, YWGT, and YMAX
C    YOLD      D   IN   Array of old values
C    YWGT      D   OUT  Array of scaled values new
C    YMAX      D   IN   Array of maximum values old
C                  OUT  Array of maximum values new
C    THREL     D   IN   Relative threshold value
C    THABS     D   IN   Absolute threshold value
C
C* Local variables:
C
C    YUSER     D   V    User defined array of maximum values
C
C* Type declaration
C
      INTEGER I, LOUT, MAXODE, N
C
      DOUBLE PRECISION DABS, DMAX1, EPMACH, ONE, THABS, THREL, U, Y,
     $     YMAX, YOLD, YUSER, YWGT, ZERO
C
      CHARACTER MODE*8
C
C* Constants:
C
C    EPMACH    D   K    Relative machine precision
C    LOUT      I   K    Output is written on logical unit LOUT
C    MAXODE    I   K    Maximal number of first-order ODEs
C    ONE       D   K    1.0
C    ZERO      D   K    0.0
C
CCRY  (adapted to Cray X-MP)
CCRY  PARAMETER ( EPMACH = 7.106D-15     ,
C
CIBM  (adapted to Siemens 7.865, IBM 370-compatible)
CIBM  PARAMETER ( EPMACH = 2.22D-16      ,
C
CSUN  (adapted to sun)
      PARAMETER ( EPMACH = 0.1085D-18    ,
     $            LOUT   = 6             ,
     $            MAXODE = 1024          ,
     $            ONE    = 1.0  D0       ,
     $            ZERO   = 0.0  D0       )
C
      DIMENSION Y(N), YOLD(N), YWGT(N), YMAX(N), YUSER(MAXODE)
      SAVE YUSER
      IF (MODE .EQ.          'INITIAL '         ) THEN
C                             --------
         DO 100 I=1,N
            YUSER(I) = DABS (YMAX(I))
            U = DABS (Y(I))
            IF (U .LT. EPMACH) U = ONE
            YMAX(I) = DMAX1 (U, YUSER(I), THABS)
            YWGT(I) = YMAX(I)
 100     CONTINUE
C
      ELSE IF (MODE .EQ.     'INTERNAL'         ) THEN
C                             --------
         DO 200 I=1,N
            YWGT(I) = DMAX1 (YMAX(I)*THREL, DABS(Y(I)),
     $                       DABS(YOLD(I)), YUSER(I), THABS)
 200     CONTINUE
C
      ELSE IF (MODE .EQ.     'ACCEPTED'         ) THEN
C                             --------
         DO 300 I=1,N
            YMAX(I) = DMAX1 (YMAX(I), DABS(Y(I)))
 300     CONTINUE
C
      ELSE
         WRITE (LOUT, '(//,A,/)')
     $      ' D1SCAL    - ERROR -   Illegal mode'
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION D1ERRN(Y, N, YWGT)
C* Title:
C
C  Scaled root mean square error
C
C
C* Parameters:
C
C    Y         D   IN   Array of values Y(1),...,Y(N)
C    N         I   IN   Length of Vectors Y and YWGT
C    YWGT      D   IN   Array of scaled values
C
C* Type declaration
C
      INTEGER I, N
C
      DOUBLE PRECISION DBLE, DSQRT, SUM, Y, YWGT, ZERO
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
         SUM = SUM + (Y(I) / YWGT(I)) * (Y(I) / YWGT(I))
 100  CONTINUE
C
      D1ERRN = DSQRT(SUM / DBLE(N))
      RETURN
      END
