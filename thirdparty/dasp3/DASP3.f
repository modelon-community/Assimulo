      SUBROUTINE DASP3 (DYDT, DZDT, OUTPDA, T,TEND,WSY,WSZ,N,M,TOL,
     *                  ABSREL,WGHT,EPS,
     *                  A,W,SLU,IPS,EQ,IND,LFLAG)
C       
C     AUTHOR G SODERLIND, DEPT OF NUMERICAL ANALYSIS,
C     ROYAL INSTITUTE OF TECHNOLOGY, STOCKHOLM
C     MAJOR REVISION G SODERLIND 1980-09-12
C     DOUBLE PRECISION VERSION:       1980-10-22
C     ------------------------------------------------------------------
C     
C     THIS SUBROUTINE IS USED FOR THE NUMERICAL INTEGRATION OF
C     THE PARTITIONED SYSTEM OF ORDINARY DIFFERENTIAL EQUATIONS
C     
C     DY/DT = F(T,Y,Z)        (N EQUATIONS)
C     
C     EPS*DZ/DT = G(T,Y,Z)        (M EQUATIONS)
C     
C     IT IS ASSUMED THAT THE FIRST SYSTEM IS NON-STIFF AND THAT
C     THE STIFFNESS OF THE SECOND SYSTEM IS DUE TO THE PARAMETER
C     EPS, A DIAGONAL "MATRIX" WITH SMALL ENTRIES. ONE OR MORE
C     DIAGONAL ENTRIES OF EPS MAY BE SET EQUAL TO ZERO IF A
C     DIFFERENTIAL-ALGEBRAIC SYSTEM IS TO BE SOLVED. THE SUB-
C     ROUTINE DASP3 REQUIRES THREE SUBROUTINES NAMED
C     
C     OUTPDA(T,Y,Z,N,M,JSTOP)
C     DYDT(T,Y,Z,F,N,M)
C     DZDT(T,Y,Z,G,N,M)
C     
C     FOR OUTPUT AND EVALUATIONS OF THE FUNCTIONS F AND G.
C     THE DERIVATIVES MUST BE COMPUTED TO DOUBLE PRECISION.
C     THE PARAMETERS TO THE SUBROUTINE DASP3 ARE:
C     
C     T       INDEPENDENT VARIABLE. ON CALL TO DASP3 IT SHOULD
C             BE SET TO THE INITIAL TIME. DOUBLE PRECISION,
C     
C     TEND    END OF INTEGRATION INTERVAL. TEND<T IS OK.
C             DOUBLE PRECISION.
C     
C     WSY     WORKSPACE VECTOR OF AT LEAST 10*N COMPONENTS FOR
C             THE NON-STIFF SYSTEM. INITIAL VALUES IN WSY(I),
C             I=1(1)N. DOUBLE PRECISION.
C     
C     WSZ     WORKSPACE VECTOR OF AT LEAST 9*M COMPONENTS FOR
C             THE STIFF SYSTEM. INITIAL VALUES IN WSZ(I),
C             I=1(1)M. DOUBLE PRECISION.
C     
C     N       NUMBER OF NON-STIFF EQUATIONS (DEP. VARIABLES)
C     
C     M       NUMBER OF STIFF AND/OR ALGEBRAIC EQUATIONS (DEP
C             VARIABLES). NOTE: M>=1.
C     
C     TOL     ERROR TOLERANCE PARAMETER. THE RELATIVE LOCAL
C             ERROR PER STEP WILL BE APPROXIMATELY EQUAL
C             TO TOL IN EACH STEP.
C     
C     ABSREL  AN ARRAY OF N+M DEPENDENT VARIABLE NORMALIZING
C             FACTORS USED FOR THE ESTIMATION OF LOCAL ERRORS.
C             ABSOLUTE ERROR IS ESTIMATED FOR COMPONENTS
C             SATISFYING ABS(Y(I)/ABSREL(I))<1, OTHERWISE THE
C             RELATIVE ERROR IS ESTIMATED.    (SIMILARLY FOR Z)
C             NOTE: FOR EVERY I,      ABSREL(I)>0.
C     
C     WGHT    AN ARRAY OF N+M ERROR WEIGHTS WHICH CAN BE USED
C             TO SET AN INDIVIDUAL ERROR TOLERANCE FOR EACH
C             DEPENDENT VARIABLE: THE MAXIMUM NORM OF WGHT SHOULD
C             EQUAL ONE. NOTE:   O<WGHT(I)<=1.
C     
C     ERR     AN M-DIMENSIONAL VECTOR CONTAINING THE M DIAGONAL
C             ELEMENTS OF THE "MATRIX" EPS. NOTE: EPS(I)>=0;
C     
C     A,W     STORAGE FOR TWO M*M MATRICES.
C     
C     SLU     A 2*M-VECTOR USED FOR LU-DECOMPOSITION OF THE
C             ITERATION MATRIX AND FOR VARIOUS OTHER TASKS.
C     
C     IPS     AN M-VECTOR OF INTEGERS USED FOR PARTIAL
C             PIVOTING IN LU-DECOMPOSITION.
C     
C     EQ      AN M-DIMENSIONAL (LOGICAL) ARRAY.
C     
C     IND     A 2*M-DIMENSIONAL VECTOR OF INTEGERS.
C     
C     LFLAG   ERROR FLAG USED IF THE INTEGRATION HAS TO BE
C             TERMINATED BEFORE TEND IS REACHED. LFLAG=0 IF
C             EVERYTHING IS OK. THE TERMINATION CODES HAVE
C             THE FOLLOWING MEANINGS:
C     
C      =1     ZERO ROW FOUND IN DECOMP (SINGULAR MATRIX)
C      =2     ZERO DIVIDE IN SOLVE
C      =3     STEP-SIZE TOO SMALL (REQUESTED ACCURACY MAY
C             NOT BE ATTAINABLE)
C      =4     TOO MANY FAILURES TO PASS THE ERROR TEST (DYDT
C             AND/OR DZDT MAY BE TOO ROUGH)
C      =5     NO CONVERGENCE IN INIVAL ITERATIONS (ONLY FOR
C             DIFFERENTIAL-ALGEBRAIC SYSTEMS)
C      =6     N<0 OR M<1
C     
C-----------------------------------------------------------------------
C  
      EXTERNAL DYDT, DZDT, OUTPDA
      DIMENSION WGHT(1),EPS(1),SLU(1),IPS(1),IND(1)
      DIMENSION A(M,1),W(M,1),ABSREL(1)
      DIMENSION BC(3),HH(3),HCOF(3)
      DOUBLE PRECISION T,TEND,RUNIT,DTEMP
      DOUBLE PRECISION WSY(1),WSZ(1)
      DOUBLE PRECISION SAVE(2)
      LOGICAL CHANGE,NEWJAC,INCORD,SKIP
      LOGICAL EQ(1)
      COMMON /COUNTS/ NYDER,NZDER,NSTEP,NREJ,MATEST,NREST,LUFACT
      COMMON /STEPS/ HH
      COMMON /JCIND/ IYN1,IZAUX,IZDIF,IZPC
      COMMON /ROFF/ RUNIT,SQRUNT
      EQUIVALENCE (H,HH(1))
      DATA BC/ 1.0,
     A       .666666667,
     B       .545454545 /
C     
C       INITIALIZING SECTION
C     
C      print *, 'Tol', tol
      IF(N.GE.0 .AND. M.GE.1) GOTO 3
      LFLAG = 6
      RETURN
    3 CONTINUE
      NREJ = 0
      NYDER = 0
      NZDER = 0
      NSTEP = 0
      NREST = 0
      LFLAG = 0
      JSTOP = 0
      MATEST = 0
      LUFACT = 0
      NOHMIN = 0
      NOSKIP = 0
C
      IYD1 = N
      IYF1 = IYD1+N
      IYF2 = IYF1+N
      IYF3 = IYF2+N
      IYF4 = IYF3+N
      IYN1 = IYF4+N
      IYK  = IYN1+N
      IYKD = IYK+N
      IYERR = IYKD+N
C
      IZD1 = M
      IZD2 = IZD1+M
      IZD3 = IZD2+M
      IZPC = IZD3+M
      IZDP = IZPC+M
      IZDIF = IZDP+M
      IZERR = IZDIF+M
      IZAUX = IZERR+M
C
C       CALCULATE ROUND-OFF UNIT
C
      RUNIT = 1.0D0
   10 RUNIT = .5D0*RUNIT
      IF(1D0+RUNIT.GT.1D0) GOTO 10
      RUNIT = 2D0*RUNIT
      SQRUNT = DSQRT(RUNIT)
C
      NALG = 0
      TOL = ABS(TOL)
      EMIN = 1.0
      HH(2) = 0.
      HH(3) = 0.
      DO 1 I=1,M
         ABSREL(I) = ABS(ABSREL(I))
         WGHT(I) = ABS(WGHT(I))
         IND(I) = I
         IND(I+M) = I+1
         EQ(I) = .FALSE.
         U = ABS(EPS(I))
         IF(U.NE.0.0) GOTO 2
         EQ(I) = .TRUE.
         NALG = NALG+1
    2    IF(U.LT.EMIN) EMIN = U
    1    CONTINUE
      IF(NALG.EQ.0) GOTO 9
         CALL INIVAL(DZDT, T,WSY,WSZ,ABSREL,N,M,A,NALG,EQ,SLU,IPS,TOL,
     *               LFLAG)
         IF(LFLAG.NE.0) RETURN
    9 TEMP = AMAX1(10.0,1E7*TOL)
      HMIN = TEMP*RUNIT*(TEND-T)
      IF(N) 20,30,20
   20 NYDER = 1
      CALL DYDT(T,WSY,WSZ,WSY(IYD1+1),N,M)
   30 STP = 1.0
      DO 40 I=1,M
         IF(EPS(I).NE.0.0) STP = AMAX1(STP,ABS(1.0/EPS(I)))
   40    CONTINUE
      CALL OUTPDA(T,WSY,WSZ,N,M,JSTOP)
      ITRY = 1
 1000 H = TEND-T
      NZDER = NZDER+1
      CALL DZDT(T,WSY,WSZ,WSZ(IZAUX+1),N,M)
      TEMP = ANORM(WSZ(IZAUX+1),WSZ,ABSREL(N+1),M)
      I = ITRY/4
      U = 1E-3/1E2**I
      HTRY = U/(STP*TEMP+1.0/ABS(H))
      HTRY = AMAX1(HTRY,ABS(1E2*HMIN))
C
C       CHECK DIRECTION OF INTEGRATION
C
      IF(H.LT.0.0) HTRY = -HTRY
      H = HTRY
      IF(NALG.EQ.M) GOTO 165
      IF (N) 710,730,710
C
C       FORWARD EULER STEP ON NONSTIFF VARIABLES
C
  710 DO 720 I=1,N
  720    WSY(I+IYF1) = WSY(I)+DBLE(H)*WSY(I+IYD1)
C
C       FORWARD EULER STEP, ON STIFF VARIABLES
C
  730 DO 740 I=1,M
         WSZ(I+IZD1) = WSZ(I)
C
C       SKIP IF EPS(I)=0
C
      IF(EPS(I).EQ.0.0) GOTO 740
      WSZ(I+IZD1) = WSZ(I)+DBLE(H/EPS(I))*WSZ(I+IZAUX)
  740 CONTINUE
C
C       FORWARD EULER STEP ON STIFF VARIABLES BACK TO T=TO
C
      NZDER = NZDER+1
      CALL DZDT(T+DBLE(H),WSY(IYF1+1),WSZ(IZD1+1),WSZ(IZD2+1),N,M)
      DO 750 I=1,M
         DTEMP = 0.0
         IF(EPS(I).NE.0.0) DTEMP = DBLE(H/EPS(I))*WSZ(I+IZD2)
         WSZ(I+IZD2) = WSZ(I+IZD1)-DTEMP
  750    CONTINUE
C
C       ESTIMATE NORM OF LOCAL ERROR ON FIRST STEP
C
      DO 760 I=1,M
         DTEMP = .5D0*DBLE(WGHT(I+N))
  760    WSZ(I+IZERR) = DTEMP*(WSZ(I+IZD2)-WSZ(I))
      TESTFF = ANORM(WSZ(IZERR+1),WSZ,ABSREL(N+1),M)
C
C      COMPUTE NEW STARTING STEP
C
      H = 0.5*H*SQRT(TOL/(TESTFF+1E-5*TOL))
      TEMP = TEND-T
      IF(H/TEMP.GT.1.0) H = 1E-2*TEMP
      IF(H/HMIN.LT.1.0) H = HMIN
  165 HHALF = .5*H
      JSTART = 1
      KORD = 1
      BETA = BC(1)
      INCORD = .FALSE.
C
C       INITIALIZE JMS (JACOBIAN MATRIX STATUS)
C       JMS = 1: JACOBIAN ESTIMATED ON CURRENT STEP
C             0: OLD JACOBTAN IN USE
C            -1: JACOBIAN NOT USED (FUNCTIONAL ITERATIONS) 
C
      JMS = 0
      NEWJAC = .TRUE.
      IF(EMIN.EQ.0.0) GOTO 45
         NEWJAC = .FALSE.
         JMS = -1
   45 IDBL = 2
      LCH = 0
C
C     CLEAR WORKSPACE VECTORS (SKIP OLD DATA ON RESTART)
C
      I0 = IZD1+1
      I1 = IZAUX
      DO 50 I=I0,I1
   50    WSZ(I) = 0D0
      HCOF(1) = H
      IF(N) 60,80,60
   60 I0 = IYF1+1
      I1 = IYERR+N
      DO 70 I=I0,I1
   70    WSY(I) = 0D0
   80 CONTINUE
      DO 100 I=1,M
         U = EPS(I)
         IF(U.NE.0.0) WSZ(I+IZD1) = WSZ(I+IZAUX)/DBLE(U)
  100    CONTINUE
      IF(EMIN.EQ.0.0) GOTO 1100
      IF(JMS.EQ.-1) GOTO 2000
C
C-----------------------------------------------------------------------       
C       NEW STEP ENTRY
C----------------------------------------------------------------------- 
C
 1100 IF(N) 130,1400,130
C
C-----------------------------------------------------------------------
C       PERFORM ONE FULL RK4 STEP ON THE NON-STIFF COMPONENTS
C-----------------------------------------------------------------------                     
C
C       COMPUTE YN+K0/6 AND YN+K0/2
C
  130 DO 140 I=1,N
          DTEMP = WSY(I+IYD1)*DBLE(HHALF)
          WSY(I+IYN1) = DTEMP
  140    WSY(I+IYK) = WSY(I)+DTEMP
C
C       COMPUTE YN+K0/6+K1/3+K2/3 AND YN+K2
C
      CALL PREPOL(.5,3,WSZ(IZPC+1),WSZ,M,1)
      DO 150 K=1,2
         DTEMP = T+DBLE(HHALF)
         CALL DYDT(DTEMP,WSY(IYK+1),WSZ(IZPC+1),WSY(IYKD+1),N,M)
         DO 150 I=1,N
            DTEMP = WSY(I+IYKD)
            WSY(I+IYN1) = WSY(I+IYN1)+DBLE(H)*DTEMP
  150       WSY(I+IYK) = WSY(I)+DBLE(HHALF)*K*DTEMP
C
C       COMPUTE YN1=YN+K0/6+K1/3+K2/3+K3/6
C
      CALL PREPOL(1.0,3,WSZ(IZPC+1),WSZ,M,1)
      CALL DYDT(T+DBLE(H),WSY(IYK+1),WSZ(IZPC+1),WSY(IYKD+1),N,M)
      DO 160 I=1,N
         WSY(I+IYN1) = (WSY(I+IYN1)+DBLE(HHALF)*WSY(I+IYKD))/3D0
  160    WSY(I+IYN1) = WSY(I)+WSY(I+IYN1)
      NYDER = NYDER+3
C
C-----------------------------------------------------------------------     
C       PREDICT Z FOR STIFF SUBSYSTEM
C-----------------------------------------------------------------------     
C
 1400 CALL PREPOL(1.0,KORD,WSZ(IZPC+1),WSZ,M,1)
      K = KORD-1
      CC = H*BETA
C
      CALL PDERIV(1.0,K,WSZ(IZDP+1),WSZ(IZD1+1),M,1)
      DO 170 I=1,M
         WSZ(I+IZDP) = DBLE(CC)*WSZ(I+IZDP)-WSZ(I+IZPC)
  170    WSZ(I+IZERR) = 0D0
C
      IF(.NOT.NEWJAC) GOTO 1600
C
C-----------------------------------------------------------------------
C       NEW DIFFERENCE APPROXIMATION TO JACOBIAN MATRIX
C-----------------------------------------------------------------------     
C
      MATEST = MATEST+1
      JMS = 1
      CALL JACEST(DZDT,T,H,A,N,M,WSY,WSZ,ABSREL,SLU,IND)
      IF(MATEST.EQ.1) CALL SPAPAT(M,A,IND,EQ,SLU)
      IF(EMIN.EQ.0.) GOTO 2000
      RHOJAC = ABS(CC)*CTRACT(A,M,N,WSZ,ABSREL,SLU,EPS)
      IF(RHOJAC.GT.0.25) GOTO 2000
      NQEST = 0
      Q = 100.0
      NEWJAC = .FALSE.
      JMS = -1
C
C-----------------------------------------------------------------------
C       ENTER CORRECTOR LOOP
C-----------------------------------------------------------------------
C
 1600 ICMAX = 3
      IF (JMS.EQ.-1) ICMAX = 4
      DO 200 ICORR=1,ICMAX
         RES = 1.0
         CALL DZDT(T+DBLE(H),WSY(IYN1+1),WSZ(IZPC+1),WSZ(IZAUX+1),N,M)
         NZDER = NZDER+1
C
C        TEST IF ITERATION MATRIX IS IN USE. NEWTON ITERATION
C        WILL BE USED IF JMS >= 0
C
         IF(JMS) 205,215,215
  205    DO 210 I=1,M
            WSZ(I+IZDIF) = DBLE(CC/EPS(I))*WSZ(I+IZAUX)-WSZ(I+IZDP)
  210       WSZ(I+IZAUX) = WSZ(I+IZPC)-WSZ(I+IZDIF)
         GOTO 260
C
C        FORM RESIDUAL VECTOR AND STORE IN SINGLE PRECISION VECTOR SLU.
C
  215    DO 220 I=1,M
            CDA = CC
            IF(EPS(I).EQ.0.0) CDA = 1.0
            DTEMP = DBLE(EPS(I))*(WSZ(I+IZPC)+WSZ(I+IZDP))
     *                     -DBLE(CDA)*WSZ(I+IZAUX)
            SLU(I) = DTEMP
  220       WSZ(I+IZDIF) = DTEMP
C
C        COMPUTE NORM OF RESIDUAL VECTOR
C
         RES = ANORM(WSZ(IZDIF+1),WSZ(IZPC+1),ABSREL(N+1),M)
C
C        SOLVE THE LINEAR SYSTEM W*X-D, WHERE D TS THE RESIDUAL
C        VECTOR AND W TS THE LU DECOMPOSITION OF THE MATRIX
C        DIAG(EPS(I))-H*BETA*J.
C
         IF(M-1) 240,230,240
  230    WSZ(9) = WSZ(7)*DBLE(W(1,1))
         GOTO 260
  240    CALL SOLVE(M,W,SLU,SLU(M+1),IPS)
C
C        MOVE SOLUTION VECTOR (CONTAINED IN SLU(M+*)) TO WSZ(*+IZAUX).
C        COMPUTE NORM OF SOLUTION VECTOR
C
          DO 250 I=1,M
  250        WSZ(I+IZAUX) = DBLE(SLU(I+M))
  260     CNORM = ANORM(WSZ(IZAUX+1),WSZ(IZPC+1),ABSREL(N+1),M)
C
C        OMIT ESTIMATION OF CONVERGENCE RATE FIRST TIME
C
          IF(ICORR.EQ.1) GOTO 280
             RHO = CNORM/CN1
          IF(RHO-1.0) 270,1700,1700
  270     Q = RHO/(1 .0-RHO)
          NQEST = 10
  280     CN1 = CNORM
C
C        COMPUTE NEW ITERATE AND ADD CORRECTION TO ERROR VECTOR.
C        PREPARE FOR CONVERGENCE TEST.
C
          SUP = 0.0
          SKIP = .TRUE.
          DO 290 I=1,M
             DTEMP = WSZ(I+IZAUX)
             TEMP = DTEMP
             DTEMP = WSZ(I+IZPC)-DTEMP
             IF(JMS.EQ.-1) DTEMP = WSZ(I+IZDIF)
             SKIP = SKIP .AND. DTEMP .EQ. WSZ(I+IZPC)
             WSZ(I+IZPC) = DTEMP
             WSZ(I+IZERR) = WSZ(I+IZERR)+TEMP
             DEN = AMAX1(ABSREL(I+N),ABS(SNGL(DTEMP)))
             CONV = WGHT(I+N)*ABS(Q*TEMP)/DEN
             IF(CONV.GT.SUP) SUP = CONV
  290        CONTINUE
C
C       CONVERGENCE TEST (SKIP IF THE CORRECTION WAS TOO SMALL)
C
         IF(SKIP .OR. RES.LE.1E-6*TOL) GOTO 1800
         IF(SUP.LT.BND.AND.NQEST.GT.0) GOTO 1800
  200    CONTINUE
C
C-----------------------------------------------------------------------         
C       THE CORRECTOR ITERATION FAILED TO CONVERGE. IF THE JACOBIAN
C       MATRIX WAS ESTIMATED ON THE CURRENT STEP THE STEP-SIZE WILL
C       BE REDUCED, OTHERWISE A NEW DIFFERENCE APPROXIMATION TO THE
C       MATRIX WILL BE MADE.
C----------------------------------------------------------------------- 
C       TEST MATRIX STATUS
C
 1700 IF(JMS.EQ.1) GOTO 300
      HNEW = H
      NEWJAC = .TRUE.
      GOTO 1400
C
  300 CHANGE = .TRUE.
      HNEW = .25*H
      JMS = 0
      GOTO 405
C
C----------------------------------------------------------------------- 
C     COMPUTE ERROR ESTIMATED AND TEST IF THE STEP WAS OK.
C-----------------------------------------------------------------------
C
 1800 IF(SKIP .AND. ICORR.EQ.1) NOSKIP = NOSKIP+1
      IF(NOSKIP.GE.100) GOTO 2100
      DO 310 I=1,M
  310    WSZ(I+IZERR) = ERCON*WGHT(I+N)*WSZ(I+IZERR)
      TESTFF = ANORM(WSZ(IZERR+1),WSZ(IZPC+1),ABSREL(N+1),M)/TOL
      TESTFF = AMAX1(TESTFF,1.0E-8)
      AZS = TESTFF**(-ERCON)*.95
      AYS = AZS
      IF(JSTART.EQ.0.OR.JSTART.GE.4) GOTO 320
C
C       CAREFUL STEPSIZE CONTROL FOR FIRST STEPS
C
      AG = AMIN1(.9*AZS,2.0)
      TEST = TESTFF
      GOTO 390
  320 IF(N) 340,330,340
  330 AYS = AZS
      TEST = TESTFF
      GOTO 380
  340 CALL PREPOL(1.0,5,WSY(IYERR+1),WSY,N,2)
      DO 370 I=1,N
         DTEMP = WSY(I+IYN1)
  370    WSY(I+IYERR) = DBLE(WGHT(I))*(WSY(I+IYERR)-DTEMP)
      TESTRK = ANORM(WSY(IYERR+1),WSY(IYN1+1),ABSREL,N)/TOL
      TESTRK = AMAX1(TESTRK,1.0E-8)
      AYS = TESTRK**(-.20)*.95
      TEST = AMAX1(TESTFF,TESTRK)
C
C       COMPUTE MAXIMUM POSSIBLE STEP SIZE
C
  380 AMAX = 10.0
      IF(EMIN.EQ.0.0) AMAX = 20.0
      AG = AMIN1(AZS,AYS,AMAX) 
C       
C        ERROR TEST
C
  390 HNEW = H
      IF(TEST.GT.1.5) GOTO 400
         IF(TEST.LT.1.0) GOTO 410
         AG = .95*AG
         GOTO 430
C
C       THE STEP WAS REJECTED DUE TO LARGE ERROR EXCESS
C
  400 HNEW = AG*H
      IF(JSTART.NE.0) IDBL = KORD+1
      CHANGE = .TRUE.
  405 ITRY = ITRY+1
      NREJ = NREJ+1
      IF(MOD(ITRY,4).NE.0) GOTO 500
      NREST = NREST+1
      IF(NREST.LT.25) GOTO 1000
      GOTO 2105
C
C       TEST IF STEP-SIZE INCREASE IS POSSIBLE
C
  410 IF(JSTART.EQ.0) GOTO 420
         JSTART = JSTART+1
         IF(IDBL.NE.0) GOTO 1900
         KORD = KORD+1
         IF(KORD.EQ.2) GOTO 415
            JSTART = 0
            AG = AMIN1(10.0,AYS,2.0*AZS)
            IF(AG.LT.1.2) AG = 1.0
            GOTO 430
  415    INCORD = .TRUE.
         CHANGE = .TRUE.
         IDBL = 4
         GOTO 1900
  420 IF(IDBL.GT.0.OR.AG.LT.1.2) GOTO 1900 
C
C        INCREASE STEP-SIZE (REDUCE THE STEP-SIZE IF THERE WAS
C        A SMALL ERROR EXCESS)
C
  430 HNEW = AG*H
      CHANGE = .TRUE.
      IDBL = 4
C
C----------------------------------------------------------------------- 
C     THE STEP WAS SUCCESSFUL - UPDATE DIVIDED DIFFERENCES &C.
C-----------------------------------------------------------------------     
C
 1900 T = T+DBLE(H)
      ITRY = 1
      IF(N) 440,470,440
  440 DTEMP = HCOF(2)
      IF(JSTART.EQ.2) DTEMP = 1.0
      DO 450 I=1,N
         WSY(I+IYK) = (WSY(I+IYN1)-WSY(I))/DBLE(H)
         WSY(I+IYKD) = (WSY(I+IYK)-WSY(I+IYD1))/DBLE(H)
  450    WSY(I) = WSY(I+IYN1)
C
C       EVALUATE THE DERIVATIVE AT THE BEGINNING OF THE NEXT STEP
C
      NYDER = NYDER+1
      CALL DYDT(T,WSY,WSZ(IZPC+1),WSY(IYD1+1),N,M)
      DO 460 I=1,N
         SAVE(1) = (WSY(I+IYKD)-WSY(I+IYF1))/DTEMP
         WSY(I+IYF1) = (WSY(I+IYD1)-WSY(I+IYK))/DBLE(H)
         SAVE(2) = (SAVE(1)-WSY(I+IYF2))/DTEMP
         WSY(I+IYF2) = (WSY(I+IYF1)-WSY(I+IYKD))/DBLE(H)
         WSY(I+IYF3) = (WSY(I+IYF2)-SAVE(1))/DTEMP
  460    WSY(I+IYF4) = (WSY(I+IYF3)-SAVE(2))/DTEMP
  470 KK = KORD
      IF(JSTART.EQ.0) GOTO 445
      KK = KORD+JSTART-2
      KK = MIN0(KK,3)
  445 DO 480 I=1,M
         IPL = 2
         INEW = IZD1
         IOLD = 0
         SAVE(1) = WSZ(I)
         WSZ(I) = WSZ(I+IZPC)
         DO 480 K=1,KK
            SAVE(IPL) = WSZ(I+INEW)
            IPL = 3-IPL
            WSZ(I+INEW) = (WSZ(I+IOLD)-SAVE(IPL))/DBLE(HCOF(K))
            IOLD = INEW
  480    INEW = INEW + M 
C
C       COMPUTE STEP-SIZE FOR KORD = 2
C
         IF(.NOT.INCORD) GOTO 7800
         INCORD = .FALSE.
         SUP = 0.0
         DO 7600 I=1,M
            TEMP = WSZ(I)
            DEN = AMAX1(ABSREL(I+N),ABS(TEMP))
            TEMP = WSZ(I+IZD3)
            TEMP = ABS(WGHT(I+N)*TEMP/DEN)
            IF(TEMP.GT.SUP) SUP = TEMP
 7600       CONTINUE
      TEMP = TOL/(2.0*SUP+1.0E-5*TOL)
      AZS = .5*TEMP**.33
      AG = AMIN1(20.0,AZS)
      AG = AMAX1(0.1,AG)
      HNEW = AG*H
 7800 CONTINUE
C
C       UPDATE COUNTS
C
      IDBL = IDBL-1
      NSTEP = NSTEP+1
      IF(JMS.EQ.1) JMS = 0
      NQEST = NQEST-1
      LCH = LCH+1
C
C       OUTPUT
C
      CALL OUTPDA(T,WSY,WSZ,N,M,JSTOP)
C
C       TEST FOR TERMINATION FORCED BY USER
C
      IF(JSTOP.EQ.-1) RETURN
C
C       LIMIT STEPSIZE USING FUNCTION HMAX
C
      TEMP = ABS(HMAX(T,TEND,WSY,WSZ,N,M))
      U = TEND-T
      IF(JSTOP.EQ.1) GOTO 485
      IF(U.LT.0.0) TEMP = -TEMP
      IF(HNEW/TEMP.LE.1.0) GOTO 485
         HNEW = TEMP
         CHANGE = .TRUE.
  485 U = U/HNEW
      IF(JSTOP.EQ.1 .AND. U.LT.1E-6) GOTO 2110
      IF(U-1.0) 455,465,475
  455 HNEW = TEND-T
      CHANGE = .TRUE.
  465 JSTOP = 1
  475 CONTINUE
      IF(.NOT.CHANGE.AND.LCH.GE.3) GOTO 1100
C
C        UPDATE STEP-SIZE DEPENDENT ARRAYS
C
      HH(3) = HH(2)
      HH(2) = HH(1)
      IF(HNEW/HMIN.GT.1.0) GOTO 500
         NOHMIN = NOHMIN+1
         IF(NOHMIN.EQ.100) GOTO 2100
         HNEW = HMIN
  500 HH(1) = HNEW
      HCOF(1) = HNEW
      HHALF = .5*HNEW
      HCOF(2) = HCOF(1)+HH(2)
      HCOF(3) = HCOF(2)+HH(3)
      IF(.NOT.CHANGE) GOTO 1100
      LCH = 0
      IF(EMIN.EQ.0 .OR. JMS.EQ.-1) GOTO 2000
      IF(RHOJAC*H/HH(2).LT.0.25) JMS = -1
C     
C-----------------------------------------------------------------------
C       FORM DIAG(EPS(I))-H*BETA*J MATRIX AND DECOMPOSE INTO LU
C       FORM, COMPUTE ERROR TOLERANCE &C.
C-----------------------------------------------------------------------
C     
 2000 CONTINUE
      BND = .5*TOL/FLOAT(KORD+2)
      ERCON = 1.0/FLOAT(KORD+1)
      CHANGE = .FALSE.
      BETA = BC(KORD)
      IF(JMS.EQ.-1) GOTO 610
      CC = H*BETA
      DO 560 J=1,M
         DO 560 I=1,M
            CDA = CC
            IF(EPS(I).EQ.0.0) CDA = 1.0
  560       W(I,J) = -CDA*A(I,J)
      DO 570 I=1,M
  570    W(I,I) = W(I,I)+EPS(I)
      IF(M-1) 590,580,590
  580 IF(W(1,1).EQ.0.0) GOTO 585
      W(1,1) = 1.0/W(1,1)
      GOTO 610
  590 CALL DECOMP(M,W,SLU,IPS,LFLAG)
      LUFACT = LUFACT+1
      IF(LFLAG) 2110,610,2110
C
C        FORCE ESTIMATION OF RHO
C
  610 Q = 100.0
      NQEST = 0
      IF(.NOT.NEWJAC) GOTO 1100
      NEWJAC = .FALSE.
      GOTO 1600 
C
C        RETURN
C
  585 LFLAG = 2
      RETURN
 2100 LFLAG = 3
      RETURN
 2105 LFLAG = 4
      RETURN
 2110 CONTINUE
      RETURN
      END
C
C       * * *     END DASP3 * * *
C
