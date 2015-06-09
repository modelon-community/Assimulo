C=======================================================================
      SUBROUTINE INIVAL(DZDT, T,WSY,WSZ,SWPT,N,M,A,NALG,EQ,SLU,IPS,TOL,
     *                  LFLAG)
C=======================================================================
C
C       SUBROUTINE INIVAL CALCULATES INITIAL VALUES OF "ALGEBRAIC
C       STATES" CORRESPONDING TO THE GIVEN SET OF INITIAL VALUES
C       OF THE DIFFERENTIAL EQUATIONS. THE USER IS SUPPOSED TO
C       PROVIDE SOME STARTING VECTOR FOR THE ITERATIONS. IF NO
C       SUCH VECTOR IS KNOWN, THE STARTING VECTOR SHOULD BE SET
C       TO ZERO. NOTE THAT INIVAL MAY FAIL TO FIND THE INITIAL
C       VALUES IF THE STARTING VECTOR IS UNKNOWN AND THE PROBLEM
C       IS STRONGLY NONLINEAR.  (LFLAG=5 ON RETURN)
C             
C       AUTHOR G SODERLIND 1980-09-01
C       DOUBLE PRECISION VERSION:       1980-10-22
C  
      EXTERNAL DZDT  
      DIMENSION SLU(1),IPS(1),SWPT(1),A(M,1)
      DOUBLE PRECISION WSY(1),WSZ(1),T,RUNIT,DTEMP
      DOUBLE PRECISION CZIT,CZOLD,CJG,ZJ,ZINC
      COMMON /JCIND/ IYN1,IZAUX,IZDIF,IZPC
      COMMON /COUNTS/ NYDER,NZDER,NSTEP,NREJ,MATEST,NREST,LUFACT
      COMMON /ROFF/ RUNIT,SQRUNT
      LOGICAL EQ(1)
      print *, 'inival', A(1,1)
      DAMP = .7
      IZERR = 7*M
      TH = 0.
C             
C      MOVE STARTING VECTOR
C             
      DO 10 I=1,M
   10    WSZ(I+IZERR) = WSZ(I)
C             
C      ENTER MAIN CONTINUATION LOOP
C             
      DO 170 IMB=1,10
C             
C      FIRST ORDER DAMPED PREDICTOR STEP
C             
         DO 20 I=1,M
            DTEMP = WSZ(I)
            WSZ(I) = DTEMP+DBLE(DAMP)*(DTEMP-WSZ(I+IZERR))
   20       WSZ(I+IZERR) = DTEMP
C             
C      WSZ(*) NOW CONTAINS PREDICTED VECTOR.
C      WSZ(*+IZERR) CONTAINS SOLUTION VECTOR CALCULATED IN
C      PREVIOUS CONTINUATION STEP.
C      COMPUTE PARAMETERS FOR THIS CONTINUATION STEP.
C             
         TH = TH+.1
         TH1 = 1.-TH
         CZOLD = DBLE(DAMP*TH1)
         CZIT = 1D0-CZOLD
         CJG = DBLE(DAMP*TH)
         NIX = -1
   30    NIX =-NIX+1
         Q = 1.0E6
C     
C      ESTIMATE JACOBIAN AT PREDICTED POINT
C             
         NZDER = NZDER+NALG+1
         CALL DZDT(T,WSY,WSZ,WSZ(IZAUX+1),N,M)
         K = 0
         DO 50 J=1 ,M
            IF(.NOT.EQ(J)) GOTO 50
            K = K+1
            ZJ = WSZ(J)
            TEMP = ZJ
            ZINC = SQRUNT*AMAX1(SWPT(I+N),ABS(TEMP))
            WSZ(J) = ZJ+ZINC
            CALL DZDT(T,WSY,WSZ,WSZ(IZDIF+1),N,M)
            WSZ(J)  = ZJ
            KK = 0
C      
C      STORE ENTRIES OF JACOBIAN IN CONSECUTIVE MEMORY
C      LOCATIONS
C               
            DO 40 I=1,M
               IF(.NOT.EQ(I)) GOTO 40
               KK = KK+1
               A(KK,K) = (WSZ(I+IZDIF)-WSZ(I+IZAUX))/ZINC
   40          CONTINUE
   50       CONTINUE
         IF(NALG.EQ.1) GOTO 60
C             
C      LU FACTORIZATION OF JACOBIAN MATRIX
C             
         CALL DECOMP(NALG,A,SLU,IPS,LFLAG)
         print *, 'decomp lflag=', lflag
         IF(LFLAG.NE.0) RETURN
C             
C      ENTER CORRECTOR LOOP.
C      UP TO FIVE NEWTON ITERATIONS ARE USED FOR A GIVEN THETA
C             
   60    DO 150 ICORR=1,5
            NZDER = NZDER+1
C             
C      COMPUTE RESIDUAL
C             
            CALL DZDT(T,WSY,WSZ,WSZ(IZAUX+1),N,M)
C             
C      MOVE RESIDUALS OF ALGEBRAIC EQUATIONS TO CONSECUTIVE
C      MEMORY LOCATIONS IN SINGLE PRECISION VECTOR SLU.
C             
            K=0
            DO 70 I=1,M
               IF(.NOT.EQ(I)) GOTO 70
               K = K+1
               SLU(K) = WSZ(I+IZAUX)
   70          CONTINUE
            IF(NALG.GT.1) GOTO 90
            TEMP =  A(1,1)
            IF(TEMP.NE.0.)  GOTO 80
            LFLAG = 1
            print *, 'here'
            RETURN
   80       SLU(2) =  SLU(1)/TEMP
            GOTO 100
   90       CALL SOLVE(NALG,A,SLU,sLU(NALG+1),IPS)
C             
C      COMPUTE NEW ITERATE
C             
  100       K = 0
            DO 110 I=1,M
               WSZ(I+IZDIF)  = WSZ(I)
               IF(.NOT.EQ(I))  GOTO 110
               K = K+1
               DTEMP = CJG*DBLE(SLU(NALG+K))
               WSZ(I) = CZIT*WSZ(I)-DTEMP+CZOLD*WSZ(I+IZERR)
  110          CONTINUE
C             
C      COMPUTE NORM OF CORRECTION
C             
            DO 120 I=1,M
  120          WSZ(I+IZAUX) = DABS(WSZ(I)-WSZ(I+IZDIF))
            CNORM  = ANORM(WSZ(IZAUX+1),WSZ,SWPT(N+1),M)
C                 
C      ESTIMATE CONVERGENCE RATE
C             
            IF(ICORR.EQ.1)  GOTO 130
                  RHO = CNORM/CN1
                  IF(RHO.GE.1.0) GOTO 160
                  Q = RHO/(1 .0-RHO)
  130          CN1 = CNORM
C             
C      CHECK FOR CONVERGENCE
C             
            SUP = Q*CNORM
            IF(SUP.GT.TOL) GOTO 150
C             
C      CONVERGED. INITIAL VALUES FOUND IF THETA=1 (IMB=10).
C             
            IF(IMB.EQ.10) RETURN
C                 
C      ELSE NEXT CONTINUATION STEP (END OF CONTINUATION LOOP).
C             
            GOTO 170
C             
C      NOT YET CONVERGED. NEXT ITERATION (SAME THETA).
C             
  150       CONTINUE
C             
C      NOT CONVERGED AFTER FIVE ITERATIONS. ESTTMATE NEW
C      JACOBIAN AND TRY ANOTHER FIVE ITERATIONS.
C             
         IF(NIX.EQ.0) GOTO 30
C             
C      ELSE DIVERGED. NO INITIAL VALUES FOUND WITH THIS
C      STARTING VECTOR.
C             
  160    LFLAG = 5
         RETURN
  170    CONTINUE
      RETURN
      END
C             
C             * * *   END INIVAL      * * *
C 
