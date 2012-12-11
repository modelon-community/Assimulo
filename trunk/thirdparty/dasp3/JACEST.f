C=======================================================================             
      SUBROUTINE JACEST(DZDT,T,H,A,N,M,WSY,WSZ,ABSREL,SLU,IND)
C=======================================================================            
      EXTERNAL DZDT
      DIMENSION A(M,1),ABSREL(1),SLU(1),IND(1)
      DOUBLE PRECISION WSY(1),WSZ(1),T,RUNIT
      COMMON /JCIND/ IYN1,IZAUX,IZDIF,IZPC
      COMMON /COUNTS/ NYDER,NZDER,NSTEP,NREJ,MATEST,NREST,LUFACT
      COMMON /ROFF/ RUNIT,SQRUNT
C             
C      SUBROUTINE JACEST ESTIMATES THE JACOBIAN MATRIX
C      OF THE SYSTEM BY NUMERICAL DIFFERENCING,
C             
C      REVISED G SODERLIND 1980-05-29
C      DOUBLE PRECISION VERSION: 1980-10-22
C             
      NZDER = NZDER+1
      CALL DZDT(T+DBLE(H),WSY(IYN1+1),WSZ(IZPC+1),WSZ(IZAUX+1),N,M)
      K0 = 1
      DO 30 II=1,M
         K = IND(II+M)
         KM1 = K-1
         DO 10 JJ=K0,KM1
            J = IND(JJ)
            TEMP = WSZ(J+IZPC)
C               *** The next line was changed: I --> II  (Claus Fuhrer)
            SLU(J) = SQRUNT*AMAX1(ABSREL(II+N),ABS(TEMP))
   10       WSZ(J+IZPC) = WSZ(J+IZPC)+DBLE(SLU(J))
C             
         NZDER = NZDER+1
         CALL DZDT(T+DBLE(H),WSY(IYN1+1),WSZ(IZPC+1),WSZ(IZDIF+1),N,M)
C             
         DO 20 JJ=K0,KM1
            J = IND(JJ)
            WSZ(J+IZPC) = WSZ(J+IZPC)-DBLE(SLU(J))
            DO 20 I=1,M
               IF(A(I,J).EQ.0.0 .AND. MATEST.NE.1) GOTO 20
               A(I,J)  = (WSZ(I+IZDIF)-WSZ(I+IZAUX))/DBLE(SLU(J))
   20          CONTINUE
         IF(KM1.EQ.M) RETURN
   30    K0 = K
      RETURN
      END
C             
C      * * * END JACEST      * * *
C             

         

                    

          



                 
                 
                 
                 
