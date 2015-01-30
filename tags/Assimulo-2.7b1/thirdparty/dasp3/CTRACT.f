C     ==================================================================
      FUNCTION CTRACT(A,M,N,WSZ,ABSREL,SLU,EPS)
C     ==================================================================
      DIMENSION A(M,1),ABSREL(1),SLU(1),EPS(1)
C             
C      CTRACT COMPUTES WEIGHTED MAXIMUM NORM OF A.
C      AUTHOR G SODERLIND 1980-09-01
C      DOUBLE PRECISION VERSION: 1980-10-22
C             
      DOUBLE PRECISION WSZ(1)
      COMMON /JCIND/ IYN1,IZAUX,IZDIF,IZPC
C             
C      STORE WEIGHTS IN ARRAY SLU
C             
      DO 10 I=1,M
            TEMP = WSZ(I+IZPC)
   10       SLU(I) = AMAX1(ABSREL(I+N),ABS(TEMP))
C             
C      COMPUTE WEIGHTED NORM OF A
C             
       CNORM = 0.
       DO 30 I=1,M
          ROWSUM = 0.
          DO 20 J=1,N
   20        ROWSUM = ROWSUM+ABS(A(I,J)*SLU(J))
          ROWSUM = ROWSUM/SLU(I)/ABS(EPS(I))
          IF(ROWSUM.GT.CNORM) CNORM = ROWSUM
   30     CONTINUE
       CTRACT = CNORM
       RETURN
       END
C             
C       * * *   END CTRACT      * * *
C    
