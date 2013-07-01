C     ==================================================================
      FUNCTION ANORM(X,REL,ABSREL,N)
C     ==================================================================
      DIMENSION ABSREL(1)
      DOUBLE PRECISION X(1),REL(1)
C            
C      ANORM COMPUTES A WEIGHTED MAX NORM OF X. THE WEIGHTS
C      ARE CALCULATED FROM ARRAYS REL AND ABSREL.
C      AUTHOR ,G SODERLIND 1980-10-22
C
      SUP = 0.0
      DO 1 I=1,N
         TEMP = REL(I)
         DEN = AMAX1(ABSREL(I),ABS(TEMP))
         TEMP = X(I)
         SUP = AMAX1(SUP,ABS(TEMP/DEN))
    1    CONTINUE
      ANORM = SUP
      RETURN
      END
C
C      * * *   END ANORM * * *
C
