C=======================================================================
      SUBROUTINE SOLVE(NN,UL,B,X,IPS)
C=======================================================================
      DIMENSION UL(NN,1),B(1),X(1),IPS(1)

C      SUBROUTINE SOLVE WAS ORIGINALLY PUBLISHED TN "COMPUTER
C      SOLUTION OF LINEAR ALGEBRAIC SYSTEMS" BY C FORSYTHE AND
C      C B MOLER, PRENTICE-HALL 1967.
C      ADAPTED TO THIS PACKAGE BY C SODERLIND 1976-11-30 
C
      N = NN
      X(1) = B(1)/UL(1,1)
      IF(N.LT.2) RETURN
      NP1 = N+1
      IP = IPS(1)
      X(1) = B(IP)
      DO 2 I=2,N
         IP = IPS(I)
         IM1 = I-1
         SUM = 0.0
         DO 1 J=1,IM1
    1       SUM = SUM+UL(IP,J)*X(J)
    2    X(I) = B(IP)-SUM
C
      IP = IPS(N)
      X(N) = X(N)/UL(IP,N)
      DO 4 IBACK=2,N
           I = NP1-IBACK
           IP = IPS(I)
           IP1 = 1+I
           SUM = 0.0
           DO 3 J=IP1,N
    3         SUM = SUM+UL(IP,J)*X(J)
    4      X(I) = (X(I)-SUM)/UL(IP,I)
      RETURN
      END
C
C      * * * END SOLVE * * *
C
