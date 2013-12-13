C=======================================================================
      SUBROUTINE DECOMP(NN,UL,SLU,IPS,LFLAG)
C=======================================================================
      DIMENSION UL(NN,1),SLU(1),IPS(1)
C             
C             
C      SUBROUTINE DECOMP WAS ORIGINALLY PUBLISHED IN 'COMPUTER
C      SOLUTION OF LINEAR ALGEBRAIC SYSTEMS" BY G FORSYTHE AND
C      C B MOLER,      PRENTICE-HALL 1967.
C      ADAPTED TO THIS PACKAGE BY G SODERLIND 1976-11-30
C          
      N = NN
      DO 5 I=1 ,N
           IPS(I) = I
           ROWNRM = 0.0
           DO 2 J=1,N
              U = ABS(UL(I,J))
              IF(U.GT.ROWNRM) ROWNRM = U
    2         CONTINUE
           IF(ROWNRM) 5,4,5
    4      LFLAG = 1
           RETURN
    5      SLU(I) = 1.0/ROWNRM
C                 
C             GAUSS ELIM + PART PIT
C             
      NM1 = N-1
      IF(NM1.LT.1) GOTO    18
      DO 17 K=1,NM1
         BIG = 0.0
         DO 11 I=K,N
         IP = IPS(I)
         SIZE = ABS(UL(IP,K))*SLU(IP)
         IF(SIZE.LE.BIG) GOTO    11
            BIG = SIZE
            IDXPIV = I
   11    CONTINUE
         IF(BIG.EQ.0.0)  GOTO 20
         IF(IDXPIV-K) 14,15,14
   14    J = IPS(K)
         IPS(K)  = IPS(IDXPIV)
         IPS(IDXPIV) = J
   15    KP = IPS(K)
         PIVOT = UL(KP,K)
         KP1 = K+1
         DO 17 I=KP1,N
            IP = IPS(I)
            IF(UL(IP,K).EQ.0.0) GOTO 17
            EM = -UL(IP,K)/PIVOT
            UL(IP,K) = -EM
            DO 16 J=KP1,N
   16          UL(IP,J) = UL(IP,J)+EM*UL(KP,J)
   17       CONTINUE
   18 CONTINUE
C             
      KP = IPS(N)
      IF(UL(KP,N).NE.0.0) RETURN
   20 LFLAG = 2
      RETURN
      END
C             
C      * * *   END DECOMP      * * *
C 
