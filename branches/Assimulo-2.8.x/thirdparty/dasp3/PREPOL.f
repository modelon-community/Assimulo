C=======================================================================             
      SUBROUTINE PREPOL(THETA,DEGREE,P,Z,M,ITYPE)
C=======================================================================
      INTEGER DEGREE, DEG
      COMMON /STEPS/ HH
      DIMENSION HH(3),HCOF(3)
      DOUBLE PRECISION P(1),Z(1),U
C             
C      THIS SUBROUTINE COMPUTES THE INTERPOLATION POLYNOMIAL
C      AT TIME T+THETA*H. ITYPE=1 FOR NEWTON INTERPOLATION
C      AND ITYPE=2 FOR HERMITE INTERPOLATION (ONLY USED FOR
C      THE NONSTIFF SUBSYSTEM). DIVIDED DIFFERENCES (INPUT)
C      ARE STORED IN THE VECTOR Z. OUTPUT IS CONTAINED TN P.
C             
C      REVISED G SODERLIND 1980-05-29
C      DOUBLE PRECISION VERSION: 1980-10-22
C   
      DEG = MAX0(DEGREE,0)
      HCOF(1) = THETA*HH(1)
      IF(DEG-1) 10,4,10
   10 DO 1 I=2,3
    1    HCOF(I) = HCOF(I-1)+HH(I)
    4 INEW = DEG*M
      DO 2 I=1,M
    2      P(I) = Z(I+INEW)
C             
      DO 3 K=1,DEG
           KK = DEG-K
           IF(ITYPE.GE.2) KK = KK/2
           U = HCOF(KK+1)
           INEW = INEW-M
           DO 3 I=1,M
    3         P(I) = U*P(I)+Z(I+INEW)
      RETURN
      END
C             
C      * * *   END PREPOL      * * *
C   
