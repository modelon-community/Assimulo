      SUBROUTINE DEC (N, NDIM, A, IP, IER)
! VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
!-----------------------------------------------------------------------
!  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
!  INPUT..
!     N = ORDER OF MATRIX.
!     NDIM = DECLARED DIMENSION OF ARRAY  A .
!     A = MATRIX TO BE TRIANGULARIZED.
!  OUTPUT..
!    A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
!    A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
!    IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
!    IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
!    IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
!          SINGULAR AT STAGE K.
! USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
! DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
! IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
!
! REFERENCE..
!    C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!    C.A.C.M. 15 (1972), P. 274.
!-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (ABS(A(I,K)) .GT. ABS(A(M,K))) M = I  
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,N
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,N
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
!----------------------- END OF SUBROUTINE DEC -------------------------
      END
!
!
      SUBROUTINE SOL (N, NDIM, A, B, IP)
!VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
!-----------------------------------------------------------------------
! SOLUTION OF LINEAR SYSTEM, A*X = B .
! INPUT..
!   N = ORDER OF MATRIX.
!   NDIM = DECLARED DIMENSION OF ARRAY  A .
!   A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
!   B = RIGHT HAND SIDE VECTOR.
!   IP = PIVOT VECTOR OBTAINED FROM DEC.
! DO NOT USE IF DEC HAS SET IER .NE. 0.
! OUTPUT..
!   B = SOLUTION VECTOR, X .
!----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        DO 10 I = KP1,N
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
!---------------------- END OF SUBROUTINE SOL -------------------------
      END
!
!
      SUBROUTINE DECH (N, NDIM, A, LB, IP, IER)
!VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J,LB,NA
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
!----------------------------------------------------------------------
! MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A HESSENBERG
! MATRIX WITH LOWER BANDWIDTH LB
! INPUT..
!    N = ORDER OF MATRIX A.
!    NDIM = DECLARED DIMENSION OF ARRAY  A .
!    A = MATRIX TO BE TRIANGULARIZED.
!    LB = LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED, LB.GE.1).
! OUTPUT..
!    A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
!    A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
!    IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
!    IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
!    IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
!          SINGULAR AT STAGE K.
! USE  SOLH  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
! DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
! IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
!
! REFERENCE..
!    THIS IS A SLIGHT MODIFICATION OF
!    C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!    C.A.C.M. 15 (1972), P. 274.
!----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        NA = MIN0(N,LB+K)
        DO 10 I = KP1,NA
          IF (ABS(A(I,K)) .GT. ABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        T = A(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(K,K)
        A(K,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = KP1,NA
 30       A(I,K) = -A(I,K)*T
        DO 50 J = KP1,N
          T = A(M,J)
          A(M,J) = A(K,J)
          A(K,J) = T
          IF (T .EQ. 0.D0) GO TO 45
          DO 40 I = KP1,NA
 40         A(I,J) = A(I,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(N,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
!---------------------- END OF SUBROUTINE DECH ------------------------
      END
!
!
      SUBROUTINE SOLH (N, NDIM, A, LB, B, IP)
!VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1,LB,NA
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
!----------------------------------------------------------------------
! SOLUTION OF LINEAR SYSTEM, A*X = B .
! INPUT..
!   N = ORDER OF MATRIX A.
!   NDIM = DECLARED DIMENSION OF ARRAY  A .
!   A = TRIANGULARIZED MATRIX OBTAINED FROM DECH.
!   LB = LOWER BANDWIDTH OF A.
!   B = RIGHT HAND SIDE VECTOR.
!   IP = PIVOT VECTOR OBTAINED FROM DEC.
! DO NOT USE IF DECH HAS SET IER .NE. 0.
! OUTPUT..
!   B = SOLUTION VECTOR, X .
!----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        NA = MIN0(N,LB+K)
        DO 10 I = KP1,NA
 10       B(I) = B(I) + A(I,K)*T
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        B(K) = B(K)/A(K,K)
        T = -B(K)
        DO 30 I = 1,KM1
 30       B(I) = B(I) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(1,1)
      RETURN
!---------------------- END OF SUBROUTINE SOLH ------------------------
      END
!
      SUBROUTINE DECC (N, NDIM, AR, AI, IP, IER)
!VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
!----------------------------------------------------------------------
! MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION
! ------ MODIFICATION FOR COMPLEX MATRICES --------
! INPUT..
!    N = ORDER OF MATRIX.
!    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI .
!    (AR, AI) = MATRIX TO BE TRIANGULARIZED.
! OUTPUT..
!    AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART.
!    AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART.
!    AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
!                                                   REAL PART.
!    AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
!                                                   IMAGINARY PART.
!    IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
!    IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
!    IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
!          SINGULAR AT STAGE K.
! USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
! IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
!
! REFERENCE..
!    C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!    C.A.C.M. 15 (1972), P. 274.
!----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (ABS(AR(I,K))+ABS(AI(I,K)) .GT. &
                ABS(AR(M,K))+ABS(AI(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        TR = AR(M,K)
        TI = AI(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        AR(M,K) = AR(K,K)
        AI(M,K) = AI(K,K)
        AR(K,K) = TR
        AI(K,K) = TI
 20     CONTINUE
        IF (ABS(TR)+ABS(TI) .EQ. 0.D0) GO TO 80
        DEN=TR*TR+TI*TI
        TR=TR/DEN
        TI=-TI/DEN
        DO 30 I = KP1,N
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          AR(I,K)=-PRODR
          AI(I,K)=-PRODI
 30       CONTINUE
        DO 50 J = KP1,N
          TR = AR(M,J)
          TI = AI(M,J)
          AR(M,J) = AR(K,J)
          AI(M,J) = AI(K,J)
          AR(K,J) = TR
          AI(K,J) = TI
          IF (ABS(TR)+ABS(TI) .EQ. 0.D0) GO TO 48
          IF (TI .EQ. 0.D0) THEN
            DO 40 I = KP1,N
            PRODR=AR(I,K)*TR
            PRODI=AI(I,K)*TR
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 40         CONTINUE
            GO TO 48
          END IF
          IF (TR .EQ. 0.D0) THEN
            DO 45 I = KP1,N
            PRODR=-AI(I,K)*TI
            PRODI=AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 45         CONTINUE
            GO TO 48
          END IF
          DO 47 I = KP1,N
            PRODR=AR(I,K)*TR-AI(I,K)*TI
            PRODI=AI(I,K)*TR+AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 47         CONTINUE
 48       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (ABS(AR(N,N))+ABS(AI(N,N)) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
!---------------------- END OF SUBROUTINE DECC ------------------------
      END
!
!
      SUBROUTINE SOLC (N, NDIM, AR, AI, BR, BI, IP)
!VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
!----------------------------------------------------------------------
! SOLUTION OF LINEAR SYSTEM, A*X = B .
! INPUT..
!   N = ORDER OF MATRIX.
!   NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI.
!   (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
!   (BR,BI) = RIGHT HAND SIDE VECTOR.
!   IP = PIVOT VECTOR OBTAINED FROM DEC.
! DO NOT USE IF DEC HAS SET IER .NE. 0.
! OUTPUT..
!   (BR,BI) = SOLUTION VECTOR, X .
!----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        TR = BR(M)
        TI = BI(M)
        BR(M) = BR(K)
        BI(M) = BI(K)
        BR(K) = TR
        BI(K) = TI
        DO 10 I = KP1,N
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 10       CONTINUE
 20     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        DEN=AR(K,K)*AR(K,K)+AI(K,K)*AI(K,K)
        PRODR=BR(K)*AR(K,K)+BI(K)*AI(K,K)
        PRODI=BI(K)*AR(K,K)-BR(K)*AI(K,K)
        BR(K)=PRODR/DEN
        BI(K)=PRODI/DEN
        TR = -BR(K)
        TI = -BI(K)
        DO 30 I = 1,KM1
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 30       CONTINUE
 40     CONTINUE
 50     CONTINUE
        DEN=AR(1,1)*AR(1,1)+AI(1,1)*AI(1,1)
        PRODR=BR(1)*AR(1,1)+BI(1)*AI(1,1)
        PRODI=BI(1)*AR(1,1)-BR(1)*AI(1,1)
        BR(1)=PRODR/DEN
        BI(1)=PRODI/DEN
      RETURN
!---------------------- END OF SUBROUTINE SOLC ------------------------
      END  
!
!
      SUBROUTINE DECHC (N, NDIM, AR, AI, LB, IP, IER)
!VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
!----------------------------------------------------------------------
! MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION
! ------ MODIFICATION FOR COMPLEX MATRICES --------
! INPUT..
!    N = ORDER OF MATRIX.
!    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI .
!    (AR, AI) = MATRIX TO BE TRIANGULARIZED.
! OUTPUT..
!    AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART.
!    AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART.
!    AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
!                                                   REAL PART.
!    AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
!                                                   IMAGINARY PART.
!    LB = LOWER BANDWIDTH OF A (DIAGONAL NOT COUNTED), LB.GE.1.
!    IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
!    IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
!    IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
!          SINGULAR AT STAGE K.
! USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
! IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
!
! REFERENCE..
!    C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!    C.A.C.M. 15 (1972), P. 274.
!----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (LB .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K 
        NA = MIN0(N,LB+K)
        DO 10 I = KP1,NA
          IF (ABS(AR(I,K))+ABS(AI(I,K)) .GT. &
                ABS(AR(M,K))+ABS(AI(M,K))) M = I
 10     CONTINUE
        IP(K) = M
        TR = AR(M,K)
        TI = AI(M,K)
        IF (M .EQ. K) GO TO 20
        IP(N) = -IP(N)
        AR(M,K) = AR(K,K)
        AI(M,K) = AI(K,K)
        AR(K,K) = TR
        AI(K,K) = TI
 20     CONTINUE
        IF (ABS(TR)+ABS(TI) .EQ. 0.D0) GO TO 80
        DEN=TR*TR+TI*TI
        TR=TR/DEN
        TI=-TI/DEN
        DO 30 I = KP1,NA
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          AR(I,K)=-PRODR
          AI(I,K)=-PRODI
 30       CONTINUE
        DO 50 J = KP1,N
          TR = AR(M,J)
          TI = AI(M,J)
          AR(M,J) = AR(K,J)
          AI(M,J) = AI(K,J)
          AR(K,J) = TR
          AI(K,J) = TI
          IF (ABS(TR)+ABS(TI) .EQ. 0.D0) GO TO 48
          IF (TI .EQ. 0.D0) THEN
            DO 40 I = KP1,NA
            PRODR=AR(I,K)*TR
            PRODI=AI(I,K)*TR
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 40         CONTINUE
            GO TO 48
          END IF
          IF (TR .EQ. 0.D0) THEN
            DO 45 I = KP1,NA
            PRODR=-AI(I,K)*TI
            PRODI=AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 45         CONTINUE
            GO TO 48
          END IF
          DO 47 I = KP1,NA
            PRODR=AR(I,K)*TR-AI(I,K)*TI
            PRODI=AI(I,K)*TR+AR(I,K)*TI
            AR(I,J) = AR(I,J) + PRODR
            AI(I,J) = AI(I,J) + PRODI
 47         CONTINUE
 48       CONTINUE
 50       CONTINUE
 60     CONTINUE
 70   K = N
      IF (ABS(AR(N,N))+ABS(AI(N,N)) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
!---------------------- END OF SUBROUTINE DECHC -----------------------
      END
!
!
      SUBROUTINE SOLHC (N, NDIM, AR, AI, LB, BR, BI, IP)
!VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
!----------------------------------------------------------------------
! SOLUTION OF LINEAR SYSTEM, A*X = B .
! INPUT..
!   N = ORDER OF MATRIX.
!   NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI.
!   (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
!   (BR,BI) = RIGHT HAND SIDE VECTOR.
!   LB = LOWER BANDWIDTH OF A.
!   IP = PIVOT VECTOR OBTAINED FROM DEC.
! DO NOT USE IF DEC HAS SET IER .NE. 0.
! OUTPUT..
!   (BR,BI) = SOLUTION VECTOR, X .
!----------------------------------------------------------------------
      IF (N .EQ. 1) GO TO 50
      NM1 = N - 1
      IF (LB .EQ. 0) GO TO 25
      DO 20 K = 1,NM1
        KP1 = K + 1
        M = IP(K)
        TR = BR(M)
        TI = BI(M)
        BR(M) = BR(K)
        BI(M) = BI(K)
        BR(K) = TR
        BI(K) = TI
        DO 10 I = KP1,MIN0(N,LB+K)
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 10       CONTINUE
 20     CONTINUE
 25     CONTINUE
      DO 40 KB = 1,NM1
        KM1 = N - KB
        K = KM1 + 1
        DEN=AR(K,K)*AR(K,K)+AI(K,K)*AI(K,K)
        PRODR=BR(K)*AR(K,K)+BI(K)*AI(K,K)
        PRODI=BI(K)*AR(K,K)-BR(K)*AI(K,K)
        BR(K)=PRODR/DEN
        BI(K)=PRODI/DEN
        TR = -BR(K)
        TI = -BI(K)
        DO 30 I = 1,KM1
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(I) = BR(I) + PRODR
          BI(I) = BI(I) + PRODI
 30       CONTINUE
 40     CONTINUE
 50     CONTINUE
        DEN=AR(1,1)*AR(1,1)+AI(1,1)*AI(1,1)
        PRODR=BR(1)*AR(1,1)+BI(1)*AI(1,1)
        PRODI=BI(1)*AR(1,1)-BR(1)*AI(1,1)
        BR(1)=PRODR/DEN
        BI(1)=PRODI/DEN
      RETURN
!---------------------- END OF SUBROUTINE SOLHC -----------------------
      END  
!
      SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER)
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
!----------------------------------------------------------------------
! MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
! MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
! INPUT..
!    N       ORDER OF THE ORIGINAL MATRIX A.
!    NDIM    DECLARED DIMENSION OF ARRAY  A.
!    A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS  
!               OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
!               THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS 
!               ML+1 THROUGH 2*ML+MU+1 OF  A.
!    ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!    MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
! OUTPUT..
!    A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND 
!               THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.  
!    IP      INDEX VECTOR OF PIVOT INDICES.
!    IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
!    IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
!               SINGULAR AT STAGE K.
! USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
! DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
! IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
!
! REFERENCE..
!    THIS IS A MODIFICATION OF
!    C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!    C.A.C.M. 15 (1972), P. 274.
!----------------------------------------------------------------------
      IER = 0
      IP(N) = 1 
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
  5   A(I,J) = 0.D0
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (ABS(A(I,K)) .GT. ABS(A(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        T = A(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        A(M,K) = A(MD,K)
        A(MD,K) = T
 20     CONTINUE
        IF (T .EQ. 0.D0) GO TO 80
        T = 1.D0/T
        DO 30 I = MD1,MDL
 30       A(I,K) = -A(I,K)*T 
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          T = A(M,J) 
          IF (M .EQ. MM) GO TO 35
          A(M,J) = A(MM,J)
          A(MM,J) = T
 35       CONTINUE
          IF (T .EQ. 0.D0) GO TO 45
          JK = J - K
          DO 40 I = MD1,MDL
            IJK = I - JK
 40         A(IJK,J) = A(IJK,J) + A(I,K)*T
 45       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (A(MD,N) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
!---------------------- END OF SUBROUTINE DECB ------------------------
      END
!
!
      SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP)
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
!----------------------------------------------------------------------
! SOLUTION OF LINEAR SYSTEM, A*X = B .
! INPUT..
!   N      ORDER OF MATRIX A.
!   NDIM   DECLARED DIMENSION OF ARRAY  A .
!   A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
!   ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!   MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!   B      RIGHT HAND SIDE VECTOR.
!   IP     PIVOT VECTOR OBTAINED FROM DECB.
! DO NOT USE IF DECB HAS SET IER .NE. 0.
! OUTPUT..
!   B      SOLUTION VECTOR, X .
!----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        T = B(M)
        B(M) = B(K)
        B(K) = T
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
 10       B(IMD) = B(IMD) + A(I,K)*T
 20     CONTINUE
 25   CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        B(K) = B(K)/A(MD,K)
        T = -B(K) 
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
 30       B(IMD) = B(IMD) + A(I,K)*T
 40     CONTINUE
 50   B(1) = B(1)/A(MD,1)
      RETURN
!---------------------- END OF SUBROUTINE SOLB ------------------------
      END
!
      SUBROUTINE DECBC (N, NDIM, AR, AI, ML, MU, IP, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
!----------------------------------------------------------------------
! MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED COMPLEX
! MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
! INPUT..
!    N       ORDER OF THE ORIGINAL MATRIX A.
!    NDIM    DECLARED DIMENSION OF ARRAY  A.
!    AR, AI     CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS  
!               OF THE MATRIX ARE STORED IN THE COLUMNS OF  AR (REAL
!               PART) AND AI (IMAGINARY PART)  AND
!               THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS 
!               ML+1 THROUGH 2*ML+MU+1 OF  AR AND AI.
!    ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!    MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
! OUTPUT..
!    AR, AI  AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND 
!               THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.  
!    IP      INDEX VECTOR OF PIVOT INDICES.
!    IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
!    IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
!               SINGULAR AT STAGE K.
! USE  SOLBC  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
! DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
! IF IP(N)=O, A IS SINGULAR, SOLBC WILL DIVIDE BY ZERO.
!
! REFERENCE..
!    THIS IS A MODIFICATION OF
!    C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
!    C.A.C.M. 15 (1972), P. 274.
!----------------------------------------------------------------------
      IER = 0
      IP(N) = 1 
      MD = ML + MU + 1
      MD1 = MD + 1
      JU = 0
      IF (ML .EQ. 0) GO TO 70
      IF (N .EQ. 1) GO TO 70
      IF (N .LT. MU+2) GO TO 7
      DO 5 J = MU+2,N
      DO 5 I = 1,ML
      AR(I,J) = 0.D0
      AI(I,J) = 0.D0
  5   CONTINUE
  7   NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = MD
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IF (ABS(AR(I,K))+ABS(AI(I,K)) .GT. &
                ABS(AR(M,K))+ABS(AI(M,K))) M = I
 10     CONTINUE
        IP(K) = M + K - MD
        TR = AR(M,K)
        TI = AI(M,K)
        IF (M .EQ. MD) GO TO 20
        IP(N) = -IP(N)
        AR(M,K) = AR(MD,K)
        AI(M,K) = AI(MD,K)
        AR(MD,K) = TR
        AI(MD,K) = TI
 20     IF (ABS(TR)+ABS(TI) .EQ. 0.D0) GO TO 80
        DEN=TR*TR+TI*TI
        TR=TR/DEN
        TI=-TI/DEN
        DO 30 I = MD1,MDL
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          AR(I,K)=-PRODR
          AI(I,K)=-PRODI
 30       CONTINUE
        JU = MIN0(MAX0(JU,MU+IP(K)),N)
        MM = MD
        IF (JU .LT. KP1) GO TO 55
        DO 50 J = KP1,JU
          M = M - 1
          MM = MM - 1
          TR = AR(M,J)
          TI = AI(M,J)
          IF (M .EQ. MM) GO TO 35
          AR(M,J) = AR(MM,J)
          AI(M,J) = AI(MM,J)
          AR(MM,J) = TR
          AI(MM,J) = TI
 35       CONTINUE
          IF (ABS(TR)+ABS(TI) .EQ. 0.D0) GO TO 48
          JK = J - K
          IF (TI .EQ. 0.D0) THEN
            DO 40 I = MD1,MDL
            IJK = I - JK
            PRODR=AR(I,K)*TR
            PRODI=AI(I,K)*TR
            AR(IJK,J) = AR(IJK,J) + PRODR
            AI(IJK,J) = AI(IJK,J) + PRODI
 40         CONTINUE
            GO TO 48
          END IF
          IF (TR .EQ. 0.D0) THEN
            DO 45 I = MD1,MDL
            IJK = I - JK
            PRODR=-AI(I,K)*TI
            PRODI=AR(I,K)*TI
            AR(IJK,J) = AR(IJK,J) + PRODR
            AI(IJK,J) = AI(IJK,J) + PRODI
 45         CONTINUE
            GO TO 48
          END IF
          DO 47 I = MD1,MDL
            IJK = I - JK
            PRODR=AR(I,K)*TR-AI(I,K)*TI
            PRODI=AI(I,K)*TR+AR(I,K)*TI
            AR(IJK,J) = AR(IJK,J) + PRODR
            AI(IJK,J) = AI(IJK,J) + PRODI
 47         CONTINUE
 48       CONTINUE
 50       CONTINUE
 55     CONTINUE
 60     CONTINUE
 70   K = N
      IF (ABS(AR(MD,N))+ABS(AI(MD,N)) .EQ. 0.D0) GO TO 80
      RETURN
 80   IER = K
      IP(N) = 0
      RETURN
!---------------------- END OF SUBROUTINE DECBC ------------------------
      END
!
!
      SUBROUTINE SOLBC (N, NDIM, AR, AI, ML, MU, BR, BI, IP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
!----------------------------------------------------------------------
! SOLUTION OF LINEAR SYSTEM, A*X = B ,
!                 VERSION BANDED AND COMPLEX-DOUBLE PRECISION.
! INPUT..
!   N      ORDER OF MATRIX A.
!   NDIM   DECLARED DIMENSION OF ARRAY  A .
!   AR, AI TRIANGULARIZED MATRIX OBTAINED FROM DECB (REAL AND IMAG. PART).
!   ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!   MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
!   BR, BI RIGHT HAND SIDE VECTOR (REAL AND IMAG. PART).
!   IP     PIVOT VECTOR OBTAINED FROM DECBC.
! DO NOT USE IF DECB HAS SET IER .NE. 0.
! OUTPUT..
!   BR, BI SOLUTION VECTOR, X (REAL AND IMAG. PART).
!----------------------------------------------------------------------
      MD = ML + MU + 1
      MD1 = MD + 1
      MDM = MD - 1
      NM1 = N - 1
      IF (ML .EQ. 0) GO TO 25
      IF (N .EQ. 1) GO TO 50
      DO 20 K = 1,NM1
        M = IP(K)
        TR = BR(M)
        TI = BI(M)
        BR(M) = BR(K)
        BI(M) = BI(K)
        BR(K) = TR
        BI(K) = TI
        MDL = MIN(ML,N-K) + MD
        DO 10 I = MD1,MDL
          IMD = I + K - MD
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(IMD) = BR(IMD) + PRODR
          BI(IMD) = BI(IMD) + PRODI
 10     CONTINUE
 20     CONTINUE
 25     CONTINUE
      DO 40 KB = 1,NM1
        K = N + 1 - KB
        DEN=AR(MD,K)*AR(MD,K)+AI(MD,K)*AI(MD,K)
        PRODR=BR(K)*AR(MD,K)+BI(K)*AI(MD,K)
        PRODI=BI(K)*AR(MD,K)-BR(K)*AI(MD,K)
        BR(K)=PRODR/DEN
        BI(K)=PRODI/DEN
        TR = -BR(K)
        TI = -BI(K)
        KMD = MD - K
        LM = MAX0(1,KMD+1)
        DO 30 I = LM,MDM
          IMD = I - KMD
          PRODR=AR(I,K)*TR-AI(I,K)*TI
          PRODI=AI(I,K)*TR+AR(I,K)*TI
          BR(IMD) = BR(IMD) + PRODR
          BI(IMD) = BI(IMD) + PRODI
 30       CONTINUE
 40     CONTINUE
        DEN=AR(MD,1)*AR(MD,1)+AI(MD,1)*AI(MD,1)
        PRODR=BR(1)*AR(MD,1)+BI(1)*AI(MD,1)
        PRODI=BI(1)*AR(MD,1)-BR(1)*AI(MD,1)
        BR(1)=PRODR/DEN
        BI(1)=PRODI/DEN
 50   CONTINUE
      RETURN
!---------------------- END OF SUBROUTINE SOLBC ------------------------
      END
!
!
      subroutine elmhes(nm,n,low,igh,a,int)
!
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real*8 a(nm,n)
      real*8 x,y
!    real*8 abs
      integer int(igh)
!
!    this subroutine is a translation of the algol procedure elmhes,
!    num. math. 12, 349-368(1968) by martin and wilkinson.
!    handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
!
!    given a real general matrix, this subroutine
!    reduces a submatrix situated in rows and columns
!    low through igh to upper hessenberg form by
!    stabilized elementary similarity transformations.
!
!    on input:
!
!     nm must be set to the row dimension of two-dimensional
!       array parameters as declared in the calling program
!       dimension statement;
!
!     n is the order of the matrix;
!
!     low and igh are integers determined by the balancing
!       subroutine  balanc.      if  balan! has not been used,
!       set low=1, igh=n;
!
!     a contains the input matrix.
!
!    on output:
!
!     a contains the hessenberg matrix.  the multipliers
!       which were used in the reduction are stored in the
!       remaining triangle under the hessenberg matrix;
!
!     int contains information on the rows and columns
!       interchanged in the reduction.
!       only elements low through igh are used.
!
!    questions and comments should be directed to b. s. garbow,
!    applied mathematics division, argonne national laboratory
!
!    ------------------------------------------------------------------
!
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
!
      do 180 m = kp1, la
       mm1 = m - 1
       x = 0.0d0
       i = m
!
       do 100 j = m, igh
          if (abs(a(j,mm1)) .le. abs(x)) go to 100
          x = a(j,mm1)
          i = j
  100   continue
!
       int(m) = i
       if (i .eq. m) go to 130
!   :::::::::: interchange rows and columns of a ::::::::::
       do 110 j = mm1, n
          y = a(i,j)
          a(i,j) = a(m,j)
          a(m,j) = y
  110   continue
!
       do 120 j = 1, igh
          y = a(j,i)
          a(j,i) = a(j,m)
          a(j,m) = y
  120   continue
!   :::::::::: end interchange ::::::::::
  130   if (x .eq. 0.0d0) go to 180
       mp1 = m + 1

       do 160 i = mp1, igh
          y = a(i,mm1)
          if (y .eq. 0.0d0) go to 160
          y = y / x
          a(i,mm1) = y
!
          do 140 j = m, n
  140      a(i,j) = a(i,j) - y * a(m,j)
!
          do 150 j = 1, igh
  150      a(j,m) = a(j,m) + y * a(j,i)
!
  160   continue
!
  180 continue
!
  200 return
!   :::::::::: last card of elmhes ::::::::::
      end

