C=======================================================================
      SUBROUTINE PDERIV(THETA,DEG,PP,Y,N,ITYPE)
C=======================================================================
      DIMENSION HH(3)
      DOUBLE PRECISION Y(1),PP(1),U(4)
      DOUBLE PRECISION C1,C2,C3,SU,PR 
C
C      THIS SUBROUTINE COMPUTES THE DERIVATIVE OF AN INTERPOLATION
C      POLYNOMIAL, ITYPE=1 FOR NEWTON INTERPOLATION AND ITYPE=2
C      FOR HERMITE INTERPOLATION (ONLY USED FOR NONSTIFF SUBSYSTEM).
C      PP CONTAINS THE DERIVATIVE AT TIME T+THETA*H, DIVIDED
C      DIFFERENCES (INPUT) ARE STORED IN THE VECTOR Y.
C
C      REVISED G SODERLIND 1980-05-29
C      DOUBLE PRECISION VERSION: 1980-10-22
C
      INTEGER DEG
      COMMON /STEPS/ HH
      K = MAX0(0,DEG)
      C1 = DBLE(THETA*HH(1))
      C2 = C1+DBLE(HH(2))
      C3 = C2+DBLE(HH(3))
      PR = C1*C2
      SU = C1 + C2
      IF(ITYPE.GT.1) GOTO 10
         K = MIN0(2,K)
         U(1) = SU
         U(2) = PR+C3*SU
         GOTO 20
   10 K = MIN0(4,K)
      U(1) = 2D0*C1
      U(2) = C1*(C2+SU)
      U(3) = 2D0*PR*SU
      U(4) = PR*(2D0*C3*SU+PR)
   20 DO 40 I=1,N
         SU=0D0
         IF(K.EQ.0) GOTO 40
         INEW = N
         DO 30 J=1,K
            SU = SU+U(J)*Y(I+INEW)
   30       INEW = INEW+N
   40    PP(I) = Y(I)+SU
      RETURN
      END
C
C      * * * END PDERIV * * *
C
