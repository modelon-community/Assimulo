C ******************************************
C     VERSION OF JUNE 22, 2000      
C ******************************************
C
      SUBROUTINE DECOMR(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,FAC1,E1,LDE1,IP1,IER,IJOB,CALHES,IPHES)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, PARAMETER :: DP=kind(1D0)
      REAL(kind=DP), dimension(LDJAC,N), intent(in) :: FJAC 
      REAL(kind=DP), dimension(LDMAS,NM1), intent(in) :: FMAS
      REAL(kind=DP), dimension(LDE1,NM1), intent(out) :: E1
      INTEGER, dimension(NM1), intent(out) :: IP1
      INTEGER, dimension(N), intent(out) :: IPHES

      LOGICAL CALHES
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO  I=1,N
            E1(I,J)=-FJAC(I,J)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=-FJAC(I,JM1)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
 45   MM=M1/M2
      DO J=1,M2
         DO I=1,NM1
            SUM=0.D0
            DO K=0,MM-1
               SUM=(SUM+FJAC(I,J+K*M2))/FAC1
            END DO
            E1(I,J)=E1(I,J)-SUM
         END DO
      END DO
      CALL DEC (NM1,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,J)
         END DO
         E1(MDIAG,J)=E1(MDIAG,J)+FAC1
      END DO
      CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,JM1)
         END DO
         E1(MDIAG,J)=E1(MDIAG,J)+FAC1
      END DO
  46  MM=M1/M2
      DO J=1,M2
         DO I=1,MBJAC
            SUM=0.D0
            DO K=0,MM-1
               SUM=(SUM+FJAC(I,J+K*M2))/FAC1
            END DO
            E1(I+MLE,J)=E1(I+MLE,J)-SUM
         END DO
      END DO
      CALL DECB (NM1,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E1(I,J)=-FJAC(I,J)
         END DO
         DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS)
            E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J)
         END DO
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=-FJAC(I,JM1)
         END DO
         DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS)
            E1(I,J)=E1(I,J)+FAC1*FMAS(I-J+MBDIAG,J)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,J)
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J)
         END DO
      END DO
      CALL DECB (N,LDE1,E1,MLE,MUE,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  14  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E1(I+MLE,J)=-FJAC(I,JM1)
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            E1(IB,J)=E1(IB,J)+FAC1*FMAS(I,J)
         END DO
      END DO
      GOTO 46
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,J)
         END DO
      END DO
      CALL DEC (N,LDE1,E1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E1(I,J)=FMAS(I,J)*FAC1-FJAC(I,JM1)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      IF (CALHES) CALL ELMHES (LDJAC,N,1,N,FJAC,IPHES) 
      CALHES=.FALSE.
      DO J=1,N-1
         J1=J+1
         E1(J1,J)=-FJAC(J1,J)
      END DO
      DO J=1,N
         DO I=1,J
            E1(I,J)=-FJAC(I,J)
         END DO
         E1(J,J)=E1(J,J)+FAC1
      END DO
      CALL DECH(N,LDE1,E1,1,IP1,IER)
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE DECOMR
C
C ***********************************************************
C
      SUBROUTINE DECOMC(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &            M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,IP2,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, PARAMETER :: DP=kind(1D0)
      REAL(kind=DP), dimension(LDJAC,N), intent(in) :: FJAC
      REAL(kind=DP), dimension(LDMAS,NM1), intent(in) :: FMAS
      REAL(kind=DP), dimension(LDE1,NM1), intent(out) :: E2R
      REAL(kind=DP), dimension(LDE1,NM1), intent(out) :: E2I
      INTEGER, dimension(NM1), intent(out) :: IP2

      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            E2R(I,J)=-FJAC(I,J)
            E2I(I,J)=0.D0
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
      CALL DECC (N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=-FJAC(I,JM1)
            E2I(I,J)=0.D0
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
  45  MM=M1/M2
      ABNO=ALPHN**2+BETAN**2
      ALP=ALPHN/ABNO
      BET=BETAN/ABNO
      DO J=1,M2
         DO I=1,NM1
            SUMR=0.D0
            SUMI=0.D0
            DO K=0,MM-1
               SUMS=SUMR+FJAC(I,J+K*M2)
               SUMR=SUMS*ALP+SUMI*BET
               SUMI=SUMI*ALP-SUMS*BET
            END DO
            E2R(I,J)=E2R(I,J)-SUMR
            E2I(I,J)=E2I(I,J)-SUMI
         END DO
      END DO
      CALL DECC (NM1,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            IMLE=I+MLE
            E2R(IMLE,J)=-FJAC(I,J)
            E2I(IMLE,J)=0.D0
         END DO
         E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN
         E2I(MDIAG,J)=BETAN
      END DO
      CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E2R(I+MLE,J)=-FJAC(I,JM1)
            E2I(I+MLE,J)=0.D0
         END DO
         E2R(MDIAG,J)=E2R(MDIAG,J)+ALPHN
         E2I(MDIAG,J)=E2I(MDIAG,J)+BETAN
      END DO
  46  MM=M1/M2
      ABNO=ALPHN**2+BETAN**2
      ALP=ALPHN/ABNO
      BET=BETAN/ABNO
      DO J=1,M2
         DO I=1,MBJAC
            SUMR=0.D0
            SUMI=0.D0
            DO K=0,MM-1
               SUMS=SUMR+FJAC(I,J+K*M2)
               SUMR=SUMS*ALP+SUMI*BET
               SUMI=SUMI*ALP-SUMS*BET
            END DO
            IMLE=I+MLE
            E2R(IMLE,J)=E2R(IMLE,J)-SUMR
            E2I(IMLE,J)=E2I(IMLE,J)-SUMI
         END DO
      END DO
      CALL DECBC (NM1,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO  J=1,N
         DO  I=1,N
            E2R(I,J)=-FJAC(I,J)
            E2I(I,J)=0.D0
         END DO
      END DO
      DO J=1,N
         DO I=MAX(1,J-MUMAS),MIN(N,J+MLMAS)
            BB=FMAS(I-J+MBDIAG,J)
            E2R(I,J)=E2R(I,J)+ALPHN*BB
            E2I(I,J)=BETAN*BB
         END DO
      END DO
      CALL DECC(N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=-FJAC(I,JM1)
            E2I(I,J)=0.D0
         END DO
         DO I=MAX(1,J-MUMAS),MIN(NM1,J+MLMAS)
            FFMA=FMAS(I-J+MBDIAG,J)
            E2R(I,J)=E2R(I,J)+ALPHN*FFMA
            E2I(I,J)=E2I(I,J)+BETAN*FFMA
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO J=1,N
         DO I=1,MBJAC
            IMLE=I+MLE
            E2R(IMLE,J)=-FJAC(I,J)
            E2I(IMLE,J)=0.D0
         END DO
         DO I=MAX(1,MUMAS+2-J),MIN(MBB,MUMAS+1-J+N)
            IB=I+MDIFF
            BB=FMAS(I,J)
            E2R(IB,J)=E2R(IB,J)+ALPHN*BB
            E2I(IB,J)=BETAN*BB
         END DO
      END DO
      CALL DECBC (N,LDE1,E2R,E2I,MLE,MUE,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  14  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,MBJAC
            E2R(I+MLE,J)=-FJAC(I,JM1)
            E2I(I+MLE,J)=0.D0
         END DO
         DO I=1,MBB
            IB=I+MDIFF
            FFMA=FMAS(I,J)
            E2R(IB,J)=E2R(IB,J)+ALPHN*FFMA
            E2I(IB,J)=E2I(IB,J)+BETAN*FFMA
         END DO
      END DO
      GOTO 46
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO J=1,N
         DO I=1,N
            BB=FMAS(I,J)
            E2R(I,J)=BB*ALPHN-FJAC(I,J)
            E2I(I,J)=BB*BETAN
         END DO
      END DO
      CALL DECC(N,LDE1,E2R,E2I,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO J=1,NM1
         JM1=J+M1
         DO I=1,NM1
            E2R(I,J)=ALPHN*FMAS(I,J)-FJAC(I,JM1)
            E2I(I,J)=BETAN*FMAS(I,J)
         END DO
      END DO
      GOTO 45
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO J=1,N-1
         J1=J+1
         E2R(J1,J)=-FJAC(J1,J)
         E2I(J1,J)=0.D0
      END DO
      DO J=1,N
         DO I=1,J
            E2I(I,J)=0.D0
            E2R(I,J)=-FJAC(I,J)
         END DO
         E2R(J,J)=E2R(J,J)+ALPHN
         E2I(J,J)=BETAN
      END DO
      CALL DECHC(N,LDE1,E2R,E2I,1,IP2,IER)
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE DECOMC
C
C ***********************************************************
C
C
      SUBROUTINE SLVRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &          M1,M2,NM1,FAC1,ALPHN,BETAN,E1,E2R,E2I,LDE1,Z1,Z2,Z3,
     &          F1,F2,F3,CONT,IP1,IP2,IPHES,IER,IJOB)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, PARAMETER :: DP=kind(1D0)
      REAL(kind=DP), dimension(LDJAC,N), intent(in) :: FJAC
      REAL(kind=DP), dimension(LDMAS,NM1), intent(in) :: FMAS
      REAL(kind=DP), dimension(LDE1,NM1), intent(in) :: E1
      REAL(kind=DP), dimension(LDE1,NM1), intent(in) :: E2R
      REAL(kind=DP), dimension(LDE1,NM1), intent(in) :: E2I
      INTEGER, dimension(NM1), intent(in) :: IP1
      INTEGER, dimension(NM1), intent(in) :: IP2
      INTEGER, dimension(N), intent(in) :: IPHES
      REAL(kind=DP), dimension(N), intent(inout) :: Z1,Z2,Z3,F1,F2,F3

      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,13,15), IJOB
C
C -----------------------------------------------------------
C
   1  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC (N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  11  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
 48   ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=1,NM1
               IM1=I+M1
               Z1(IM1)=Z1(IM1)+FJAC(I,JKM)*SUM1
               Z2(IM1)=Z2(IM1)+FJAC(I,JKM)*SUM2
               Z3(IM1)=Z3(IM1)+FJAC(I,JKM)*SUM3
            END DO
         END DO
      END DO
      CALL SOL (NM1,LDE1,E1,Z1(M1+1),IP1)
      CALL SOLC (NM1,LDE1,E2R,E2I,Z2(M1+1),Z3(M1+1),IP2)
 49   CONTINUE
      DO I=M1,1,-1
         MPI=M2+I
         Z1(I)=(Z1(I)+Z1(MPI))/FAC1
         Z2I=Z2(I)+Z2(MPI)
         Z3I=Z3(I)+Z3(MPI)
         Z3(I)=(Z3I*ALPHN-Z2I*BETAN)/ABNO
         Z2(I)=(Z2I*ALPHN+Z3I*BETAN)/ABNO
      END DO
      RETURN
C
C -----------------------------------------------------------
C
   2  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      CALL SOLBC (N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  12  CONTINUE
C ---  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
  45  ABNO=ALPHN**2+BETAN**2
      MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUM2=0.D0
         SUM3=0.D0
         DO K=MM-1,0,-1
            JKM=J+K*M2
            SUM1=(Z1(JKM)+SUM1)/FAC1
            SUMH=(Z2(JKM)+SUM2)/ABNO
            SUM3=(Z3(JKM)+SUM3)/ABNO
            SUM2=SUMH*ALPHN+SUM3*BETAN
            SUM3=SUM3*ALPHN-SUMH*BETAN
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               FFJA=FJAC(I+MUJAC+1-J,JKM)
               Z1(IM1)=Z1(IM1)+FFJA*SUM1
               Z2(IM1)=Z2(IM1)+FFJA*SUM2
               Z3(IM1)=Z3(IM1)+FFJA*SUM3
            END DO
         END DO
      END DO
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,Z1(M1+1),IP1)
      CALL SOLBC (NM1,LDE1,E2R,E2I,MLE,MUE,Z2(M1+1),Z3(M1+1),IP2)
      GOTO 49
C
C -----------------------------------------------------------
C
   3  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  13  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         J1B=MAX(1,I-MLMAS)
         J2B=MIN(NM1,I+MUMAS)
         DO J=J1B,J2B
            JM1=J+M1
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(JM1)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      IF (IJOB.EQ.14) GOTO 45
      GOTO 48
C
C -----------------------------------------------------------
C
   4  CONTINUE
C ---  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            BB=FMAS(I-J+MBDIAG,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOLB (N,LDE1,E1,MLE,MUE,Z1,IP1)
      CALL SOLBC(N,LDE1,E2R,E2I,MLE,MUE,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
   5  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=1,N
            BB=FMAS(I,J)
            S1=S1-BB*F1(J)
            S2=S2-BB*F2(J)
            S3=S3-BB*F3(J)
         END DO
         Z1(I)=Z1(I)+S1*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      CALL SOL (N,LDE1,E1,Z1,IP1)
      CALL SOLC(N,LDE1,E2R,E2I,Z2,Z3,IP2)
      RETURN
C
C -----------------------------------------------------------
C
  15  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO I=1,NM1
         IM1=I+M1
         S1=0.0D0
         S2=0.0D0
         S3=0.0D0
         DO J=1,NM1
            JM1=J+M1
            BB=FMAS(I,J)
            S1=S1-BB*F1(JM1)
            S2=S2-BB*F2(JM1)
            S3=S3-BB*F3(JM1)
         END DO
         Z1(IM1)=Z1(IM1)+S1*FAC1
         Z2(IM1)=Z2(IM1)+S2*ALPHN-S3*BETAN
         Z3(IM1)=Z3(IM1)+S3*ALPHN+S2*BETAN
      END DO
      GOTO 48
C
C -----------------------------------------------------------
C
   6  CONTINUE
C ---  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ---  THIS OPTION IS NOT PROVIDED
      RETURN
C
C -----------------------------------------------------------
C
   7  CONTINUE
C ---  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N
         S2=-F2(I)
         S3=-F3(I)
         Z1(I)=Z1(I)-F1(I)*FAC1
         Z2(I)=Z2(I)+S2*ALPHN-S3*BETAN
         Z3(I)=Z3(I)+S3*ALPHN+S2*BETAN
      END DO
      DO MM=N-2,1,-1
          MP=N-MM
          MP1=MP-1
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 746
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE 
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 746      CONTINUE
          DO I=MP+1,N 
             E1IMP=FJAC(I,MP1)
             Z1(I)=Z1(I)-E1IMP*Z1(MP)
             Z2(I)=Z2(I)-E1IMP*Z2(MP)
             Z3(I)=Z3(I)-E1IMP*Z3(MP)
          END DO
       END DO
       CALL SOLH(N,LDE1,E1,1,Z1,IP1)
       CALL SOLHC(N,LDE1,E2R,E2I,1,Z2,Z3,IP2)
       DO MM=1,N-2
          MP=N-MM
          MP1=MP-1
          DO I=MP+1,N 
             E1IMP=FJAC(I,MP1)
             Z1(I)=Z1(I)+E1IMP*Z1(MP)
             Z2(I)=Z2(I)+E1IMP*Z2(MP)
             Z3(I)=Z3(I)+E1IMP*Z3(MP)
          END DO
          I=IPHES(MP)
          IF (I.EQ.MP) GOTO 750
          ZSAFE=Z1(MP)
          Z1(MP)=Z1(I)
          Z1(I)=ZSAFE
          ZSAFE=Z2(MP)
          Z2(MP)=Z2(I)
          Z2(I)=ZSAFE 
          ZSAFE=Z3(MP)
          Z3(MP)=Z3(I)
          Z3(I)=ZSAFE
 750      CONTINUE
      END DO
      RETURN
C
C -----------------------------------------------------------
C
  55  CONTINUE
      RETURN
      END
C
C     END OF SUBROUTINE SLVRAD
C
C ***********************************************************
C

      SUBROUTINE ESTRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS,
     &                  H,G0,DD1,DD2,DD3,CL1,CL3,CQ1,CQ2,CQ3,CERLQ,
     &                  FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,E1,LDE1,ALPHA,
     &                  Z1,Z2,Z3,CONT,F1,F2,F3,IP1,IPHES,SCAL,ERR,CERR,
     &                  FIRST,REJECT,FAC1,ARGLAG,PHI,RPAR,IPAR,
     &                  IOUT,PAST,IPAST,NRDS,JEFLAG,IEFLAG,LRPAST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER, PARAMETER :: DP=kind(1D0)
      REAL(kind=DP), dimension(LDJAC,N), intent(in) :: FJAC
      REAL(kind=DP), dimension(LDMAS,NM1), intent(in) :: FMAS
      REAL(kind=DP), dimension(LDE1,NM1), intent(in) :: E1
      INTEGER, dimension(NM1), intent(in) :: IP1
      INTEGER, dimension(N), intent(in) :: IPHES
      REAL(kind=DP), dimension(N), intent(in)  :: Z1,Z2,Z3
      REAL(kind=DP), dimension(N), intent(in)  :: Y0,Y
      REAL(kind=DP), dimension(N) :: F1,F2,F3
      REAL(kind=DP), dimension(N) :: CONT,SCAL
      REAL(kind=DP), dimension(LRPAST), intent(in)  :: PAST
      REAL(kind=DP), dimension(:), allocatable :: W1,W2,Q1,Q2
      INTEGER, dimension(NRDS), intent(in) :: IPAST
      REAL(kind=DP), dimension(1), intent(in)  :: RPAR
      INTEGER, dimension(1), intent(in) :: IPAR
      INTEGER :: LRPAST

      LOGICAL FIRST,REJECT,LEFT
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C --- COMMON/BPCOM/BPP,ILBP,LEFT
      EXTERNAL ARGLAG,PHI
C      WRITE (6,*) FIRST,REJECT,LEFT
      ALLOCATE (W1(N),W2(N),Q1(N),Q2(N))
      HEE1=DD1/H
      HEE2=DD2/H
      HEE3=DD3/H
      GOTO (1,2,3,4,5,6,7,55,55,55,11,12,13,14,15), IJOB
C
   1  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX
      DO  I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      F3=G0*H*CONT
      IF (ALPHA.NE.0.D0) THEN
       W2=W1/(G0*H)
       Q2=Q1/(G0*H)
       CALL SOL (N,LDE1,E1,W2,IP1) 
       CALL SOL (N,LDE1,E1,Q2,IP1) 
      END IF
      CALL SOL (N,LDE1,E1,CONT,IP1) 
      GOTO 77
C
  11  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      W2=W1
      Q2=Q1
  48  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUML=0.D0
         SUMQ=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            SUML=(W1(J+K*M2)+SUML)/FAC1
            SUMQ=(Q1(J+K*M2)+SUMQ)/FAC1
            DO I=1,NM1
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
               W2(IM1)=W1(IM1)+FJAC(I,J+K*M2)*SUML
               Q2(IM1)=Q1(IM1)+FJAC(I,J+K*M2)*SUMQ
            END DO
         END DO
      END DO
      F3=G0*H*CONT
      IF (ALPHA.NE.0.D0) THEN
       W2=W2/(G0*H)
       Q2=Q2/(G0*H)
       CALL SOL (NM1,LDE1,E1,W2(M1+1),IP1) 
       CALL SOL (NM1,LDE1,E1,Q2(M1+1),IP1) 
      END IF
      CALL SOL (NM1,LDE1,E1,CONT(M1+1),IP1) 
      DO I=M1,1,-1
        IF (ALPHA.NE.0.D0) THEN
         W2(I)=(W2(I)+W2(M2+I))/FAC1
         Q2(I)=(Q2(I)+Q2(M2+I))/FAC1
        END IF
        CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   2  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX
      DO I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      F3=G0*H*CONT
      IF (ALPHA.NE.0.D0) THEN
       W2=W1/(G0*H)
       Q2=Q1/(G0*H)
       CALL SOLB (N,LDE1,E1,MLE,MUE,W2,IP1)
       CALL SOLB (N,LDE1,E1,MLE,MUE,Q2,IP1)
      END IF
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  12  CONTINUE
C ------  B=IDENTITY, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      W2=W1
      Q2=Q1
      F3=G0*H*CONT
  45  MM=M1/M2
      DO J=1,M2
         SUM1=0.D0
         SUML1=0.D0
         SUMQ1=0.D0
         DO K=MM-1,0,-1
            SUM1=(CONT(J+K*M2)+SUM1)/FAC1
            SUML1=(W2(J+K*M2)+SUM1)/FAC1
            SUMQ1=(Q2(J+K*M2)+SUM1)/FAC1
            DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
               IM1=I+M1
               CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
               W2(IM1)=W2(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUML1
               Q2(IM1)=Q2(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUMQ1
            END DO
         END DO
      END DO
      IF (ALPHA.NE.0.D0) THEN
       W2=W2/(G0*H)
       Q2=Q2/(G0*H)
       CALL SOLB (NM1,LDE1,E1,MLE,MUE,W2(M1+1),IP1)
       CALL SOLB (NM1,LDE1,E1,MLE,MUE,Q2(M1+1),IP1)
       DO I=M1,1,-1
          W2(I)=(W2(I)+W2(M2+I))/FAC1
          Q2(I)=(Q2(I)+Q2(M2+I))/FAC1
       END DO
      END IF
      CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
      DO I=M1,1,-1
         CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
      END DO
      GOTO 77
C
   3  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         SUML=0.D0
         SUMQ=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J)
            SUML=SUML+FMAS(I-J+MBDIAG,J)*W1(J)
            SUMQ=SUMQ+FMAS(I-J+MBDIAG,J)*Q1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
         W2(I)=SUML
         Q2(I)=SUMQ
      END DO
      F3=G0*H*CONT
      IF (ALPHA.NE.0.D0) THEN
       W2=W2/(G0*H)
       Q2=Q2/(G0*H)
       CALL SOL (N,LDE1,E1,W2,IP1) 
       CALL SOL (N,LDE1,E1,Q2,IP1) 
      END IF
      CALL SOL (N,LDE1,E1,CONT,IP1) 
      GOTO 77
C
  13  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
         W2(I)=CL1*Z1(I)+CL3*Z3(I)
         Q2(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         SUML=0.D0
         SUMQ=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1)
            SUML=SUML+FMAS(I-J+MBDIAG,J)*W1(J+M1)
            SUMQ=SUMQ+FMAS(I-J+MBDIAG,J)*Q1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         W2(IM1)=SUML
         Q2(IM1)=SUMQ
         CONT(IM1)=SUM+Y0(IM1)
      END DO
      F3=G0*H*CONT
      GOTO 48
C
   4  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         SUML=0.D0
         SUMQ=0.D0
         DO J=MAX(1,I-MLMAS),MIN(N,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J)
            SUML=SUML+FMAS(I-J+MBDIAG,J)*W1(J)
            SUMQ=SUMQ+FMAS(I-J+MBDIAG,J)*Q1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
         W2(I)=SUML
         Q2(I)=SUMQ
      END DO
      IF (ALPHA.NE.0.D0) THEN
       W2=W2/(G0*H)
       Q2=Q2/(G0*H)
       CALL SOLB (N,LDE1,E1,MLE,MUE,W2,IP1)
       CALL SOLB (N,LDE1,E1,MLE,MUE,Q2,IP1)
      END IF
      F3=G0*H*CONT
      CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
      GOTO 77
C
  14  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A BANDED MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
         W2(I)=CL1*Z1(I)+CL3*Z3(I)
         Q2(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         SUML=0.D0
         SUMQ=0.D0
         DO J=MAX(1,I-MLMAS),MIN(NM1,I+MUMAS)
            SUM=SUM+FMAS(I-J+MBDIAG,J)*F1(J+M1)
            SUML=SUML+FMAS(I-J+MBDIAG,J)*W1(J+M1)
            SUMQ=SUMQ+FMAS(I-J+MBDIAG,J)*Q1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
         W2(IM1)=SUML
         Q2(IM1)=SUMQ
      END DO
      F3=G0*H*CONT
      GOTO 45
C
   5  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A FULL MATRIX
      DO I=1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      DO I=1,N
         SUM=0.D0
         SUML=0.D0
         SUMQ=0.D0
         DO J=1,N
            SUM=SUM+FMAS(I,J)*F1(J)
            SUML=SUML+FMAS(I,J)*W1(J)
            SUMQ=SUMQ+FMAS(I,J)*Q1(J)
         END DO
         F2(I)=SUM
         CONT(I)=SUM+Y0(I)
         W2(I)=SUML
         Q2(I)=SUMQ
      END DO
      F3=G0*H*CONT
      IF (ALPHA.NE.0.D0) THEN
       W2=W2/(G0*H)
       Q2=Q2/(G0*H)
       CALL SOL (N,LDE1,E1,W2,IP1) 
       CALL SOL (N,LDE1,E1,Q2,IP1) 
      END IF
      CALL SOL (N,LDE1,E1,CONT,IP1) 
      GOTO 77
C
  15  CONTINUE
C ------  B IS A BANDED MATRIX, JACOBIAN A FULL MATRIX, SECOND ORDER
      DO I=1,M1
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
         W2(I)=CL1*Z1(I)+CL3*Z3(I)
         Q2(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      DO I=M1+1,N
         F1(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      DO I=1,NM1
         SUM=0.D0
         SUML=0.D0
         SUMQ=0.D0
         DO J=1,NM1
            SUM=SUM+FMAS(I,J)*F1(J+M1)
            SUML=SUML+FMAS(I,J)*W1(J+M1)
            SUMQ=SUMQ+FMAS(I,J)*Q1(J+M1)
         END DO
         IM1=I+M1
         F2(IM1)=SUM
         CONT(IM1)=SUM+Y0(IM1)
         W2(IM1)=SUM
         Q2(IM1)=SUM
      END DO
      F3=G0*H*CONT
      GOTO 48
C
   6  CONTINUE
C ------  B IS A FULL MATRIX, JACOBIAN A BANDED MATRIX
C ------  THIS OPTION IS NOT PROVIDED
      DEALLOCATE (W1,W2,Q1,Q2)
      RETURN
C
   7  CONTINUE
C ------  B=IDENTITY, JACOBIAN A FULL MATRIX, HESSENBERG-OPTION
      DO I=1,N 
         F2(I)=HEE1*Z1(I)+HEE2*Z2(I)+HEE3*Z3(I)
         CONT(I)=F2(I)+Y0(I)
         W1(I)=CL1*Z1(I)+CL3*Z3(I)
         Q1(I)=CQ1*Z1(I)+CQ2*Z2(I)+CQ3*Z3(I)
      END DO
      W2=W1
      Q2=Q1
      F3=G0*H*CONT
      DO MM=N-2,1,-1
         MP=N-MM
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 310
         ZSAFE=CONT(MP)
         ZSAFEL=W2(MP)
         ZSAFEQ=Q2(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
         W2(MP)=W2(I)
         W2(I)=ZSAFEL
         Q2(MP)=Q2(I)
         Q2(I)=ZSAFEQ
 310     CONTINUE
         DO I=MP+1,N 
            CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
            W2(I)=W2(I)-FJAC(I,MP-1)*W2(MP)
            Q2(I)=Q2(I)-FJAC(I,MP-1)*Q2(MP)
         END DO
      END DO
      CALL SOLH(N,LDE1,E1,1,CONT,IP1)
      CALL SOLH(N,LDE1,E1,1,W2,IP1)
      CALL SOLH(N,LDE1,E1,1,Q2,IP1)
      DO MM=1,N-2
         MP=N-MM
         DO I=MP+1,N 
            CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
            W2(I)=W2(I)+FJAC(I,MP-1)*W2(MP)
            Q2(I)=Q2(I)+FJAC(I,MP-1)*Q2(MP)
         END DO
         I=IPHES(MP)
         IF (I.EQ.MP) GOTO 440
         ZSAFE=CONT(MP)
         CONT(MP)=CONT(I)
         CONT(I)=ZSAFE
         ZSAFEL=W2(MP)
         W2(MP)=W2(I)
         W2(I)=ZSAFEL
         ZSAFEQ=Q2(MP)
         Q2(MP)=Q2(I)
         Q2(I)=ZSAFE
 440     CONTINUE
      END DO
C
C --------------------------------------
C
  77  CONTINUE
C ********************
C --- ERROR ESTIMATION
C ********************
CCC
CCC   STIMA DELL'ERRORE QUANDO SI USA UN'INTERPOLAZIONE LINEARE
CCC   DA INTERPRETARE COME P_LIN(0)-Y0
CCC   Errore lin. standard (non stiff)
      ERRB=0.D0
      DO  I=1,N
         ERRB=ERRB+(W1(I)/SCAL(I))**2
      END DO
      ERRLB=MAX(SQRT(ERRB/N),1.D-10)
CCC   Errore scalato (stiff)
      IF (ALPHA.NE.0.D0) THEN
       ERR=0.D0
       DO  I=1,N
          ERR=ERR+(W2(I)/SCAL(I))**2
       END DO
       ERRL=MAX(SQRT(ERR/N),1.D-10)
      ELSE
       ERRL=0.D0
      END IF
CCC
CCC   STIMA DELL'ERRORE QUANDO SI USA UN'INTERPOLAZIONE QUADRATA
CCC   DA INTEROPRETARE COME P_QUAD(0)-Y0
CCC   Errore quad. standard (non stiff)
      ERRB=0.D0
      DO  I=1,N
         ERRB=ERRB+(Q1(I)/SCAL(I))**2
      END DO
      ERRQB=MAX(SQRT(ERRB/N),1.D-10)
CCC   Errore scalato (stiff)
      IF (ALPHA.NE.0.D0) THEN
       ERR=0.D0
       DO  I=1,N
          ERR=ERR+(Q2(I)/SCAL(I))**2
       END DO
       ERRQ=MAX(SQRT(ERR/N),1.D-10)
      ELSE
       ERRQ=0.D0
      END IF
      
CCC   CONTINUOUS ERROR EST.
      CERRB=ERRQB*(ERRQB/SQRT(ERRQB*ERRQB+
     &             CERLQ*CERLQ*MIN(ERRLB,ERRQB/CERLQ)**2))
      IF (ALPHA.NE.0.D0) THEN
       CERR=ERRQ*(ERRQ/SQRT(ERRQ*ERRQ+
     &            CERLQ*CERLQ*MIN(ERRL,ERRQ/CERLQ)**2))
      ELSE
       CERR=0.D0
      END IF
CCC   Se ne prende una c.lin
      CERR=ALPHA*CERR+(1.D0-ALPHA)*CERRB

CCC   Errore standard (non stiff)
      ERRB=0.D0
      DO  I=1,N
         ERRB=ERRB+(F3(I)/SCAL(I))**2
      END DO
      ERRB=MAX(SQRT(ERRB/N),1.D-10)
CCC   Errore scalato (stiff)
      ERR=0.D0
      DO  I=1,N
         ERR=ERR+(CONT(I)/SCAL(I))**2
      END DO
      ERR=MAX(SQRT(ERR/N),1.D-10)
CCC   Se ne prende il min
      ERR=MIN(ERRB,ERR)
CCC
CCC   Si evita la seconda pre-moltiplicazione per (I-hgJ)^-1
CCC   se  ------------\/------------
      IF (ERR.LT.1.D0.OR.JEFLAG.GT.0) THEN 
        DEALLOCATE (W1,W2,Q1,Q2)
        RETURN
      ELSE IF (FIRST.OR.REJECT) THEN
CCC
        DO I=1,N
           CONT(I)=Y(I)+CONT(I)
        END DO
C ---
C ---   VOGLIAMO VALUTARE LA FUNZIONE IN X^+
CD      XURO=1.D-15
CD      IF (X.NE.0.D0) THEN
CD       XX = X*(1.D0+XURO)
CD      ELSE
CD       XX = XURO
CD      ENDIF
        XX=X
C ---   VALUTA IN  XX
C        WRITE(6,*) 'DECDEL FCN',N,XX,NRDS, SIZE(PAST), SIZE(IPAST), 
C     &   SIZE(CONT), SIZE(F1), SIZE(RPAR), SIZE(IPAR)
        CALL FCN(N,XX,CONT,F1,ARGLAG,PHI,RPAR,IPAR,
     &           PAST,IPAST,NRDS,LRPAST)
C        WRITE(6,*) 'END DECDEL FCN'
CW      WRITE (6,*) 'Seconda premoltiplicazione per la stima di errore '
CW      WRITE (6,*) 'X, ERR = ',X,ERR
C ---
        NFCN=NFCN+1
        DO I=1,N
           CONT(I)=F1(I)+F2(I)
        END DO
        GOTO (31,32,31,32,31,32,33,55,55,55,41,42,41,42,41), IJOB
C ----- FULL MATRIX OPTION
  31    CONTINUE
        CALL SOL(N,LDE1,E1,CONT,IP1) 
        GOTO 88
C ----- FULL MATRIX OPTION, SECOND ORDER
 41     CONTINUE
        DO J=1,M2
           SUM1=0.D0
           DO K=MM-1,0,-1
              SUM1=(CONT(J+K*M2)+SUM1)/FAC1
              DO I=1,NM1
                 IM1=I+M1
                 CONT(IM1)=CONT(IM1)+FJAC(I,J+K*M2)*SUM1
              END DO
           END DO
        END DO
        CALL SOL(NM1,LDE1,E1,CONT(M1+1),IP1) 
        DO I=M1,1,-1
           CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
        END DO
        GOTO 88
C ----- BANDED MATRIX OPTION
 32     CONTINUE
        CALL SOLB (N,LDE1,E1,MLE,MUE,CONT,IP1)
        GOTO 88
C ----- BANDED MATRIX OPTION, SECOND ORDER
 42     CONTINUE
        DO J=1,M2
           SUM1=0.D0
           DO K=MM-1,0,-1
              SUM1=(CONT(J+K*M2)+SUM1)/FAC1
              DO I=MAX(1,J-MUJAC),MIN(NM1,J+MLJAC)
                 IM1=I+M1
                 CONT(IM1)=CONT(IM1)+FJAC(I+MUJAC+1-J,J+K*M2)*SUM1
              END DO
           END DO
        END DO
        CALL SOLB (NM1,LDE1,E1,MLE,MUE,CONT(M1+1),IP1)
        DO I=M1,1,-1
           CONT(I)=(CONT(I)+CONT(M2+I))/FAC1
        END DO
        GOTO 88
C ----- HESSENBERG MATRIX OPTION
  33    CONTINUE
        DO MM=N-2,1,-1
           MP=N-MM
           I=IPHES(MP)
           IF (I.EQ.MP) GOTO 510
           ZSAFE=CONT(MP)
           CONT(MP)=CONT(I)
           CONT(I)=ZSAFE
 510       CONTINUE
           DO I=MP+1,N 
              CONT(I)=CONT(I)-FJAC(I,MP-1)*CONT(MP)
           END DO
        END DO
        CALL SOLH(N,LDE1,E1,1,CONT,IP1)
        DO MM=1,N-2
           MP=N-MM
           DO I=MP+1,N 
              CONT(I)=CONT(I)+FJAC(I,MP-1)*CONT(MP)
           END DO
           I=IPHES(MP)
           IF (I.EQ.MP) GOTO 640
           ZSAFE=CONT(MP)
           CONT(MP)=CONT(I)
           CONT(I)=ZSAFE
 640       CONTINUE
        END DO
C -----------------------------------
   88   CONTINUE
        SERR=ERR
        ERR=0.D0 
        DO I=1,N
           ERR=ERR+(CONT(I)/SCAL(I))**2
        END DO
        ERR=MAX(SQRT(ERR/N),1.D-10)
        ERR=MIN(SERR,ERR)
      END IF
      DEALLOCATE (W1,W2,Q1,Q2)
      RETURN
C -----------------------------------------------------------
  55  CONTINUE
      DEALLOCATE (W1,W2,Q1,Q2)
      RETURN
      END
C
C     END OF SUBROUTINE ESTRAD
C
C ***********************************************************
C
C
