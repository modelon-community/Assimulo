*********************************************************
C
      DOUBLE PRECISION FUNCTION CONTR5(I,N,X,CONT,XSOL,HSOL) 
C ----------------------------------------------------------
C     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN
C     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X.
C     IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR
C     THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAR5).
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER, PARAMETER :: DP=kind(1D0)
      REAL(kind=DP), dimension(4*N), intent(in) ::  CONT
C --- REQUIRED CONSTANTS
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2

      N2=2*N
      N3=3*N

      S=(X-XSOL)/HSOL

      CONTR5=CONT(I)+S*(CONT(I+N)+(S-C2M1)*(CONT(I+N2)
     &     +(S-C1M1)*CONT(I+N3)))

      RETURN
      END
C
C     END OF FUNCTION CONTR5
C
C ***********************************************************

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
     &                  IOUT,PAST,IPAST,NRDS,JEFLAG,IEFLAG,LPAST)
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
      REAL(kind=DP), dimension(LPAST), intent(in)  :: PAST
      REAL(kind=DP), dimension(:), allocatable :: W1,W2,Q1,Q2
      INTEGER, dimension(1), intent(in) :: IPAST
      REAL(kind=DP), dimension(1), intent(in)  :: RPAR
      INTEGER, dimension(1), intent(in) :: IPAR
      INTEGER :: LPAST

      LOGICAL FIRST,REJECT,LEFT
      COMMON/LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG
C --- COMMON/BPCOM/BPP,ILBP,LEFT
      EXTERNAL ARGLAG,PHI

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
        CALL FCN(N,XX,CONT,F1,ARGLAG,PHI,RPAR,IPAR,
     &           PAST,IPAST,NRDS,LPAST)
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
      SUBROUTINE DEC (N, NDIM, A, IP, IER)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAY  A .
C     A = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
C     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
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
C----------------------- END OF SUBROUTINE DEC -------------------------
      END
C
C
      SUBROUTINE SOL (N, NDIM, A, B, IP)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAY  A .
C    A = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    B = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    B = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
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
C----------------------- END OF SUBROUTINE SOL -------------------------
      END
c
c
      SUBROUTINE DECH (N, NDIM, A, LB, IP, IER)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J,LB,NA
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A HESSENBERG
C  MATRIX WITH LOWER BANDWIDTH LB
C  INPUT..
C     N = ORDER OF MATRIX A.
C     NDIM = DECLARED DIMENSION OF ARRAY  A .
C     A = MATRIX TO BE TRIANGULARIZED.
C     LB = LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED, LB.GE.1).
C  OUTPUT..
C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U .
C     A(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOLH  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N).
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A SLIGHT MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
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
C----------------------- END OF SUBROUTINE DECH ------------------------
      END
C
C
      SUBROUTINE SOLH (N, NDIM, A, LB, B, IP)
C VERSION REAL DOUBLE PRECISION
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1,LB,NA
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX A.
C    NDIM = DECLARED DIMENSION OF ARRAY  A .
C    A = TRIANGULARIZED MATRIX OBTAINED FROM DECH.
C    LB = LOWER BANDWIDTH OF A.
C    B = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DECH HAS SET IER .NE. 0.
C  OUTPUT..
C    B = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
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
C----------------------- END OF SUBROUTINE SOLH ------------------------
      END
C
      SUBROUTINE DECC (N, NDIM, AR, AI, IP, IER)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION
C  ------ MODIFICATION FOR COMPLEX MATRICES --------
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI .
C     (AR, AI) = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART.
C     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART.
C     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    REAL PART.
C     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    IMAGINARY PART.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
      IER = 0
      IP(N) = 1
      IF (N .EQ. 1) GO TO 70
      NM1 = N - 1
      DO 60 K = 1,NM1
        KP1 = K + 1
        M = K
        DO 10 I = KP1,N
          IF (ABS(AR(I,K))+ABS(AI(I,K)) .GT.
     &          ABS(AR(M,K))+ABS(AI(M,K))) M = I
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
C----------------------- END OF SUBROUTINE DECC ------------------------
      END
C
C
      SUBROUTINE SOLC (N, NDIM, AR, AI, BR, BI, IP)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI.
C    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    (BR,BI) = RIGHT HAND SIDE VECTOR.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    (BR,BI) = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
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
C----------------------- END OF SUBROUTINE SOLC ------------------------
      END  
C
C
      SUBROUTINE DECHC (N, NDIM, AR, AI, LB, IP, IER)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NDIM,IP,IER,NM1,K,KP1,M,I,J
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION
C  ------ MODIFICATION FOR COMPLEX MATRICES --------
C  INPUT..
C     N = ORDER OF MATRIX.
C     NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI .
C     (AR, AI) = MATRIX TO BE TRIANGULARIZED.
C  OUTPUT..
C     AR(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; REAL PART.
C     AI(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U ; IMAGINARY PART.
C     AR(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    REAL PART.
C     AI(I,J), I.GT.J = MULTIPLIERS = LOWER TRIANGULAR FACTOR, I - L.
C                                                    IMAGINARY PART.
C     LB = LOWER BANDWIDTH OF A (DIAGONAL NOT COUNTED), LB.GE.1.
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW.
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER = 0 IF MATRIX A IS NONSINGULAR, OR K IF FOUND TO BE
C           SINGULAR AT STAGE K.
C  USE  SOL  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  IF IP(N)=O, A IS SINGULAR, SOL WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
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
          IF (ABS(AR(I,K))+ABS(AI(I,K)) .GT.
     &          ABS(AR(M,K))+ABS(AI(M,K))) M = I
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
C----------------------- END OF SUBROUTINE DECHC -----------------------
      END
C
C
      SUBROUTINE SOLHC (N, NDIM, AR, AI, LB, BR, BI, IP)
C VERSION COMPLEX DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NDIM,IP,NM1,K,KP1,M,I,KB,KM1
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N = ORDER OF MATRIX.
C    NDIM = DECLARED DIMENSION OF ARRAYS  AR AND AI.
C    (AR,AI) = TRIANGULARIZED MATRIX OBTAINED FROM DEC.
C    (BR,BI) = RIGHT HAND SIDE VECTOR.
C    LB = LOWER BANDWIDTH OF A.
C    IP = PIVOT VECTOR OBTAINED FROM DEC.
C  DO NOT USE IF DEC HAS SET IER .NE. 0.
C  OUTPUT..
C    (BR,BI) = SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
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
C----------------------- END OF SUBROUTINE SOLHC -----------------------
      END  
C
      SUBROUTINE DECB (N, NDIM, A, ML, MU, IP, IER)
      DOUBLE PRECISION A,T
      DIMENSION A(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED
C  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
C  INPUT..
C     N       ORDER OF THE ORIGINAL MATRIX A.
C     NDIM    DECLARED DIMENSION OF ARRAY  A.
C     A       CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS  
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  A  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS 
C                ML+1 THROUGH 2*ML+MU+1 OF  A.
C     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C  OUTPUT..
C     A       AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND 
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.  
C     IP      INDEX VECTOR OF PIVOT INDICES.
C     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
C                SINGULAR AT STAGE K.
C  USE  SOLB  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
C  IF IP(N)=O, A IS SINGULAR, SOLB WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
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
C----------------------- END OF SUBROUTINE DECB ------------------------
      END
C
C
      SUBROUTINE SOLB (N, NDIM, A, ML, MU, B, IP)
      DOUBLE PRECISION A,B,T
      DIMENSION A(NDIM,N), B(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B .
C  INPUT..
C    N      ORDER OF MATRIX A.
C    NDIM   DECLARED DIMENSION OF ARRAY  A .
C    A      TRIANGULARIZED MATRIX OBTAINED FROM DECB.
C    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    B      RIGHT HAND SIDE VECTOR.
C    IP     PIVOT VECTOR OBTAINED FROM DECB.
C  DO NOT USE IF DECB HAS SET IER .NE. 0.
C  OUTPUT..
C    B      SOLUTION VECTOR, X .
C-----------------------------------------------------------------------
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
C----------------------- END OF SUBROUTINE SOLB ------------------------
      END
C
      SUBROUTINE DECBC (N, NDIM, AR, AI, ML, MU, IP, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AR(NDIM,N), AI(NDIM,N), IP(N)
C-----------------------------------------------------------------------
C  MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION OF A BANDED COMPLEX
C  MATRIX WITH LOWER BANDWIDTH ML AND UPPER BANDWIDTH MU
C  INPUT..
C     N       ORDER OF THE ORIGINAL MATRIX A.
C     NDIM    DECLARED DIMENSION OF ARRAY  A.
C     AR, AI     CONTAINS THE MATRIX IN BAND STORAGE.   THE COLUMNS  
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  AR (REAL
C                PART) AND AI (IMAGINARY PART)  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS 
C                ML+1 THROUGH 2*ML+MU+1 OF  AR AND AI.
C     ML      LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C     MU      UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C  OUTPUT..
C     AR, AI  AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND 
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.  
C     IP      INDEX VECTOR OF PIVOT INDICES.
C     IP(N)   (-1)**(NUMBER OF INTERCHANGES) OR O .
C     IER     = 0 IF MATRIX A IS NONSINGULAR, OR  = K IF FOUND TO BE
C                SINGULAR AT STAGE K.
C  USE  SOLBC  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C  DETERM(A) = IP(N)*A(MD,1)*A(MD,2)*...*A(MD,N)  WITH MD=ML+MU+1.
C  IF IP(N)=O, A IS SINGULAR, SOLBC WILL DIVIDE BY ZERO.
C
C  REFERENCE..
C     THIS IS A MODIFICATION OF
C     C. B. MOLER, ALGORITHM 423, LINEAR EQUATION SOLVER,
C     C.A.C.M. 15 (1972), P. 274.
C-----------------------------------------------------------------------
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
          IF (ABS(AR(I,K))+ABS(AI(I,K)) .GT.
     &          ABS(AR(M,K))+ABS(AI(M,K))) M = I
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
C----------------------- END OF SUBROUTINE DECBC ------------------------
      END
C
C
      SUBROUTINE SOLBC (N, NDIM, AR, AI, ML, MU, BR, BI, IP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION AR(NDIM,N), AI(NDIM,N), BR(N), BI(N), IP(N)
C-----------------------------------------------------------------------
C  SOLUTION OF LINEAR SYSTEM, A*X = B ,
C                  VERSION BANDED AND COMPLEX-DOUBLE PRECISION.
C  INPUT..
C    N      ORDER OF MATRIX A.
C    NDIM   DECLARED DIMENSION OF ARRAY  A .
C    AR, AI TRIANGULARIZED MATRIX OBTAINED FROM DECB (REAL AND IMAG. PART).
C    ML     LOWER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    MU     UPPER BANDWIDTH OF A (DIAGONAL IS NOT COUNTED).
C    BR, BI RIGHT HAND SIDE VECTOR (REAL AND IMAG. PART).
C    IP     PIVOT VECTOR OBTAINED FROM DECBC.
C  DO NOT USE IF DECB HAS SET IER .NE. 0.
C  OUTPUT..
C    BR, BI SOLUTION VECTOR, X (REAL AND IMAG. PART).
C-----------------------------------------------------------------------
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
C----------------------- END OF SUBROUTINE SOLBC ------------------------
      END
c
C
      subroutine elmhes(nm,n,low,igh,a,int)
C
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      real*8 a(nm,n)
      real*8 x,y
C     real*8 abs
      integer int(igh)
C
C     this subroutine is a translation of the algol procedure elmhes,
C     num. math. 12, 349-368(1968) by martin and wilkinson.
C     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
C
C     given a real general matrix, this subroutine
C     reduces a submatrix situated in rows and columns
C     low through igh to upper hessenberg form by
C     stabilized elementary similarity transformations.
C
C     on input:
C
C      nm must be set to the row dimension of two-dimensional
C        array parameters as declared in the calling program
C        dimension statement;
C
C      n is the order of the matrix;
C
C      low and igh are integers determined by the balancing
C        subroutine  balanc.      if  balanc  has not been used,
C        set low=1, igh=n;
C
C      a contains the input matrix.
C
C     on output:
C
C      a contains the hessenberg matrix.  the multipliers
C        which were used in the reduction are stored in the
C        remaining triangle under the hessenberg matrix;
C
C      int contains information on the rows and columns
C        interchanged in the reduction.
C        only elements low through igh are used.
C
C     questions and comments should be directed to b. s. garbow,
C     applied mathematics division, argonne national laboratory
C
C     ------------------------------------------------------------------
C
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
C
      do 180 m = kp1, la
       mm1 = m - 1
       x = 0.0d0
       i = m
C
       do 100 j = m, igh
          if (abs(a(j,mm1)) .le. abs(x)) go to 100
          x = a(j,mm1)
          i = j
  100   continue
C
       int(m) = i
       if (i .eq. m) go to 130
C    :::::::::: interchange rows and columns of a ::::::::::
       do 110 j = mm1, n
          y = a(i,j)
          a(i,j) = a(m,j)
          a(m,j) = y
  110   continue
C
       do 120 j = 1, igh
          y = a(j,i)
          a(j,i) = a(j,m)
          a(j,m) = y
  120   continue
C    :::::::::: end interchange ::::::::::
  130   if (x .eq. 0.0d0) go to 180
       mp1 = m + 1
C
       do 160 i = mp1, igh
          y = a(i,mm1)
          if (y .eq. 0.0d0) go to 160
          y = y / x
          a(i,mm1) = y
C
          do 140 j = m, n
  140      a(i,j) = a(i,j) - y * a(m,j)
C
          do 150 j = 1, igh
  150      a(j,m) = a(j,m) + y * a(j,i)
C
  160   continue
C
  180 continue
C
  200 return
C    :::::::::: last card of elmhes ::::::::::
      end

C ***********************************************************
C
      DOUBLE PRECISION FUNCTION DONTR5(I,N,X,CONT,XSOL,HSOL) 
C ----------------------------------------------------------
C     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN
C     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION DERIVATIVE AT .
C     X. IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR
C     THE LAST SUCCESSFULLY COMPUTED STEP (BY DELAY5).
C ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER, PARAMETER :: DP=kind(1D0)
      REAL(kind=DP), dimension(1), intent(in) ::  CONT
C --- REQUIRED CONSTANTS
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2

      N2=2*N
      N3=3*N

      S=(X-XSOL)/HSOL

        DONTR5=(CONT(I+N)+
     &   (S-C2M1)*(CONT(I+N2)+CONT(I+N3)*(S-C1M1))+
     &    S*(CONT(I+N2)+CONT(I+N3)*(2.*S-C1M1-C2M1)))/HSOL
     
      RETURN
      END
C
C     END OF FUNCTION DONTR5
C
C ***********************************************************

C ---------------------------------------------------------- 
C     NUMERICAL SOLUTION OF A STIFF DIFFERENTIAL  
C     (OR DIFFERENTIAL ALGEBRAIC) SYSTEM OF FIRST 0RDER  
C     DELAY DIFFERENTIAL EQUATIONS   
C                     M*Y'(X)=F(X,Y(X),Y(X-A),...). 
C     THE SYSTEM CAN BE (LINEARLY) IMPLICIT (MASS-MATRIX M .NE. I) 
C     OR EXPLICIT (M=I). 
C     NOTE: THIS FORM ALSO ALLOWS TO SOLVE NEUTRAL DIFFERENTIAL PROBLEMS 
C 
C     NOTE: THIS VERSION ALLOWS ARBITRARILY LARGE STEPSIZES 
C     (POSSIBLY LARGER THAN THE DELAYS) 
C 
C     THE METHOD USED IS AN IMPLICIT RUNGE-KUTTA METHOD (3 STAGE  
C     RADAU IIA) OF ORDER 5 WITH STEP SIZE CONTROL AND CONTINUOUS  
C     EXTENSION OF ORDER 3  (C.F. SECTION IV.8 OF (HW)) 
C 
C     AUTHORS: N. GUGLIELMI(*) AND E. HAIRER($)  
C          (*) UNIVERSITA` DELL'AQUILA, DIP. DI MATEMATICA 
C              VIA VETOIO (COPPITO), 67010 L'AQUILA, ITALY 
C          ($) UNIVERSITE DE GENEVE, DEPT. DE MATHEMATIQUES 
C              CH-1211 GENEVE 24, SWITZERLAND  
C                        ---------------------- 
C              E-MAIL ADRESSES:   
C                                         guglielm@univaq.it 
C                                     Ernst.Hairer@math.unige.ch 
C      
C     THIS PROGRAM EXTENDS THE CODE RADAU5 (BY E. HAIRER AND G. WANNER) 
C     TO THE CASE OF DELAY DIFFERENTIAL EQUATIONS.  
C     DETAILS ABOUT RADAU5 CAN BE FOUND IN THE BOOK: 
C     (HW)  E. HAIRER AND G. WANNER, SOLVING ORDINARY DIFFERENTIAL 
C           EQUATIONS II. STIFF AND DIFFERENTIAL-ALGEBRAIC PROBLEMS. 
C           SPRINGER SERIES IN COMPUTATIONAL MATHEMATICS 14, 
C           SPRINGER-VERLAG 1991, SECOND EDITION 1996. 
C     DETAILS ABOUT RADAR5 CAN BE FOUND IN THE PAPERS:
C     (GH1)
C     (GH2)       
C
C     VERSION 2.1 OF JULY 21, 2005 
C 
C ---------------------------------------------------------- 
      MODULE IP_ARRAY 
C        THIS VECTOR IPOSV HAS THE DIMENSION OF THE MAXIMUM NUMBER OF 
C        ALLOWED DELAYS; IF A LARGER NUMBER OF RETARDED ARGUMENT IS  
C        REQUIRED CHANGE THE DIMENSION TO THE DESIRED VALUE AND RECOMPILE 
C        HERE THE DIMENSION IS SET TO 100 
         INTEGER, dimension(100) :: IPOSV 
      END MODULE IP_ARRAY 
C 
      SUBROUTINE RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H, 
     &                  RTOL,ATOL,ITOL, 
     &                  JAC,IJAC,MLJAC,MUJAC, 
     &                  JACLAG,NLAGS,NJACL, 
     &                  IMAS,SOLOUT,IOUT, 
     &                  WORK,IWORK,RPAR,IPAR,IDID, 
     &                  GRID,IPAST,MAS,MLMAS,MUMAS,
     &                  LIPAST, LGRID, LPAST, PAST)
C ---------------------------------------------------------- 
C     INPUT PARAMETERS   
C --------------------   
C     N           DIMENSION OF THE SYSTEM  
C 
C     FCN         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE RIGHT- 
C                 HAND-SIDE OF THE DELAY EQUATION, E.G., 
C                    SUBROUTINE FCN(N,X,Y,F,RPAR,IPAR,...) 
C                    DOUBLE PRECISION X,Y(N),F(N) 
C                    EXTERNAL PHI 
C                    F(1)=G1(X,Y(*),YLAGR5(*,X-TAU(X,Y(*))),PHI,...)) 
C                    F(2)=G2(X,Y(*),YLAGR5(*,X-TAU(X,Y(*))),PHI,...)) 
C                    ETC. 
C                    (*) MEANS ALL POSSIBLE COMPONENTS 
C                 FOR AN EXPLICATION OF YLAGR5 SEE BELOW. 
C                 DO NOT USE YLAGR5(I,X-0.D0,PHI,RPAR,IPAR,...) ! 
C                 Note: 
C                 THE INITIAL FUNCTION HAS TO BE SUPPLIED BY: 
C                    FUNCTION PHI(I,X,RPAR,IPAR) 
C                    DOUBLE PRECISION PHI,X 
C                 WHERE I IS THE COMPONENT AND X THE ARGUMENT 
C                 RPAR, IPAR (SEE BELOW) 
C 
C     X           INITIAL X-VALUE 
C 
C     Y(N)        INITIAL VALUES FOR Y (MAY BE DIFFERENT FROM PHI (I,X), 
C                 IN THIS CASE IT IS HIGHLY RECOMMENDED TO SET IWORK(13) 
C                 AND GRID(1),..., (SEE BELOW) 
C 
C     XEND        FINAL X-VALUE (XEND-X HAS TO BE POSITIVE) 
C 
C     H           INITIAL STEP SIZE GUESS; 
C                 FOR STIFF EQUATIONS WITH INITIAL TRANSIENT,  
C                 H=1.D0/(NORM OF F'), USUALLY 1.D-3 OR 1.D-5, IS GOOD. 
C                 THIS CHOICE IS NOT VERY IMPORTANT, THE STEP SIZE IS 
C                 QUICKLY ADAPTED. (IF H=0.D0, THE CODE PUTS H=1.D-6). 
C 
C     RTOL,ATOL   RELATIVE AND ABSOLUTE ERROR TOLERANCES. THEY 
C                 CAN BE BOTH SCALARS OR ELSE BOTH VECTORS OF LENGTH N. 
C 
C     ITOL        SWITCH FOR RTOL AND ATOL: 
C                   ITOL=0: BOTH RTOL AND ATOL ARE SCALARS. 
C                     THE CODE KEEPS, ROUGHLY, THE LOCAL ERROR OF 
C                     Y(I) OVER RTOL*ABS(Y(I))+ATOL 
C                   ITOL=1: BOTH RTOL AND ATOL ARE VECTORS. 
C                     THE CODE KEEPS THE LOCAL ERROR OF Y(I) OVER 
C                     RTOL(I)*ABS(Y(I))+ATOL(I). 
C 
C     JAC         NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES 
C                 THE PARTIAL DERIVATIVES OF F(X,Y) WITH RESPECT TO Y 
C                 (THIS ROUTINE IS ONLY CALLED IF IJAC=1;  
C                 THE USER HAS TO SUPPLY A DUMMY SUBROUTINE  
C                 IN THE CASE IJAC=0). 
C                 FOR IJAC=1, THIS SUBROUTINE MUST HAVE THE FORM 
C                    SUBROUTINE JAC(N,X,Y,DFY,LDFY,RPAR,IPAR,...) 
C                    DOUBLE PRECISION X,Y(N),DFY(LDFY,N) 
C                    DFY(1,1)= ... 
C                 LDFY, THE COLUMN-LENGTH OF THE ARRAY, IS 
C                 FURNISHED BY THE CALLING PROGRAM. 
C                 IF (MLJAC.EQ.N) THE JACOBIAN IS SUPPOSED TO 
C                    BE FULL AND THE PARTIAL DERIVATIVES ARE 
C                    STORED IN DFY AS 
C                       DFY(I,J) = PARTIAL F(I) / PARTIAL Y(J) 
C                 ELSE, THE JACOBIAN IS TAKEN AS BANDED AND 
C                    THE PARTIAL DERIVATIVES ARE STORED 
C                    DIAGONAL-WISE AS 
C                       DFY(I-J+MUJAC+1,J) = PARTIAL F(I) / PARTIAL Y(J). 
C 
C     IJAC        SWITCH FOR THE COMPUTATION OF THE JACOBIAN: 
C                    IJAC=0: JACOBIAN IS COMPUTED INTERNALLY BY FINITE 
C                       DIFFERENCES, SUBROUTINE "JAC" IS NEVER CALLED. 
C                    IJAC=1: JACOBIAN IS SUPPLIED BY SUBROUTINE JAC. 
C 
C     MLJAC       SWITCH FOR THE BANDED STRUCTURE OF THE JACOBIAN: 
C                    MLJAC=N: JACOBIAN IS A FULL MATRIX. THE LINEAR 
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. 
C                    0<=MLJAC<N: MLJAC IS THE LOWER BANDWITH OF JACOBIAN  
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW 
C                       THE MAIN DIAGONAL). 
C 
C     MUJAC       UPPER BANDWITH OF JACOBIAN  MATRIX (>= NUMBER OF NON- 
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). 
C                 DOES NOT NEED TO BE DEFINED IF MLJAC=N. 
C 
C     JACLAG      NAME (EXTERNAL) OF THE SUBROUTINE WHICH COMPUTES 
C                 THE PARTIAL DERIVATIVES OF F(X,Y,YLAG) WITH RESPECT TO  
C                 YLAG(*) (YLAG DENOTE THE DELAYED VARIABLES) 
C 
C     NLAGS       DENOTES THE NUMBER OF DELAY ARGUMENTS.  
C                 THIS PARAMETER IS OF INTEREST FOR THE COMPUTATION OF THE 
C                 JACOBIAN. 
C                 TO BE SET = 0 IF ONE DOES WANT TO COMPUTE THE TRADITIONAL 
C                 JACOBIAN;  
C                 TO BE SET = NUMBER OF DISTINCT DELAY ARGUMENTS 
C                 IF ONE WANTS TO CORRECT THE STANDARD JACOBIAN (THROUGH 
C                 THE SUBROUTINE JACLAG) WHEN ADVANCED ARGUMENTS ARE USED. 
C 
C     NJACL       NUMBER OF TERMS IN THE JACOBIAN W.R.T. 
C                 RETARDED COMPONENTS (WHICH IS THOUGHT AS A SPARSE MATRIX). 
C 
C     ----   MAS,IMAS,MLMAS, AND MUMAS HAVE ANALOG MEANINGS      ----- 
C     ----   FOR THE "MASS MATRIX" (THE MATRIX "M" OF SECTION IV.8): - 
C 
C     MAS         NAME (EXTERNAL) OF SUBROUTINE COMPUTING THE MASS- 
C                 MATRIX M. 
C                 IF IMAS=0, THIS MATRIX IS ASSUMED TO BE THE IDENTITY 
C                 MATRIX AND NEEDS NOT TO BE DEFINED; 
C                 THE USER HAS TO SUPPLY A DUMMY SUBROUTINE IN THIS CASE. 
C                 IF IMAS=1, THE SUBROUTINE MAS IS OF THE FORM 
C                    SUBROUTINE MAS(N,AM,LMAS,RPAR,IPAR) 
C                    DOUBLE PRECISION AM(LMAS,N) 
C                    AM(1,1)= .... 
C                    IF (MLMAS.EQ.N) THE MASS-MATRIX IS STORED 
C                    AS FULL MATRIX LIKE 
C                         AM(I,J) = M(I,J) 
C                    ELSE, THE MATRIX IS TAKEN AS BANDED AND STORED 
C                    DIAGONAL-WISE AS 
C                         AM(I-J+MUMAS+1,J) = M(I,J). 
C 
C     IMAS       GIVES INFORMATION ON THE MASS-MATRIX: 
C                    IMAS=0: M IS SUPPOSED TO BE THE IDENTITY 
C                       MATRIX, MAS IS NEVER CALLED. 
C                    IMAS=1: MASS-MATRIX  IS SUPPLIED. 
C 
C     MLMAS       SWITCH FOR THE BANDED STRUCTURE OF THE MASS-MATRIX: 
C                    MLMAS=N: THE FULL MATRIX CASE. THE LINEAR 
C                       ALGEBRA IS DONE BY FULL-MATRIX GAUSS-ELIMINATION. 
C                    0<=MLMAS<N: MLMAS IS THE LOWER BANDWITH OF THE 
C                       MATRIX (>= NUMBER OF NON-ZERO DIAGONALS BELOW 
C                       THE MAIN DIAGONAL). 
C                 MLMAS IS SUPPOSED TO BE <= MLJAC. 
C 
C     MUMAS       UPPER BANDWITH OF MASS-MATRIX (>= NUMBER OF NON- 
C                 ZERO DIAGONALS ABOVE THE MAIN DIAGONAL). 
C                 DOES NOT NEED TO BE DEFINED IF MLMAS=N. 
C                 MUMAS IS SUPPOSED TO BE <= MUJAC. 
C 
C     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE 
C                 NUMERICAL SOLUTION DURING INTEGRATION.  
C                 IF IOUT=1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP. 
C                 THE USER HAS TO SUPPLY A DUMMY SUBROUTINE IF IOUT=0.  
C                 IT MUST HAVE THE FORM 
C                    SUBROUTINE SOLOUT (NR,XOLD,X,HSOL,Y,CONT,LRC,N, 
C                                       RPAR,IPAR,IRTRN) 
C                    DOUBLE PRECISION X,Y(N),CONT(LRC) 
C                    ....   
C                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH 
C                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS 
C                    THE FIRST GRID-POINT). 
C                 "XOLD" IS THE PRECEEDING GRID-POINT. 
C                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN 
C                    IS SET <0, RADAR5 RETURNS TO THE CALLING PROGRAM. 
C            
C          -----  CONTINUOUS OUTPUT: ----- 
C                 DURING CALLS TO "SOLOUT" AS WELL AS TO "FCN", A 
C                 CONTINUOUS SOLUTION IS AVAILABLE THROUGH HTHE FUNCTION 
C                        >>>   YLAGR5(I,S,PHI,RPAR,IPAR,...)   <<< 
C                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH 
C                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE S 
C                 HAS TO LIE IN AN INTERVAL WHERE THE NUMERICAL SOLUTION 
C                 IS ALREADY COMPUTED. IT DEPENDS ON THE SIZE OF LRPAST 
C                 (SEE BELOW) HOW FAR BACK THE SOLUTION IS AVAILABLE. 
C 
C     IOUT        SWITCH FOR CALLING THE SUBROUTINE SOLOUT: 
C                    IOUT=0: SUBROUTINE IS NEVER CALLED 
C                    IOUT=1: SUBROUTINE IS AVAILABLE FOR OUTPUT. 
C 
C     WORK        ARRAY OF STATE VARIABLES OF REAL TYPE FOR EXECUTION. 
C                 WORK(1), WORK(2),.., WORK(20) SERVE AS PARAMETERS 
C                 FOR THE CODE. FOR STANDARD USE OF THE CODE 
C                 WORK(1),..,WORK(20) MUST BE SET TO ZERO BEFORE 
C 
C     IWORK       INTEGER WORKING SPACE OF LENGTH "LIWORK". 
C                 IWORK(1),IWORK(2),...,IWORK(20) SERVE AS PARAMETERS 
C                 FOR THE CODE. FOR STANDARD USE, SET IWORK(1),.., 
C                 IWORK(20) TO ZERO BEFORE CALLING. 
C 
C     GRID        CONTAINS PRESCRIBED GRID POINTS, WHICH THE 
C                 INTEGRATION METHOD HAS TO TAKE AS GRID-POINTS 
C                 NORMALLY, IF GRID(1) > X, THEN ONE HAS 
C                 X < GRID(1) < GRID(2) < ... < GRID(NGRID) <= XEND 
C                 IN SOME CASES IF THERE ARE DISCONTINUITIES IN THE  
C                 INITIAL FUNCTIONS THEY ARE ALSO SET IN THE GRID  
C                 VECTOR; THEN X < GRID(J) < GRID(J+1) ... < XEND 
C                 WHERE J ADDRESSES THE FIRST ASCISSA IN GRID > X 
C 
C     LGRID       DECLARED LENGTH OF GRID VECTOR, 
C                 GRID(LGRID), 
C                 WHICH MUST BE DECLARED IN THE CALLING PROGRAM. 
C                 "LGRID" MUST BE AT LEAST 
C                              NGRID + 1 
C 
C     RPAR, IPAR  REAL AND INTEGER PARAMETERS (OR PARAMETER ARRAYS) WHICH   
C                 CAN BE USED FOR COMMUNICATION BETWEEN YOUR CALLING 
C                 PROGRAM AND THE FCN, JAC, MAS, SOLOUT SUBROUTINES.  
C 
C ---------------------------------------------------------------------- 
C  
C     SOPHISTICATED SETTING OF PARAMETERS 
C     ----------------------------------- 
C              SEVERAL PARAMETERS OF THE CODE ARE TUNED TO MAKE IT WORK  
C              WELL. THEY MAY BE DEFINED BY SETTING WORK(1),... 
C              AS WELL AS IWORK(1),... DIFFERENT FROM ZERO. 
C              FOR ZERO INPUT, THE CODE CHOOSES DEFAULT VALUES: 
C 
C    IWORK(1)  IF IWORK(1).NE.0, THE CODE TRANSFORMS THE JACOBIAN 
C              MATRIX TO HESSENBERG FORM. THIS IS PARTICULARLY 
C              ADVANTAGEOUS FOR LARGE SYSTEMS WITH FULL JACOBIAN. 
C              IT DOES NOT WORK FOR BANDED JACOBIAN (MLJAC<N) 
C              AND NOT FOR IMPLICIT SYSTEMS (IMAS=1). 
C 
C    IWORK(2)  THIS IS THE MAXIMAL NUMBER OF ALLOWED STEPS. 
C              THE DEFAULT VALUE (FOR IWORK(2)=0) IS 100000. 
C 
C    IWORK(3)  THE MAXIMUM NUMBER OF NEWTON ITERATIONS FOR THE 
C              SOLUTION OF THE IMPLICIT SYSTEM IN EACH STEP. 
C              THE DEFAULT VALUE (FOR IWORK(3)=0) IS 7. 
C 
C    IWORK(4)  IF IWORK(4).EQ.0 THE EXTRAPOLATED COLLOCATION SOLUTION 
C              IS TAKEN AS STARTING VALUE FOR NEWTON'S METHOD. 
C              IF IWORK(4).NE.0 ZERO STARTING VALUES ARE USED. 
C              THE LATTER IS RECOMMENDED IF NEWTON'S METHOD HAS 
C              DIFFICULTIES WITH CONVERGENCE (THIS IS SEEN IN THE CASE WHEN 
C              NSTEP IS LARGER THAN NACCPT + NREJCT; SEE OUTPUT PARAM.). 
C              DEFAULT IS IWORK(4)=0. 
C 
C       THE FOLLOWING 3 PARAMETERS ARE IMPORTANT FOR 
C       DELAY DIFFERENTIAL-ALGEBRAIC SYSTEMS OF INDEX > 1. 
C       THE FUNCTION-SUBROUTINE SHOULD BE WRITTEN SUCH THAT 
C       THE INDEX 1,2,3 VARIABLES APPEAR IN THIS ORDER.  
C       IN ESTIMATING THE ERROR THE INDEX 2 VARIABLES ARE 
C       MULTIPLIED BY H, THE INDEX 3 VARIABLES BY H**2. 
C 
C    IWORK(5)  DIMENSION OF THE INDEX 1 VARIABLES (MUST BE > 0). FOR  
C              DDE'S THIS EQUALS THE DIMENSION OF THE SYSTEM. 
C              DEFAULT IWORK(5)=N. 
C 
C    IWORK(6)  DIMENSION OF THE INDEX 2 VARIABLES. DEFAULT IWORK(6)=0. 
C 
C    IWORK(7)  DIMENSION OF THE INDEX 3 VARIABLES. DEFAULT IWORK(7)=0. 
C 
C    IWORK(8)  SWITCH FOR STEP SIZE STRATEGY 
C              IF IWORK(8).EQ.1  MODIFIED PREDICTIVE CONTROLLER  
C              (GUSTAFSSON) 
C              IF IWORK(8).EQ.2  CLASSICAL STEP SIZE CONTROL 
C              THE DEFAULT VALUE (FOR IWORK(8)=0) IS IWORK(8)=1. 
C              THE CHOICE IWORK(8).EQ.1 SEEMS TO PRODUCE SAFER RESULTS; 
C              FOR SIMPLE PROBLEMS, THE CHOICE IWORK(8).EQ.2 PRODUCES 
C              OFTEN SLIGHTLY FASTER RUNS 
C 
C       IF THE DIFFERENTIAL SYSTEM HAS THE SPECIAL STRUCTURE THAT 
C            Y(I)' = Y(I+M2)   FOR  I=1,...,M1, 
C       WITH M1 A MULTIPLE OF M2, A SUBSTANTIAL GAIN IN COMPUTERTIME 
C       CAN BE ACHIEVED BY SETTING THE PARAMETERS IWORK(9) AND IWORK(10). 
C       E.G., FOR SECOND ORDER SYSTEMS P'=V, V'=G(P,V), WHERE P AND V ARE  
C       VECTORS OF DIMENSION N/2, ONE HAS TO PUT M1=M2=N/2. 
C       FOR M1>0 SOME OF THE INPUT PARAMETERS HAVE DIFFERENT MEANINGS: 
C       - JAC: ONLY THE ELEMENTS OF THE NON-TRIVIAL PART OF THE 
C              JACOBIAN HAVE TO BE STORED 
C              IF (MLJAC.EQ.N-M1) THE JACOBIAN IS SUPPOSED TO BE FULL 
C                 DFY(I,J) = PARTIAL F(I+M1) / PARTIAL Y(J) 
C                FOR I=1,N-M1 AND J=1,N. 
C              ELSE, THE JACOBIAN IS BANDED ( M1 = M2 * MM ) 
C                 DFY(I-J+MUJAC+1,J+K*M2) = PARTIAL F(I+M1) / PARTIAL Y(J+K*M2) 
C                FOR I=1,MLJAC+MUJAC+1 AND J=1,M2 AND K=0,MM. 
C       - MLJAC: MLJAC=N-M1: IF THE NON-TRIVIAL PART OF THE JACOBIAN IS FULL 
C                0<=MLJAC<N-M1: IF THE (MM+1) SUBMATRICES (FOR K=0,MM) 
C                     PARTIAL F(I+M1) / PARTIAL Y(J+K*M2),  I,J=1,M2 
C                    ARE BANDED, MLJAC IS THE MAXIMAL LOWER BANDWIDTH 
C                    OF THESE MM+1 SUBMATRICES 
C       - MUJAC: MAXIMAL UPPER BANDWIDTH OF THESE MM+1 SUBMATRICES; 
C                DOES NOT NEED TO BE DEFINED IF MLJAC=N-M1 
C       - MAS: IF IMAS=0 THIS MATRIX IS ASSUMED TO BE THE IDENTITY AND 
C              DOES NOT NEED TO BE DEFINED.  
C              THE USER HAS TO SUPPLY A DUMMY SUBROUTINE IN THIS CASE. 
C              IT IS ASSUMED THAT ONLY THE ELEMENTS OF RIGHT LOWER BLOCK OF 
C              DIMENSION N-M1 DIFFER FROM THAT OF THE IDENTITY MATRIX. 
C              IF (MLMAS.EQ.N-M1) THIS SUBMATRIX IS SUPPOSED TO BE FULL 
C                 AM(I,J) = M(I+M1,J+M1)     FOR I=1,N-M1 AND J=1,N-M1. 
C              ELSE, THE MASS MATRIX IS BANDED 
C                 AM(I-J+MUMAS+1,J) = M(I+M1,J+M1) 
C       - MLMAS: MLMAS=N-M1: IF THE NON-TRIVIAL PART OF M IS FULL 
C                0<=MLMAS<N-M1: LOWER BANDWIDTH OF THE MASS MATRIX 
C       - MUMAS: UPPER BANDWIDTH OF THE MASS MATRIX 
C                DOES NOT NEED TO BE DEFINED IF MLMAS=N-M1 
C 
C    IWORK(9)  THE VALUE OF M1.  DEFAULT M1=0. 
C 
C    IWORK(10) THE VALUE OF M2.  DEFAULT M2=M1. 
C 
C 
C 
C    IWORK(11) SELECT THE TYPE OF ERROR CONTROL:  
C              -1: FOR A PURE CONTROL OF THE DENSE OUTPUT  
C                  (MAKES USE OF A QUADRATIC AND A LINEAR INTERPOLATING 
C                   POLYNOMIALS);  
C               1: FOR A MIXED CONTROL OF DENSE OUTPUT AND DISCRETE OUTPUT 
C               2: FOR A PURE CONTROL OF THE DISCRETE OUTPUT 
C                  (ERROR CONTROL PROVIDED BY THE SUBROUTINE ESTRAD). 
C               3: FOR A SIMPLER MIXED CONTROL OF DENSE OUTPUT AND   
C                  DISCRETE OUTPUT 
C               DEFAULT VALUE IWORK(11)=2. 
C 
C    IWORK(12) = MXST = (ON INPUT)  
C              DECLARED NUMBER OF STEPS STORED IN THE ``PAST VECTOR'',  
C              PAST(LRPAST), 
C              WHICH MUST BE DECLARED IN THE CALLING PROGRAM. 
C              "MXST" MUST BE SUFFICIENTLY LARGE. IF THE DENSE 
C              OUTPUT OF MXST BACK STEPS HAS TO BE STORED,  
C              THE DIMENSION OF PAST MUST BE  
C                       LRPAST=MXST*(4*NRDENS+2) 
C              WHERE NRDENS=IWORK(15) (SEE BELOW). 
C 
C    IWORK(13) = NGRID = (ON INPUT) 
C              NUMBER OF PRESCRIBED POINTS IN THE 
C              INTEGRATION INTERVAL WHICH HAVE TO BE GRID-POINTS 
C              IN THE INTEGRATION. USUALLY, AT THESE POINTS THE 
C              SOLUTION OR ONE OF ITS DERIVATIVE HAS A DISCONTINUITY. 
C              DEFINE THESE POINTS IN GRID(1),...,GRID(NGRID) 
C              DEFAULT VALUE:  IWORK(13)=0 
C 
C    IWORK(14) = SELECTOR FOR FULL ITERATION (2) OR SIMPLIFIED  
C              ITERATION (1) (TAKING INTO ACCOUNT POSSIBLE  
C              ADVANCED ARGUMENTS BUT PRESERVING TENSOR STRUCTURE  
C              OF THE JACOBIAN. 
C              DEFAULT VALUE:  IWORK(14)=1 
C 
C    IWORK(15) = NRDENS = (ON INPUT)  
C              NUMBER OF COMPONENTS, FOR WHICH DENSE OUTPUT 
C              IS REQUIRED (EITHER BY "SOLOUT" OR BY "FCN"); 
C              DEFAULT VALUE (FOR IWORK(15)=0) IS IWORK(15)=N; 
C              FOR   0 < NRDENS < N   THE COMPONENTS (FOR WHICH DENSE 
C              OUTPUT IS REQUIRED) HAVE TO BE SPECIFIED IN 
C              IPAST(1),...,IPAST(NRDENS); 
C              FOR  NRDENS=N  THIS IS DONE BY THE CODE. 
C    IWORK(16) = NDIMN = (ON INPUT)
C              OPTION VALID FOR NEUTRAL PROBLEMS
C              NUMBER OF DERIVATIVE COMPONENTS (Z) OF THE NEUTRAL PROBLEM 
C              EXCLUDED TRUE SOLUTION COMPONENTS
C 
C ---------- 
C 
C    WORK(1)   UROUND, THE ROUNDING UNIT, DEFAULT 1.D-16. 
C 
C    WORK(2)   THE SAFETY FACTOR IN STEP SIZE PREDICTION, 
C              DEFAULT 0.9D0. 
C 
C    WORK(3)   DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; 
C              INCREASE WORK(3), TO 0.1 SAY, WHEN JACOBIAN EVALUATIONS 
C              ARE COSTLY. FOR SMALL SYSTEMS WORK(3) SHOULD BE SMALLER  
C              (0.001D0, SAY). NEGATIV WORK(3) FORCES THE CODE TO 
C              COMPUTE THE JACOBIAN AFTER EVERY ACCEPTED STEP.      
C              DEFAULT 0.001D0. 
C 
C    WORK(4)   STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. 
C              SMALLER VALUES OF WORK(4) MAKE THE CODE SLOWER, BUT SAFER. 
C              DEFAULT MIN(0.03D0,RTOL(1)**0.5D0) 
C 
C    WORK(5) AND WORK(6) : IF WORK(5) < HNEW/HOLD < WORK(6), THEN THE 
C              STEP SIZE IS NOT CHANGED. THIS SAVES, TOGETHER WITH A 
C              LARGE WORK(3), LU-DECOMPOSITIONS AND COMPUTING TIME FOR 
C              LARGE SYSTEMS. FOR SMALL SYSTEMS ONE MAY HAVE 
C              WORK(5)=1.D0, WORK(6)=1.2D0, FOR LARGE FULL SYSTEMS 
C              WORK(5)=0.99D0, WORK(6)=2.D0 MIGHT BE GOOD. 
C              DEFAULTS WORK(5)=1.D0, WORK(6)=1.2D0 . 
C 
C    WORK(7)   MAXIMAL STEP SIZE, DEFAULT XEND-X. 
C 
C    WORK(8), WORK(9)   PARAMETERS FOR STEP SIZE SELECTION 
C              THE NEW STEP SIZE IS CHOSEN SUBJECT TO THE RESTRICTION 
C                 WORK(8) <= HNEW/HOLD <= WORK(9) 
C              DEFAULT VALUES: WORK(8)=0.2D0, WORK(9)=8.D0 
C 
C    WORK(10)  PARAMETER FOR CONTROLLING THE ERROR CONTROL OF DENSE 
C              OUTPUT (0 <= WORK(10) <= 1). (0: STRONG CONTROL, 1: WEAKER) 
C              SUGGESTED VALUES: 
C              FOR PROBLEMS WITH `ALMOST DISCONTINUOUS' SOLUTIONS  
C              (LIKE SHOCKS):  WORK(10)=0.D0 
C              FOR PROBLEMS WITH FAIRLY SMOOTH SOLUTION:  WORK(10)=1.D0 
C              FOR INTERMEDIATE PROBLEMS:  WORK(10)=1.D-M (M=1,2,3,..,) 
C              DEFAULT VALUE: WORK(10)=0.D0  
C    WORK(11)  PARAMETER FOR CONTROLLING THE SEARCH OF BREAKING POINTS:
C              IF THE ERROR INCREASES OF A FACTOR LARGER THAN WORK(11)
C              FROM A STEP TO THE SUBSEQUENT, THE ROUTINE SEARCHING BREAKING
C              POINTS ACTIVATES.
C              DEFAULT VALUE: WORK(11)=5.D0
C----------------------------------------------------------------------- 
C 
C     OUTPUT PARAMETERS  
C     -----------------  
C     X           X-VALUE FOR WHICH THE SOLUTION HAS BEEN COMPUTED 
C                 (AFTER SUCCESSFUL RETURN X=XEND). 
C 
C     Y(N)        NUMERICAL SOLUTION AT X 
C  
C     H           PREDICTED STEP SIZE OF THE LAST ACCEPTED STEP 
C 
C     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN: 
C                   IDID= 1  COMPUTATION SUCCESSFUL, 
C                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT) 
C                   IDID=-1  INPUT IS NOT CONSISTENT, 
C                   IDID=-2  LARGER NMAX IS NEEDED, 
C                   IDID=-3  STEP SIZE BECOMES TOO SMALL, 
C                   IDID=-4  MATRIX IS REPEATEDLY SINGULAR. 
C                   IDID=-5  COMPUTATION INTERRUPTED BY YLAGR5.    
C                   IDID=-6  THE EQUATION USES ADVANCED ARGUMENTS 
C 
C   IWORK(13)  NFULL   NUMBER OF FULL NEWTON ITERATIONS                    
C   IWORK(14)  NFCN    NUMBER OF FUNCTION EVALUATIONS (THOSE FOR NUMERICAL 
C                      EVALUATION OF THE JACOBIAN ARE NOT COUNTED)   
C   IWORK(15)  NJAC    NUMBER OF JACOBIAN EVALUATIONS (EITHER ANALYTICALLY 
C                      OR NUMERICALLY) 
C   IWORK(16)  NSTEP   NUMBER OF COMPUTED STEPS 
C   IWORK(17)  NACCPT  NUMBER OF ACCEPTED STEPS 
C   IWORK(18)  NREJCT  NUMBER OF REJECTED STEPS (DUE TO ERROR TEST), 
C                      (STEP REJECTIONS IN THE FIRST STEP ARE NOT COUNTED) 
C   IWORK(19)  NDEC    NUMBER OF LU-DECOMPOSITIONS OF BOTH MATRICES 
C   IWORK(20)  NSOL    NUMBER OF FORWARD-BACKWARD SUBSTITUTIONS, OF BOTH 
C                      SYSTEMS; THE NSTEP FORWARD-BACKWARD SUBSTITUTIONS, 
C                      NEEDED FOR STEP SIZE SELECTION, ARE NOT COUNTED 
C----------------------------------------------------------------------- 
C *** *** *** *** *** *** *** *** *** *** *** *** *** 
C          DECLARATIONS  
C *** *** *** *** *** *** *** *** *** *** *** *** *** 
      USE IP_ARRAY 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
      REAL(kind=DP), dimension(N), intent(inout) ::  
     &                             Y 
      REAL(kind=DP), dimension(1), intent(inout) ::  
     &                             WORK 
      REAL(kind=DP), dimension(1), intent(inout) ::  
     &                             ATOL,RTOL 
      INTEGER, dimension(1), intent(inout) :: IWORK 
      REAL(kind=DP), dimension(1), intent(inout) :: GRID 
      INTEGER, dimension(1), intent(inout) :: IPAST 
      REAL(kind=DP), dimension(1), intent(in) :: RPAR 
      INTEGER, dimension(1), intent(in) :: IPAR
      INTEGER :: LIPAST
      INTEGER :: LGRID
      INTEGER :: LPAST
      REAL(kind=DP), dimension(LPAST),intent(inout) :: PAST
      
      LOGICAL IMPLCT,NEUTRAL,JBAND,ARRET,STARTN,PRED 
      LOGICAL FLAGS,FLAGN 
 
      EXTERNAL FCN,PHI,ARGLAG,JAC,JACLAG,MAS,SOLOUT 
C ----> COMMON BLOCKS <---- 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
 
C *** *** *** *** *** *** *** 
C        SETTING THE PARAMETERS  
C *** *** *** *** *** *** *** 
C      IF (NLAGS.GT.0) THEN 
C       ALLOCATE (IPOSV(NLAGS)) 
C      ELSE 
C       ALLOCATE (IPOSV(5)) 
C      END IF 
       NN=N 
 
       NFCN=0 
       NJAC=0 
       NSTEP=0 
       NACCPT=0 
       NREJCT=0 
       NDEC=0 
       NSOL=0 
       ARRET=.FALSE. 
       FLAGS=.FALSE. 
       FLAGN=.FALSE. 
 
       IF (IOUT.EQ.1) WRITE (6,*) 'STARTING INTEGRATION...' 
C        
C ------> OPERATIONS RELEVANT TO THE DELAY DEPENDENCE <------ 
C 
C -------- ERROR CONTROL 
      IF (IWORK(11).EQ.0) THEN 
       IEFLAG=2 
      ELSE 
       IEFLAG=IWORK(11) 
      END IF 
      IF (IEFLAG.EQ.2) WORK(10)=1.D0 
C -------- NGRID   NUMBER OF PRESCRIBED GRID-POINTS 
      NGRID=IWORK(13) 
      IF (NGRID.LT.0) NGRID=0 
      IF (IOUT.EQ.1) WRITE(6,*)  
     &           'NUMBER OF PRESCRIBED GRID POINTS: ',NGRID 
C ------- NRDENS   NUMBER OF DENSE OUTPUT COMPONENTS 
      NRDENS=IWORK(15)
C ------- NDIMN   NUMBER OF COMPONENTS OF A NEUTRAL PROBLEM
      IF (IMAS.EQ.2) THEN
         IF (IWORK(16).EQ.0) THEN
          WRITE(6,*) 'NUMBER OF Y COMPONENTS HAS TO BE SPECIFIED'
          ARRET=.TRUE.
         END IF
         NDIMN=IWORK(16)
        ELSE 
         NDIMN=N
        END IF           
C ------- LIPAST   DIMENSION OF VECTOR IPAST 
      LIPAST=NRDENS+1  
      IF(NRDENS.LT.0.OR.NRDENS.GT.N) THEN 
         IF (IOUT.GT.0) WRITE(6,*) 
     &           ' CURIOUS INPUT IWORK(15)=',IWORK(15) 
         ARRET=.TRUE. 
      ELSE IF (NRDENS.EQ.0) THEN 
            NRDS=N 
      ELSE 
            NRDS=NRDENS 
      END IF 
      IF (NRDS.EQ.N) THEN 
            DO 16 I=1,NRDS 
  16           IPAST(I)=I 
      END IF 
      IF (IOUT.EQ.1) WRITE(6,*) 'NUMBER OF DELAYED COMPONENTS: ',NRDS 
C ------- LRPAST   DIMENSION OF VECTOR PAST 
      MXST=IWORK(12) 
C ------- CONTROL OF LENGTH OF PAST  ------- 
      IF(MXST.LT.1)THEN 
         IF (IOUT.GT.0) WRITE(6,*) 
     & ' INSUFFICIENT STORAGE FOR PAST, MIN. LRPAST=',1 
         ARRET=.TRUE. 
      END IF 
C ------- DIM. of PAST  -------- 
      IDIF=4*NRDS+2 
      LRPAST=MXST*IDIF 
C -------------------------------------------------  
C ------- CONTROL OF SIMPLE NEWTON ITERATION  ------- 
      ISWJL=IWORK(14) 
      IF (ISWJL.EQ.0) ISWJL=1 
 
C -------- UROUND : SMALLEST NUMBER SATISFYING 1.0D0+UROUND>1.0D0   
      IF (WORK(1).EQ.0.0D0) THEN 
         UROUND=1.0D-16 
      ELSE 
         UROUND=WORK(1) 
         IF (UROUND.LE.1.0D-19.OR.UROUND.GE.1.0D0) THEN 
            WRITE(6,*)' COEFFICIENTS HAVE 20 DIGITS, UROUND=',WORK(1) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
C -------> CHECK AND CHANGE THE TOLERANCES <------ 
      EXPM=2.0D0/3.0D0 
      IF (ITOL.EQ.0) THEN 
          IF (ATOL(1).LE.0.D0.OR.RTOL(1).LE.10.D0*UROUND) THEN 
              WRITE (6,*) ' TOLERANCES ARE TOO SMALL' 
              ARRET=.TRUE. 
          ELSE 
              QUOT=ATOL(1)/RTOL(1) 
              RTOL(1)=0.1D0*RTOL(1)**EXPM 
              ATOL(1)=RTOL(1)*QUOT 
          END IF 
      ELSE 
          DO I=1,N 
          IF (ATOL(I).LE.0.D0.OR.RTOL(I).LE.10.D0*UROUND) THEN 
              WRITE (6,*) ' TOLERANCES(',I,') ARE TOO SMALL' 
              ARRET=.TRUE. 
          ELSE 
              QUOT=ATOL(I)/RTOL(I) 
              RTOL(I)=0.1D0*RTOL(I)**EXPM 
              ATOL(I)=RTOL(I)*QUOT 
          END IF 
          END DO 
      END IF 
 
C -------> NMAX : THE MAXIMAL NUMBER OF STEPS <------- 
      IF (IWORK(2).EQ.0) THEN 
         NMAX=100000 
      ELSE 
         NMAX=IWORK(2) 
         IF (NMAX.LE.0) THEN 
            WRITE(6,*)' WRONG INPUT IWORK(2)=',IWORK(2) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
C -------> NIT :  MAXIMAL NUMBER OF NEWTON ITERATIONS <------- 
      IF (IWORK(3).EQ.0) THEN 
         NIT=7 
      ELSE 
         NIT=IWORK(3) 
         IF (NIT.LE.0) THEN 
            WRITE(6,*)' CURIOUS INPUT IWORK(3)=',IWORK(3) 
            ARRET=.TRUE. 
         END IF 
      END IF 
C -------- STARTN : SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS 
      IF(IWORK(4).EQ.0)THEN 
         STARTN=.FALSE. 
      ELSE 
         STARTN=.TRUE. 
      END IF 
 
C -------> PARAMETERS (IF ANY) FOR DIFFERENTIAL-ALGEBRAIC COMPONENTS <------- 
      NIND1=IWORK(5) 
      NIND2=IWORK(6) 
      NIND3=IWORK(7) 
      IF (NIND1.EQ.0) NIND1=N 
      IF (NIND1+NIND2+NIND3.NE.N) THEN 
       WRITE(6,*)' CURIOUS INPUT FOR IWORK(5,6,7)=',NIND1,NIND2,NIND3 
       ARRET=.TRUE. 
      END IF 
 
C -------> PRED   STEP SIZE CONTROL <------- 
      IF(IWORK(8).LE.1)THEN 
         PRED=.TRUE. 
      ELSE 
         PRED=.FALSE. 
      END IF 
 
C -------> PARAMETER FOR SECOND ORDER EQUATIONS <------- 
      M1=IWORK(9) 
      M2=IWORK(10) 
      NM1=N-M1 
      IF (M1.EQ.0) M2=N 
      IF (M2.EQ.0) M2=M1 
      IF (M1.LT.0.OR.M2.LT.0.OR.M1+M2.GT.N) THEN 
       WRITE(6,*)' CURIOUS INPUT FOR IWORK(9,10)=',M1,M2 
       ARRET=.TRUE. 
      END IF 
 
C -------> SAFE  :  SAFETY FACTOR IN STEP SIZE PREDICTION <------- 
      IF (WORK(2).EQ.0.0D0) THEN 
         SAFE=0.9D0 
      ELSE 
         SAFE=WORK(2) 
         IF (SAFE.LE.0.001D0.OR.SAFE.GE.1.0D0) THEN 
            WRITE(6,*)' CURIOUS INPUT FOR WORK(2)=',WORK(2) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
C ------> THET : DETERMINES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED; 
      IF (WORK(3).EQ.0.D0) THEN 
         THET=0.001D0 
      ELSE 
         THET=WORK(3) 
         IF (THET.GE.1.0D0) THEN 
            WRITE(6,*)' CURIOUS INPUT FOR WORK(3)=',WORK(3) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
C ---> FNEWT : STOPPING CRITERION FOR NEWTON'S METHOD, USUALLY CHOSEN <1. <--- 
      TOLST=RTOL(1) 
      IF (WORK(4).EQ.0.D0) THEN 
         FNEWT=MAX(10*UROUND/TOLST,MIN(0.03D0,TOLST**0.5D0)) 
      ELSE 
         FNEWT=WORK(4) 
         IF (FNEWT.LE.UROUND/TOLST) THEN 
            WRITE(6,*)' CURIOUS INPUT FOR WORK(4)=',WORK(4) 
            ARRET=.TRUE. 
         END IF 
      END IF 
 
C ---> QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST. <--- 
      IF (WORK(5).EQ.0.D0) THEN 
         QUOT1=1.D0 
      ELSE 
         QUOT1=WORK(5) 
      END IF 
      IF (WORK(6).EQ.0.D0) THEN 
         QUOT2=1.2D0 
      ELSE 
         QUOT2=WORK(6) 
      END IF 
      IF (QUOT1.GT.1.0D0.OR.QUOT2.LT.1.0D0) THEN 
         WRITE(6,*)' CURIOUS INPUT FOR WORK(5,6)=',QUOT1,QUOT2 
         ARRET=.TRUE. 
      END IF 
C -------------------------------------------------------  
 
C ---->    GRID WITH DISCONTINUITIES  <---- 
      XURO=100*UROUND*ABS(XEND) 
      IF (NGRID.GT.0) THEN 
         IF (GRID(NGRID)-XEND.GE.XURO) THEN 
            IF(IOUT.GT.0) WRITE(6,*) 
     &          ' GRID(NGRID) HAS TO BE <= XEND' 
            ARRET=.TRUE. 
         ENDIF 
         IF (ABS(GRID(NGRID)-XEND).GE.XURO) NGRID=NGRID+1 
      ELSE 
         NGRID=NGRID+1 
      END IF 
      GRID(NGRID)=XEND 
C -------------------------------------------------------  
 
C -------> MAXIMAL STEP SIZE <------- 
      IF (WORK(7).EQ.0.D0) THEN 
         DO I=1,NGRID 
          IF (GRID(I).GT.X) THEN  
           IGRID=I 
           GO TO 2 
          END IF 
         END DO 
 2       CONTINUE 
         HMAX=GRID(IGRID)-X 
      ELSE 
         HMAX=WORK(7) 
      END IF  
 
C ------->  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION <------- 
      IF(WORK(8).EQ.0.D0)THEN 
         FACL=5.D0 
      ELSE 
         FACL=1.D0/WORK(8) 
      END IF 
      IF(WORK(9).EQ.0.D0)THEN 
         FACR=1.D0/8.0D0 
      ELSE 
         FACR=1.D0/WORK(9) 
      END IF 
      IF (FACL.LT.1.0D0.OR.FACR.GT.1.0D0) THEN 
            WRITE(6,*)' CURIOUS INPUT WORK(8,9)=',WORK(8),WORK(9) 
            ARRET=.TRUE. 
      END IF 
C ------->  PARAMETER FOR THE CONTROL OF DENSE OUTPUT <------- 
      ALPHA=WORK(10) 
      IF (ALPHA.LT.0.D0.OR.ALPHA.GT.1.D0) THEN 
            WRITE(6,*)' CURIOUS INPUT WORK(10)=',WORK(10) 
            ARRET=.TRUE. 
      END IF 
C ------->   PARAMETER FOR CONTROLLING THE SEARCH OF BP <-------
      TCKBP=WORK(11)
        IF (TCKBP.LE.0.D0) THEN
              TCKBP=5.D0
        END IF               
   
C *** *** *** *** *** *** *** *** *** *** *** *** *** 
C         COMPUTATION OF ARRAY ENTRIES 
C *** *** *** *** *** *** *** *** *** *** *** *** *** 
C ---- IMPLICIT, BANDED OR NOT ? 
      IMPLCT=IMAS.NE.0 
       NEUTRAL=IMAS.EQ.2
      JBAND=MLJAC.LT.NM1 
C -------- COMPUTATION OF THE ROW-DIMENSIONS OF THE 2-ARRAYS --- 
C -- JACOBIAN  AND  MATRICES E1, E2 
      IF (JBAND) THEN 
         LDJAC=MLJAC+MUJAC+1 
         LDE1=MLJAC+LDJAC 
      ELSE 
         MLJAC=NM1 
         MUJAC=NM1 
         LDJAC=NM1 
         LDE1=NM1 
      END IF 
C -- MASS MATRIX 
      IF (IMPLCT) THEN 
          IF (MLMAS.NE.NM1) THEN 
              LDMAS=MLMAS+MUMAS+1 
              IF (JBAND) THEN 
                 IJOB=4 
              ELSE 
                 IJOB=3 
              END IF 
          ELSE 
              LDMAS=NM1 
              IJOB=5 
          END IF 
C ------ BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF "JAC" 
          IF (MLMAS.GT.MLJAC.OR.MUMAS.GT.MUJAC) THEN 
             WRITE (6,*) 'BANDWITH OF "MAS" NOT SMALLER THAN BANDWITH OF 
     & "JAC"' 
            ARRET=.TRUE. 
          END IF 
      ELSE 
          LDMAS=0 
          IF (JBAND) THEN 
             IJOB=2 
          ELSE 
             IJOB=1 
             IF (N.GT.2.AND.IWORK(1).NE.0) IJOB=7 
          END IF 
      END IF 
      LDMAS2=MAX(1,LDMAS) 
C ------ HESSENBERG OPTION ONLY FOR EXPLICIT EQU. WITH FULL JACOBIAN 
      IF ((IMPLCT.OR.JBAND).AND.IJOB.EQ.7) THEN 
         WRITE(6,*)' HESSENBERG OPTION ONLY FOR EXPLICIT EQUATIONS WITH  
     &FULL JACOBIAN' 
         ARRET=.TRUE. 
      END IF 
 
C ------ WHEN A FAIL HAS OCCURED, WE RETURN WITH IDID=-1 
      IF (ARRET) THEN 
         IDID=-1 
         RETURN 
      END IF 
 
C     NUMERICAL KERNEL
      WRITE(6,*) 'INTEGRATION...'
C -------- CALL TO CORE INTEGRATOR ------------ 
      CALL RADCOR(N,X,Y,XEND,H,FCN,PHI,ARGLAG,RTOL,ATOL,ITOL, 
     &   JAC,IJAC,MLJAC,MUJAC,JACLAG,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID, 
     &   NMAX,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,IJOB,STARTN, 
     &   NIND1,NIND2,NIND3,PRED,FACL,FACR,M1,M2,NM1, 
     &   IMPLCT,NEUTRAL,NDIMN,JBAND,LDJAC,LDE1,LDMAS2, 
     &   NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,NFULL,RPAR,IPAR, 
     &   IPAST,GRID,NRDS,NLAGS,NJACL, 
     &   NGRID,IEFLAG,WORK(7),TCKBP,ALPHA,ISWJL, LRPAST, PAST) 
      IWORK(13)=NFULL 
      IWORK(14)=NFCN 
      IWORK(15)=NJAC 
      IWORK(16)=NSTEP 
      IWORK(17)=NACCPT 
      IWORK(18)=NREJCT 
      IWORK(19)=NDEC 
      IWORK(20)=NSOL 
C -------- RESTORE TOLERANCES 
      EXPM=1.0D0/EXPM 
      IF (ITOL.EQ.0) THEN 
              QUOT=ATOL(1)/RTOL(1) 
              RTOL(1)=(10.0D0*RTOL(1))**EXPM 
              ATOL(1)=RTOL(1)*QUOT 
      ELSE 
          DO I=1,N 
              QUOT=ATOL(I)/RTOL(I) 
              RTOL(I)=(10.0D0*RTOL(I))**EXPM 
              ATOL(I)=RTOL(I)*QUOT 
          END DO 
      END IF 
C ----------- RETURN ----------- 
C     DEALLOCATE (IPOSV) 
      RETURN 
      END 
C 
C     END OF SUBROUTINE RADAR5 
C 
C *********************************************************** 
C 
 
      SUBROUTINE RADCOR(N,X,Y,XEND,H,FCN,PHI,ARGLAG,RTOL,ATOL,ITOL, 
     &   JAC,IJAC,MLJAC,MUJAC,JACLAG,MAS,MLMAS,MUMAS,SOLOUT,IOUT,IDID, 
     &   NMAX,SAFE,THET,FNEWT,QUOT1,QUOT2,NIT,IJOB,STARTN, 
     &   NIND1,NIND2,NIND3,PRED,FACL,FACR,M1,M2,NM1, 
     &   IMPLCT,NEUTRAL,NDIMN,BANDED,LDJAC,LDE1,LDMAS, 
     &   NFCN,NJAC,NSTEP,NACCPT,NREJCT,NDEC,NSOL,NFULL,RPAR,IPAR, 
     &   IPAST,GRID,NRDS,NLAGS,NJACL, 
     &   NGRID,IEFLAG,WORK7,TCKBP,ALPHA,ISWJL, LPAST, PAST) 
C ---------------------------------------------------------- 
C     CORE INTEGRATOR FOR RADAR5 
C     PARAMETERS SAME AS IN RADAR5 WITH WORKSPACE ADDED  
C ----------------------------------------------------------  
C         DECLARATIONS  
C ----------------------------------------------------------  
C     use definitions 
      USE IP_ARRAY 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
      REAL(kind=DP), dimension(1), intent(inout) :: Y 
      REAL(kind=DP), dimension(:), allocatable ::  
     &                              Z1,Z2,Z3,Y0,SCAL,F1,F2,F3,CONT 
      REAL(kind=DP), dimension(:), allocatable :: BPV,UCONT 
      REAL(kind=DP), dimension(:,:), allocatable ::  
     &                              FJAC,FJACS,FMAS,E1,E2R,E2I 
      REAL(kind=DP), dimension(1), intent(inout) ::  
     &                              ATOL,RTOL 
      REAL(kind=DP), dimension(1), intent(in) ::  
     &                              RPAR 
      INTEGER, dimension(1), intent(in) ::  
     &                              IPAR,IPAST 
      REAL(kind=DP), dimension(1), intent(in) ::  
     &                              GRID 
      REAL(kind=DP), dimension(:,:), allocatable :: 
     &                              FJACL,XLAG 
      REAL(kind=DP), dimension(:), allocatable :: 
     &                              FJACLAG,ZL  
      INTEGER, dimension(:,:), allocatable :: ICOUN 
      INTEGER, dimension(:), allocatable ::  
     &                              IVL,IVE,IVC,ILS 
      INTEGER, dimension(:), allocatable ::  
     &                              IP1,IP2,IPHES,IPJ
      INTEGER :: LPAST
      REAL(kind=DP), dimension(LPAST), intent(inout) :: PAST
 
      LOGICAL FLAGS,FLAGN,FLAGUS 
      LOGICAL QUADR 
      LOGICAL BPC,BPD,BPDMEM,LEFT 
      LOGICAL REPEAT
C ----> COMMON BLOCKS <---- 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
      COMMON /BPLOG/FIRST,LAST,REJECT,BPD 
      COMMON /BPCOM/BPP,ILBP,LEFT 
      COMMON /LINAL/MLE,MUE,MBJAC,MBB,MDIAG,MDIFF,MBDIAG 
 
C ---- 
      LOGICAL REJECT,FIRST,IMPLCT,NEUTRAL,PROJECT,BANDED,CALJAC 
      LOGICAL STARTN,CALHES,CALJACL,CALLAG 
      LOGICAL INDEX1,INDEX2,INDEX3,LAST,PRED 
      EXTERNAL FCN,PHI,ARGLAG 
C *** *** *** *** *** *** *** 
C  INITIALISATIONS 
C *** *** *** *** *** *** *** 
       
      ALLOCATE (Z1(N),Z2(N),Z3(N),Y0(N),SCAL(N),F1(N),F2(N),F3(N)) 
      ALLOCATE(BPV(10000)) 
      ALLOCATE (FJAC(LDJAC,N),ZL(3*N)) 
      IF (IMPLCT) ALLOCATE(FMAS(LDMAS,NM1)) 
      ALLOCATE (IP1(NM1),IP2(NM1),IPHES(NM1)) 
      ALLOCATE (E1(LDE1,NM1),E2R(LDE1,NM1),E2I(LDE1,NM1)) 
C      ALLOCATE (PAST(MXST*IDIF)) 
      IF (NLAGS.GT.0) THEN 
       ALLOCATE (FJACS(LDJAC,N),FJACLAG(NJACL)) 
       ALLOCATE (IVL(NJACL),IVE(NJACL),IVC(NJACL), 
     &           ILS(2*NLAGS+NJACL),ICOUN(3,NLAGS)) 
       IF (ISWJL.NE.1) THEN 
        ALLOCATE (IPJ(3*N),FJACL(3*N,3*N)) 
       END IF 
       ALLOCATE (XLAG(3,NLAGS)) 
      END IF 
  
C     AMPLITUDE OF CONT 
      LRC=4*N 
      ALLOCATE (CONT(LRC)) 
      ALLOCATE (UCONT(LRC+2)) 
C --- 
       OPEN(8,FILE='radar5.log')
       REWIND 8

C -------------------------------------------------  
      BPC=.FALSE. 
      BPD=.FALSE. 
      QUADR=.FALSE. 
       
C --- INITIAL PREPARATIONS 
      CALLAG=.FALSE. 
      IACT=1 
      IPOS=1 
      DO I=1,NLAGS 
       IPOSV(I)=1 
      END DO 
      X0B=X 
      DO I=1,NGRID 
       IF (GRID(I).GT.X0B) THEN  
        IGRID=I 
        GO TO 2 
       END IF 
      END DO 
 2    XEND=GRID(IGRID) 
      IBP=1 
      BPV(1)=X0B 
      BTOL=1.D1*UROUND 
      KMAX=10
      IMANT=0
       
C --- GUSTAFFSON TECHNIQUE AFTER BREAKING POINTS IS NOT APPLIED        
      FLAGUS=.FALSE. 
      ERRACC=1.D0
 
      IRTRN=2 
      CALL FCN(N,X,Y,Y0,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS,LPAST)
      IRTRN=1 
 
C     TOLMIN
      IF (ITOL.EQ.0) THEN
         RTOLM=RTOL(1)
        ELSE
       RTOLM=RTOL(1)
       DO I=2,N
        IF (RTOL(I).LT.RTOLM) RTOLM=RTOL(I)
         END DO
        END IF
 
C -------- CHECK THE INDEX OF THE PROBLEM -----  
      INDEX1=NIND1.NE.0 
      INDEX2=NIND2.NE.0 
      INDEX3=NIND3.NE.0 
C ------- COMPUTE MASS MATRIX FOR IMPLICIT CASE ---------- 
      IF (IMPLCT) CALL MAS(NM1,FMAS,LDMAS,RPAR,IPAR) 
C ---------> REQUIRED CONSTANTS <--------- 
      SQ6=DSQRT(6.D0) 
      C1=(4.D0-SQ6)/10.D0 
      C2=(4.D0+SQ6)/10.D0 
      C1M1=C1-1.D0 
      C2M1=C2-1.D0 
      C1MC2=C1-C2 
      CQ1=(2.D0+3.D0*SQ6)/6.D0 
      CQ2=(2.D0-3.D0*SQ6)/6.D0 
      CQ3=1.D0/3.D0 
      CL1=10.D0/(6.D0+SQ6)                   
      CL2=0.D0                 
      CL3=(-4.D0+SQ6)/(6.D0+SQ6)  
      CERS=5.D-1 
      CERC=5.D-1 
      CERLQ=1.D-2 
      THRS=100.D0 
      DD1=-(13.D0+7.D0*SQ6)/3.D0 
      DD2=(-13.D0+7.D0*SQ6)/3.D0 
      DD3=-1.D0/3.D0 
      U1=(6.D0+81.D0**(1.D0/3.D0)-9.D0**(1.D0/3.D0))/30.D0 
      ALPH=(12.D0-81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))/60.D0 
      BETA=(81.D0**(1.D0/3.D0)+9.D0**(1.D0/3.D0))*DSQRT(3.D0)/60.D0 
      CNO=ALPH**2+BETA**2 
      U1=1.0D0/U1 
      ALPH=ALPH/CNO 
      BETA=BETA/CNO 
      T11=9.1232394870892942792D-02 
      T12=-0.14125529502095420843D0 
      T13=-3.0029194105147424492D-02 
      T21=0.24171793270710701896D0 
      T22=0.20412935229379993199D0 
      T23=0.38294211275726193779D0 
      T31=0.96604818261509293619D0 
      TI11=4.3255798900631553510D0 
      TI12=0.33919925181580986954D0 
      TI13=0.54177053993587487119D0 
      TI21=-4.1787185915519047273D0 
      TI22=-0.32768282076106238708D0 
      TI23=0.47662355450055045196D0 
      TI31=-0.50287263494578687595D0 
      TI32=2.5719269498556054292D0 
      TI33=-0.59603920482822492497D0 
C 
C     INVERSE OF A 
      IF (NLAGS.GT.0.AND.ISWJL.NE.1) THEN 
       AI11= 3.22474487139158904909864D0 
       AI12= 1.16784008469040549492404D0 
       AI13=-0.25319726474218082618594D0 
       AI21=-3.56784008469040549492404D0 
       AI22= 0.77525512860841095090136D0 
       AI23= 1.05319726474218082618594D0 
       AI31= 5.53197264742180826185942D0 
       AI32=-7.53197264742180826185942D0 
       AI33= 5.00000000000000000000000D0 
      END IF 
C 
      IF (M1.GT.0) IJOB=IJOB+10 
      HMAXN=MIN(HMAX,XEND-X)  
      IF (H.LE.10.D0*UROUND) H=1.0D-6  
      H=MIN(H,HMAXN) 
      HOLD=H 
      REJECT=.FALSE. 
      FIRST=.TRUE. 
      LAST=.FALSE. 
        NITER=0
      IF ((X+H*1.0001D0-XEND).GE.0.D0) THEN 
         H=XEND-X 
         LAST=.TRUE. 
      END IF 
C ---  INITIALIZATION FOR THE ARRAY PAST    
       DO 3 I=0,MXST-1 
          PAST(1+IDIF*I)=X  
   3   CONTINUE 
       IPA=(MXST-1)*IDIF+1 
       DO J=1,NRDS 
          K=IPAST(J) 
          PAST(J+IPA)=Y(K) 
          PAST(J+1*NRDS+IPA)=0.D0 
          PAST(J+2*NRDS+IPA)=0.D0 
          PAST(J+3*NRDS+IPA)=0.D0 
       ENDDO 
          PAST(IPA+IDIF-1)=H  
C ---  END OF THE INITIALIZATION      
      FACCON=1.D0 
      CFAC=SAFE*(1+2*NIT) 
      NSING=0 
      XOLD=X 
      IF (IOUT.NE.0) THEN 
          IRTRN=1 
          NRSOL=1 
          XOSOL=XOLD 
          XSOL=X 
          CONT(1:N)=Y(1:N) 
          NSOLU=N 
          HSOL=HOLD
          CALL SOLOUT(NRSOL,XOSOL,XSOL,HSOL,Y,CONT,LRC,NSOLU, 
     &                RPAR,IPAR,IRTRN) 
          IF (IRTRN.LT.0) GOTO 179 
      END IF 
      MLE=MLJAC 
      MUE=MUJAC 
      MBJAC=MLJAC+MUJAC+1 
      MBB=MLMAS+MUMAS+1 
      MDIAG=MLE+MUE+1 
      MDIFF=MLE+MUE-MUMAS 
      MBDIAG=MUMAS+1 
      N2=2*N 
      N3=3*N 
      IF (ITOL.EQ.0) THEN 
          DO I=1,N 
             SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I)) 
          END DO 
      ELSE 
          DO I=1,N 
             SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I)) 
          END DO 
      END IF 
      HHFAC=H 
      NFCN=NFCN+1 
 
      NFULL=0 
C -------------------------- 
C --- BASIC INTEGRATION STEP   
C -------------------------- 
  10  CONTINUE 
C *** *** *** *** *** *** *** 
C  COMPUTATION OF THE JACOBIAN 
C *** *** *** *** *** *** *** 
C ----------------------- 
      FLAGS=.FALSE. 
      FLAGN=.FALSE. 
C ----------------------- 
      ALOPT=0.D0 
      NJAC=NJAC+1 
        IF (BPD) THEN
         BPDMEM=.TRUE.
         BPD=.FALSE.
      END IF
      IF (IJAC.EQ.0) THEN 
C --- COMPUTE JACOBIAN MATRIX NUMERICALLY 
         IF (BANDED) THEN 
C --- JACOBIAN IS BANDED 
            MUJACP=MUJAC+1 
            MD=MIN(MBJAC,M2) 
            DO MM=1,M1/M2+1 
               DO K=1,MD 
                  J=K+(MM-1)*M2 
 12               F1(J)=Y(J) 
                  F2(J)=DSQRT(UROUND*MAX(1.D-5,ABS(Y(J)))) 
                  Y(J)=Y(J)+F2(J) 
                  J=J+MD 
                  IF (J.LE.MM*M2) GOTO 12  
                  CALL FCN(N,X,Y,CONT,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,
     &                     NRDS,LPAST) 
                  J=K+(MM-1)*M2 
                  J1=K 
                  LBEG=MAX(1,J1-MUJAC)+M1 
 14               LEND=MIN(M2,J1+MLJAC)+M1 
                  Y(J)=F1(J) 
                  MUJACJ=MUJACP-J1-M1 
                  DO L=LBEG,LEND 
                     FJAC(L+MUJACJ,J)=(CONT(L)-Y0(L))/F2(J)  
                  END DO 
                  J=J+MD 
                  J1=J1+MD 
                  LBEG=LEND+1 
                  IF (J.LE.MM*M2) GOTO 14 
               END DO 
            END DO 
         ELSE 
C --- JACOBIAN IS FULL 
            DO I=1,N 
               YSAFE=Y(I) 
                 DELT=DSQRT(UROUND*MAX(1.D-5,ABS(YSAFE))) 
                         Y(I)=YSAFE+DELT 
                 CALL FCN(N,X,Y,CONT,ARGLAG,PHI, 
     1                      RPAR,IPAR,PAST,IPAST,NRDS,LPAST)
                 DO J=M1+1,N 
                 FJAC(J-M1,I)=(CONT(J)-Y0(J))/DELT 
               END DO 
               Y(I)=YSAFE 
            END DO
         END IF 
      ELSE 
C --- COMPUTE JACOBIAN MATRIX ANALYTICALLY
         CALL JAC(N,X,Y,FJAC,LDJAC,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS,
     &            LPAST)
      END IF
        IF (BPDMEM) THEN
         BPDMEM=.FALSE.
         BPD=.TRUE.
      END IF 
      CALJAC=.TRUE. 
      CALHES=.TRUE. 
      JLFLAG=0 
C --- SAVE FJAC 
      IF (NLAGS.GT.0) FJACS=FJAC 
 
C ------------------------------------------------ 
C --- GLOBAL ITERATION BEGINS HERE
C ------------------------------------------------ 
  20  CONTINUE 
 
C *** *** *** *** *** *** *** 
C  STARTING VALUES FOR NEWTON ITERATION 
C *** *** *** *** *** *** *** 
    
      IF (FIRST.OR.STARTN) THEN 
            DO I=1,N 
             Z1(I)=0.D0 
             Z2(I)=0.D0 
             Z3(I)=0.D0 
             F1(I)=0.D0 
             F2(I)=0.D0 
             F3(I)=0.D0 
C 
             A1=Y(I) 
             ZL(I)=A1 
             ZL(I+N)=A1 
             ZL(I+N2)=A1 
            END DO 
      ELSE 
         C3Q=H/HOLD 
         C1Q=C1*C3Q 
         C2Q=C2*C3Q 
         DO I=1,N 
            A1=Y(I) 
            AK1=CONT(I+N) 
            AK2=CONT(I+N2) 
            AK3=CONT(I+N3) 
            Z1I=C1Q*(AK1+(C1Q-C2M1)*(AK2+(C1Q-C1M1)*AK3)) 
            Z2I=C2Q*(AK1+(C2Q-C2M1)*(AK2+(C2Q-C1M1)*AK3)) 
            Z3I=C3Q*(AK1+(C3Q-C2M1)*(AK2+(C3Q-C1M1)*AK3)) 
            Z1(I)=Z1I 
            Z2(I)=Z2I 
            Z3(I)=Z3I 
C 
            ZL(I)=A1+Z1I 
            F1(I)=TI11*Z1I+TI12*Z2I+TI13*Z3I 
            F2(I)=TI21*Z1I+TI22*Z2I+TI23*Z3I 
            F3(I)=TI31*Z1I+TI32*Z2I+TI33*Z3I 
         END DO 
            IF (NLAGS.GT.0) THEN 
             DO I=1,N 
              A1=Y(I) 
              ZL(I+N)=A1+Z2(I) 
              ZL(I+N2)=A1+Z3(I) 
             END DO 
            END IF 
      END IF 
C --- 
      IF (JLFLAG.EQ.0) THEN 
       IF (NLAGS.EQ.0) THEN 
        IF (.NOT.(REJECT.OR.LAST)) THEN  
         IF (IMANT.EQ.1) GOTO 30 
        END IF 
        GOTO 22  
       END IF 
C --------------------- 
C      INITIALIZATION   
C 
C ---  LAGS CYCLE --- C 
C --- 
       CALJACL=.FALSE. 
       X1=X+C1*H 
       X2=X+C2*H 
       X3=X+H 
       DO IL=1,NLAGS 
        XLAG(1,IL)=0.D0 
        XLAG(2,IL)=0.D0 
        XLAG(3,IL)=0.D0 
        ICOUN(1,IL)=0 
        ICOUN(2,IL)=0 
        ICOUN(3,IL)=0 
       END DO 
C ---  LOOP ON LAG TERMS 
       DO IL=1,NLAGS 
C ---   DELAYED ARGUMENTS ARE COMPUTED 
        XLAG(1,IL)=ARGLAG(IL,X1,ZL,RPAR,IPAR,PHI,PAST,IPAST,NRDS,LPAST) 
           IF (XLAG(1,IL).GT.X) ICOUN(1,IL)=1 
        XLAG(2,IL)=ARGLAG(IL,X2,ZL(N+1),RPAR,IPAR,PHI,PAST,IPAST,NRDS,
     &          LPAST)
           IF (XLAG(2,IL).GT.X) ICOUN(2,IL)=1 
        XLAG(3,IL)=ARGLAG(IL,X3,ZL(N2+1),RPAR,IPAR,PHI,PAST,IPAST,NRDS,
     &          LPAST) 
           IF (XLAG(3,IL).GT.X) ICOUN(3,IL)=1 
        IF (ICOUN(1,IL)+ICOUN(2,IL)+ICOUN(3,IL).GE.1) CALJACL=.TRUE.
       END DO 
 
       IF (CALJACL) THEN 
        CALL JACLAG(N,X,Y,FJACLAG,ARGLAG,PHI,IVE,IVC,IVL, 
     &                  RPAR,IPAR,PAST,IPAST,NRDS,LPAST) 
        IF (.NOT.CALLAG) THEN 
         CALLAG=.TRUE. 
C --     ORDERING STEP 
         LL=2*NLAGS+1 
         DO L=1,NLAGS 
          NL=0 
          DO I=1,NJACL 
           IF (IVL(I).EQ.L) THEN 
            ILS(LL)=I 
            NL=NL+1 
            LL=LL+1 
           END IF 
          END DO 
          ILS(2*L-1)=NL 
          ILS(2*L)=LL-1 
         END DO 
        END IF 
C 
        DO IL=1,NLAGS 
          DO IS=1,3 
           SELECT CASE (IS) 
            CASE (1) 
             XACT=X1 
            CASE (2) 
             XACT=X2 
            CASE (3) 
             XACT=X3 
           END SELECT  
           IF (XLAG(IS,IL).GT.XACT) THEN 
            IF (IOUT.EQ.2) WRITE (6,*) 
     &          ' WARNING!: ADVANCED ARGUMENTS ARE USED AT X= ',XACT 
            XLAG(IS,IL)=XACT 
           END IF 
          END DO 
 
C ------ JACOBIAN MAINTAINS THE TENSOR STRUCTURE
C ------ UPDATING CONDITION
          ALOPT=0.D0 
          IF (ICOUN(1,IL).EQ.1) THEN 
             S1=DIM(XLAG(1,IL),X)/H 
             ALOPT=(-1.D0+S1)*S1*(-13.D0-7.D0*Sqrt(6.D0) + 
     &        5.D0*(2.D0+3.D0*Sqrt(6.D0))*S1) 
          END IF 
          IF (ICOUN(2,IL).EQ.1) THEN 
             S2=DIM(XLAG(2,IL),X)/H 
             ALOPT=ALOPT-(-1+S2)*S2* 
     &       (13.D0-7.D0*Sqrt(6.D0)+5.D0*(-2.D0+3.D0*Sqrt(6.D0))*S2) 
          END IF 
          IF (ICOUN(3,IL).EQ.1) THEN 
             S3=DIM(XLAG(3,IL),X)/H 
             ALOPT=ALOPT+S3*(1.D0 - 8.D0*S3 + 10.D0*S3**2) 
          END IF 
          ALOPT=ALOPT/9.D0         
          
C         OPTIMAL COEFFICIENT (W.R.T. FROBENIUS NORM)
C         JACLAG ~= ALOPT*I 
C 
C         ACTIVATES IF ALOPT DIFFERENT FROM ZERO
          IF (ABS(ALOPT).GE.1.D-8) THEN 
           NL =ILS(2*IL-1) 
           ILE=ILS(2*IL) 
           DO K=1,NL 
            KK=ILS(ILE-K+1) 
            IK=IVE(KK) 
            JK=IVC(KK) 
            FJAC(IK,JK)=FJAC(IK,JK)+ALOPT*FJACLAG(KK)
C                       ----S                   
           END DO 
          ELSE 
           CALJACL=.FALSE. 
          END IF 
        END DO 
C 
       ELSE 
C ---   FACTORIZATION IS CONSERVED  
       END IF 
       IF (.NOT.(REJECT.OR.LAST.OR.CALJACL)) THEN  
        IF (IMANT.EQ.1) GOTO 30 
       END IF 
       GOTO 22  
C --- 
      ELSE IF (JLFLAG.EQ.1) THEN 
C ---  THE FULL ITERATION MAKES USE OF THE SAME STEPSIZE OF THE SIMPLIFIED
       GOTO 23 
C --- 
      END IF 
 
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C - - THE DIFFERENCE IN THE SOLUTION OF THE SYSTEM STARTS HERE 
C ------------------------------- - - - - - - - - - - - - - - - 
C --- SIMPLE NEWTON ITERATION 
C ------------------------------- 
  22  CONTINUE 
C 
C --- COMPUTE THE MATRICES E1 AND E2 AND THEIR DECOMPOSITIONS 
      FAC1=U1/H 
      ALPHN=ALPH/H 
      BETAN=BETA/H 
C --- RK JACOBIAN FACTORIZATION 
      CALL DECOMR(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS, 
     &            M1,M2,NM1,FAC1,E1,LDE1,IP1,IER,IJOB,CALHES,IPHES) 
      IF (IER.NE.0) THEN 
          GOTO 78 
      END IF 
      CALL DECOMC(N,FJAC,LDJAC,FMAS,LDMAS,MLMAS,MUMAS, 
     &            M1,M2,NM1,ALPHN,BETAN,E2R,E2I,LDE1,IP2,IER,IJOB) 
      IF (IER.NE.0) THEN 
          GOTO 78 
      END IF 
      NDEC=NDEC+1 
C --- 
  30  CONTINUE 
C --- UPDATE NSTEP 
      NSTEP=NSTEP+1 
      IF (NSTEP.GT.NMAX) GOTO 178 
      IF (0.1D0*H.LE.ABS(X)*UROUND) GOTO 177 
          IF (INDEX2) THEN 
             DO I=NIND1+1,NIND1+NIND2 
                SCAL(I)=SCAL(I)/HHFAC 
             END DO 
          END IF 
          IF (INDEX3) THEN 
             DO I=NIND1+NIND2+1,NIND1+NIND2+NIND3 
                SCAL(I)=SCAL(I)/(HHFAC*HHFAC) 
             END DO 
          END IF 
      XPH=X+H 
 
C *** *** *** *** *** *** *** 
C  LOOP FOR THE SIMPLE NEWTON ITERATION 
C *** *** *** *** *** *** *** 
C              -----------      
            NEWT=0 
            IMANT=0 
C ----------------------- 
            FLAGS=.FALSE. 
            FLAGN=.FALSE. 
C ----------------------- 
            FACCON=MAX(FACCON,UROUND)**0.8D0 
            THETA=ABS(THET) 
C ------------------------------------------------------- 
C           REFERENCE POINT FOR THE SIMPLE AZIONE SIMPLE 
C ------------------------------------------------------- 
  40        CONTINUE  
C ------------------------- 
            IF (FLAGS) THEN 
             FLAGN=.TRUE. 
C ---        CALLS SUBROUTINE LAGR5 
C ------------------------------------------------------------- 
C ---------- DYNAMIC UPDATE OF INTERPOLANT (in PAST).  
C ------------------------------------------------------------- 
             DO J=1,NRDS 
                I=IPAST(J) 
                PAST(J+IACT)=Y(I)+Z3(I) 
                 Z2I=Z2(I) 
                 Z1I=Z1(I) 
                A1=(Z2I-Z3(I))/C2M1 
                PAST(J+1*NRDS+IACT)=A1 
                 AK=(Z1I-Z2I)/C1MC2 
                 ACONT3=Z1I/C1 
                 ACONT3=(AK-ACONT3)/C2 
                A2=(AK-A1)/C1M1  
                PAST(J+2*NRDS+IACT)=A2 
                PAST(J+3*NRDS+IACT)=A2-ACONT3 
             ENDDO 
C ---        UPDATE DI PAST 
             PAST(IACT)=X 
C            LEFT ENDPOINT 
             PAST(IACT+IDIF-1)=H 
C            STEPSIZE 
            END IF 
C ---------------- 
            IF (NEWT.GE.NIT) THEN 
             INREJ=2 
             GOTO 421 
            END IF 
C ----------------------------------- 
C ---     COMPUTE THE RIGHT-HAND SIDE 
            DO I=1,N 
             A1=Y(I) 
             ZL(I)=A1+Z1(I) 
             ZL(I+N)=A1+Z2(I) 
             ZL(I+N2)=A1+Z3(I) 
            END DO 
C           COMPUTATION OF STAGE VALUES
            CALL FCN(N,X+C1*H,ZL,Z1,ARGLAG,PHI,RPAR,IPAR,PAST,
     &               IPAST,NRDS,LPAST) 
            CALL FCN(N,X+C2*H,ZL(N+1),Z2,ARGLAG,PHI,RPAR,IPAR,PAST,
     &               IPAST,NRDS,LPAST)
            IF (BPD) THEN
C ---------------------------------------------------------------------- 
C ---         A BREAKING POINT HAS BEEN DETECTED
C ---------------------------------------------------------------- 
              LEFT=.FALSE. 
              DBP=ARGLAG(ILBP,X+C2*H,ZL(N+1),RPAR,IPAR,PHI, 
     &                   PAST,IPAST,NRDS,LPAST) 
              IF (DBP.LT.BPP) LEFT=.TRUE. 
            ELSE IF (LAST) THEN 
C ---         DISC. FLAG 
              XX=(X+H)*(1.D0-BTOL) 
              DO I=1,N 
                 ZL(N2+I)=(1.D0-BTOL)*ZL(N2+I)+BTOL*ZL(N+I) 
              END DO 
              CALL FCN(N,XX,ZL(N2+1),Z3,ARGLAG,PHI,RPAR,IPAR,PAST, 
     &                 IPAST,NRDS,LPAST) 
              GO TO 42 
            END IF 
            CALL FCN(N,X+H,ZL(N2+1),Z3,ARGLAG,PHI,RPAR,IPAR,PAST, 
     &               IPAST,NRDS,LPAST) 
 42         CONTINUE
    
            NFCN=NFCN+3 
C ----------------------------------- 
C ---     SOLVE THE LINEAR SYSTEMS 
           DO I=1,N 
              A1=Z1(I) 
              A2=Z2(I) 
              A3=Z3(I) 
              Z1(I)=TI11*A1+TI12*A2+TI13*A3 
              Z2(I)=TI21*A1+TI22*A2+TI23*A3 
              Z3(I)=TI31*A1+TI32*A2+TI33*A3 
           END DO 
            
        CALL SLVRAD(N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS, 
     &          M1,M2,NM1,FAC1,ALPHN,BETAN,E1,E2R,E2I,LDE1,Z1,Z2,Z3, 
     &          F1,F2,F3,CONT,IP1,IP2,IPHES,IER,IJOB) 
            NSOL=NSOL+1 
            NEWT=NEWT+1 
C ---       NORM OF DY 
                  DYNO=0.D0 
              DO I=1,N
            END DO
                  DO I=1,N 
               DENOM=SCAL(I) 
               DYNO=DYNO+(Z1(I)/DENOM)**2+(Z2(I)/DENOM)**2 
     &          +(Z3(I)/DENOM)**2 
            END DO 
            DYNO=DSQRT(DYNO/N3) 
C -------------------------------------------------------- 
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE 
C --------------------------------------------------------
            IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN 
                THQ=DYNO/DYNOLD 
                IF (NEWT.EQ.2) THEN 
                   THETA=THQ 
                ELSE 
                   THETA=SQRT(THQ*THQOLD) 
                END IF 
                THQOLD=THQ 
                INREJ=0 
C --- 1 
                IF (THETA.LT.0.99D0) THEN 
                    FACCON=THETA/(1.0D0-THETA) 
                    DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT 
                    IF (DYTH.GE.1.0D0) INREJ=1 
C ----------------  SLOW CONVERGENCE --- 
                ELSE 
                    INREJ=2 
C ----------------  DIVERGENCE --- 
                END IF 
C --- 1 
            ELSE 
C ---           NEWT=1 
                INREJ=0 
            END IF 
C ---------------------------------------- 
C ------------------ THE STEP IS REPEATED   
C ---------------------------------------- 
C ---- 2 
 421        IF (INREJ.GT.0) THEN 
C 
C ----- 3 
             IF (.NOT.(BPC.OR.FIRST)) THEN    
              IF (BPD) THEN
C ---          BP IS WRONG                      
               IBP=IBP-1
               BPD=.FALSE.
              END IF
              HP=H*0.99D0
              CALL BPDTCT(N,X,HP,Y,ARGLAG,RPAR,IPAR,UCONT,GRID,NLAGS, 
     &                    FIRST,LAST,XEND,IGRID,BPV,IBP,ILBP,BPP,BPD, 
     &                    KMAX,PHI,PAST,IPAST,NRDS,LPAST) 
            
              BPC=.TRUE. 
C ------ 4 
              IF (BPD) THEN 
               H=HP 
               LAST=.TRUE.
              ELSE IF (INREJ.EQ.1) THEN 
               QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH)) 
               HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT)) 
               H=HHFAC*H 
               LAST=.FALSE. 
C
              ELSE IF (INREJ.EQ.2) THEN 
               H=H*0.55D0 
               HHFAC=0.55D0 
               LAST=.FALSE. 
C
              END IF 
C ------ 4 
             ELSE   
C ----- 3 
              IF (BPD) THEN  
C              BP WRONG
               IBP=IBP-1 
               BPD=.FALSE. 
               LAST=.FALSE. 
              END IF 
C ------------------------------------------- 
              IF (.NOT.CALJACL.OR.ISWJL.EQ.1) THEN  
C ---          THERE ARE NOT SMALL DELAYS 
               IF (REJECT.AND.CALJAC) THEN  
                H=H*0.12D0 
                HHFAC=0.12D0 
               ELSE IF (INREJ.EQ.1) THEN 
                QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH)) 
                HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT)) 
                H=HHFAC*H 
                LAST=.FALSE. 
               ELSE IF (INREJ.EQ.2) THEN 
                H=H*0.55D0 
                HHFAC=0.55D0 
                LAST=.FALSE. 
               END IF 
              ELSE 
C ---          CALJACL IS TRUE: FULL ITEARATION IS DONE          
C ---          THE STEPSIZE DOES NOT CHANGET
               JLFLAG=1 
              END IF 
             END IF 
C ----- 3 
             REJECT=.TRUE. 
             IF (CALJAC) GOTO 20 
             GOTO 10 
            END IF 
C ---- 2 
C -------------------------------------------------------- 
            DYNOLD=MAX(DYNO,UROUND)
            DO I=1,N 
               F1I=F1(I)+Z1(I) 
               F2I=F2(I)+Z2(I) 
               F3I=F3(I)+Z3(I) 
               F1(I)=F1I 
               F2(I)=F2I 
               F3(I)=F3I 
               Z1(I)=T11*F1I+T12*F2I+T13*F3I 
               Z2(I)=T21*F1I+T22*F2I+T23*F3I 
               Z3I=T31*F1I+    F2I 
               Z3(I)=Z3I 
C ---          APPROX DELLA SOLUZIONE 
               ZL(I+N2)=Y(I)+Z3I 
            END DO 
C -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
            IF (NEWT.EQ.1.OR.FACCON*DYNO.GT.FNEWT) THEN 
C ---        NEWTON PROCEDE
             GOTO 40 
            END IF 

C ----------------------------------------------------------------
C ---       ITERATIVE CORRECTION OF THE BREAKING POINT
C ----------------------------------------------------------------
            IF (BPD) THEN
             HNEWT=H
             NITER = NITER+1
             IF (NITER.GT.NIT) THEN  
C              BP WRONG                 
               IBP=IBP-1 
               BPD=.FALSE. 
               LAST=.FALSE. 
               H=H*0.55D0 
               HHFAC=0.55D0
               REJECT=.TRUE.
               NITER=0 
               IF (CALJAC) GOTO 20 
               GOTO 10   
             END IF  
                    
             CALL BPACC(N,X,H,Y,ARGLAG,RPAR,IPAR,Z1,Z2,Z3,FIRST,
     &                  BPV,IBP,ILBP,BPP,KMAX,PHI,PAST,IPAST,NRDS,LPAST) 
             IF (ABS(H-HNEWT)/HNEWT.GE.MAX(BTOL,RTOLM*1.D-2)) THEN
              GOTO 20
C             REF POINT
             ELSE
                H=HNEWT
                NITER=0
             END IF
            END IF  
C ----------------------------------------------------------------


C *** *** *** *** *** *** *** *** *** *** *** 
C END LOOP 
C *** *** *** *** *** *** *** *** *** *** *** 
      GOTO 55 
C *** *** *** *** *** 
C 
 
C --- FULL NEWTON ITERATION 
  23  CONTINUE 
      NFULL =NFULL +1 
C --- 
 
C ---------------------------------------------------- 
C --- ALTERNATIVE FULL NEWTON JACOB. INTEGRATION STEP   
C ---------------------------------------------------- 
C 
      DO I=1,N 
       DO J=1,N 
        FJACL(I+N,J)=0.D0 
        FJACL(I+2*N,J)=0.D0 
        FJACL(I+2*N,J+N)=0.D0 
        FJACL(I,J+N)=0.D0 
        FJACL(I,J+2*N)=0.D0 
        FJACL(I+N,J+2*N)=0.D0 
        FIJ=FJACS(I,J) 
        FJACL(I,J)=FIJ 
        FJACL(I+N,J+N)=FIJ 
        FJACL(I+2*N,J+2*N)=FIJ 
       END DO 
      END DO 
 
      QUADR=FIRST 
      IF (.NOT.QUADR) THEN 
       DL1=C1*(C1-C2)*(C1-1.D0) 
       DL2=C2*(C2-C1)*(C2-1.D0) 
       DL3=(1.D0-C1)*(1.D0-C2) 
      ELSE 
       DL1=(C1-C2)*(C1-1.D0) 
       DL2=(C2-C1)*(C2-1.D0) 
       DL3=(1.D0-C1)*(1.D0-C2) 
      ENDIF 
      DO IL=1,NLAGS 
       DO I=1,3 
            XL=XLAG(I,IL) 
             IF (XL.GT.X) THEN 
CCC 
CCC            DERIVATIVES OF THE COLLOCATION POLYNOMIAL WRT STAGES 
CCC            d u/d Y_k = L_k(xlag_i) 
CCC                                    
               IF (.NOT.QUADR) THEN 
                DCOLI1=((XL-X)/H)*((XL-X2)/H)*((XL-X3)/H)/DL1 
                DCOLI2=((XL-X)/H)*((XL-X1)/H)*((XL-X3)/H)/DL2 
                DCOLI3=((XL-X)/H)*((XL-X1)/H)*((XL-X2)/H)/DL3 
               ELSE 
                DCOLI1=((XL-X2)/H)*((XL-X3)/H)/DL1 
                DCOLI2=((XL-X1)/H)*((XL-X3)/H)/DL2 
                DCOLI3=((XL-X1)/H)*((XL-X2)/H)/DL3 
               END IF 
 
C -----------> JACOBIAN IS UPDATED 
               NL=ILS(2*IL-1) 
               ILE=ILS(2*IL) 
 
C -----------> FULL JACOBIAN MATRIX FJACL 
               DO K=1,NL 
                 KK=ILS(ILE-K+1) 
                 IK=IVE(KK) 
                 JK=IVC(KK) 
                 FJLK=FJACLAG(KK) 
                 IKI=IK+(I-1)*N 
C 
                  FJACL(IKI,JK)=FJACL(IKI,JK)+FJLK*DCOLI1 
                  FJACL(IKI,JK+N)=FJACL(IKI,JK+N)+FJLK*DCOLI2    
                  FJACL(IKI,JK+2*N)=FJACL(IKI,JK+2*N)+FJLK*DCOLI3 
               END DO 
             END IF 
       END DO            
CCC 
CCC --> NLAGS 
      END DO 
CCC <-- 
 
      AI11H   =-AI11/H 
      AI12H   =-AI12/H 
      AI13H   =-AI13/H 
      AI21H   =-AI21/H 
      AI22H   =-AI22/H 
      AI23H   =-AI23/H 
      AI31H   =-AI31/H 
      AI32H   =-AI32/H 
      AI33H   =-AI33/H 
 
C --- FJACL 
      IF (IMPLCT) THEN 
       DO I1=1,N 
         DO J1=1,N 
           FJACL(I1,J1)=FJACL(I1,J1)+AI11H*FMAS(I1,J1) 
           FJACL(I1,J1+N)=FJACL(I1,J1+N)+AI12H*FMAS(I1,J1) 
           FJACL(I1,J1+2*N)=FJACL(I1,J1+2*N)+AI13H*FMAS(I1,J1) 
 
           FJACL(I1+N,J1)=FJACL(I1+N,J1)+AI21H*FMAS(I1,J1) 
           FJACL(I1+N,J1+N)=FJACL(I1+N,J1+N)+AI22H*FMAS(I1,J1) 
           FJACL(I1+N,J1+2*N)=FJACL(I1+N,J1+2*N)+AI23H*FMAS(I1,J1) 
 
           FJACL(I1+2*N,J1)=FJACL(I1+2*N,J1)+AI31H*FMAS(I1,J1) 
           FJACL(I1+2*N,J1+N)=FJACL(I1+2*N,J1+N)+AI32H*FMAS(I1,J1) 
           FJACL(I1+2*N,J1+2*N)=FJACL(I1+2*N,J1+2*N)+AI33H*FMAS(I1,J1) 
         END DO 
       END DO 
      ELSE 
C ---  EXPLICIT CASE 
       DO I1=1,N 
         FJACL(I1,I1)=FJACL(I1,I1)+AI11H 
         FJACL(I1,I1+N)=FJACL(I1,I1+N)+AI12H 
         FJACL(I1,I1+2*N)=FJACL(I1,I1+2*N)+AI13H 
 
         FJACL(I1+N,I1)=FJACL(I1+N,I1)+AI21H 
         FJACL(I1+N,I1+N)=FJACL(I1+N,I1+N)+AI22H 
         FJACL(I1+N,I1+2*N)=FJACL(I1+N,I1+2*N)+AI23H 
 
         FJACL(I1+2*N,I1)=FJACL(I1+2*N,I1)+AI31H 
         FJACL(I1+2*N,I1+N)=FJACL(I1+2*N,I1+N)+AI32H 
         FJACL(I1+2*N,I1+2*N)=FJACL(I1+2*N,I1+2*N)+AI33H 
       END DO 
      END IF 
       
 
  
C ----- FACTORIZATION OF THE FULL JACOBIAN 
      CALL DEC(3*N,3*LDJAC,FJACL,IPJ,IER)  
      IF (IER.NE.0) THEN 
       GOTO 78 
      END IF 
CCC --->                                
      NDEC=NDEC+1 
  33  CONTINUE 
C --- 
C --- R: EVERY STEP STARTS WITH A SIMPLE ITERATION
      NSTEP=NSTEP+1 
      IF (NSTEP.GT.NMAX) GOTO 178 
      IF (0.1D0*H.LE.ABS(X)*UROUND) GOTO 177 
C --- 
      XPH=X+H 
 
C *** *** *** *** *** *** *** 
C  LOOP FOR NEWTON ITERATION 
C *** *** *** *** *** *** *** 
C              -----------      
            NEWT=0 
C ----------------------- 
            FLAGS=.FALSE. 
            FLAGN=.FALSE. 
C ----------------------- 
            FACCON=MAX(FACCON,UROUND)**0.8D0 
            THETA=ABS(THET) 
CCC --- --- --- --- --- --- --- --- --- --- --- --- --- 
CCC         REFERENCE POINT FOR FULL ITERATION          
CCC --- --- --- --- --- --- --- --- --- --- --- --- --- 
  43        CONTINUE 
C --- --- --- --- --- --- --- --- --- --- --- --- 
            IF (FLAGS) THEN 
             FLAGN=.TRUE. 
C ***************** 
CCC ---      DYNAMIC UPDATE OF THE CURRENT INTETPOLANT (in PAST)        
             DO J=1,NRDS 
                I=IPAST(J) 
                PAST(J+IACT)=Y(I)+Z3(I) 
                 Z2I=Z2(I) 
                 Z1I=Z1(I) 
                PAST(J+1*NRDS+IACT)=(Z2I-Z3(I))/C2M1 
                 AK=(Z1I-Z2I)/C1MC2 
                 ACONT3=Z1I/C1 
                 ACONT3=(AK-ACONT3)/C2 
                PAST(J+2*NRDS+IACT)=(AK-PAST(J+1*NRDS+IACT))/C1M1 
                PAST(J+3*NRDS+IACT)=PAST(J+2*NRDS+IACT)-ACONT3 
             ENDDO 
CCC          UPDATE 
             PAST(IACT)=X 
CCC          LEFT ENDPOINT OF CURRENT INTERVAL 
             PAST(IACT+IDIF-1)=H 
CCC          USED STEPSIZE                                                 
            END IF 
C --- --- --- --- --- --- --- --- --- --- --- --- 
            IF (NEWT.GE.NIT) THEN  
             INREJ=2 
             GOTO 431 
C ---------> UNEXPECTED STEP-REJECTION 
            END IF 
C ---     COMPUTE THE RIGHT-HAND SIDE 
            CONT(1:N)=Y(1:N)+Z1(1:N) 
            CALL FCN(N,X+C1*H,CONT,F1,ARGLAG,PHI,RPAR,IPAR,PAST,
     &               IPAST,NRDS,LPAST) 
            CONT(1:N)=Y(1:N)+Z2(1:N) 
            CALL FCN(N,X+C2*H,CONT,F2,ARGLAG,PHI,RPAR,IPAR,PAST,
     &               IPAST,NRDS,LPAST)
            CONT(1:N)=Y(1:N)+Z3(1:N) 
            CALL FCN(N,XPH,CONT,F3,ARGLAG,PHI,RPAR,IPAR,PAST,
     &               IPAST,NRDS,LPAST)
            NFCN=NFCN+3 
C 
CCC --->    RHS COMPUTATION 
            IF (IMPLCT) THEN 
             ZL=0.D0 
             DO I1=1,N 
               DO J1=1,N 
                 ZL(I1)=ZL(I1)+AI11H*FMAS(I1,J1)*Z1(J1) 
                 ZL(I1)=ZL(I1)+AI12H*FMAS(I1,J1)*Z2(J1) 
                 ZL(I1)=ZL(I1)+AI13H*FMAS(I1,J1)*Z3(J1) 
 
                 ZL(I1+N)=ZL(I1+N)+AI21H*FMAS(I1,J1)*Z1(J1) 
                 ZL(I1+N)=ZL(I1+N)+AI22H*FMAS(I1,J1)*Z2(J1) 
                 ZL(I1+N)=ZL(I1+N)+AI23H*FMAS(I1,J1)*Z3(J1) 
 
                 ZL(I1+2*N)=ZL(I1+2*N)+AI31H*FMAS(I1,J1)*Z1(J1) 
                 ZL(I1+2*N)=ZL(I1+2*N)+AI32H*FMAS(I1,J1)*Z2(J1) 
                 ZL(I1+2*N)=ZL(I1+2*N)+AI33H*FMAS(I1,J1)*Z3(J1) 
               END DO 
             END DO 
            ELSE 
             ZL(1:N)=AI11H*Z1(1:N)+AI12H*Z2(1:N)+AI13H*Z3(1:N) 
             ZL(N+1:2*N)=AI21H*Z1(1:N)+AI22H*Z2(1:N)+AI23H*Z3(1:N) 
             ZL(2*N+1:3*N)=AI31H*Z1(1:N)+AI32H*Z2(1:N)+AI33H*Z3(1:N) 
            END IF 
C 
            ZL(1:N)=-ZL(1:N)-F1(1:N) 
            ZL(N+1:2*N)=-ZL(N+1:2*N)-F2(1:N) 
            ZL(2*N+1:3*N)=-ZL(2*N+1:3*N)-F3(1:N) 
               
C --------> SOLVE THE LINEAR SYSTEMS 
            CALL SOL(3*N,3*LDJAC,FJACL,ZL,IPJ) 
            NSOL=NSOL+1 
            NEWT=NEWT+1 
            DYNO=0.D0 
            DO I=1,N 
               DENOM=SCAL(I) 
               DYNO=DYNO+(ZL(I)/DENOM)**2+(ZL(I+N)/DENOM)**2 
     &          +(ZL(I+2*N)/DENOM)**2 
            END DO 
            DYNO=DSQRT(DYNO/N3) 
C -------------------------------------------------------- 
C ---     BAD CONVERGENCE OR NUMBER OF ITERATIONS TO LARGE 
C -------------------------------------------------------- 
              IF (NEWT.GT.1.AND.NEWT.LT.NIT) THEN 
                THQ=DYNO/DYNOLD 
                IF (NEWT.EQ.2) THEN 
                   THETA=THQ 
                ELSE 
                   THETA=SQRT(THQ*THQOLD) 
                END IF 
                THQOLD=THQ 
                INREJ=0 
C --- 1 
                IF (THETA.LT.0.99D0) THEN 
                    FACCON=THETA/(1.0D0-THETA) 
                    DYTH=FACCON*DYNO*THETA**(NIT-1-NEWT)/FNEWT 
                    IF (DYTH.GE.1.0D0) INREJ=1 
C ----------------  CONVERGENZA LENTA --- 
                ELSE 
                    INREJ=2 
C ----------------  DIVERGENZA --- 
                END IF 
C --- 1 
            ELSE 
                INREJ=0 
            END IF 
C ---------------------------------------- 
C ------------------ THE STEP IS REPEATED  
C ---------------------------------------- 
 431        IF (INREJ.GT.0) THEN 
             REJECT=.TRUE. 
             LAST=.FALSE.
             JLFLAG=0 
C -------------------------------------------- 
C ---        TURN TO A SIMPLE ITERATION         
             IF (INREJ.EQ.1) THEN 
                QNEWT=DMAX1(1.0D-4,DMIN1(20.0D0,DYTH)) 
                HHFAC=.8D0*QNEWT**(-1.0D0/(4.0D0+NIT-1-NEWT)) 
                H=HHFAC*H 
             ELSE IF (INREJ.EQ.2) THEN 
                H=H*0.55D0 
                HHFAC=0.55D0 
             END IF 
 
             IF (CALJAC) GOTO 20 
             GOTO 10 
            END IF 
C -------------------------------------------------------- 
            DYNOLD=MAX(DYNO,UROUND) 
C --        UPDATE OF Z VALUES  
            Z1(1:N)=Z1(1:N)+ZL(1:N) 
            Z2(1:N)=Z2(1:N)+ZL(N+1:2*N) 
            Z3(1:N)=Z3(1:N)+ZL(2*N+1:3*N) 
C -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
            IF (FACCON*DYNO.GT.FNEWT) THEN 
             GOTO 43 
            END IF 
C - - - - - END FULL NEWTON ITERATION
C ----------------------------------------------------------------- 
C -----------------------------------------------------------------
C          RK EQUATIONS SUCCESFULLY SOLVED 
C -----------------------------------------------------------------
C ----------------------------------------------------------------- 
C *** *** *** *** *** *** *** *** *** *** *** 
C END LOOP 
C *** *** *** *** *** *** *** *** *** *** *** 
 
C -----------------------------------------------------------------
  55  CONTINUE 
C ******************** 
C --- ERROR ESTIMATION   
C ******************** 
C 
CCC    ERROR ESTIMATES                                        
       CALL ESTRAD (N,FJAC,LDJAC,MLJAC,MUJAC,FMAS,LDMAS,MLMAS,MUMAS, 
     &              H,U1,DD1,DD2,DD3,CL1,CL3,CQ1,CQ2,CQ3,CERLQ, 
     &              FCN,NFCN,Y0,Y,IJOB,X,M1,M2,NM1,E1,LDE1,ALPHA, 
     &              Z1,Z2,Z3,CONT,F1,F2,F3,IP1,IPHES,SCAL,SERR,CERR, 
     &              FIRST,REJECT,FAC1,ARGLAG,PHI,RPAR,IPAR, 
     &              IOUT,PAST,IPAST,NRDS,JLFLAG,IEFLAG,LPAST) 
 
      FAC=MIN(SAFE,CFAC/(NEWT+2*NIT)) 
 
      IF (FIRST) THEN 
C ------------------------------------------------------------ 
C ---  AFTER A GRID OR BREAKING POINT            
C --------------------------------------------- 
       ERR=SERR
C ---  WE REQUIRE .2<=HNEW/H<=8. 
       QUOT=MAX(FACR,MIN(FACL,ERR**0.25D0/FAC)) 
       HNEW=H/QUOT 
      ELSE 
C   ---------------------------------------------------------- 
       IF (IEFLAG.EQ.-1) THEN 
        CERR2=H*CERR 
C ---    
        ERR=CERR2 
       ELSE IF (IEFLAG.EQ.1) THEN 
C -----    
        CERR2=H*CERR 
        ERR=CERS*SERR+CERC*CERR2 
       ELSE IF (IEFLAG.EQ.2) THEN 
C ----- STANDARD DISCRETE ERROR 
        ERR=SERR 
       ELSE IF (IEFLAG.EQ.3) THEN 
C ----- NOT AVAILABLE AT THE MOMENT 
        CERR2=H*CERR 
        ERR=CERS*SERR+CERC*CERR2 
       END IF 
C ------------------------------------ 
C ---  COMPUTATION OF HNEW 
C ---  AS PREVIOUSLY COMPUTED: 
C ------------------------------------ 
C ---  WE REQUIRE .2<=HNEW/H<=8. 
C --- 
C ---  LINEAR COMBINATION OF BOTH ERROR COMPONENTS
       QUOT=MAX(FACR,MIN(FACL,ERR**.25D0/FAC)) 
 
       HNEW=H/QUOT 
      END IF 
C --- 
C *** *** *** *** *** *** *** 
C  DOES THE ERROR INCREASE MUCH? 
C *** *** *** *** *** *** *** 
      REPEAT=.FALSE.
      IF (ERR.GE.1.D0.OR.((ERR/ERRACC.GE.TCKBP).AND.(.NOT.BPD))) THEN 
C ---   KIND OF STEPSIZE
        IF (BPD) THEN
C ---     BP IS WRONG! REPEAT
          IBP=IBP-1
          BPD=.FALSE.
          REPEAT=.TRUE.
          IF (ERR.LT.1.D0) HNEW=H*0.55D0
        ELSE
C ---     LOOK FOR A BP        
          HP=H*0.99D0
          CALL BPDTCT(N,X,HP,Y,ARGLAG,RPAR,IPAR,UCONT,GRID,NLAGS, 
     &                FIRST,LAST,XEND,IGRID,BPV,IBP,ILBP,BPP,BPD,
     &                KMAX,PHI,PAST,IPAST,NRDS,LPAST) 
          IF (BPD) REPEAT=.TRUE.
        END IF
      END IF 
      BPC=.FALSE. 
C ---
C *** *** *** *** *** *** *** 
C  IS THE ERROR SMALL ENOUGH ? 
C *** *** *** *** *** *** *** 
C     IF (ERR.LT.1.D0) THEN 
 56   IF (.NOT.REPEAT.AND.ERR.LT.1.D0) THEN
C -------------------- 
C --- STEP IS ACCEPTED   
C -------------------- 
         FLAGN=.FALSE. 
         CALJACL=.FALSE. 
C 
C ---    UPDATE NUMBER OF ACCEPTED STEPS              
         NACCPT=NACCPT+1 
         IF (PRED) THEN 
C --------> PREDICTIVE CONTROLLER OF GUSTAFSSON 
            IF (NACCPT.GT.1) THEN 
             IF (FLAGUS) THEN 
               FACGUS=(HACC/H)*(ERR**2/ERRACC)**0.25D0/SAFE 
               FACGUS=MAX(FACR,MIN(FACL,FACGUS)) 
               QUOT=MAX(QUOT,FACGUS) 
               HNEW=H/QUOT 
             ELSE 
               FLAGUS=.TRUE. 
             END IF 
            END IF 
            HACC=H 
         END IF 
         ERRACC=MAX(1.0D-2,ERR)
C        ERRACC=ERR 
         XOLD=X 
         HOLD=H 
         X=XPH  

C ---    AGGIORNAMENTO DELLA SOLUZIONE 
         DO I=1,N 
            Z3I=Z3(I) 
            YI=Y(I)+Z3I   
            Y(I)=YI 
C ------------------------------------------------------------------------- 
C --------- UPDATE THE STAGES AND SOLUTION 
C ------------------------------------------------------------------------- 
C --- 
            CONT(I)=YI 
C ---    
              Z2I=Z2(I) 
              Z1I=Z1(I) 
C ------------------------------------------------------------------- 
C ---         INVERSE DEVIDED DIFFERENCES
C ------------------------------------------------------------------- 
              A1=(Z2I-Z3I)/C2M1 
              CONT(I+N)=A1 
              AK=(Z1I-Z2I)/C1MC2 
              ACONT3=Z1I/C1 
              ACONT3=(AK-ACONT3)/C2 
              A2=(AK-CONT(I+N))/C1M1 
              CONT(I+N2)=A2 
            IF (.NOT.FIRST) THEN 
              CONT(I+N3)=A2-ACONT3 
            ELSE 
C ---         QUADRATIC APPROXIMATION
              CONT(I+N3)=0.D0                        
C ---         INVECE DI:  
C ---         CONT(I+N3)=CONT(I+N2)-ACONT3 
            END IF 
         END DO 
C ---------------------------------------- 
C ---    SAVE LAST ACCEPTED STEP INFORMATION   
         DO I=1,LRC 
          UCONT(I)=CONT(I) 
         END DO 
         UCONT(LRC+1)=X 
         UCONT(LRC+2)=H 
C ---    FOR POSSIBLE SEARCH OF BREAKING POINTS     
C ------------------------------------------------- 
 
C ------------------------------------------------------------------ 
C ---    STEP IS ACCEPTED> DENSE OUTPUT IS STORED IN PAST        
C ---    
         DO J=1,NRDS 
            I=IPAST(J) 
            PAST(J+IACT)=CONT(I) 
            PAST(J+1*NRDS+IACT)=CONT(I+N) 
            PAST(J+2*NRDS+IACT)=CONT(I+N2) 
            IF (.NOT.FIRST) THEN 
              PAST(J+3*NRDS+IACT)=CONT(I+N3) 
            ELSE 
C ---       QUADRATIC APPROXIMATION                        
              PAST(J+3*NRDS+IACT)=0.D0 
            END IF 
         ENDDO 
C ---------> AGGIORNAMENTO DI PAST <--------- 
         PAST(IACT)=XOLD 
C ---    
         IACT=IACT+IDIF 
C ---    POINTER TO NEXT STEP 
         PAST(IACT-1)=H 
C ---    
         IF (IACT+IDIF-1.GT.MXST*IDIF) IACT=1  
C ---    CONTROL ON THE MEMORY DIMENSION
C -------------------------------------------------------------- 

C ----------------------------------------------------------------- 
C ---    COMPUTATION AT BP FOR NEXT STEP
         IF (LAST) THEN 
C ---     LAST HAS TO BE RE/DEFINED
          IF (BPD) THEN 
           IGRID=IGRID-1
          END IF 
C --- 
          IF (.NOT.IMPLCT.OR.NEUTRAL) THEN  ! EXPLICIT PROBLEM
            HE=DMAX1(H/1.D4,10.D0*UROUND) 
C -------------------- 
            LEFT=.TRUE. 
C -------------------- 
C --- 
C ---       EULER STEP
            CALL FCN(N,X,Y,F2,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS,
     &               LPAST) 
            IF (NEUTRAL) THEN
                DO I=1,N-NDIMN 
               Z2(I)=Y(I)+HE*F2(I) 
              END DO
                    DO I=1,NDIMN
                 Z2(N-NDIMN+I)=F2(IPAST(NRDS+I))+HE*F2(N-NDIMN+I)
                END DO
            ELSE
                    DO I=1,N 
               Z2(I)=Y(I)+HE*F2(I) 
              END DO
                  END IF  
C -------------------- 
            LEFT=.FALSE. 
C -------------------- 
            CALL FCN(N,X,Y,F3,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS,
     &               LPAST)
            IF (NEUTRAL) THEN
                DO I=1,N-NDIMN 
               Z3(I)=Y(I)+HE*F3(I) 
              END DO
                    DO I=1,NDIMN
                 Z3(N-NDIMN+I)=F3(IPAST(NRDS+I))+HE*F3(N-NDIMN+I)
                END DO
            ELSE
                    DO I=1,N 
               Z3(I)=Y(I)+HE*F3(I) 
              END DO 
            END IF
            XL =ARGLAG(ILBP,X+HE,Z2,RPAR,IPAR,PHI,PAST,IPAST,NRDS,LPAST) 
            XLR=ARGLAG(ILBP,X+HE,Z3,RPAR,IPAR,PHI,PAST,IPAST,NRDS,LPAST)
            IF (XL.GE.BPP.AND.XLR.GE.BPP) THEN 
             LEFT=.FALSE. 
            ELSE IF (XL.LT.BPP.AND.XLR.LT.BPP) THEN 
             LEFT=.TRUE. 
            ELSE 
             IF (IOUT.EQ.1) THEN 
                    IF (XL.GT.BPP) THEN
                         WRITE (6,*) 
     &        ' WARNING!: SOLUTION DOES NOT EXIST AT X= ',X
              ELSE
                 WRITE (6,*) 
     &        ' WARNING!: SOLUTION IS  NOT UNIQUE AT X= ',X
                END IF
              GO TO 980
C             RETURN 
             END IF
            END IF
C ---       PROJECTION FOR DERIVATIVE COMPONENTS OF NEUTRAL EXPLICIT PROBLEMS
            PROJECT=.TRUE.
            IF (NEUTRAL.AND.PROJECT) THEN
               IF (LEFT) THEN
                    DO J=1,NDIMN
                 Y(N-NDIMN+J)=F2(IPAST(NRDS+J))
                END DO
            ELSE
                DO J=1,NDIMN
                 Y(N-NDIMN+J)=F3(IPAST(NRDS+J))
                END DO
               END IF
            END IF
          ELSE ! GENERAL IMPLICIT
              LEFT=.TRUE.
          END IF
C --- 
          BPD=.FALSE. 
         END IF  
C ----------------------------------------------------------------- 
C 
C --- 
         IF (ITOL.EQ.0) THEN 
             DO I=1,N 
                SCAL(I)=ATOL(1)+RTOL(1)*ABS(Y(I)) 
             END DO 
         ELSE 
             DO I=1,N 
                SCAL(I)=ATOL(I)+RTOL(I)*ABS(Y(I)) 
             END DO 
         END IF 
         IF (IOUT.NE.0) THEN 
C --- 
             NRSOL=NACCPT+1 
             XSOL=X 
             XOSOL=XOLD 
             NSOLU=N 
             HSOL=HOLD 
C --- 
             CALL SOLOUT(NRSOL,XOSOL,XSOL,HSOL,Y,CONT,LRC,NSOLU, 
     &                   RPAR,IPAR,IRTRN)
             IF (IRTRN.LT.0) GOTO 179 
         END IF 
         CALJAC=.FALSE. 
C ----------------------------- 
C ---    FIRST IS RESTORED 
         FIRST=.FALSE. 
         IF (LAST) FIRST=.TRUE. 
C ---    COMPUTATION OF Y0 
         CALL FCN(N,X,Y,Y0,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS,LPAST)
C           
         NFCN=NFCN+1 
         FIRST=.FALSE.  
C ----------------------------- 
         
C ------ FINAL POINT                               
         IF (LAST) THEN 
 45         CONTINUE
            IF (IGRID.EQ.NGRID) THEN 
               H=HOPT 
               IDID=1 
C ---          END OF COMPUTATION 
               GOTO 980 
            ELSE 
               IGRID=IGRID+1 
               LAST=.FALSE. 
C              
C              LEFT=.FALSE. 
               FIRST=.TRUE. 
               XEND=GRID(IGRID)
               IF (ABS(XEND-X).LE.(H*1.D-2)) THEN
                  IGRID=IGRID+1
                  GO TO 45 
               END IF 
               FLAGUS=.FALSE. 
               IF (WORK7.EQ.0.D0) THEN 
                  HMAXN=XEND-X  
                  HMAX=HMAXN 
               END IF 
               HNEW=MIN(HNEW,H)
            END IF 
         END IF 
 
         HNEW=MIN(HNEW,HMAXN) 
         HOPT=MIN(H,HNEW) 
         IF (REJECT) HNEW=MIN(HNEW,H)  
         REJECT=.FALSE. 
         IF ((X+HNEW/QUOT1-XEND).GE.0.D0) THEN 
            H=XEND-X 
            IF (H.LT.0.D0) THEN 
             WRITE(6,*) 'ERROR!: NEGATIVE STEPSIZE! AT ' 
             WRITE(6,*) 'X > XEND = ',X,XEND 
             STOP 
            END IF 
            LAST=.TRUE. 
         ELSE 
C ----------------------------------------------------- 
C --------  IN ORDER TO AVOID VERY SMALL END-STEPSIZES:    
            IF ((X+1.8D0*HNEW-XEND).GT.0.D0) THEN  
              H=(XEND-X)*0.55D0 
            ELSE 
              QT=HNEW/H  
              HHFAC=H 
C ------------------ 
              IMANT=0 
C ---         STEP IS MAINTAINED     
              IF (THETA.LE.THET.AND.QT.GE.QUOT1.AND.QT.LE.QUOT2) THEN 
               IF ((JLFLAG.EQ.0).AND.(.NOT.FIRST)) IMANT=1 
              ELSE  
               H=HNEW  
              END IF 
            END IF 
         END IF 
         HHFAC=H 
C ------------------------------------------- 
C ---    SIMPLE ITERATION FIRST                 
         JLFLAG=0 
C ---    
         IF ((.NOT.FIRST).AND.(THETA.LE.THET)) GOTO 20 
         GOTO 10 
C -------------------- 
C --- STEP IS ACCEPTED   
C -------------------- 
      ELSE 
C -------------------- 
C --- STEP IS REJECTED   
C -------------------- 
C --- 
         IF (BPD) THEN  
          IF (.NOT.FIRST) THEN 
           HNEW=HP 
           LAST=.TRUE. 
          END IF 
         ELSE 
          LAST=.FALSE. 
         END IF 
  
         FLAGN=.FALSE. 
C ------------------------------- 
         IF (IRTRN.LT.0) GOTO 179 
 
         REJECT=.TRUE. 
C --- 
         IF (FIRST) THEN 
             H=H*0.12D0 
             HHFAC=0.12D0 
         ELSE  
             HHFAC=HNEW/H 
             H=HNEW  
         END IF 
         IF (NACCPT.GE.1) NREJCT=NREJCT+1 
C ---> 
         JLFLAG=0  
         IF (CALJAC) GOTO 20 
         GOTO 10 
C -------------------- 
C --- STEP IS REJECTED   
C -------------------- 
      END IF 
C ---------------------------------- 
C --- END OF ERROR CONTROL PROCEDURE 
C ---------------------------------- 
C ----------------------------- 
C --- UNEXPECTED STEP-REJECTION 
  78  CONTINUE 
 
      REJECT=.TRUE. 
      FLAGN=.FALSE. 
C --- 
      LAST=.FALSE. 
C --- 
      IF (IER.NE.0) THEN 
          NSING=NSING+1 
          IF (NSING.GE.5) GOTO 176 
      END IF 
C --- 
      H=H*0.55D0  
      HHFAC=0.55D0 
C --- 
      IF (CALJAC) GOTO 20 
      IF (IRTRN.LT.0) GOTO 175 
      GOTO 10 
C --- FAIL EXIT 
 175  CONTINUE 
      IDID=-5 
      GOTO 980 
 176  CONTINUE 
      WRITE(6,979)X    
      WRITE(6,*) ' MATRIX IS REPEATEDLY SINGULAR, IER=',IER 
      IDID=-4 
      GOTO 980 
 177  CONTINUE 
      WRITE(6,979)X    
      WRITE(6,*) ' STEP SIZE TOO SMALL, H=',H 
      IDID=-3 
      GOTO 980 
 178  CONTINUE 
      WRITE(6,979)X    
      WRITE(6,*) ' MORE THAN NMAX =',NMAX,'STEPS ARE NEEDED'  
      IDID=-2 
      GOTO 980 
C --- EXIT CAUSED BY SOLOUT 
 179  CONTINUE 
      WRITE(6,979)X 
 979  FORMAT(' EXIT OF RADAR5 AT X=',E18.4)  
      IDID=2 
 
C --- RETURN LABEL 
 980  CONTINUE 
      WRITE(6,*) IBP-1,' COMPUTED BREAKING POINTS: '
      WRITE(8,*) 'BREAKING POINTS: '
      DO I=1,IBP
        WRITE(8,*) BPV(I)
      END DO
      WRITE(8,*) ' -------------- '
      CLOSE(8)

C --- DEALLOCATION OF THE MEMORY 
      DEALLOCATE (Z1,Z2,Z3,Y0,SCAL,F1,F2,F3) 
      DEALLOCATE(BPV) 
      DEALLOCATE(FJAC,ZL) 
      IF (IMPLCT) DEALLOCATE(FMAS) 
      DEALLOCATE (IP1,IP2,IPHES) 
      DEALLOCATE (E1,E2R,E2I) 
C      DEALLOCATE (PAST) 
      IF (NLAGS.GT.0) THEN 
       DEALLOCATE (FJACS,FJACLAG) 
       DEALLOCATE (IVL,IVE,IVC,ILS,ICOUN) 
       IF (ISWJL.NE.1) DEALLOCATE (IPJ,FJACL,XLAG) 
      END IF 
      DEALLOCATE (CONT,UCONT) 
      RETURN 
      END 
C 
C     END OF SUBROUTINE RADCOR 
C 
 
 
 
C *********************************************************** 
C 
      SUBROUTINE LAGR5(IL,X,Y,ARGLAG,PAST,THETA,IPOS,RPAR,IPAR,
     &                 PHI,IPAST,NRDS,LPAST)
C ---------------------------------------------------------- 
C     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION 
C     WITH THE OUTPUT-SUBROUTINE FOR RADAR5. IT PROVIDES THE 
C     POSITION OF THE DENSE OUTPUT AT THE IL-TH DELAY. 
C ---------------------------------------------------------- 
      USE IP_ARRAY 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
      REAL(kind=DP), dimension(1), intent(in) :: Y
      REAL(kind=DP), dimension(LPAST), intent(in) :: PAST 
      REAL(kind=DP), dimension(1), intent(in) :: RPAR 
      INTEGER, dimension(1), intent(in) :: IPAR 
      INTEGER, dimension(1), intent(in) :: IPAST 
      EXTERNAL PHI
      INTEGER :: LPAST

      LOGICAL FLAGS,FLAGN,FIRST,LAST,REJECT,BPD,LEFT 
C --- COMMON BLOCKS 
      COMMON /BPLOG/FIRST,LAST,REJECT,BPD 
      COMMON /BPCOM/BPP,ILBP,LEFT 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
C 
C --- COMPUTE DEVIATED ARGUMENT FOR IL-TH DELAY 
      XLAG=ARGLAG(IL,X,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS,LPAST)

C --- INITIAL PHASE 
      THETA=XLAG 
      IPOS=-1 

C     EPSILON GIVES THE OVERLAPPING
C     MIN VALUE FOR THE SUPER-POSITION NEIGHBOURHOOD OF A BP

      COMPAR=UROUND*MAX(ABS(XLAG),ABS(X0B))
      EPSACT=10.D0*COMPAR
        IF (IACT.GT.1) EPSACT=DMAX1(PAST(IACT-1)*1.D-2,EPSACT)

      IF (XLAG.LE.X0B) THEN 
C ---     DEVIATING ARGUMENT ON THE INITIAL SEGMENT         
        IF (.NOT.((IL.EQ.ILBP).AND.(BPD.OR.FIRST))) THEN 
            IF (XLAG-X0B.LT.0.D0) THEN
              RETURN
            ELSE
              IPOS=1
              THETA=-1.D0 
            END IF
          ELSE
         IF (ABS(XLAG-X0B).LE.EPSACT) THEN
            IF (LEFT) THEN
                IPOS=-1 
                THETA=XLAG
          ELSE
                IPOS=1
                THETA=(XLAG-(PAST(IPOS)+PAST(IPOS+IDIF-1)))/
     1             PAST(IPOS+IDIF-1) 
          END IF
         ELSE IF (ABS(XLAG-BPP).LE.EPSACT) THEN 
            IPOS=-1 
          IF (LEFT) THEN 
           IF (XLAG.GT.BPP) THEN 
            IF (BPP.GT.0.D0) THEN 
              THETA=BPP*(1.D0-100*UROUND) 
            ELSE 
              THETA=BPP*(1.D0+100*UROUND) 
            END IF 
           END IF 
          ELSE 
           IF (XLAG.LT.BPP) THEN 
            IF (BPP.GT.0) THEN 
              THETA=BPP*(1.D0+100*UROUND) 
            ELSE 
              THETA=BPP*(1.D0-100*UROUND) 
            END IF 
           END IF 
          END IF 
         END IF
          END IF    
          RETURN
        END IF 

C --- COMPUTE THE POSITION OF XLAG 
      IPA = IACT+IDIF 
      IF (IPA.GT.(MXST-1)*IDIF+1) IPA=1 
      IF (XLAG-PAST(IPA).LT.0.D0) THEN
         WRITE (6,*) ' MEMORY FULL, MXST = ',MXST 
         IRTRN=-1
           STOP 
C        RETURN 
      END IF 

      INEXT=IACT-IDIF 
      IF (INEXT.LT.1) INEXT=(MXST-1)*IDIF+1 
      XRIGHT=PAST(INEXT)+PAST(INEXT+IDIF-1) 

C -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
C --- INITIALIZE POSITION INSIDE THE MEMORY 
      IPOS=IPOSV(IL) 
C --- STEP AND DELAYS 
      IF (XLAG-XRIGHT.GT.0.D0) THEN 
        IF (.NOT.FLAGN) THEN 
           IPOS=IACT-IDIF 
           IF (IPOS.LT.1) IPOS=(MXST-1)*IDIF+1 
        ELSE 
           IPOS=IACT 
        END IF 
        FLAGS=.TRUE. 
C ------------------------ 
C ----- COMPUTE THETA (>0)     
C ------------------------ 
        H=PAST(IPOS+IDIF-1) 
        THETA=(XLAG-(PAST(IPOS)+H))/H 
C -----    EXTRAPOLATION USE OF THE COLLOCATION POLYNOMIAL            
C ----- 
C -----     
        RETURN                       
      ELSE 
   1    CONTINUE 
        IF (XLAG-PAST(IPOS).LE.0.D0) THEN 
           IPOS=IPOS-IDIF 
           IF (IPOS.LT.1) IPOS=(MXST-1)*IDIF+1 
           GOTO 1 
        END IF 
   2    CONTINUE 
        INEXT=IPOS+IDIF 
        IF (INEXT.GT.(MXST-1)*IDIF+1) INEXT=1 
        IF (XLAG.GT.PAST(INEXT).AND.INEXT.NE.IACT) THEN 
           IPOS=INEXT 
           GOTO 2 
        END IF
          IF (IPOS.EQ.1) THEN 
            IPREV=(MXST-1)*IDIF+1 
        ELSE 
            IPREV=IPOS-IDIF 
        END IF
C --- 
C ---   IN CORRESPONDENCE OF BREAKING POINTS                
C --- 
 3      CONTINUE 
        IF (.NOT.((IL.EQ.ILBP).AND.(BPD.OR.FIRST))) GOTO 10
 
        IF (BPP.EQ.X0B) THEN
            IF (LEFT) THEN
              IPOS=-1
              THETA=XLAG
              RETURN
          ELSE
              IF (IPOS.EQ.-1) IPOS=1
               GO TO 10
              END IF
        END IF

          IPOSB=0  
          IF (ABS(BPP-PAST(IPOS)).LE.10.D0*UROUND) THEN
           IPOSB=IPOS
          ELSE IF (ABS(BPP-PAST(INEXT)).LE.10.D0*UROUND) THEN
             IPOSB=INEXT
C          ELSE IF (ABS(BPP-PAST(IPREV)).LE.10.D0*UROUND) THEN
C             IPOSB=IPREV
          END IF

        IF (IPOSB.EQ.0) THEN
           GO TO 10
          END IF

        IF (IPOSB.EQ.1) THEN
            EPSILON=(PAST(IPOSB+IDIF)-PAST(IPOSB))
          ELSE IF (IPOSB.EQ.(MXST-1)*IDIF+1) THEN
            EPSILON=(PAST(IPOSB)-PAST(IPOSB-IDIF))
          ELSE
            EPSILON=DMIN1(PAST(IPOSB+IDIF)-PAST(IPOSB),
     1                  PAST(IPOSB)-PAST(IPOSB-IDIF))
          END IF
          EPSILON=DMAX1(EPSILON*1.D-2,EPSACT)
        
        IF (ABS(XLAG-BPP).GT.EPSILON) GOTO 10 

        IF (IPOSB.EQ.1) THEN
            IF (LEFT) THEN 
                IPOS=-1 
C                  IF (BPP.GT.0.D0) THEN 
C                  THETA=BPP*(1.D0-100*UROUND) 
C               ELSE 
C                  THETA=BPP*(1.D0+100*UROUND) 
C               END IF
C               SE PERO' IL DATO INIZIALE E' ESTENDIBILE
                THETA=XLAG
                RETURN
            ELSE 
                IPOS=1 
                GO TO 10
            END IF
        END IF
         
        IF (LEFT) THEN           
C ---     PREVIOUS INTERVAL HAS TO BE SELECTED
          IPOS=IPOSB-IDIF
        ELSE 
C ---     NEXT INTERVAL HAS TO BE SELECTED
          IPOS=IPOSB
        END IF 
      
C --------------------------------------------------------- 

C ----- COMPUTE THETA (<0): SITUAZIONE PIU' TIPICA     
  10    THETA=(XLAG-(PAST(IPOS)+PAST(IPOS+IDIF-1)))/PAST(IPOS+IDIF-1) 
C ----- REM: THETA IS NEGATIVE                                  
        END IF 
C --- UPDATE POSITION INSIDE THE MEMORY 
      IPOSV(IL)=IPOS 

      RETURN 
      END 
C 
C     END OF FUNCTION LAGR5 
C 
C -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

C *********************************************************** 
C 
      FUNCTION YLAGR5(IC,THETA,IPOS,PHI,RPAR,IPAR,PAST,IPAST,NRDS,LPAST) 
C ---------------------------------------------------------- 
C     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION 
C     WITH THE SUBROUTINE LAGR5. IT PROVIDES AN APPROXIMATION  
C     TO THE IC-TH COMPONENT OF THE SOLUTION AT XLAG. 
C ---------------------------------------------------------- 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
C --- REAL(kind=DP) PHI,YLAGR5
      REAL(kind=DP), dimension(1), intent(in) :: RPAR 
      INTEGER, dimension(1), intent(in) :: IPAR 
C --- 
      REAL(kind=DP), dimension(LPAST), intent(in) :: PAST 
      INTEGER, dimension(1), intent(in) :: IPAST 
      INTEGER :: LPAST
      LOGICAL FLAGS,FLAGN 
C---- COMMON BLOCKS 
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
 
C 
C --- INITIAL PHASE 
      IF (IPOS.EQ.-1) THEN 
            YLAGR5=PHI(IC,THETA,RPAR,IPAR,LPAST)   
C ---       CALL PHI(IC,THETA,YLAGR5,RPAR,IPAR)
            RETURN 
      END IF 
C --- 
C --- COMPUTE PLACE OF IC-TH COMPONENT  
      I=0  
      DO 5 J=1,NRDS  
      IF (IPAST(J).EQ.IC) I=J 
   5  CONTINUE 
      IF (I.EQ.0) THEN 
         WRITE (6,*) ' NO DENSE OUTPUT AVAILABLE FOR COMP.',IC 
         RETURN 
      END IF   
C ----- COMPUTE DESIRED APPROXIMATION 
      I=I+IPOS 
      YLAGR5=PAST(I)+THETA*(PAST(NRDS+I)+(THETA-C2M1)*(PAST(2*NRDS+I) 
     &            +(THETA-C1M1)*(PAST(3*NRDS+I))))
      RETURN 
      END 
C 
C     END OF FUNCTION YLAGR5 
C 
 
C *********************************************************** 
C 
      FUNCTION DLAGR5(IC,THETA,IPOS,PHI,RPAR,IPAR,PAST,IPAST,NRDS,LPAST) 
C ---------------------------------------------------------- 
C     THIS FUNCTION CAN BE USED FOR CONTINUOUS OUTPUT IN CONNECTION 
C     WITH THE SUBROUTINE LAGR5. IT PROVIDES AN APPROXIMATION  
C     TO THE IC-TH COMPONENT OF THE SOLUTION DERIVATIVE AT XLAG. 
C ---------------------------------------------------------- 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
C --- REAL(kind=DP) PHI,DLAGR5
      REAL(kind=DP), dimension(1), intent(in) :: RPAR 
      INTEGER, dimension(1), intent(in) :: IPAR 
CCC  
      REAL(kind=DP), dimension(LPAST), intent(in) :: PAST 
      INTEGER, dimension(1), intent(in) :: IPAST 
      INTEGER :: LPAST
      LOGICAL FLAGS,FLAGN 
C---- COMMON BLOCKS 
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
 
C 
C --- INITIAL PHASE 
      IF (IPOS.EQ.-1) THEN 
C           DLAGR5=DPHI(IC,THETA,RPAR,IPAR)   
            DLAGR5=0.D0 
            RETURN 
      END IF 
C --- 
C --- COMPUTE PLACE OF IC-TH COMPONENT  
      I=0  
      DO 5 J=1,NRDS  
      IF (IPAST(J).EQ.IC) I=J 
   5  CONTINUE 
      IF (I.EQ.0) THEN 
         WRITE (6,*) ' NO DENSE OUTPUT AVAILABLE FOR COMP.',IC 
         RETURN 
      END IF   
C ----- COMPUTE DESIRED APPROXIMATION 
      H=PAST(IPOS+IDIF-1) 
      I=I+IPOS 
      DLAGR5=PAST(NRDS+I)+(THETA-C2M1)*(PAST(2*NRDS+I) 
     &     +(THETA-C1M1)*PAST(3*NRDS+I)) 
     &     +THETA*(PAST(2*NRDS+I)+(2.D0*THETA-C2M1-C1M1)*PAST(3*NRDS+I))
      DLAGR5=DLAGR5/H 
      RETURN 
      END 
C 
C     END OF FUNCTION DLAGR5 
C 
 
C 
      SUBROUTINE BPDTCT(N,X,H,Y,ARGLAG,RPAR,IPAR,UCONT,GRID,NLAGS, 
     &                  FIRST,LAST,XEND,IGRID,BPV,IBP,ILBP,BPP,BPD,
     &                  KMAX,PHI,PAST,IPAST,NRDS,LPAST) 
                       
C ---------------------------------------------------------- 
C     THIS SUBROUTINE CAN BE USED FOR DETECTING BREAKING POINTS  
C     WITH THE OUTPUT-SUBROUTINE FOR RADAR5. IT PROVIDES AN 
C     APPROXIMATION TO THE IC-TH COMPONENT OF THE SOLUTION AT X. 
C ---------------------------------------------------------- 
        IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
        INTEGER, PARAMETER :: DP=kind(1D0) 
        REAL(kind=DP), dimension(N) :: Y 
        REAL(kind=DP), dimension(LPAST), intent(in) :: PAST 
        INTEGER, dimension(1), intent(in) :: IPAST 
        REAL(kind=DP), dimension(:), allocatable :: YADV 
        REAL(kind=DP), dimension(1) :: UCONT 
        REAL(kind=DP), dimension(1), intent(inout) :: GRID 
        REAL(kind=DP), dimension(1) :: BPV 
        INTEGER, dimension(1) :: IPAR 
        REAL(kind=DP), dimension(1) :: RPAR 
        EXTERNAL PHI
        INTEGER :: LPAST
        LOGICAL FIRST,LAST,BPD,FLAGS,FLAGN 
C----   COMMON BLOCKS 
        COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
        COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
C 
C ---     
        IF (FIRST) RETURN
        BPD=.FALSE.  
        LRC=4*N 
        EPSILON=1.D-10
        ALLOCATE(YADV(N)) 
        COMPAR=UROUND*MAX(ABS(X),ABS(X+H)) 
        
        XLAST=UCONT(LRC+1) 
        HLAST=UCONT(LRC+2) 
        DO IL=1,NLAGS 
         ALS = ARGLAG(IL,X,Y,RPAR,IPAR,PHI,PAST,IPAST,NRDS,LPAST) 
C -----  DEVIATING ARGUMENT AT X 
C -----  EXTRAPOLATION OF THE COLLOCATION POLYNOMIAL       
         DO IC=1,N 
          YADV(IC)=CONTR5(IC,N,X+H,UCONT,XLAST,HLAST) 
         END DO 
         ALD = ARGLAG(IL,X+H,YADV,RPAR,IPAR,PHI,PAST,IPAST,NRDS,LPAST) 
C -----  DEVIATING ARGUMENT AT X+H        
 
         IF (ABS(ALS-ALD).LE.COMPAR) GO TO 33
         DO L=1,IGRID-1 
          BPP=GRID(L) 
          IF ((ALS-BPP)*(ALD-BPP).LT.COMPAR) THEN 
           BPD=.TRUE. 
C          BREAKING POINT! 
           GO TO 33 
          END IF 
         END DO    
         DO L=IBP,1,-1 
          BPP=BPV(L) 
          IF ((ALS-BPP)*(ALD-BPP).LT.COMPAR) THEN 
C          BREAKING POINT! 
           BPD=.TRUE. 
           GO TO 33 
          END IF 
         END DO    
 33      CONTINUE 
         IF (BPD) THEN
C --- 
           THLIM=1.D0 
           THRIGH=THLIM 
C -------------------------------------------------- 
           THLEFT=0.D0
C ---      
           DO K=1,KMAX   
            THNEW = THLEFT - (ALS-BPP)*(THRIGH-THLEFT)/(ALD-ALS)  
C ---       TEST DI CONVERGENZA 
            IF (ABS(THRIGH-THNEW).LE.EPSILON.OR. 
     &          ABS(THLEFT-THNEW).LE.EPSILON) GOTO 36 
            XA=X+THNEW*H 
            DO IC=1,N 
             YADV(IC)=CONTR5(IC,N,XA,UCONT,XLAST,HLAST) 
            END DO 
            ALN = ARGLAG(IL,XA,YADV,RPAR,IPAR,PHI,PAST,IPAST,NRDS,LPAST) 
            IF ((ALS-BPP)*(ALN-BPP).LE.0.D0) THEN 
             ALD=ALN 
             THRIGH=THNEW 
            ELSE 
             ALS=ALN 
             THLEFT=THNEW 
            END IF 
C --- 
           END DO 
 36        CONTINUE 
C ---      BP FOUND!      
           IF ((THNEW.GT.COMPAR).AND.(THNEW.LT.THLIM)) THEN 
            HP=THNEW*H 
C ---
            IF (HP.LE.1.D2*COMPAR) THEN 
             BPD=.FALSE. 
             GO TO 37 
            END IF 
C ---
            IBP=IBP+1 
            BPV(IBP)=X+HP 
            H=HP 
            ILBP=IL 
            GOTO 37 
           ELSE 
C ---       BP ALREADY PRESENT
            BPD=.FALSE. 

           END IF  
         END IF 
        END DO 
 
 37     CONTINUE 
 
        DEALLOCATE(YADV) 

        RETURN 
        END 
 

      SUBROUTINE BPACC(N,X,H,Y,ARGLAG,RPAR,IPAR,Z1,Z2,Z3,FIRST,
     &                 BPV,IBP,ILBP,BPP,KMAX,PHI,PAST,IPAST,NRDS,LPAST)
                       
C ---------------------------------------------------------- 
C     THIS SUBROUTINE CAN BE USED FOR APPROXIMATING BREAKING POINTS  
C     IN TANDEM WITH THE SIMPLIFIED NEWTON ITERATION.. 
C ---------------------------------------------------------- 
        IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
        INTEGER, PARAMETER :: DP=kind(1D0) 
        REAL(kind=DP), dimension(N) :: Y,Z1,Z2,Z3 
        REAL(kind=DP), dimension(LPAST), intent(in) :: PAST
        INTEGER, dimension(1), intent(in) :: IPAST
        REAL(kind=DP), dimension(:), allocatable :: YCONT,YAPP
        REAL(kind=DP), dimension(1) :: BPV 
        INTEGER, dimension(1) :: IPAR 
        REAL(kind=DP), dimension(1) :: RPAR 
        EXTERNAL PHI
        INTEGER :: LPAST
        LOGICAL FIRST 
C----   COMMON BLOCKS 
        COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2 
        COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
C 
C ---     
        ALLOCATE(YCONT(4*N),YAPP(N))
        EPSILON=UROUND*1.D3
C ---   DYNAMIC UPDATE
        DO I=1,N 
            Z3I=Z3(I) 
            YI=Y(I)+Z3I   
            YCONT(I)=YI 
            Z2I=Z2(I) 
            Z1I=Z1(I) 
            A1=(Z2I-Z3I)/C2M1 
            YCONT(I+N)=A1 
            AK=(Z1I-Z2I)/C1MC2 
            ACONT3=Z1I/C1 
            ACONT3=(AK-ACONT3)/C2 
            A2=(AK-YCONT(I+N))/C1M1 
            YCONT(I+2*N)=A2 
            IF (.NOT.FIRST) THEN 
              YCONT(I+3*N)=A2-ACONT3 
            ELSE 
              YCONT(I+3*N)=0.D0 
            END IF 
        END DO 

C ---   INITIAL VALUES FOR THE COMPUTATION
        XSOL=X+H
        HSOL=H
        THLEFT=0.9D0
        THRIGH=1.0D0
        XL=X+THLEFT*H
        DO I=1,N 
            YAPP(I)=CONTR5(I,N,XL,YCONT,XSOL,HSOL) 
        END DO 
        ALS = ARGLAG(ILBP,XL,YAPP,RPAR,IPAR,PHI,PAST,IPAST,NRDS,LPAST)
C ---
        XR=X+THRIGH*H
        ALD = ARGLAG(ILBP,XR,YCONT,RPAR,IPAR,PHI,PAST,IPAST,NRDS,LPAST)
        DO K=1,KMAX 
            THNEW = THRIGH - (ALD-BPP)*(THRIGH-THLEFT)/(ALD-ALS)  
            THLEFT= THRIGH
C ---       TEST DI CONVERGENZA 
            IF (ABS(THNEW-THRIGH).LE.EPSILON) GOTO 36 
            IF ((THNEW.LE.0.5D0).OR.(THNEW.GE.1.5D0)) THEN 
              DEALLOCATE(YCONT,YAPP)
              RETURN
            END IF
            THRIGH=THNEW
            XAP=X+THRIGH*H 
            ALS=ALD
            DO I=1,N 
             YAPP(I)=CONTR5(I,N,XAP,YCONT,XSOL,HSOL) 
            END DO 
            ALD = ARGLAG(ILBP,XAP,YAPP,RPAR,IPAR,PHI,PAST,IPAST,NRDS,
     &                  LPAST) 
            IF (ABS(ALD-ALS).LE.EPSILON) GOTO 36 
        END DO 

 36     CONTINUE
C ---   BP FOUND 
        H=MIN(THRIGH,THLEFT)*H
        BPV(IBP)=X+H 

        DEALLOCATE(YCONT,YAPP)
 
        RETURN 
        END 

