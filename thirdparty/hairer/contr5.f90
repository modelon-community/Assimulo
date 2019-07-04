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
      REAL(kind=DP), dimension(1), intent(in) ::  CONT
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

