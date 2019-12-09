! ***********************************************************
!
      DOUBLE PRECISION FUNCTION DONTR5(I,N,X,CONT,XSOL,HSOL) 
! ----------------------------------------------------------
!     THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN
!     APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION DERIVATIVE AT .
!     X. IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR
!     THE LAST SUCCESSFULLY COMPUTED STEP (BY DELAY5).
! ----------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER, PARAMETER :: DP=kind(1D0)
      REAL(kind=DP), dimension(1), intent(in) ::  CONT
! --- REQUIRED CONSTANTS
      COMMON /CONSTN/C1,C2,C1M1,C2M1,C1MC2

      N2=2*N
      N3=3*N

      S=(X-XSOL)/HSOL

      DONTR5=(CONT(I+N)+ &
         (S-C2M1)*(CONT(I+N2)+CONT(I+N3)*(S-C1M1))+ &
          S*(CONT(I+N2)+CONT(I+N3)*(2.*S-C1M1-C2M1)))/HSOL
     
      RETURN
      END
!
!     END OF FUNCTION DONTR5
!
! ***********************************************************

