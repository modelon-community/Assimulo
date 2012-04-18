
      SUBROUTINE ASSIMULO_RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H, 
     &                  RTOL,ATOL,ITOL, 
     &                  JAC,IJAC,MLJAC,MUJAC, 
     &                  JACLAG,NLAGS,NJACL, 
     &                  IMAS,SOLOUT,IOUT, 
     &                  WORK,IWORK,RPAR,IPAR,IDID, 
     &                  GRID,IPAST,MAS,MLMAS,MUMAS,
     &                  LIPAST, LGRID,
     &                  PAST,LRPAST)
      USE IP_ARRAY 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      INTEGER, PARAMETER :: DP=kind(1D0) 
      REAL(kind=DP), dimension(N), intent(inout) ::  
     &                             Y 
      REAL(kind=DP), dimension(30), intent(inout) ::  
     &                             WORK 
      REAL(kind=DP), dimension(1), intent(inout) ::  
     &                             ATOL,RTOL 
      INTEGER, dimension(30), intent(inout) :: IWORK 
      REAL(kind=DP), dimension(LGRID), intent(inout) :: GRID 
      INTEGER, dimension(LIPAST), intent(inout) :: IPAST 
      REAL(kind=DP), dimension(1), intent(in) :: RPAR 
      INTEGER, dimension(1), intent(in) :: IPAR 
      
      REAL(kind=DP), dimension(LRPAST), intent(inout) :: PAST

      INTEGER, intent(in) :: LIPAST
      INTEGER, intent(in) :: LGRID
      INTEGER, intent(in) :: LRPAST
      
      LOGICAL FLAGS, FLAGN
      EXTERNAL FCN,PHI,ARGLAG,JAC,JACLAG,MAS,SOLOUT 
C ----> COMMON BLOCKS <---- 
      COMMON /POSITS/X0B,UROUND,HMAX,IACT,IRTRN,IDIF,MXST,FLAGS,FLAGN 
      
      
      CALL  RADAR5(N,FCN,PHI,ARGLAG,X,Y,XEND,H, 
     &                  RTOL,ATOL,ITOL, 
     &                  JAC,IJAC,MLJAC,MUJAC, 
     &                  JACLAG,NLAGS,NJACL, 
     &                  IMAS,SOLOUT,IOUT, 
     &                  WORK,IWORK,RPAR,IPAR,IDID, 
     &                  GRID,IPAST,MAS,MLMAS,MUMAS,
     &                  PAST)
     
      END
