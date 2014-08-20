      SUBROUTINE ZIBSEC (CPTIM, IFAIL)
C
C*********************************************************************
C  Set CPTIM to cpu time in seconds.
C  This routine is machine dependent.
C*********************************************************************
C
C  Output parameters
C    CPTIM    REAL        Cpu time in seconds
C    IFAIL    INTEGER     Errorcode
C 
C  CSUN code works on sun and MacOSX (g77).
C  CXLF code works on MacOSX (xlf)
C  CVIS code works with Visual Fortran (Windows)
C
C*********************************************************************
C
      REAL CPTIM
      INTEGER IFAIL
C
CSUN      REAL RTIME(2)
CXLF      REAL(4) etime_
CXLF      TYPE TB_TYPE
CXLF      SEQUENCE
CXLF      REAL(4) USRTIME
CXLF      REAL(4) SYSTIME
CXLF      END TYPE
CXLF      TYPE (TB_TYPE) ETIME_STRUCT
C
      CPTIM = 0
      IFAIL = 0
C
CSUN      CPTIM = ETIME(RTIME)
CSUN      CPTIM = RTIME(1)
C
CXLF      CPTIM = etime_(ETIME_STRUCT)
C
CVIS      CALL CPU_TIME(CPTIM)
C
      RETURN
      END
