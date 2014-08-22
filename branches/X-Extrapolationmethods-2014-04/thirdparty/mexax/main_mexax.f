      PROGRAM PLFDR
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
C
C* Begin Prologue PLFDR
C
C ---------------------------------------------------------------------
C
C* Purpose           Driver testing MEXX with
C                    example pendulum (PL)
C                    using default options i.e. full (F) matrix mode
C
C* File              plfdr.f
C* Version           1.1
C* Latest Change     96/07/19 (1.5)
C* Modification history
C    1.1             Correct workspace formula.
C                    Avoid free format because of compiler dependency.
C* Library           CodeLib
C* Code              Fortran 77
C                    Double Precision
C* Environment       Standard version for FORTRAN77 environments on
C                    PCs, workstations, and hosts
C
C* Subroutines supplied
C
C     PLFIN          Provide initial values
C     PLFF           Specify the problem
C     DUOUT          Dummy routine (instaed of solution output)
C     DUDOUT         Dummy routine (instaed of dense output)
C     DUSWI          Dummy routine (instead of switching functions)
C      
C  Problem dimensions for a maximum of 30 bodies.
C  The actual number of bodies will be requested from standard input.
C
      PARAMETER (
     $     NBODY  = 30,
     $     NPDIM  = 3*NBODY,
     $     NVDIM  = NPDIM,
     $     NUDIM  = 0,
     $     NGDIM  = 2*NBODY,
     $     NLDIM  = NGDIM,
     $     NPVUD  = NPDIM*2 + NUDIM,
     $     NZGDIM = 8*NBODY - 4
     $     )
C      INTEGER     NL1, NU1, NG1, NGL
C      PARAMETER ( NL1 = MAX (1,NLDIM) ,
C     $            NU1 = MAX (1,NUDIM) ,
C     $            NG1 = MAX (1,NGDIM) ,
C     $            NGL = MAX (1,NGDIM,NLDIM) )
C
C  Workspace dimension
C
      PARAMETER (
     $     LIWK   = NPDIM + 4*NVDIM + NLDIM + NUDIM + 60,
     $     LRWK   = (NVDIM+NLDIM+1)**2 + NPDIM*(NGDIM+NLDIM+1+18) +
     $               NVDIM*(NVDIM+45) + 28*(NLDIM+1) + 18*(NUDIM+1) + 
     $               NGDIM+1 + 50
     $     )
C      PARAMETER (
C     $     LIWK   = NPDIM + 4*NVDIM + NLDIM + NUDIM + 60,
C     $     LRWK   = (NVDIM+NL1)**2 + NPDIM*(NGL+18) +
C     $               NVDIM*(NVDIM+45) + 28*NL1 + 18*NU1 + NG1 + 50
C     $     )
C
      EXTERNAL PLFIN, PLFF, DUOUT, DUDOUT, DUSWI
C
      INTEGER MXJOB(150), IWK(LIWK)
      DOUBLE PRECISION RWK(LRWK)
      DIMENSION P(NPDIM), V(NVDIM), RTOL(NPVUD), ATOL(NPVUD)
      DIMENSION A(NVDIM), RLAM(NLDIM+1), U(NUDIM+1)
C
      COMMON / PEND / PMASS(NBODY), PLEN(NBODY), NPEND
C
      CHARACTER*30 PROTXT
      DATA PROTXT /'Pendulum'/
C
      CHARACTER*67 CMODUL
      DATA CMODUL /'@(#) plfdr.f 1.5 from 96/07/19 (mexx 1.2)
     $     '/
C
      DATA MXJOB /150*0/, RWK /LRWK*0.0D0/, IWK /LIWK*0/
C
C -- Machine constants
      LUSIN = 5
      LUSOUT = 6
C
C -- Start message
      WRITE (LUSOUT,9000) ' '
      WRITE (LUSOUT,9000)
     $     ' ===================================================='
      WRITE (LUSOUT,9000) CMODUL
      WRITE (LUSOUT,9000)
     $     ' ===================================================='
      WRITE (LUSOUT,9000) ' '
C
C -- Initialize problem
      NP = NPDIM
      NG = NGDIM
      NU = NUDIM
      NV=NP
      NL=0
      NBLK=0
      NMRC=0
      NZGMAX=0
      NZFMAX=0
      NSWMAX=0
C
      CALL PLFIN (NP, NV, NL, NG, NU,
     $     T0, P, V, U, A, RLAM, PROTXT)
C
C  Set options for integration
C
      T=T0
      TEND = 1.D0
      H = 1.D-3
      ITOL = 0
      TOL = 1.D-3
      RTOL(1) = TOL
      ATOL(1) = TOL
C
      IDID = 0
C
      CALL MEXX (NP, NV, NL, NG, NU,
     $     PLFF, T, TEND, P, V, U, A, RLAM,
     $     ITOL, RTOL, ATOL,
     $     H, MXJOB, IDID, LIWK, IWK, LRWK, RWK,
     $     DUOUT, DUDOUT, DUSWI)
C     
      IF (IDID .EQ. 0) THEN
         WRITE (LUSOUT,9000) ' o.k.'
      ELSE
         WRITE (LUSOUT,9100) ' Error return code from MEXX:', IDID
      ENDIF
C
C
      WRITE (LUSOUT,9000) ' '
      WRITE (LUSOUT,9200) ' Position values at the end point:', T
      WRITE (LUSOUT,9000) ' '
      WRITE (LUSOUT,9000) ' P(I),I=1,2,..'
      WRITE (LUSOUT,9000) ' '
      WRITE (LUSOUT,9210) (P(I),I=1,NP)
C
C
      WRITE (LUSOUT,9100) ' '
      WRITE (LUSOUT,9100) 'C==================================='
      WRITE (LUSOUT,9100) ' '
      WRITE (LUSOUT,9200) ' TOL       ', RTOL(1)
      WRITE (LUSOUT,9100) ' '
      WRITE (LUSOUT,9100) ' Statistics:  NSTEP, NACCPT, NREJCT'
      WRITE (LUSOUT,9100) '           ', MXJOB(51), MXJOB(52), MXJOB(53)
      WRITE (LUSOUT,9100) ' '
      WRITE (LUSOUT,9100) '             NFPROB,   NFCN,   NGCN'
      WRITE (LUSOUT,9100) '           ', MXJOB(54), MXJOB(55), MXJOB(56)
      WRITE (LUSOUT,9100) ' '
      WRITE (LUSOUT,9100) '              NMULT,   NDEC,   NSOL'
      WRITE (LUSOUT,9100) '           ', MXJOB(57), MXJOB(58), MXJOB(59)
      WRITE (LUSOUT,9100) ' '
      WRITE (LUSOUT,9100) 'C==================================='
C      
C
      WRITE (LUSOUT,9000) ' '
      WRITE (LUSOUT,9000) ' End of example ', PROTXT
      WRITE (LUSOUT,9000) ' Bye' 
C
      STOP
C
 9000 FORMAT (2A)
 9100 FORMAT (A, 3I8)
 9200 FORMAT (A, 1PD8.1)
 9210 FORMAT ((4(1PD19.10)))   
C
C  end of program PLFDR
C
      END
C
      SUBROUTINE PLFIN (NP, NV, NL, NG, NU,
     $     T0, P, V, U, A, RLAM, PROTXT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER NP, NV, NL, NG, NU
      DIMENSION P(*), V(*), U(*), A(*), RLAM(*)
      CHARACTER*(*) PROTXT
C
C **********************************************************************
C
C   Pendulum (variable number of bodies < 31)
C
C
C   Parameter:  masses                  pmass(i) 
C               lengths                 plen(i)  
C
C
C **********************************************************************
C
C* File              plfin.f
C* Version           1.1
C* Latest Change     96/07/19 (1.4)
C* Modification history
C    1.1             Correct initial values of RLAM.
C                    New free format starts text output with ' '.
C
C -------------------------------------------------------
C
      PARAMETER (NBMAX = 30)
      COMMON / PEND / PMASS(NBMAX), PLEN(NBMAX), NPEND
C
      SAVE
C
      CHARACTER*67 CMODUL
      DATA CMODUL /'@(#) plfin.f 1.4 from 96/07/19 (mexx 1.2)
     $     '/
C
      PROTXT = 'Pendulum (default options)'
C
C  Machine constants (logical units stdin and stdout)
      LUSIN = 5
      LUSOUT = 6
C
      WRITE (LUSOUT,*) 'Start example ', PROTXT
      WRITE (LUSOUT,*) ' '
C
C  Dimensions
C  Get problem parameters
C
      WRITE (LUSOUT,*) 'Enter number of bodies (default 8):'
C      READ (LUSIN,*, ERR=1000, END=1000) NPEND
C      GOTO 1010
 1000 NPEND = 8
 1010 WRITE (LUSOUT,*) 'Given:', NPEND
      IF (NPEND .LT. 1 .OR. NPEND .GT. NBMAX) THEN
         WRITE (LUSOUT,*)
     $        'Number of bodies is out of range (default value chosen)'
         NPEND = 8
      ENDIF
      WRITE (LUSOUT,*) ' '
C
      NP = 3*NPEND
      NV = NP
      NG = 2*NPEND
      NL = NG
      NU = 0
C
C  Model parameters
      T0 = 0.D0
      PMASSH = 1.0D0/DBLE(NPEND)
      PLENH = 0.5D0/DBLE(NPEND)
      DO 1020 I=1,NPEND
         PMASS(I) = PMASSH
         PLEN(I) = PLENH  
 1020 CONTINUE
C     
C  Initial values
      PIHALF = ASIN(1.D0)
      P(1) = PLEN(1)
      P(2) = 0.D0 
      P(3) = PIHALF  
C     
      V(1) = 0.D0
      V(2) = 0.D0
      V(3) = 0.D0
C     
      RLAM(1) = 0.D0
      RLAM(2) = 0.D0
C     
      DO 1030 I=2,NPEND
         L = (I-1)*3 + 1
         P(L) = P(L-3) + PLEN(I-1) + PLEN(I)
         P(L+1) = P(2) 
         P(L+2) = P(3)
         V(L) = 0.D0
         V(L+1) = 0.D0
         V(L+2) = 0.D0
         L = (I-1)*2 + 1
         RLAM(L) = 0.D0
         RLAM(L+1) = 0.D0
 1030 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PLFF (NP, NV, NL, NG, NU, LDG,
     $     T, P, V, U, RLAM,
     $     AM, GP, F, PDOT, UDOT, G, GI, FL,
     $     QFLAG, IFAIL)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  Input parameters:
      INTEGER NP, NV, NL, NG, NU, LDG
      LOGICAL QFLAG(9)
      DOUBLE PRECISION T, P(NP), V(NV), U(NU), RLAM(NL)
C  Output parameters:
      DOUBLE PRECISION AM(NV,NV), GP(LDG,*), F(NV), PDOT(NP), UDOT(NU),
     $     G(NG), GI(NG), FL(NV,NL)
C
C
C* File              plff.f
C* Version           1.1
C* Latest Change     96/05/03 (1.5)
C* Modification history
C    1.1             New free format starts text output with ' '.
C
C      
      PARAMETER (NBMAX = 30)
      COMMON / PEND / PMASS(NBMAX), PLEN(NBMAX), NPEND
      SAVE
C
      DATA GRAV/10./
C
C  Machine constant (logical unit stdout)
      LUSOUT = 6
C
      IFAIL = 0
C
      IF (QFLAG(1)) THEN
         DO 1010 J=1,NV
            DO 1000 I=1,NV
               AM(I,J) = 0.D0
 1000       CONTINUE
 1010    CONTINUE
         DO 1020 I=1,NPEND
            L = 3*(I-1) + 1
            AM(L  ,L  ) = PMASS(I)
            AM(L+1,L+1) = PMASS(I)
            AM(L+2,L+2) = PLEN(I)**2*PMASS(I)/3.D0
 1020    CONTINUE
      ENDIF
C
      IF (QFLAG(2) .OR. QFLAG(3)) THEN
         DO 1040 J=1,NV
            DO 1030 I=1,NG
               GP(I,J) = 0.D0
 1030       CONTINUE
 1040    CONTINUE               
         DO 1050 I=1,NPEND
            IROW = (I-1)*2 + 1
            ICOL = (I-1)*3 + 1
            GP(IROW, ICOL) = -1.D0
            GP(IROW, ICOL+2) = PLEN(I)*COS(P(ICOL+2))
            GP(IROW+1, ICOL+1) = -1.D0
            GP(IROW+1, ICOL+2) = PLEN(I)*SIN(P(ICOL+2))
            IF (I .GT. 1) THEN
               GP(IROW, ICOL-3) = 1.D0
               GP(IROW, ICOL-1) = PLEN(I-1)*COS(P(ICOL-1))
               GP(IROW+1, ICOL-2) = 1.D0
               GP(IROW+1, ICOL-1) = PLEN(I-1)*SIN(P(ICOL-1))
            ENDIF
 1050    CONTINUE
      ENDIF         
      IF (QFLAG(4)) THEN
         DO 1060 I=1,NPEND
            L = 3*(I-1) + 1
            F(L)= 0.D0
            F(L+1) = -GRAV*PMASS(I)  
            F(L+2) = 0.D0
 1060    CONTINUE
      ENDIF
      IF (QFLAG(5)) THEN
         DO 1070 I=1,NV
            PDOT(I) = V(I)
 1070    CONTINUE
      ENDIF
      IF (QFLAG(7)) THEN
         G(1) = -P(1) + PLEN(1)*SIN(P(3))
         G(2) = -P(2) - PLEN(1)*COS(P(3))
         DO 1080 I=2,NPEND
            L2 = 2*(I-1) + 1
            L3 = 3*(I-1) + 1
            G(L2) = P(L3-3) + PLEN(I-1)*SIN(P(L3-1))
     $           -P(L3) + PLEN(I)*SIN(P(L3+2))
            G(L2+1) = P(L3-2) - PLEN(I-1)*COS(P(L3-1))
     $           -P(L3+1) - PLEN(I)*COS(P(L3+2))
 1080    CONTINUE
      ENDIF
      IF (QFLAG(8)) THEN
         DO 1090 I=1,NG
            GI(I) = 0.D0
 1090    CONTINUE  
      ENDIF
      IF (QFLAG(9)) THEN
         DO 1110 J=1,NL
            DO 1100 I=1,NV
               FL(I,J) = 0.D0
 1100       CONTINUE
 1110    CONTINUE               
         WRITE (LUSOUT,*) 'PLFF (FL not available)'
      ENDIF      
C
      RETURN
      END
C
      SUBROUTINE DUOUT (NP, NV, NU, NL, T, P, V, U, A, RL, INFO, IRTRN)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER INFO(50), NP, NV, NU
      INTEGER IRTRN
      DOUBLE PRECISION T
      DIMENSION P(NP), V(NV), U(NU), A(NV), RL(NL)
C
C* Begin Prologue DUOUT
C
C* Purpose           Dummy routine instead of solution output for MEXX
C
C* File              duout.f
C* Version           1.0
C* Latest Change     96/03/05 (1.2)
C
C  
C* Input arguments:
C
C   NP        Int   Length of the position vector P
C
C   NV        Int   Length of the velocity vector V
C
C   NU        Int   Length of the vector of dynamic variables U
C
C   NA        Int   Length of the acceleration vector A
C
C   NL        Int   Length of the vector of Lagrange multipliers
C
C   T         Dble  Current time
C
C   P(NP)     Dble  The position vector. 
C
C   V(NV)     Dble  The velocity vector. 
C
C   U(NU)     Dble  The dynamic variable vector. 
C
C   A(NV      Dble  the dense output acceleration vector. 
C
C   RL(NL)    Dble  the dense output Lagrange multiplier vector. 
C  
C   INFO(50)  Int   Info array
C                   INFO(1) .. INFO(40) see MXJOB(1) .. MXJOB(40) 
C                                           in MEXX
C                   INFO(41): Type of call: 
C                             = 0: the first call
C                             = 1: an intermediate call
C                             = 2: the last call at TFIN
C                             = 3: at a root of switching function
C                                  and integration continues
C                             = 4: at a root of switching function
C                                  and integration will terminate
C                                  (this is the last call)
C
C                   INFO(42): number of current integration step
C
C                   INFO(43): interpolation indicator
C                             = 0: current t is an integration point
C                             = 1: the solution values are interpolated
C                                  values
C
C                   INFO(44): call counter
C   
C   Do NOT alter any of the input arguments
C
C* Output argument:
C
C   IRTRN     Int  Stop indicator
C                  If IRTRN .NE. 0 on return, MEXX will terminate
C 
C
      LUSOUT = 6
C
      WRITE (LUSOUT,*) ' +++ Warning in DUOUT'
      WRITE (LUSOUT,*)
     $     '     This dummy routine should never have been called.'
C
      IRTRN = 0
C
      RETURN
      END
C
      SUBROUTINE DUDOUT (NDP, NDV, NDU, NDA, NDL, T,
     $                   DP, DV, DU, DA, DL, INDP,
     $                   INDV, INDU, INDA, INDL, INFO, IFAIL)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER NDP, NDV, NDU, NDA, NDL
      DOUBLE PRECISION T, DP(*), DV(*), DU(*), DA(*), DL(*)
      INTEGER INDP(*), INDV(*), INDU(*), INDA(*), INDL(*), INFO(50),
     $     IFAIL 
C
C
C* File              dudout.f
C* Version           1.0
C* Latest Change     96/03/05 (1.2)
C      
C* Input arguments:
C
C   NDP       Int   Number of dense output components of position 
C                   vector P. Dimension of arrays DP and INDP. The
C                   current value is either the user prescibed value
C                   IWK(51) (if MDIND=INFO(33)=1 was set) or NP, the
C                   dimension of the position vector P (set by MEXX if
C                   INFO(33)=0 was set).
C
C   NDV       Int   as NDP but now for the velocity vector V.
C
C   NDU       Int   as NDP but now for the vector of dynamic 
C                   variables U.
C
C   NDA       Int   as NDP but now for the acceleration vector A.
C
C   NDL       Int   as NDP but now for the vector of Lagrange 
C                   multipliers.
C
C   T         Dble  Current time
C
C   DP        Dble  the dense output position vector. 
C
C   DV        Dble  the dense output velocity vector. 
C
C   DU        Dble  the dense output dynamic variable vector. 
C
C   DA        Dble  the dense output acceleration vector. 
C
C   DL        Dble  the dense output Lagrange multiplier vector. 
C
C   INDP(NDP) Int   indicator array of length NP. List of components of
C                   the position vecor p which was selected for dense
C                   output. INDP(I)=J means that component J of the
C                   vector P is available as component I of vector DP.
C                   (I.e. DP(I)=P(INDP(I))
C                   Note, that MEXX automatically provides INDP(I)=I
C                   if INFO(33)=0 was set.
C
C   INDV(NDV) Int   as INDP but now for the velocity vector V.
C
C   INDU(NDU) Int   as INDP but now the vector of dynamic variables U.
C
C   INDU(NDA) Int   as INDP but now for the acceleration vector A.
C
C   INDL(NDL) Int   as INDP but now the vector of Lagrange multipliers.
C  
C   INFO(50)  Int   Info array
C                   INFO(1) .. INFO(40) see MXJOB(1) .. MXJOB(40) 
C                                           in MEXX
C                   INFO(41): Type of call: 
C                             = 0: the first call
C                             = 1: an intermediate call
C                             = 2: the last call at TFIN
C                             = 3: at a root of switching function
C                                  and integration continues
C                             = 4: at a root of switching function
C                                  and integration will terminate
C                                  (this is the last call)
C
C                   INFO(42): number of current integration step
C
C                   INFO(43): interpolation indicator
C                             = 0: current t is an integration point
C                             = 1: the solution values are interpolated
C                                  values
C
C                   INFO(44): call counter
C
C   IFAIL     Int   Error indicator. IFAIL=0 on input
C   
C   Do NOT alter any of the input arguments.
C
C
C* Output arguments:
C
C   IFAIL     Int  Error indicator
C                  If IFAIL.LT.0  on return, MEXX will terminate
C 
C----------------------------------------------------------------------
C
C
      LUSOUT = 6
C
      WRITE (LUSOUT,*) ' +++ Warning in DUDOUT'
      WRITE (LUSOUT,*)
     $     '     This dummy routine should never have been called.'
C
      IFAIL = 0
C
      RETURN 
C
C end of subroutine DUDOUT
C
      END
C
      SUBROUTINE DUSWI (NP, NV, NL, NU, T, P, V, U, A, RLAM, NSWIF, G)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION P(NP), V(NV), U(*), A(NV), RLAM(NL), G(*)
C
C* File              duswi.f
C* Version           1.0
C* Latest Change     96/03/05 (1.2)
C      
C
      LUSOUT = 6
C
      WRITE (LUSOUT,*) ' +++ Warning in DUSWI'
      WRITE (LUSOUT,*)
     $     '     This dummy routine should never have been called.'
C
      IFAIL = 0
C
      RETURN
C
      END
