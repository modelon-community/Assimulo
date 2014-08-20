      SUBROUTINE MEXX (NP, NV, NL, NG, NU,
     $     FPROB, T, TFIN, P, V, U, A, RLAM,
     $     ITOL, RTOL, ATOL,
     $     H, MXJOB, IERR, LIWK, IWK, LRWK, RWK,
     $     SOLOUT, DENOUT, FSWIT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NP, NV, NL, NG, NU
      EXTERNAL FPROB
      DOUBLE PRECISION T, TFIN, P(NP), V(NV), U(*), A(NV), RLAM(NL)
      INTEGER ITOL
      DOUBLE PRECISION RTOL(*), ATOL(*), H
      INTEGER MXJOB(150), IERR, LIWK, IWK(LIWK), LRWK
      DOUBLE PRECISION RWK(LRWK)
      EXTERNAL SOLOUT, DENOUT, FSWIT
C
C* Begin Prologue MEXX
C
C ---------------------------------------------------------------------
C
C* Title
C
C    Numerical integration of the equations of motion of a
C    constrained mechanical system, including dry friction, 
C    nonholonomic constraints, generalized velocities and 
C    external dynamics.
C
C* Written by        Ch. Engstler, Ch. Lubich, U. Nowak, U. Poehle
C* Purpose           Solution of the equations of motion of a
C                    constrained mechanical system.
C* Method            Extrapolated half-explicit Euler method due to /1/.
C* Category          i1a2c. - Stiff and mixed algebraic-differential
C                             equations, Special Applications
C* Keywords          extrapolation, ODE, constrained mechanical system
C* File              mexx.f
C* Version           pre-release 1.2
C* Latest Change     96/07/18 (1.10)
C* Modification history
C    1.2             Correct mxroot.f: If more than one root is found
C                    report the solution values for the first root
C                    instead of the root detected last.
C    1.1             Correct and enhance documentation.
C                    Adapt to ZIB GUI.
C                    Correct workspace formula.
C* Library           CodeLib
C* Code              Fortran 77
C                    Double Precision
C* Environment       Standard version for FORTRAN77 environments on
C                    PCs, workstations, and hosts
C* Copyright     (c) Konrad-Zuse-Zentrum fuer Informationstechnik
C                    Berlin (ZIB)
C                    Takustrasse 7, D-14195 Berlin-Dahlem
C                    phone : + 49/30/84185-0
C                    fax   : + 49/30/84185-125
C* Contact           Universitaet Tuebingen
C                    Christian Lubich
C                    e-mail: lubich@na.mathematik.uni-tuebingen.de
C
C                 or ZIB, Numerical Analysis and Modelling
C* Contact           Rainald Ehrig
C                    ZIB, Numerical Analysis and Modelling
C                    phone : + 49/30/84185-282
C                    fax   : + 49/30/84185-107
C                    e-mail: ehrig@zib.de
C
C  ---------------------------------------------------------------------
C
C* Licence
C  -------
C  You may use or modify this code for your own non-commercial purposes
C  for an unlimited time.  In any case you should not deliver this code
C  without a special permission of ZIB.  In case you intend to use the
C  code commercially, we oblige you to sign an according licence
C  agreement with ZIB.
C
C
C* Warranty
C  --------
C  
C  This code has been tested up to a certain level. Defects and
C  weaknesses, which may be included in the code, do not establish any
C  warranties by ZIB. ZIB does not take over any liabilities which may
C  follow from acquisition or application of this code.
C
C
C* Software status
C  ---------------
C  
C  This code is under care of ZIB and belongs to ZIB software class I.
C
C
C  ---------------------------------------------------------------------
C
C* References
C  ----------
C
C /1/ Ch. Lubich:
C     Extrapolation integrators for constrained multibody systems.
C     Report, Univ. Innsbruck, 1990.
C
C /2/ Ch. Lubich, U. Nowak, U. Poehle, Ch. Engstler:
C     MEXX - Numerical Software for the Integration of Constrained
C     Mechanical Systems
C     Konrad-Zuse-Zentrum fuer Informationstechnik Berlin,
C     Technical Report SC 92-12 (December 1992)
C
C /3/ E. Hairer and A. Ostermann:
C     Dense output for extrapolation methods.
C     Report, Universite de Geneve, 1989
C
C /4/ E. Hairer, S. P. Norsett, G. Wanner:
C     Solving Ordinary Differential Equations I.
C     Springer-Verlag, 1987.
C
C  ---------------------------------------------------------------------
C
C* Summary
C  -------
C
C  Numerical integration of the equations of motion of a
C  constrained mechanical system, including dry friction
C  and external dynamics:
C
C  a)           p' = T(t,p) * v
C  b)     M(t,p)v' = f(t,p,v,lambda,u) - G(t,p)^T * lambda
C  c)           0  = G(t,p) * v + gI(t,p)
C  d)           u' = d(t,p,v,lambda,u)
C
C  with position constraints treated as invariants:
C
C               0  = g(t,p)
C
C
C  This is an extrapolation algorithm, following reference /1/.
C
C  Step size control and order selection are adapted from:
C  E.Hairer, S.P.Norsett, G.Wanner (ref. /4/).
C
C  Dense output is adapted from the code ELEX described by
C  E.Hairer and A. Ostermann in ref. /3/.
C
C
C* Software used:
C  
C  The numerical solution of the arising linear equations is done by
C  means of LINPACK or Harwell subroutines:
C  
C    DGEFA and DGESL (Gauss-algorithm with column-pivoting and
C    row-interchange) in the case MXJOB(2)=0 or MXJOB(2)=3,
C  
C    DSIFA and DSISL (Gauss-algorithm for symmetric matrices with
C    symmetric pivoting) in the case MXJOB(2)=1,
C  
C    DCHDC (Cholesky decomposition of a positive definite matrix), DTRSL
C    (solving triangular systems), and DQRDC (using householder
C    transformations to compute the QR-factorization) in the case
C    MXJOB(2)=2,
C  
C    MA28AD (LU-factorization of a sparse matrix), MA28BD (factorize a
C    matrix of a known sparsity pattern), and MA28CD (solve a system of
C    equations without iterative refinement) in the case MXJOB(2)=10.
C  
C    Machine dependent constants are provided by the subroutines
C    ZIBCONST.
C
C  ---------------------------------------------------------------------
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C     FPROB       Name (external) of subroutine specifying the problem
C
C              1. calling sequence for the standard case:
C                 ---------------------------------------
C                 SUBROUTINE FPROB (NP, NV, NL, NG, NU, LDG,
C                                   T, P, V, U, RLAM,
C                                   AM, GP, F, PDOT, UDOT, G, GI, FL,
C                                   QFLAG, IFAIL) 
C
C               - Input parameters:
C
C                 LOGICAL QFLAG(9)
C                         see below
C
C                 INTEGER NP, NV, NL, NG, NU, LDG
C                         NP, NV, NL, NG, NU: see above
C                         LDG: Leading dimension of GP
C                              LDG = NL  if QFLAG(2)
C                                  = NG  if QFLAG(3)
C
C                 DOUBLE PRECISION T, P(NP), V(NV), U(NU), RLAM(NL)
C                         The actual values of time, position, velocity,
C                         dynamic variables and Lagrange multipliers
C
C               - Output parameters:
C
C                  DOUBLE PRECISION AM(NV,NV), GP(LDG,*), F(NV),
C                 $       PDOT(NP), UDOT(NU), G(NG), GI(NL), FL(NV,NL)
C                  INTEGER IFAIL - error flag
C                         for IFAIL .NE. 0 MEXX terminates
C
C               - Remark on dimensions: Even if one of the above 
C                    dimensions (NG or NU) is zero, the calling program
C                    provides a minimum of 1 for each argument. Such, if
C                    e.g. NU=0, the statement 'DIMENSION U(*)' is safe 
C                    even if NU=0.
C
C                 According to the actual values of QFLAG, FPROB must
C                 return:
C
C                 if (QFLAG(1)) the mass matrix  M(t,p) in AM
C                 
C                 if (QFLAG(2)) the velocity constraint matrix G(t,p)
C                               in GP
C                 
C                 if (QFLAG(3)) the position constraint Jacobian
C                               C(t,p) := dg/dp(t,p) * T(t,p) in GP
C                 
C                 if (QFLAG(4)) the forces  f(t,p,v,lambda,u)  IN  F
C                  
C                 if (QFLAG(5)) the derivative  p'=T(t,p)*v in PDOT
C                  
C                 if (QFLAG(6)) the derivative  u'=d(t,p,v,lambda,u)
C                               in UDOT
C                
C                 if (QFLAG(7)) the constraint residual g(t,p) in G
C                 
C                 if (LFLAG(8)) the second term 'gI' of the velocity 
C                               constraint equation 'G(t,p)v+gI(t,p)=0'
C                               in GI
C
C                 if (LFLAG(9)) the Jacobian df/dlambda in FL (only 
C                               evaluated if the optional parameter
C                               MXJOB(1)=MDISC is set to 1)
C
C              2. Calling sequence if MMODE=1 and JOB<10 (see below):
C                 ---------------------------------------------------
C                 SUBROUTINE FPROB (NP, NV, NL, NG, NU, LDG,
C                $                  T, P, V, U, RLAM,
C                $                  AM, GP, F, PDOT, UDOT, G, GI, FL,
C                $                  QFLAG, NBLK, NMRC, IFAIL)
C
C               - Additional input parameters:
C
C                    NBLK - number of diagonal blocks of AM
C                           (prescribed by user in MXJOB(4))
C                    NMRC - size of diagonal blocks of AM
C                           (prescribed by user in MXJOB(5))
C
C                    Storage of AM: AM(NV,NMRC)   (NV = NBLK*NMRC)
C
C              3. Calling sequence if MMODE=1 and JOB=10 (see below):
C                 ---------------------------------------------------
C                 SUBROUTINE FPROB (NP, NV, NL, NG, NU, LDG,
C                $                  T, P, V, U, RLAM,
C                $                  AM, GP, F, PDOT, UDOT, G, GI, FL,
C                $                  QFLAG, NBLK, NMRC,
C                $                  NZGMAX, NZFMAX,
C                $                  IROWG, JCOLG, NZGACT, IPCG,
C                $                  IROWF, JCOLF, NZFACT, IPCF, IFAIL)
C
C               - Additional input parameter:
C                    NZGMAX - expected (maximum) number of non-zero
C                         elements of GP (prescribed by user in
C                         MXJOB(7))
C                    NZFMAX - expected (maximum) number of non-zero
C                         elements of FL (prescribed by user in
C                         MXJOB(8))
C
C                 Additional input/output parameters:
C
C                    IROWG(NZGACT) - row indices of elements of GP
C                    JCOLG(NZGACT) - associated column indices
C                    NZGACT - actual number of non-zero elements
C                         (NZGACT .LE. NZGMAX)
C                    IPCG - pattern change indicator for the sparse
C                         matrix GP
C
C                         On input:
C                         =0 Pattern is known from previous calls
C                            to FPROB
C                         =1 Pattern is unknown i.e. IROWG, JCOLG, and
C                            NZGACT are undefined and have to be set by
C                            FPROB (the very first time FPROB is called
C                            with IPCG=1).
C
C                         On output:
C                         =0 Pattern was not changed by FPROB during the
C                            last call
C                         =1 Pattern is new i.e. IROWG, JCOLG, and
C                            NZGACT have been changed (as IPCG=1
C                            requires an expensive refactorization of
C                            the matrix, IPCG should be set to 0 if
C                            pattern doesn't change).
C
C                    IROWF(NZFACT) - row indices of elements of GP
C                    JCOLF(NZFACT) - associated column indices
C                    NZFACT - actual number of non-zero elements
C                         (NZFACT .LE. NZFMAX)
C                    IPCF - pattern change indicator for the sparse
C                         matrix FL just as IPCG for matrix GP
C
C                 Storage of AM: AM(NV,NMRC)   (NV = NBLK*NMRC)
C                 Storage of GP and FL: GP(NZGACT) and FL(NZFACT)
C
C                 The calling program provides GP(NZGMAX),
C                 IROWG(NZGMAX), JCOLG(NZGMAX), FL(NZFMAX),
C                 IROWF(NZGMAX), and JCOLF(NZFMAX)
C
C                 The following holds for FL similarly:
C                 For the sparse matrix GP, the non-zero elements of the
C                 matrix are stored in a REAL array and two INTEGER
C                 arrays all of length NZGMAX (=MXJOB(7)). GP(i) holds
C                 the value and IROWG(i) and JCOLG(i) hold the
C                 corresponding row index and column index respectively.
C
C
C
C     SOLOUT      Name (external) of subroutine providing the numerical
C                 solution during integration. If MOUT=1 (see MXJOB(30))
C                 it is called after every successful step. (Supply
C                 at least a dummy subroutine if MOUT=0)
C
C                   SUBROUTINE SOLOUT (NP, NV, NU, NL, T, P, V, U,
C                                      A, RLAM, INFO, IRTRN)
C                   INTEGER NP, NV, NU, NL
C                   DOUBLE PRECISION T, P(NP), V(NV), U(NU), A(NV),
C                                    RLAM(NL)
C                   INTEGER INFO(50), IRTRN
C                 furnishes the solution (P,V,U,A,RLAM) at the NR-th
C                 grid-point T (the initial value is considered
C                 as the first grid-point).
C                 IRTRN  serves to interrupt the integration. If IRTRN
C                 is set negative, MEXX returns to the calling
C                 program.
C                 See the sample routine SOLOUT for details.
C
C                 
C     DENOUT      Name (external) of subroutine providing dense output
C                 of the solution by interpolation between integration
C                 points. If MDOUT>0 (see MXJOB(31)) it is called after
C                 every successful step. (Supply at least a dummy
C                 subroutine if MDOUT=0)
C
C                   SUBROUTINE DENOUT (NDP, NDV, NDU, NDA, NDL,
C                                      T, DP, DV, DU, DA, DL,
C                                      INDP, INDV, INDU, INDA, INDL
C                                      INFO, IRTRN)
C                   INTEGER NDP, NDV, NDU, NDA, NDL
C                   DOUBLE PRECISION X, DP, DV, DU, DA, DL,
C                   INTEGER INDP, INDV, INDU, INDA, INDL, INFO, IRTRN
C
C                 furnishes user selected components of the solution at
C                 the time T by interpolation from grid point TOLD.
C                 IRTRN  serves to interrupt the integration. If IRTRN
C                 is set negative, MEXX returns to the calling
C                 program.
C                 See the sample routine DENOUT for details.
C
C
C     FSWIT       Name (external) of subroutine specifying a maximum
C                 number of NSWMAX (see MXJOB(36)) functions g(y,t) for
C                 the root finder
C
C                 Calling sequence:
C                 
C                 SUBROUTINE FSWIT (NP, NV, NL, NU, T, P, V, U, A, RLAM,
C                                   NSWIF, G)
C                 Input parameters: NP, NV, NL, NU
C                   INTEGER NP, NV, NL, NU
C                   DOUBLE PRECISION T, P(NP), V(NV), U(NU), RLAM(NL)
C                     Actual time value and the actual values of time,
C                     position, velocity, dynamic variable and Lagrange
C                     multiplier 
C
C                 Output parameters:
C                   INTEGER NSWIF - Actual number of user defined
C                     equations to be checked by root finder 
C                   DOUBLE PRECISION G(NSWIF) - Actual values of the
C                     NSWIF switching functions 
C                 See the sample driver routine for details.
C                    
C
C
C
C* Input parameters (* marks input/output parameters)
C  ----------------
C  
C  These first parameters NP, NV, ... , ATOL set up the problem and the
C  basic solution requirements for the constrained mechanical system
C  to solve.
C  
C  Several options may be set via the vector MXJOB, e.g. for different
C  kinematic formulations or the extent of optional printout.
C
C  Additional input parameters for special purposes can be given in
C  IWK and RWK.
C  
C  User defined subroutines SOLOUT and DENOUT allow output of solution
C  values during the integration. SOLOUT and/or DENOUT will be called
C  only if the respective options are set.
C  
C  The last parameter FSWIT is used if the user wants to activate the
C  root finder during the integration. FSWIT will be called only if the
C  root finding option is set. 
C
C  Due to FORTRAN restrictions, at least some dummy routines have to be
C  present for entries SOLOUT, DENOUT, and FSWIT even if solution output
C  and/or root finding is not desired. 
C  
C
C  NP         Int     Dimension of the position vector  P  (NP .GT. 0)
C
C  NV         Int     Dimension of the velocity vector V
C                     (0 .LT. NV .LE. NP)
C
C  NL         Int     Number of velocity constraints
C                     (0 .LE. NL .LE. NV)
C                     Dimension of Lagrange multipliers lambda (RLAM)
C
C  NG         Int     Number of position constraints  g
C                     (0 .LE. NG .LE. NL)
C                     (if NG .EQ. 0 : MXJOB(2) = job = 0 or 
C                                                job = 1 must be used)
C
C  NU         Int     Dimension of dynamic variable  U  (NU .GE. 0)
C
C  T          Dble *  Initial time value
C  TFIN       Dble    Final time value
C
C  P(NP)      Dble *  Initial position vector
C  V(NV)      Dble *  Initial velocity vector
C  U(NU)      Dble *  Initial dynamic variable vector
C  A(NV)      Dble *  Initial acceleration vector
C  RLAM(NL)   Dble *  Initial Lagrange multipliers vector
C
C  ITOL       Int     Switch for RTOL and ATOL:
C                     ITOL=0: both RTOL and ATOL are scalars.
C                       The code keeps, roughly, the local error
C                       of P(I) below RTOL*ABS(P(I)) + ATOL, that
C                       of V(I) below RTOL*ABS(V(I)) + ATOL, and that
C                       of U(I) below RTOL*ABS(U(I)) + ATOL.
C                     ITOL=1: both RTOL and ATOL are vectors.
C                       The code keeps, roughly, the local error
C                       of P(I) below RTOL(I)*ABS(P(I)) + ATOL(I), that
C                       of V(I) below RTOL(NP+I)*ABS(V(I)) + ATOL(NP+I),
C                       and that of U(I) below
C                       RTOL(NP+NV+I)*ABS(U(I)) + ATOL(NP+NV+I) .
C
C  RTOL       Dble    Relative and absolute error tolerances. They can
C  ATOL       Dble    be both scalars or else both vectors of length
C                     NP+NV+NU .
C
C  H          Dble *  Initial step size guess
C
C  MXJOB(150) Int  *  Elements no. 1 to 50 are used as input options for
C                     MEXX performance.
C                     Zero initiation by user generates internal default
C                     setting. Description see below.
C
C
C
C* Output parameters (* marks input/output parameters)
C  -----------------
C
C  Additional output parameters for special purposes may be found in
C  IWK and RWK.
C
C  T          Dble *  Time value where the solution is computed
C                     (after successful return T=TFIN)
C
C  P(NP)      Dble *  Position vector at time T
C  V(NV)      Dble *  Velocity vector at time T
C  U(NU)      Dble *  Dynamic variables at time T
C  A(NV)      Dble *  Acceleration vector at time T
C  RLAM(NL)   Dble *  Lagrange multipliers at time T
C
C  H          Dble *  Actual step size guess
C
C  MXJOB(150) Int  *  Elements no. 51 to 150 are used as optional output
C                     from MEXX. Description see below.
C
C  IERR       Int     Reports on success upon return:
C                     IERR = 0  computation successful,
C                     IERR < 0  computation failed,
C                     IERR > 0  computation successful but stopped
C                               before TFIN.
C                      =  -2: Invalid input: TFIN .LE. T
C                      =  -3: INTEGER work space exhausted
C                      =  -4: REAL work space exhausted
C                      =  -9: Step size H .LE. machine precision
C                      = -10: Limit of integration steps reached
C                      = -11: No convergence in projection step
C                      = -15: Decomposition failed (from MEXXEX)
C                      = -16: Decomposition failed (from MXPRJC)
C                             For IERR -15 or -16:
C                             The error code of the decomposition
C                             routine can be found at MXJOB(73) and
C                             MXJOB(74) tells which decomposition
C                             routine failed.
C                      = -17: user interrupt in FPROB (from MEXXEX)
C                      = -18: user interrupt in FPROB (from MXPRJC)
C                      = -19: got negative return code from DENOUT
C                      = -20: got negative return code from MXRTC/F
C                             For IERR -17, -18, -19, or -20:
C                             The error code of the user routine
C                             can be found at MXJOB(73).
C                      = -21: invalid input: NP .LE. 0
C                      = -22: invalid input: NG .LT. 0
C                      = -23: invalid input: NU .LT. 0
C                      = -25: invalid input: NG .GT. NP
C
C                     IERR < -100 and  IERR > -151:
C                     Invalid input for position -100-IERR of 
C                     array MXJOB
C                      =-101: JOB is not 0, 1, 2, 3, or 10
C                      =-102: MMODE is not 0 or 1
C                      =-112, -114, -120: Illegal logical unit number
C                      =-122: IALAM .LT. 0 .OR. IALAM .GT. 1
C                      =-130: MOUT .LT. 0 .OR. MOUT .GT. 1
C                      =-131: MDOUT .LT. 0 .OR. MDOUT .GT. 4
C                      =-132: NDOUT .LT. 0
C                      =-133: MDIND .LT. 0 .OR. MDIND .GT. 1
C
C                     IERR < -150:
C                      =-151: Invalid input: JOB .EQ. 10 .AND.
C                             MMODE .NE. 1
C                      =-152: Invalid input: JOB .EQ. 3 .AND. NG .LE. 0
C                      =-153: Invalid input: NBLK .NE. NP/NMRC
C                      =-154: Invalid input: NZGMAX .LE. 0 .AND.
C                             JOB .EQ. 10
C                      =-155: Invalid input: JOB .EQ. 0 .AND. 
C                             MMODE .NE. 0
C
C                     IERR < -350: Invalid input for dense output
C                      =-351: NDP .LT. 0 .OR. NDP .GT. NP
C                      =-352: NDV .LT. 0 .OR. NDV .GT. NV
C                      =-353: NDA .LT. 0 .OR. NDA .GT. NV
C                      =-354: NDU .LT. 0 .OR. NDU .GT. NU
C                      =-355: NDL .LT. 0 .OR. NDL .GT. NL
C
C
C* Work space parameters
C  ---------------------
C  
C  The calling program has to provide MEXX with sufficient work space.
C  The minimum work space required depends on the dimension of the
C  position vector (NP), the number of constraints (NL), and the
C  dimension of the dynamic variable (NU) mainly. It has to be adjusted
C  according to the selected storage mode of the mass matrix (AM), the
C  choice of the linear equation solver, and other options set by the
C  user.
C
C  LIWK       Int     declared size of INTEGER work space for MEXX
C  IWK(LIWK)  Int  *  INTEGER work space for MEXX
C                     Required minimum for the default case:
C
C                       NP + 4*NV + NL + NU + 60
C
C                     Subtract NL if JOB=3
C
C                     If MGPFL=1 add
C                       12*(NV+NL) + 2*NZGMAX +2*NZFMAX 
C                                  + (2+FILLR+FILLC)*NZA
C
C                     with
C                       FILLR = 1.1 (see subroutine MXLIN)
C                       FILLC = 2.0 (see subroutine MXLIN)
C                       NZA = 2*NZGMAX + 2*NZFMAX + NBLK*NMRC**2
C
C                     If MSWIT = 1 or = 2 add 5*NSWMAX
C
C
C  LRWK       Int     Declared size of REAL work space for MEXX
C  RWK(LRWK)  Dble *  REAL work space for MEXX
C                     Required minimum for the default case:
C
C                       (NV+NL)**2 + NP*(NGL+18) + NV*(NV+45) + 28*NL
C                                  + NG + 18*NU + 50
C
C                     with NGL = MAX (1,NG,NL)
C
C                     NOTE: Even if NL=0 or/and NU=0 or/and NU=0 use a
C                           value of 1 at least for the estimation of
C                           the work space. 
C
C                     If special options are set the minimum work space
C                     required can be computed by:
C
C                     LD + LA + LG + LF + LX + LO + LR
C                        + NP*(NGL+6) + NV*9 + NL1*6 + NG1 + NU1*6 + 50
C
C                     with
C                       NGL = MAX (1,NG,NL)
C                       NL1 = MAX (1,NL)
C                       NG1 = MAX (1,NG)
C                       NU1 = MAX (1,NU)
C                 
C                       LD = space for decomposition
C                       LA = space for matrix AM
C                       LG = space for matrix G
C                       LF = space for matrix FL
C                       LX = space for extrapolation
C                       LO = space for dense output
C                       LR = space for root finder
C
C                       LD = (NV+NL)**2             , if JOB=0 (default)
C                          = (NV+NL)**2             , if JOB=1
C                          = (NV+NL)**2 + 2*(NV+NL) , if JOB=2
C                          = (NV+NL)**2             , if JOB=3
C                          = NV + NL + FILLC*NZA    , if JOB=10
C
C                            NZA = 2*NZGMAX + NV*NMRC
C                            FILLC = 2.0 (see subroutine MXLIN)
C
C                       LA = NV*NV                 , if MMODE=0
C                                                      (default)
C                          = NBLK*NMRC**2          , if MMODE=1
C
C                       LG = NV*NL                 , if MGPFL=0
C                                                      (default)
C                          = NZGMAX                , if MGPFL=1
C
C                       LX = KDMAX*(NP+3*NV+2*NL+NU)
C
C                            KDMAX = 12 , maximum number of
C                                         extrapolation steps
C                                         (2 .LE. KDMAX .LE. 12)
C                                         see MXJOB(9)
C
C                       LF = 0                     , if MDISC=0
C                                                      (default)
C                          = NV*NL                 , if MGPFL=0
C                          = NZFMAX                , if MGPFL=1
C
C                       LO = 0                     , if MDOUT=0 and
C                                                       MGLOB=0 and
C                                                       MSWIT=0
C                                                      (default)
C                          = ND*(LSAFE+KDENS+1)    , if MDOUT>0 or
C                                                       MGLOB>0 or
C                                                       MSWIT>0
C                          = ND*156  (default)
C                       with
C                          ND = NDP+NDV+NDU+NDA+NDL, if MSWIT=0 and
C                                                       MDIND=1
C                               NP+NV+NU+NV+NL     , if MSWIT>0 or
C                                                       MDIND=0
C
C                          LSAFE = NJ(1)+NJ(2)+ ... +NJ(KDMAX)
C                                  NJ: Step size sequence
C                                = 140             (default)  
C
C                          KDENS = KDMAX + 3
C                                = 15              (default)
C
C                       add 1 if MDOUT=3
C                       add NDOUT if MDOUT=4
C
C                       LR = 0                     , if MSWIT=0
C                                                       (default)
C                          = 6*NSWMAX              , if MSWIT>0
C
C
C ---------------------------------------------------------------------
C
C* Optional input parameters of MXJOB
C  ----------------------------------
C
C  Position   Default Meaning
C
C  MXJOB(1)    0      MDISC  Mode of discretization:
C                      =0   standard case, F=df/dlambda << 1
C                      =1   F is large i.e. there is substantial dry
C                           friction in the system
C
C  MXJOB(2)    0      JOB    Switch for linear algebra according to
C                           different kinematic formulations:
C                      = 0  non-symmetrical full-mode LU decomposition
C                           generally applicable 
C                           for invertible system   |AM GP^T-df/dl|
C                           matrices                |GP      0    |
C                      = 1  symmetrical full-mode LU decomposition
C                           generally applicable if MDISC=0
C                           for invertible symmetric    |AM GP^T|
C                           system matrices             |GP  0  |
C                      = 2  constraint formulation
C                           mass matrix  AM  is symmetric positive
C                           definite.
C                           (Cholesky dec. of AM, QR-dec. of GP)
C                      = 3  minimal coordinate formulation
C                           AM=|M 0|  with symm. positive def. matrix  M
C                              |0 0|    (only  M  is used)
C                           GP=(K J)  with invertible matrix  K
C                           (full mode LU-dec. of K and MM)
C                      = 10 sparse mode LU-decomposition
C                           generally applicable (as JOB=0) but:
C                           AM must be given in block diagonal mode (in
C                           FPROB) i.e. MMODE=1
C                           GP must be given in sparse storage mode (in
C                           FPROB) i.e. MGPFL=1
C                           for MDISC=1 FL must be given in sparse
C                           storage mode (in FPROB) i.e. MGPFL=1
C
C  MXJOB(3)    0      MMODE Storage mode of mass matrix AM:
C                      = 0  full mode
C                      = 1  block diagonal mode
C
C  MXJOB(4)    0      NBLK  If (MMODE .EQ. 1): number of blocks of M
C
C  MXJOB(5)    0      NMRC  If (MMODE .EQ. 1): dimension of a block of
C                           M 
C
C  MXJOB(6)    0      MGPFL Storage mode of matrices GP and FL
C                      = 0  full mode
C                      = 1  sparse mode
C
C  MXJOB(7)    0      NZGMAX If (MGPFL .EQ. 1): maximum expected number
C                           of elements of GP
C
C  MXJOB(8)    0      NZFMAX If (MGPFL .EQ. 1 .and. MDISC .EQ. 1):
C                           maximum expected number of elements of FL
C
C  MXJOB(9)   12      KDMAX Maximum number of extrapolation steps for
C                           integration  (2 .LE. KDMAX .LE. 12)
C                           KDMAX = 0  is reset to 12 internally
C
C  MXJOB(11)   0      MNGEN General print output
C                      = 0  No print output
C                      = 1  Only Error/Warning messages
C                      = 2  additionally information on MEXX
C
C  MXJOB(12)   6      LUGEN Associated logical unit
C                           LUGEN = 0  is reset to 6
C
C  MXJOB(13)   0      MNINT Integration monitor
C                      = 0  No print output
C                      = 1  Short integration monitor
C                      = 2  Note step size reductions
C                      = 3  Explain step size reductions
C                      = 4  Error estimates
C                     Additionally, for sparse linear algebra
C                      = 1  Print error analysis for MA28 routines
C                      = 2  Print warning messages for MA28
C                      = 4  Decompositions with analyze factor (MA28A)
C                      = 5  Decompositions without analyze factor
C                           (MA28B)
C  MXJOB(14)   6      LUINT  Associated logical unit
C                            LUINT = 0  is reset to 6
C
C  MXJOB(15)   0      MGLOB Write data describing global solution
C                     (coefficients for the interpolation polynomial)
C                     for selected components (see MXJOB(33))
C                      = 0  OFF
C                      = 1  ON
C  MXJOB(16)   40     LUGLO  Associated logical unit
C                            LUGLO = 0 is reset to 40
C
C  MXJOB(17)   0      MNSWI Root finder monitor
C                      = 0  No print output
C                      = 1,..,3 
C  MXJOB(18)   30     LUSWI  Associated logical unit
C                            LUSWI = 0 is reset to 6
C
C  MXJOB(19)   0      MNTIM Time monitor/statistics
C                      = 0  OFF
C                      = 1  ON
C  MXJOB(20)   6      LUTIM  Associated logical unit
C                            LUTIM = 0 is reset to 6
C
C  MXJOB(22)   0      IALAM  Computation of initial acceleration and
C                           Lagrange multipliers
C                      = 0  the values given in A(NV), RLAM(NL) are used
C                      = 1  A(NV), RLAM(NL) are determined by
C                           extrapolation 
C
C  MXJOB(30)   0      MOUT  Print mode for intermediate values of
C                           solution
C                      = 0  No output
C                      = 1  Call user provided subroutine SOLOUT after
C                           each accepted step
C  MXJOB(31)   0      MDOUT Print mode for dense output i.e.
C                           interpolated values for special components
C                           of the solution vectors
C                      = 0  No output
C                      > 0  Call user provided subroutine DENOUT
C                           according to the following rules:
C
C                      = 1  Interpolation at NDOUT points which are
C                           uniformly distributed between T0,TFIN.
C                           Additionally at T0,TFIN. Thus, DENOUT is
C                           called NDOUT+2 times (NDOUT see MXJOB(32)).
C
C                      = 2  Interpolation at NDOUT points between all
C                           integration points. Additionally at the
C                           integration points. Thus, DENOUT is called
C                           NSTEP*(NDOUT+1) + 1 times.
C
C                      = 3  Interpolation such that the maximum interval
C                           is .LE. RWK(51). If the current integration
C                           interval is .GT. RWK(51), it is uniformly
C                           subdivided such that the required condition
C                           holds with a minimum of interpolation
C                           points.
C                           
C                      = 4  Interpolation at NDOUT user prescribed
C                           output points stored in RWK(51), ...,
C                           RWK(50+NDOUT). All integration points may be
C                           ignored.  Thus, DENOUT is called NDOUT
C                           times.
C                           
C  MXJOB(32)   0      NDOUT Number of output points for dense output
C                           (see MXJOB(31)).
C
C  MXJOB(33)   0      MDIND Switch if dense output (see MXJOB(32)) or
C                           global solution (see MXJOB(15)) is desired
C                           for all or part of the solution components
C                      = 0  Use all components for dense output
C                      = 1  Lists of components selected for dense
C                           output are given by the user in IWK.
C
C  MXJOB(35)   0      MSWIT Switch if a root of the user defined
C                           function FSWIT should be found
C                      = 0  Off, FSWIT may be dummy
C                      = 1  On,  MEXX will stop if a root is found
C                      = 2  On,  MEXX will report the root and continue
C
C  MXJOB(36)   0      NSWMAX Maximum number of components in the user
C                           defined root function FSWIT
C                           (1 .LE. NSWMAX, if MSWIT .GT. 0, see
C                           MXJOB(35))
C  MXJOB(37)   1      NSUB  Number of sub-intervals that should be
C                           checked for sign change during root finding
C                           NSUB = 0  is reset to 1
C
C     Output:
C
C  MXJOB(101)    NZGACT If (MGPFL .EQ. 1): number of non-zero elements
C                   of GP encountered actually
C
C  MXJOB(102)    IPCG   If (MGPFL .EQ. 1): flag if sparse matrix GP has
C                   a new pattern
C
C  MXJOB(103)    If (MGPFL .EQ. 1): reserved (internal flag if pattern
C                   of sparse matrix GP has been analyzed at least
C                   once)
C
C  MXJOB(104)    NZFACT If (MGPFL .EQ. 1): number of non-zero elements
C                   of FL encountered actually
C
C  MXJOB(105)    IPCF   If (MGPFL .EQ. 1): flag if sparse matrix FL has
C                   a new pattern
C
C  MXJOB(106)    If (MGPFL .EQ. 1): reserved (internal flag if pattern
C                   of sparse matrix FL has been analyzed at least
C                   once)
C
C  MXJOB(41)     LRWKMI Length of minimum required REAL work space
C  MXJOB(42)     LIWKMI Length of minimum required INTEGER work space
C  MXJOB(43)     LRWKUS Length of available REAL work space
C  MXJOB(44)     LIWKUS Length of available INTEGER work space (in case
C                   of sparse linear algebra mode (MGPFL=1) the minimum
C                   required work space and the work space actually
C                   available may differ)
C
C  51 to 59      Statistics
C
C  MXJOB(51)     NSTEP  Number of all computed steps
C  MXJOB(52)     NACCPT Number of accepted steps
C  MXJOB(53)     NREJCT Number of steps, where ERR .GT. TOL
C                       (after at least one step has been accepted)
C  MXJOB(54)     NFPROB Number of calls to FPROB
C  MXJOB(55)     NFCN   Number of f-evaluations
C  MXJOB(56)     NGCN   Number of g-evaluations
C  MXJOB(57)     NMUL   Number of calls to MMULT
C  MXJOB(58)     NDEC   Number of decompositions (calls to ADEC)
C  MXJOB(59)     NSOL   Number of forward-backward substitutions (ASOL)
C
C  MXJOB(73)     IFAILS Error code encountered in any subroutine which
C                   caused the integrator to stop.
C  MXJOB(74)     ISUB   Decomposition subroutine causing an error
C                   = 10 DGEFA with KIN=0
C                   = 1  DSIFA with KIN=1
C                   = 2  DCHDC with KIN=2, MODE=0
C                   = 3  DCHDC with KIN=2, MODE=1
C                   = 4  DTRSL with KIN=2, MODE=0
C                   = 5  DTRSL with KIN=2, MODE=1
C                   = 6  DGEFA with KIN=3, decomposition of K
C                   = 7  DGEFA with KIN=3, decomposition of MM
C                   = 8  MA28AD with KIN=10
C                   = 9  MA28BD with KIN=10
C
C  80 to 98      Reserved for work space management
C
C  MXJOB(80)     LRWKM  Length of REAL work space left for linear 
C                   algebra
C  MXJOB(81)     LIWKM  Length of INTEGER work space left for linear
C                   algebra
C  MXJOB(82)     NVLDIM NVDIM + NLDIM
C  MXJOB(83)            Length of REAL work space required for
C                   discretization
C  MXJOB(84)            Length of INTEGER work space required for
C                   discretization, dense output and root finder
C  MXJOB(85)     LRWKO  Length of REAL work space required for dense
C                   output and root finder
C
C  MXJOB(86)     LRWKMI Minimum length of REAL work space required for
C                   decomposition
C  MXJOB(87)     LIWKMI Minimum length of INTEGER work space required
C                   for decomposition
C                   NOTE: in case of MGPFL=1 i.e. sparse
C                   decomposition, the work space required is dependent
C                   on fill-in parameters (see subroutine MXLIN)
C
C  MXJOB(88)     LRWKMact Length of REAL work space available for
C                   decomposition
C  MXJOB(89)     LIWKMact Length of INTEGER work space available for
C                   decomposition
C
C                For sparse decomposition only (MGPFL=1)
C                a priori settings:
C
C  MXJOB(90)     NZA    Number of non-zero elements of matrix A
C  MXJOB(91)     LIRNmin  Number of row indices for decomposed matrix,
C                   = NZA*FILLR (see subroutine MXLIN)
C  MXJOB(92)     LICNmin  Length of decomposed matrix A and array of
C                   column indices, = NZA*FILLC (see subroutine MXLIN)
C  MXJOB(93)     LIRNact  Length of work space available for row
C                   indices, LIRNact .GE. LIRNmin
C  MXJOB(94)     LICNact  Length of work space available for decomposed
C                   matrix and column indices, LICNact .GE. LICNmin
C
C                a posteriori measurements:
C                (see Harwell Subroutine Library Specification)
C
C  MXJOB(95)     MINIRN   and
C  MXJOB(96)     MINICN   Minimum values for LIRN and LICN enabling a
C                   successful run on an identical matrix
C  MXJOB(97)     IRNCP    and
C  MXJOB(98)     ICNCP    Flag if it could pay to increase or to reduce
C                   LIRN or LICN for subsequent runs
C
C
C* Optional INTEGER input/output in IWK:
C  -------------------------------------
C
C  Position           Meaning
C
C   51..55            Number of components selected for dense output or
C                     global solution if MDIND=1 (see MXJOB(33)).
C  IWK(51)    in      NDP .LE. NP  Number of components from vector P
C  IWK(52)    in      NDV .LE. NV                       from vector V
C  IWK(53)    in      NDU .LE. NU                       from vector U
C  IWK(54)    in      NDA .LE. NV                       from vector A
C  IWK(55)    in      NDL .LE. NV                       from vector RLAM
C
C   56-60             Reserved
C
C                     List of components selected for dense output or
C                     global solution (in increasing order !!!!!) if
C                     MDIND=1 (see MXJOB(33)).
C  IWK(61)    in      INDP(NP) List of components from vector P
C  IWK(LDV)   in      INDV(NV)                    from vector V
C                       with LDV = 61 + NP
C  IWK(LDU)   in      INDU(NU)                    from vector U
C                       with LDU = 61 + NP +NV
C  IWK(LDA)   in      INDA(NV)                    from vector A
C                       with LDA = 61 + NP + NV + NU
C  IWK(LDL)   in      INDL(NV)                    from vector RLAM
C                       with LDL = 61 + NP +2*NV + NU
C
C  IWK(LISW)  out     ISWIT(NSWMAX) 
C                       with LISW = 61 + NP +3*NV + NU
C                     indicates the components of FSWIT for which a root
C                     was encountered 
C                     ( ISWIT(l)=1 if component l of G equals zero ).
C
C
      PARAMETER (ZERO = 0.D0)
      INTEGER MNGEN
      LOGICAL QPRTIM
C
C --- Arrays of length KD
C
      PARAMETER (KD = 12)
      DIMENSION NJ(KD), HH(KD), WK(KD), AA(KD)
C     NJ(1).GE.2 required
C
      SAVE
C
      CHARACTER*67 CMODUL
      DATA CMODUL /'@(#) MEXX pre-release 1.2 (1.10) from 96/07/18
     $     '/
C
      DATA NJ/2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 30/,
     $     HH /KD*0.0D0/, WK /KD*0.0D0/, AA /KD*0.0D0/
      
C
      DATA IRET /0/
C
      DOUBLE PRECISION UROUND, EPMACH, SMALL
C
C  Initiation
C  ----------
C
      CALL ZIBCONST(EPMACH,SMALL)
      UROUND = 1.0D1*EPMACH
      IERR = 0
      QPRTIM = .FALSE.
      DO 1000 I=1,100
 1000 MXJOB(50+I) = 0
C
C
C-----------------------------------------------------------------
C  Check for valid input
C-----------------------------------------------------------------
C
C  Check t, NP, NV, NL, NG, .....
C---------------------
C
      IF (NP .LE. 0) THEN
         IERR = -21
         GOTO 1080
      ENDIF
      IF (NG .LT. 0) THEN
         IERR = -22
         GOTO 1080
      ENDIF
      IF (NU .LT. 0) THEN
         IERR = -23
         GOTO 1080
      ENDIF
      IF (TFIN .LE. T) THEN
         IERR = -2
      ENDIF
      IF (NG .GT. NP) THEN
         IERR = -25
         GOTO 1080
      ENDIF
C
C  Check optional input MXJOB
C----------------------------
C
C
C  Check I/O
C-----------
C
      MNGEN = MXJOB(11)
      IF (MXJOB(12) .EQ. 0) MXJOB(12) = 6
      LUGEN = MXJOB(12)
      IF (LUGEN .LE. 0 .AND. MNGEN .GT. 0) THEN
         LUGEN = 6
         IERR = -112
         GOTO 1080
      ENDIF
C
C  first printout
C
      IF (MNGEN .GE. 2) THEN
         WRITE (LUGEN,*) ' '
         WRITE (LUGEN,*) ' ***** M E X X 2 1  rev. 1.2 ',
     $        '(1.10) *****'
         WRITE (LUGEN,*) ' from 96/07/18'
         WRITE (LUGEN,*) ' '
      ENDIF
C
C  Check matrix storage/decomposition mode
C
      MDISC = MXJOB(1)
      JOB = MXJOB(2)
      MMODE = MXJOB(3)
      NBLK = MXJOB(4)
      NMRC = MXJOB(5)
      MGPFL = MXJOB(6)
      IF (MGPFL .EQ. 1) THEN
         NZGMAX = MXJOB(7)
         NZFMAX = MXJOB(8)
      ELSE
         NZGMAX = 1
         NZFMAX = 1
      ENDIF
C
      IF (JOB .LT. 0 .OR. JOB .GT. 10) THEN
         IERR = -101
         GOTO 1080
      ENDIF
      IF (JOB .GT. 3 .AND. JOB .LT. 10) THEN
         IERR = -101
         GOTO 1080
      ENDIF
      IF (MMODE .LT. 0 .OR. MMODE .GT. 1) THEN
         IERR = -102
         GOTO 1080
      ENDIF
      IF (JOB .EQ. 10 .AND. MMODE .NE. 1) THEN
         IERR = -151
         GOTO 1080
      ENDIF
      IF (JOB .EQ. 3 .AND. NG .LE. 0) THEN
         IERR = -152
         GOTO 1080
      ENDIF
      IF (MMODE .EQ. 1) THEN
         NBLKH = NV/NMRC
         IF (NBLKH .NE. NBLK) THEN
            IERR = -153
            GOTO 1080
         ENDIF
      ENDIF
      IF (JOB .EQ. 10 .AND. NZGMAX .LE. 0) THEN
         IERR = -154
         GOTO 1080
      ENDIF
C
      IF (MXJOB(9) .EQ. 0) THEN
         KDMAX = KD
      ELSE IF (MXJOB(9) .GE. 2 .AND. MXJOB(9) .LE. KD) THEN
         KDMAX = MXJOB(9)
      ELSE
         IF (MNGEN .GT. 0) THEN
            WRITE (LUGEN,*) ' MEXX - Warning:'
            WRITE (LUGEN,*) ' Illegal value of KDMAX, ', MXJOB(9),
     $           ', reset to default value ', KD, '.'
         ENDIF
         KDMAX = KD
      ENDIF
C
      MNINT = MXJOB(13)
      IF (MXJOB(14) .EQ. 0) MXJOB(14) = 6
      LUINT = MXJOB(14)
      IF (LUINT .LE. 0 .AND. MNINT .GT. 0) THEN
         LUINT = 6
         IERR = -114
         GOTO 1080
      ENDIF
      MGLOB = MXJOB(15)
      IF (MXJOB(16) .EQ. 0) MXJOB(16) = 40
      LUGLO = MXJOB(16)
      IF (LUGLO .LE. 0 .AND. MGLOB .GT. 0) THEN
         LUGLO = 40
         IERR = -116
         GOTO 1080
      ENDIF
      MNSWI = MXJOB(17)
      IF (MXJOB(18) .EQ. 0) MXJOB(18) = 6
      LUSWI = MXJOB(18)
      IF (LUSWI .LE. 0 .AND. MNSWI .GT. 0) THEN
         LUSWI = 6
         IERR = -118
         GOTO 1080
      ENDIF
      MNTIM = MXJOB(19)
      IF (MXJOB(20) .EQ. 0) MXJOB(20) = 6
      LUTIM = MXJOB(20)
      IF (LUTIM .LE. 0 .AND. MNTIM .GT. 0) THEN
         LUTIM = 6
         IERR = -120
         GOTO 1080
      ENDIF
C
      IALAM = MXJOB(22)
      IF (IALAM .LT. 0 .OR. IALAM .GT. 1) THEN
         IERR = -122
         GOTO 1080
      ENDIF
C
C------------------------------
C  Check (dense) output options
C------------------------------
C
      MOUT = MXJOB(30)
      MDOUT = MXJOB(31)
      NDOUT = MXJOB(32)
      MDIND = MXJOB(33)
C
      IF (MOUT .LT. 0 .OR. MOUT .GT. 1) THEN
         IERR = -130
         GOTO 1080
      ENDIF
C
      IF (MDOUT .LT. 0 .OR. MDOUT .GT. 4) THEN
         IERR = -131
         GOTO 1080
      ENDIF    
C
      IF (NDOUT .LT. 0) THEN
         IERR = -132
         GOTO 1080
      ENDIF    
C
      IF (MDIND .LT. 0 .OR. MDIND .GT. 1) THEN
         IERR = -133
         GOTO 1080
      ENDIF    
C
      IF (MDOUT .EQ. 0) THEN
         IF (NDOUT .GT. 0) THEN
            IF (MNGEN .GT. 0) THEN
               WRITE (LUGEN,*) ' MEXX - Warning:'
               WRITE (LUGEN,*)
     $              ' MDOUT = 0  but NDOUT > 0.  NDOUT ignored'
            ENDIF
         ENDIF
      ENDIF
      IF (MDOUT .EQ. 0 .AND. MGLOB .EQ. 0) THEN
         IF (MDIND .EQ. 1) THEN
            IF (MNGEN .GT. 0) THEN
               WRITE (LUGEN,*) ' MEXX - Warning:'
               WRITE (LUGEN,*)
     $              ' MDOUT = MGLOB = 0  but MDIND = 1.  MDIND ignored'
            ENDIF
         ENDIF
      ENDIF
      IF ((MDOUT .GT. 0 .OR. MGLOB .GT. 0) .AND. MDIND .EQ. 1) THEN
         NDP = IWK(51)
         NDV = IWK(52)
         NDU = IWK(53)
         NDA = IWK(54)
         NDL = IWK(55)
C  test for wrong numbers of d.o. components
         IF (NDP .LT. 0 .OR. NDP .GT. NP) THEN
            IERR = -351
            GOTO 1080
         ENDIF  
         IF (NDV .LT. 0 .OR. NDV .GT. NV) THEN
            IERR = -352
            GOTO 1080
         ENDIF  
         IF (NDA .LT. 0 .OR. NDA .GT. NV) THEN
            IERR = -353
            GOTO 1080
         ENDIF  
         IF (NDU .LT. 0 .OR. NDU .GT. NU) THEN
            IERR = -354
            GOTO 1080
         ENDIF  
         IF (NDL .LT. 0 .OR. NDL .GT. NL) THEN
            IERR = -355
            GOTO 1080
         ENDIF  
      ENDIF  
C
      IF (MDOUT .EQ. 3) THEN
         OUTMAX = RWK(51)
      ELSE
         OUTMAX = ZERO
      ENDIF
C
C
C------------------------------
C  Check root finder options
C------------------------------
C
      MSWIT = MXJOB(35)
      IF (MSWIT .EQ. 0) THEN
         NSWMAX = 0
      ELSE
         NSWMAX = MXJOB(36)
      ENDIF
C      NSUB = MXJOB(37)
      IF (MSWIT .GT. 0 .AND. MXJOB(37) .EQ. 0) MXJOB(37) = 1
C
C  set defaul values in iwk for call of user dense o. function
C
      IF ((MDOUT .GT. 0 .OR. MGLOB .GT. 0) .AND. MDIND .EQ. 0) THEN
         IWK(51) = NP
         IWK(52) = NV
         IWK(53) = NU
         IWK(54) = NV
         IWK(55) = NL
         INDPP=60
         INDVP = INDPP + NP
         INDUP = INDVP + NV
         INDAP = INDUP + NU
         INDLP = INDAP + NV
         DO 1010 I=1,NP
            IWK(INDPP+I)=I
 1010    CONTINUE
         DO 1020 I=1,NV
            IWK(INDVP+I)=I
 1020    CONTINUE
         DO 1030 I=1,NU
            IWK(INDUP+I)=I
 1030    CONTINUE
         DO 1040 I=1,NV
            IWK(INDAP+I)=I
 1040    CONTINUE
         DO 1050 I=1,NL
            IWK(INDLP+I)=I
 1050    CONTINUE
      ENDIF
C
C  determine type of dense output generation
C  (dense output selection by user may conflict with necissity of dense
C   output generation for all comp. for root finding)
C
      IDOFLG=0
      IF (MDOUT .GT. 0 .OR. MGLOB .GT. 0) THEN
         IF (MDIND .EQ. 1) THEN
            IDOFLG=1
         ELSE
            IDOFLG=2
         ENDIF
      ENDIF
      IF (MSWIT .GT. 0) IDOFLG=2
C
      IF (IDOFLG .EQ. 0) THEN
         NDP = 0
         NDV = 0
         NDU = 0
         NDA = 0
         NDL = 0
      ENDIF
C
      IF (IDOFLG .EQ. 2) THEN
         NDP = NP
         NDV = NV
         NDU = NU
         NDA = NV
         NDL = NL
      ENDIF
C
C  Check tolerances
C------------------
C
      IF (ITOL .EQ. 0) THEN
         IF (ATOL(1) .LE. 0.D0 .OR. RTOL(1) .LE. 10.D0*UROUND) THEN
            IERR = 11
            IF (MNGEN .GT. 0) THEN
               WRITE (LUGEN,*) ' Warning:'
               WRITE (LUGEN,*) ' A prescribed tolerance is too small'
            ENDIF
         ENDIF
         TOLMIN = ATOL(1) + RTOL(1)
      ELSE
         TOLMIN = 1.D0
         N = NP + NV + NU
         DO 1060 I=1,N
            IF (ATOL(I) .LE. 0.D0 .OR. RTOL(I) .LE. 10.D0*UROUND) THEN
               IERR = 12
               IF (MNGEN .GT. 0) THEN
                  WRITE (LUGEN,*) ' Warning:'
                  WRITE (LUGEN,*)
     $                 ' A prescribed tolerance is too small'
               ENDIF
            ENDIF
            TOLMIN = MIN (TOLMIN, ATOL(I) + RTOL(I))
 1060    CONTINUE
      ENDIF
      RWK(50) = TOLMIN
C
C  Check optional input from RWK
C-------------------------------
C
      HMAX = RWK(22)
      IF (HMAX .LT. UROUND) HMAX = TFIN - T
C
C
C  Initialize time monitor
C-------------------------
C
      QPRTIM = MNTIM .GT. 0
      IF (QPRTIM) THEN
         CALL MONINI ('MEXX-rev. 1.2 (1.10)', LUTIM)
         CALL MONSTR (IRET)
         IF (IRET .EQ. 0) THEN
            CALL MONDEF ( 1, 'PROJEC')
            CALL MONDEF ( 2, 'MEXXEX')
            CALL MONDEF ( 3, 'Dense output')
            CALL MONDEF ( 4, 'MXIPOL')
            CALL MONDEF ( 5, 'F-f-pdot-udot')
            CALL MONDEF ( 6, 'F-am-gp-gi')
            CALL MONDEF ( 7, 'F-am-dgdp-g')
            CALL MONDEF ( 8, 'F-pdot')
            CALL MONDEF ( 9, 'F-g')
            CALL MONDEF (10, 'F-f-pdot-udot-dfdl')
            CALL MONDEF (11, 'Switch function')
            CALL MONDEF (12, 'ASOL')
            CALL MONDEF (13, 'ADEC')
            IF (JOB .EQ. 2) THEN
               CALL MONDEF (14, 'DCHDC')
               CALL MONDEF (15, 'DTRSL')
               CALL MONDEF (16, 'DQRDC')
            ELSE IF (JOB .EQ. 10) THEN
               CALL MONDEF (14, 'MA28A')
               CALL MONDEF (15, 'MA28B')
            ENDIF
         ELSE
            QPRTIM = .FALSE.
         ENDIF
      ENDIF
C
C-----------------------------------------------------------------
C  Work space distribution
C-----------------------------------------------------------------
C
C
C  some dimensions
C
      NPDIM = MAX (1, NP)
      NVDIM = MAX (1, NV)
      NLDIM = MAX (1, NL)
      NGDIM = MAX (1, NG)
      NUDIM = MAX (1, NU)
      NZGDIM = NZGMAX
      NZFDIM = NZFMAX
      N6UDIM = NPDIM + 3*NVDIM + 2*NLDIM + NUDIM
      N2UDIM = NPDIM + NVDIM + NUDIM
      NVLDIM = NVDIM + NLDIM
      NTODIM = 1
C
      IF (ITOL .NE. 0) NTODIM = N2UDIM
C
      IF (MMODE .EQ. 0) THEN
         MAMDIM = NVDIM**2
         LDA = NVDIM
         MDA = NVDIM
      ELSE IF (MMODE .EQ. 1) THEN
         MAMDIM = NBLK*NMRC**2
         LDA = NBLK*NMRC
         MDA = NMRC
      ENDIF
      IF (MGPFL .EQ. 0) THEN
         LDG = MAX(NLDIM, NGDIM)
         MDG = MAX(NVDIM, NPDIM)
C              GP(NLDIM,NVDIM): velocity constraint matrix G(t,p)
C              GP(NGDIM,NPDIM): position constraint Jacobian dg/dp(t,p)
         LDF = NVDIM
         MDF = NLDIM
      ELSE
         MGPDIM = NZGDIM
         LDG = NZGDIM
         MDG = 1
         MFLDIM = NZFDIM
         LDF = NZFDIM
         MDF = 1
      ENDIF
      MGPDIM = LDG*MDG
      IF (MDISC .EQ. 0) THEN
         LDF = 1
         MDF = 0
      ENDIF
      MFLDIM = LDF*MDF
      IF (MFLDIM .NE. 0) THEN
         NFL = NV
      ELSE
         NFL = 0
      ENDIF
C
C  Pointers to vectors for discretization:
C  IWOPT(50), RWOPT(50), DOUTS(NDOUT), P1(NP), V1(NV), A1(NV), RL1(NL),
C  U1(NU), TAB(N6UDIM,KDMAX),
C
C  WF(NV), WP(NP), WV(NV), WA(NV), WR(NL),
C  F0(NV), F1(NV), WPDOT(NP), PDOT0(NP), PDOT1(NP),
C  G(NG), GI(NL), WU(NU), WUDOT(NU), UDOT0(NU), UDOT1(NU),
C  WZ(NV+NL), SCAL(N2UDIM)
C
C  If MDISC=1: FL0(NV), FL1(NV)
C
      JDOUTS = 51
      MDOUT=MXJOB(31)
      NDOUT=MXJOB(32)
      IF (MDOUT .EQ. 3) THEN
         JP1 = JDOUTS + 1
      ELSE IF (MDOUT .EQ. 4) THEN
         JP1 = JDOUTS + NDOUT
      ELSE
         JP1 = JDOUTS
      ENDIF
C
      JV1 = JP1 + NPDIM
      JA1 = JV1 + NVDIM
      JRL1 = JA1 + NVDIM
      JU1 = JRL1 + NLDIM
      JTAB = JU1 + NUDIM
      JF = JTAB + KDMAX*N6UDIM
      JP = JF + NVDIM
      JV = JP + NPDIM
      JA = JV + NVDIM
      JRL = JA + NVDIM
      JF0 = JRL + NLDIM
      JF1 = JF0 + NVDIM
      JPDOT = JF1 + NVDIM
      JPDOT0 = JPDOT + NPDIM
      JPDOT1 = JPDOT0 + NPDIM
      JG = JPDOT1 + NPDIM
      JGI = JG + NGDIM
      JU = JGI + NLDIM
      JUDOT = JU + NUDIM
      JUDOT0 = JUDOT + NUDIM
      JUDOT1 = JUDOT0 + NUDIM
      JZ = JUDOT1 + NUDIM
      JSCAL = JZ + NVLDIM
      JFL0 = JSCAL + N2UDIM
      JFL1 = JFL0 + NFL
C
C  Pointers to matrices AM(LDA,MDA), GP(LDG,MDG), FL(LDF,MDF)
C  (passed to FPROBL and MXLIN)
C
      JAM = JFL1 + NFL
      JGP = JAM + MAMDIM
      JFL = JGP + MGPDIM
C  
C  INTEGER work space for dense output INDP(NP), INDV(NV), INDU(NU),
C  INDA(NV), INDL(NV)
C
      JINDP = 61
      JINDV = JINDP + NP
      JINDU = JINDV + NV
      JINDA = JINDU + NU
      JINDL = JINDA + NV
C  
C  INTEGER work space for root finder ISWIT(NSWMAX), IRESL(NSWMAX),
C  IRESR(NSWMAX), IRES0(NSWMAX), IRES1(NSWMAX)
C
      JISW = JINDL + NV
      JRESL = JISW + NSWMAX
      JRESR = JRESL + NSWMAX
      JRES0 = JRESR + NSWMAX
      JRES1 = JRES0 + NSWMAX
C  
C  
C  REAL work space for dense output YSAFE(ND,LSAFE), DENS(ND,KDENS), and
C  YIP(ND)
C  
      IF (IDOFLG .NE. 0) THEN
         LSAFE = 0
         DO 1070 I=1,KDMAX
            LSAFE = LSAFE + NJ(I)
 1070    CONTINUE 
C
         ND = NDP + NDV + NDU + NDA + NDL
         NDH = ND
C
C        KDENS .GE. ILAM(KDOMAX)+IRHO(KDOMAX)+3
         KDENS = KDMAX + 3
C
      ELSE
C
         LSAFE = 1
         ND = 1
         NDH = 0
         KDENS = 1
      ENDIF
C
      JYSAFE = JFL + MFLDIM
      JDENS = JYSAFE + NDH*LSAFE
      JYIP = JDENS + NDH*KDENS
C
C  Pointer to work space for root finder GVALW1(NSWMAX), GVAL0(NSWMAX), 
C  GVAL1(NSWMAX), GVALW2(NSWMAX), GVALW3(NSWMAX), GVALW4(NSWMAX)
C
      JGVW1 = JYIP + NDH
      JGV0 = JGVW1 + NSWMAX
      JGV1 = JGV0 + NSWMAX
      JGVW2 = JGV1 + NSWMAX
      JGVW3 = JGVW2 + NSWMAX
      JGVW4 = JGVW3 + NSWMAX
C
C  Pointer to work space for linear algebra
C  RWK(LRWKM), IWK(LIWKM)
C
      JLINR = JGVW4 + NSWMAX
      JLINI = JRES1 + NSWMAX
C
C---------------------------------
C  check for sufficient work space
C---------------------------------
C
      MXJOB(83) = JYSAFE - 1
      MXJOB(84) = JLINI - 1
C     
      LRWKO = JLINR - JYSAFE
      LRWKM = LRWK - MXJOB(83) - LRWKO
      LIWKM = LIWK - MXJOB(84)
      MXJOB(80) = LRWKM
      MXJOB(81) = LIWKM
      MXJOB(82) = NVLDIM
      MXJOB(85) = LRWKO
C
C
      IFAIL = 0
      CALL MXLIN (NV, NL, NVLDIM, MXJOB, IFAIL)
C
      MXJOB(41) = MXJOB(83) + MXJOB(85)  + MXJOB(86)
      MXJOB(42) = MXJOB(84) + MXJOB(87)
      MXJOB(43) = MXJOB(83) + MXJOB(85)  + MXJOB(88)
      MXJOB(44) = MXJOB(84) + MXJOB(89)
C
C-----------------------------------------------------------------
C  general initial output
C-----------------------------------------------------------------
C
      IF (MNGEN .GE. 2) THEN
         WRITE (LUGEN,*) ' Required tolerance:'
         IF (ITOL .EQ. 0) THEN
            WRITE (LUGEN,*) '  Scalar RTOL = ', RTOL(1)
            WRITE (LUGEN,*) '  Scalar ATOL = ', ATOL(1)
         ENDIF
         IF (ITOL .EQ. 1) THEN
            WRITE (LUGEN,*) '  Componentwise RTOL '
            WRITE (LUGEN,*) '  Componentwise ATOL '
         ENDIF
         WRITE (LUGEN,*) ' '
         WRITE (LUGEN,*) ' Selected linear algebra mode:'
         WRITE (LUGEN,*) '  MMODE = ', MMODE
         WRITE (LUGEN,*) '  JOB = ', JOB
         WRITE (LUGEN,*) ' Selected solution output mode:'
         WRITE (LUGEN,*) '  MOUT = ', MOUT
         WRITE (LUGEN,*) '  MDOUT = ', MDOUT
         IF (MDOUT .EQ. 1 .OR. MDOUT .EQ. 2 .OR. MDOUT .EQ. 4)
     $        WRITE (LUGEN,*) '  NDOUT = ', NDOUT
         IF (MDOUT .EQ. 3)
     $        WRITE (LUGEN,*) '  OUTMAX = ', OUTMAX
         WRITE (LUGEN,*) ' '
         WRITE (LUGEN,*) ' Work space distribution:'
         WRITE (LUGEN,*) '  REAL:'
         WRITE (LUGEN,*) '   provided:                    ', LRWK
         WRITE (LUGEN,*) '   minimum required:            ', MXJOB(41)
         WRITE (LUGEN,*) '    discretization and matrices:', MXJOB(83)
         WRITE (LUGEN,*) '    matrix decomposition:       ', MXJOB(86)
         WRITE (LUGEN,*) '   actually available:          ', MXJOB(43)
         WRITE (LUGEN,*) '    discretization and matrices:', MXJOB(83)
         WRITE (LUGEN,*) '    matrix decomposition:       ', MXJOB(88)
         WRITE (LUGEN,*) '   dense output:                ', MXJOB(85)
         WRITE (LUGEN,*) '  INTEGER:'
         WRITE (LUGEN,*) '   provided:                    ', LIWK
         WRITE (LUGEN,*) '   minimum required:            ', MXJOB(42)
         WRITE (LUGEN,*) '    discretization and matrices:', MXJOB(84)
         WRITE (LUGEN,*) '    matrix decomposition:       ', MXJOB(87)
         WRITE (LUGEN,*) '   actually available:          ', MXJOB(44)
         WRITE (LUGEN,*) '    discretization and matrices:', MXJOB(84)
         WRITE (LUGEN,*) '    matrix decomposition:       ', MXJOB(89)
      ENDIF
C
      MXJOB(103) = -1
C
C
C
      IF (IFAIL .LT. 0) THEN
C        Error exit if workspace not sufficient
         IF (IFAIL .EQ. -1) IERR = -3
         IF (IFAIL .EQ. -2) IERR = -4
         GOTO 1080
      ENDIF
C
C
      CALL MEXX21 (NP, NV, NL, NG, NU, FPROB, SOLOUT, DENOUT,
     $     T, TFIN, P, V, U, A, RLAM,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM, N2UDIM, N6UDIM,
     $     NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF, KDMAX,
     $     MXJOB, H, HMAX, RTOL, ATOL, ITOL,
     $     NJ, HH, WK, AA, RWK(1), IWK(1), RWK(JDOUTS),
     $     RWK(JP1), RWK(JV1), RWK(JRL1), RWK(JA1), RWK(JU1),
     $     RWK(JTAB), RWK(JF), RWK(JP), RWK(JV),
     $     RWK(JA), RWK(JRL), RWK(JF0), RWK(JF1), RWK(JPDOT),
     $     RWK(JPDOT0), RWK(JPDOT1), RWK(JG), RWK(JGI),
     $     RWK(JU), RWK(JUDOT), RWK(JUDOT0), RWK(JUDOT1), RWK(JZ),
     $     RWK(JSCAL), RWK(JFL0), RWK(JFL1), RWK(JAM), RWK(JGP),
     $     RWK(JFL), RWK(JLINR), IWK(JLINI),
     $     LRWKM, LIWKM, IERR, UROUND, MOUT, MDOUT, NDOUT, IDOFLG, 
     $     NDP, NDV, NDU, NDA, NDL,
     $     IWK(JINDP), IWK(JINDV), IWK(JINDU), IWK(JINDA), IWK(JINDL),
     $     ND, LSAFE, KDENS, RWK(JYSAFE), RWK(JDENS), RWK(JYIP),
     $     NSWMAX, FSWIT, IWK(JISW), IWK(JRESL), IWK(JRESR), IWK(JRES0),
     $     IWK(JRES1), RWK(JGVW1), RWK(JGV0), RWK(JGV1), RWK(JGVW2),
     $     RWK(JGVW3), RWK(JGVW4))  
C
C
      IF (IERR .LT. 0) GOTO 1080
C
C
      IF (MNGEN .GE. 2) THEN
         WRITE (LUGEN,*) ' '
         WRITE (LUGEN,*) ' MEXX performed successfully'
      ENDIF
      IF (QPRTIM) THEN
         CALL MONHLT
         CALL MONPRT
      ENDIF
C
      RETURN
C
C-----------------------------------------------------------------
C  Error exits
C-----------------------------------------------------------------
C
 1080 CONTINUE
C
      IF (MNGEN .GE. 1) THEN
         WRITE (LUGEN,*) ' MEXX fails'
         WRITE (LUGEN,*) '  Error code:', IERR
C
         IF (IERR .LE. -101 .AND. IERR .GE. -250) THEN
            IERRH = -(IERR + 100)
            WRITE (LUGEN,*) '  - Invalid input for position:', IERRH
            WRITE (LUGEN,*) '  - of array MXJOB'
         ELSE IF (IERR .LE. -301 .AND. IERR .GE. -399) THEN
            IERRH = -(IERR + 100)
            WRITE (LUGEN,*) '  - Invalid input for position:', IERRH
            WRITE (LUGEN,*) '  - of array IWK'
         ENDIF
C
         IF      (IERR .EQ.   -2) THEN
            WRITE (LUGEN,*) '  - Invalid input: TFIN .LE. T '
         ELSE IF (IERR .EQ.   -3) THEN
            WRITE (LUGEN,*) '  - INTEGER work space exhausted'
            WRITE (LUGEN,*) '  - LIWK must be at least:', MXJOB(42)
         ELSE IF (IERR .EQ.   -4) THEN
            WRITE (LUGEN,*) '  - REAL work space exhausted'
            WRITE (LUGEN,*) '  - LRWK must be at least:', MXJOB(41)
         ELSE IF (IERR .EQ.   -9) THEN
            WRITE (LUGEN,*) '  - Stepsize H .LE. machine precision'
         ELSE IF (IERR .EQ.  -10) THEN
            WRITE (LUGEN,*) '  - Limit of integration steps reached'
         ELSE IF (IERR .EQ.  -11) THEN
            WRITE (LUGEN,*) '  - No convergence in projection step'
         ELSE IF (IERR .EQ.  -15) THEN
            WRITE (LUGEN,*) '  - Decomposition failed'
            WRITE (LUGEN,*) '  - Error return from ADEC in MEXXEX'
            WRITE (LUGEN,*) '  - From routine   :', MXJOB(74)
            WRITE (LUGEN,*)  ' - Error code     :', MXJOB(73)
         ELSE IF (IERR .EQ.  -16) THEN
            WRITE (LUGEN,*) '  - Decomposition failed'
            WRITE (LUGEN,*) '  - Error return from ADEC in PROJEC'
            WRITE (LUGEN,*) '  - From routine   :', MXJOB(74)
            WRITE (LUGEN,*) '  - Error code     :', MXJOB(73)
         ELSE IF (IERR .EQ.  -17) THEN
            WRITE (LUGEN,*) '  - User interrupt'
            WRITE (LUGEN,*) '  - Error return from FPROB in MEXXEX'
            WRITE (LUGEN,*) '  - Error code:', MXJOB(73)
         ELSE IF (IERR .EQ.  -18) THEN
            WRITE (LUGEN,*) '  - User interrupt'
            WRITE (LUGEN,*) '  - Error return from FPROB in PROJEC'
            WRITE (LUGEN,*) '  - Error code:', MXJOB(73)
         ELSE IF (IERR .EQ.  -19) THEN
            WRITE (LUGEN,*) '  - User interrupt'
            WRITE (LUGEN,*) '  - Got negative return code from DENOUT'
            WRITE (LUGEN,*) '  - Error code:', MXJOB(73)
         ELSE IF (IERR .EQ.  -20) THEN
            WRITE (LUGEN,*) '  - User interrupt'
            WRITE (LUGEN,*) '  - Got negative return code from MXRTC/F'
            WRITE (LUGEN,*) '  - Error code:', MXJOB(73)
         ELSE IF (IERR .EQ.  -21) THEN
            WRITE (LUGEN,*) '  - Invalid input: NP .LE. 0 '
         ELSE IF (IERR .EQ.  -22) THEN
            WRITE (LUGEN,*) '  - Invalid input: NG .LT. 0 '
         ELSE IF (IERR .EQ.  -23) THEN
            WRITE (LUGEN,*) '  - Invalid input: NU .LT. 0 '
         ELSE IF (IERR .EQ.  -25) THEN
            WRITE (LUGEN,*) '  - Invalid input: NG .GT. NP'
         ELSE IF (IERR .EQ. -101) THEN
            WRITE (LUGEN,*) '  - Invalid input:',
     $           ' JOB is not 0, 1, 2, 3, or 10'
         ELSE IF (IERR .EQ. -102) THEN
            WRITE (LUGEN,*) '  - Invalid input:',
     $           ' MMODE is not 0 or 1'
         ELSE IF (IERR .EQ. -112
     $       .OR. IERR .EQ. -114
     $       .OR. IERR .EQ. -118
     $       .OR. IERR .EQ. -120) THEN
            WRITE (LUGEN,*) '  - Invalid input:',
     $           ' illegal logical unit number set to 6'
         ELSE IF (IERR .EQ. -116) THEN
            WRITE (LUGEN,*) '  - Invalid input:',
     $           ' illegal logical unit number set to 40'
         ELSE IF (IERR .EQ. -122) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: IALAM .LT. 0 .OR. IALAM .GT. 1'
         ELSE IF (IERR .EQ. -130) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: MOUT .LT. 0 .OR. MOUT .GT. 1'
         ELSE IF (IERR .EQ. -131) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: MDOUT .LT. 0 .OR. MDOUT .GT. 4'
         ELSE IF (IERR .EQ. -132) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: NDOUT .LT. 0'
         ELSE IF (IERR .EQ. -133) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: MDIND .LT. 0 .OR. MDIND .GT. 1'
         ELSE IF (IERR .EQ. -151) THEN
            WRITE (LUGEN,*) '  - Invalid input:',
     $           ' JOB .EQ. 10 .AND. MMODE .NE. 1'
         ELSE IF (IERR .EQ. -152) THEN
            WRITE (LUGEN,*) '  - Invalid input:',
     $           ' JOB .EQ. 3 .AND. NG .LE. 0'
         ELSE IF (IERR .EQ. -153) THEN
            WRITE (LUGEN,*) '  - Invalid input: NBLK .NE. NP/NMRC'
         ELSE IF (IERR .EQ. -154) THEN
            WRITE (LUGEN,*) '  - Invalid input:',
     $           ' NZGMAX .LE. 0 .AND. JOB .EQ. 10'
         ELSE IF (IERR .EQ. -155) THEN
            WRITE (LUGEN,*) '  - Invalid input:',
     $           ' JOB .EQ. 10 .AND. MMODE .NE. 0'
         ELSE IF (IERR .EQ. -351) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: NDP .LT. 0 .OR. NDP .GT. NP'
         ELSE IF (IERR .EQ. -352) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: NDV .LT. 0 .OR. NDV .GT. NV'
         ELSE IF (IERR .EQ. -353) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: NDA .LT. 0 .OR. NDA .GT. NV'
         ELSE IF (IERR .EQ. -354) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: NDU .LT. 0 .OR. NDU .GT. NU'
         ELSE IF (IERR .EQ. -355) THEN
            WRITE (LUGEN,*)
     $           '  - Invalid input: NDL .LT. 0 .OR. NDL .GT. NL'
         ENDIF
      ENDIF
C
C
      IF (QPRTIM) THEN
         CALL MONHLT
         CALL MONPRT
      ENDIF
C
      RETURN
C
C  End of subroutine MEXX
C
      END
C
      SUBROUTINE MEXX21 (NP, NV, NL, NG, NU, FPROB, SOLOUT, DENOUT,
     $     T, TFIN, P, V, U, A, RLAM,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM, N2UDIM, N6UDIM,
     $     NVLDIM,
     $     LDA, LDG, LDF, MDA, MDG, MDF, KD,
     $     MXJOB, H, HMAX, RTOL, ATOL, ITOL,
     $     NJ, HH, WK, AA, RWOPT, IWOPT, DOUTS,
     $     P1, V1, RL1, A1, U1,
     $     TAB, WF, WP, WV, WA, WRL,
     $     F0, F1, WPDOT, PDOT0, PDOT1, G, GI,
     $     WU, WUDOT, UDOT0, UDOT1, WZ, SCAL, FL0, FL1, AM, GP, FL,
     $     RWKL, IWKL,
     $     LRWKL, LIWKL, IERR, UROUND, MOUT, MDOUT, NDOUT, IDOFLG, 
     $     NDP, NDV, NDU, NDA, NDL,
     $     INDP, INDV, INDU, INDA, INDL,
     $     ND, LSAFE, KDENS,
     $     YSAFE, DENS, YIP,
     $     NSWIF, FSWIT, ISWIT, IRESL, IRESR, IRES0, IRES1,
     $     GVALW1, GVAL0, GVAL1, GVALW2, GVALW3, GVALW4)
C ---------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER NP, NV, NL, NG, NU
      EXTERNAL FPROB, SOLOUT, DENOUT
      DOUBLE PRECISION T, TFIN, P(NPDIM), V(NVDIM), U(NUDIM), A(NVDIM),
     $     RLAM(NLDIM)
      INTEGER NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM, N2UDIM, N6UDIM,
     $    NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF, KD, MXJOB(150)
      DOUBLE PRECISION H, HMAX, RTOL(NTODIM), ATOL(NTODIM)
      INTEGER ITOL, NJ(KD)
      DOUBLE PRECISION HH(KD), WK(KD), AA(KD)
C
C  Optional parameters
      DOUBLE PRECISION RWOPT(50)
      INTEGER IWOPT(60)
C
C  Work space for discretization
      DOUBLE PRECISION DOUTS(NDOUT), P1(NPDIM), V1(NVDIM), RL1(NLDIM),
     $     A1(NVDIM), U1(NUDIM), TAB(N6UDIM,KD), WF(NVDIM), WP(NPDIM),
     $     WV(NVDIM), WA(NVDIM), WRL(NLDIM), F0(NVDIM), F1(NVDIM),
     $     WPDOT(NPDIM), PDOT0(NPDIM), PDOT1(NPDIM), G(NGDIM),
     $     GI(NLDIM), WU(NUDIM), WUDOT(NUDIM), UDOT0(NUDIM),
     $     UDOT1(NUDIM), WZ(NVLDIM), SCAL(N2UDIM), FL0(LDF), FL1(LDF),
     $     AM(LDA,MDA), GP(LDG,MDG), FL(LDF,MDF)
C
C  Work space for linear algebra
      DOUBLE PRECISION RWKL(LRWKL)
      INTEGER IWKL(LIWKL), LRWKL, LIWKL
C
C  Work space for dense output formula
      INTEGER MOUT, MDOUT, NDOUT, NDP, NDV, NDU, NDA, NDL,
     $     INDP(NDP), INDV(NDV), INDU(NDU), INDA(NDA), INDL(NDL),
     $     ND, LSAFE, KDENS
      DOUBLE PRECISION YSAFE(ND,LSAFE), DENS(ND,KDENS), YIP(ND)
C
C  Parameters and work space for root finder
      INTEGER NSWIF
      EXTERNAL FSWIT
      INTEGER ISWIT(NSWIF), IRESL(NSWIF), IRESR(NSWIF), IRES0(NSWIF),
     $     IRES1(NSWIF)
      DOUBLE PRECISION GVALW1(NSWIF), GVAL0(NSWIF), GVAL1(NSWIF),
     $     GVALW2(NSWIF), GVALW3(NSWIF), GVALW4(NSWIF)
C
C
C* File              mexx21.f
C* Version           1.1
C* Latest Change     96/03/14 (1.6)
C* Modification history
C    1.1             Call MXDDU1 before MXDOUT
C      
C
C --- Parameters for projection step
C
C     EPSREQ   Required accuracy for modified Newton iteration
C              error will be smaller than
C                     (RTOL*ABS(P(I)) + ATOL)*EPSREQ
C
C     ITMAX    Maximum number of modified Newton steps
C
C
      PARAMETER (EPSREQ = 1.D-2, ITMAX = 10)
C
C --- Maximum number of steps:
C
      PARAMETER (NMAX = 3000)
C
C ---
C
C
      INTEGER IERR, INFOS(50), INFOD(50), MNGEN
      DOUBLE PRECISION EPSEST, UROUND
      LOGICAL QPRTIM
      LOGICAL QRJECT, QFAIL, QLAST, QINIT, QSWIT
      LOGICAL QFLAG(9)
C
C ---
C
      COMMON /ERRCOM/ TOLD4, URD, FAC, FAC1, FAC2,
     $                SAFE2, SAFE3, HMAXIM, ERROLD
      COMMON /LOGCOM/ QRJECT, QFAIL, QINIT
C
      SAVE
C
C ---------------------------------------------------------------------
C
      DATA FAC3 /.9D0/, FAC5 /1.D0/, FACP /5.D0/, FACV /5.D0/,
     $     SAFE1 /.80D0/, SAFE4 /0.45D0/
C
      DATA EPSEST /0.0D0/, ILAM /0/, IRHO /0/, NITER /0/, TZERO /0.D0/
C
C
      N = NV
      NPV = NP + NV
C
      NSC = 0
      IDOCNT = 0
C
      T0=T
C
      MDISC = MXJOB(1)
      MGPFL = MXJOB(6)
      IF (MGPFL .EQ. 1) THEN
         NZGMAX = MXJOB(7)
         NZFMAX = MXJOB(8)
      ELSE
         NZGMAX = 1
         NZFMAX = 1
      ENDIF
C
      MNGEN = MXJOB(11)
      LUGEN = MXJOB(12)
      MNINT = MXJOB(13)
      LUINT = MXJOB(14)
      MGLOB = MXJOB(15)
      MNTIM = MXJOB(19)
      QPRTIM = MNTIM .GT. 0
C 
      MDIND = MXJOB(33)
C
      IALAM = MXJOB(22)
C
      TOLMIN = RWOPT(50)
C
      DO 1000 I=51,59
         MXJOB(I) = 0
 1000 CONTINUE
C
      IF (MOUT .GT. 0) THEN
         DO 1010 I=1,40
            INFOS(I) = MXJOB(I)
 1010    CONTINUE 
C
         DO 1020 I=41,50
            INFOS(I) = 0
 1020    CONTINUE
      ENDIF
C
      IF (MDOUT .GT. 0 .OR. MGLOB .GT. 0) THEN
         DO 1030 I=1,40
            INFOD(I) = MXJOB(I)
 1030    CONTINUE 
C
         DO 1040 I=41,50
            INFOD(I) = 0
 1040    CONTINUE
      ENDIF
C
C -- h-Extrapolation: work for j-th line
C
      AA(1) = NJ(1) + 4
      DO 1050 J=2, KD
         AA(J) = AA(J-1) + NJ(J)
 1050 CONTINUE
C
C --
C
      K = MAX (3, MIN (KD-2, INT(-LOG10(TOLMIN)*.6D0 +1.5D0)))
      KC = 0
      H = MIN (H, HMAX, TFIN-T)
C
      URD = UROUND
      HMAXIM = HMAX
      FAC1 = 0.1D0
      FAC2 = 4.D0
      SAFE2 = .93D0
      SAFE3 = .5D0
C
      TOL = 1.D0
      ERR = 0.D0
      WK(1) = 0.D0
      TOLD4 = TOL*SAFE1
C
      FACEND = (1.D0 + SAFE2)*0.5D0
C
      IF (ITOL .EQ. 0) THEN
         TTOL=RTOL(1)
         SCALE = ATOL(1)/RTOL(1)
         DO 1060 I=1,NP
            SCAL(I)  = MAX (ABS(P(I)), SCALE)
 1060    CONTINUE
         DO 1070 I=1,NV
            SCAL(NP+I) = MAX (ABS(V(I)), SCALE)
 1070    CONTINUE
         DO 1080 I=1,NU
            SCAL(NPV+I) = MAX (ABS(U(I)), SCALE)
 1080    CONTINUE
      ELSE
         TTOL=RTOL(1)
         DO 1090 I=1,NP
            SCALE  = ATOL(I)/RTOL(I)
            SCAL(I) = MAX (ABS(P(I)), SCALE)
            TTOL = MIN (RTOL(I),TTOL)
 1090    CONTINUE
         DO 1100 I=1,NV
            SCALE    = ATOL(NP+I)/RTOL(NP+I)
            SCAL(NP+I) = MAX (ABS(V(I)), SCALE)
            TTOL = MIN (RTOL(NP+I),TTOL)
 1100    CONTINUE
         DO 1110 I=1,NU
            SCALE      = ATOL(NPV+I)/RTOL(NPV+I)
            SCAL(NPV+I) = MAX (ABS(U(I)), SCALE)
            TTOL = MIN (RTOL(NPV+I),TTOL)
 1110    CONTINUE
      ENDIF
C
      QRJECT = .FALSE.
      QLAST = .FALSE.
      QSWIT = .FALSE.
C
C ---------------------------------------------------------------------
C
 1120 CONTINUE
      IF (MNINT .GE. 1) THEN
         IF (.NOT. QFAIL .AND. .NOT. QRJECT) THEN
            WRITE (LUINT, 9000) MXJOB(52), MXJOB(55), T, H, KC, K
 9000       FORMAT (' NS, NF, T, H, K, KOPT:', I4, I6, 2(1PD15.7), 2I3)
         ELSE
            IF (MNINT .GE. 2) THEN
               WRITE (LUINT,*)
     $              '  Restart step with H = ', H, '  KOPT = ', K
            ENDIF
         ENDIF
      ENDIF
      QFAIL = .FALSE.
C ------------------------------------------------
C --- is TFIN reached in the next step ? 
C ------------------------------------------------
      H1 = TFIN - T
      IF (H1 .LE. UROUND) GOTO 1470
      H = MIN (H, H1, HMAX)
      IF (H .LE. UROUND) GOTO 1540
      IF (H .GE. H1*FACEND) THEN
         QLAST = .TRUE.
         H = H1
      ENDIF
C ------------------------------------------------
C --- the first and the last step
C ------------------------------------------------
      IF (MXJOB(51) .EQ. 0) THEN
         IF (QPRTIM) CALL MONON (1)
C
         IFAIL = 0
         CALL MXPRJC (NP, NV, NL, NG, NU,
     $        NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $        N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $        FPROB,
     $        T, P, V, U, A, RLAM, SCAL,
     $        RTOL, ATOL, ITOL,
     $        ITMAX, NITER, EPSREQ, EPSEST,
     $        WF, WP, WV, WZ, G, GI, WPDOT, WUDOT, AM, GP, FL,
     $        MXJOB, RWKL, IWKL, LRWKL, LIWKL, IFAIL)
C 
C 
         IF (QPRTIM) CALL MONOFF (1)
         IF (IFAIL .EQ.- 9) GOTO 1600
         IF (IFAIL .NE. 0) GOTO 1570
C
C        Function evaluations: f, pdot, udot, df/dl
C
         QFLAG(1) = .FALSE.
         QFLAG(2) = .FALSE.
         QFLAG(3) = .FALSE.
         QFLAG(4) = .TRUE.
         QFLAG(5) = .TRUE.
         QFLAG(6) = .TRUE.
         QFLAG(7) = .FALSE.
         QFLAG(8) = .FALSE.
         QFLAG(9) = MDISC .EQ. 1
C
         JFAIL = 0
C
         IF (QFLAG(9)) THEN
            IMON = 10
         ELSE
            IMON=5
         ENDIF
         IF (QPRTIM) CALL MONON (IMON)
         CALL MXPROB (NP, NV, NL, NG, NU, T, P, V, U, RLAM,
     $        NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $        N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $        AM, GP, FL, F0, PDOT0, UDOT0, G, GI, QFLAG,
     $        FPROB, FPROB, FPROB,
     $        MXJOB, LIWKL, IWKL,
     $        JFAIL)
         IF (QPRTIM) CALL MONOFF (IMON)
C
         IF (JFAIL .NE. 0) GOTO 1620
C
         IF (MDISC .EQ. 1) CALL MMULTF (NV, NL, LDF, MDF, FL, FL0, RLAM,
     $        MXJOB)
      ENDIF
C
      IF (MXJOB(51) .EQ. 0 .OR. QLAST) THEN
         MXJOB(51) = MXJOB(51) + 1
         ISAFE = 0
C
         DO 1130 J=1,K
            KC = J
            IF (QPRTIM) CALL MONON (2)
            CALL MEXXEX (J, NP, NV, NL, NG, NU, FPROB,
     $           T, P, V, U, RLAM, SCAL, H,
     $           NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $           N2UDIM, N6UDIM, NVLDIM,
     $           LDA, LDG, LDF, MDA, MDG, MDF,
     $           NJ, TAB, KD, ITOL, RTOL, ATOL, ERR, HH, WK, AA,
     $           WF, WP, WV, WA, WRL, F0, WPDOT, PDOT0, G, GI,
     $           WU, WUDOT, UDOT0, WZ, FL0, FL1, AM, GP, FL,
     $           MXJOB, IFAIL, LIWKL, IWKL, LRWKL, RWKL,
     $           IDOFLG, ND, NDP, NDV, NDU, NDA, NDL,
     $           INDP, INDV, INDU, INDA, INDL, LSAFE, ISAFE, YSAFE)
            IF (QPRTIM) CALL MONOFF (2)
            IF (IFAIL .EQ. -9) GOTO 1610
            IF (QFAIL) THEN
               IF (MNINT .GE. 3) THEN
                  WRITE (LUINT,*) ' - MEXXEX: Stability check fails',
     $                 ' (MEXX21 DO 100)'
               ENDIF
               GOTO 1120
            ENDIF
            IF (J .GT. 1 .AND. ERR .LE. TOL) GOTO 1170
 1130    CONTINUE
         GOTO 1160
      ENDIF
C ------------------------------------------------
C --- basic integration step
C ------------------------------------------------
      MXJOB(51) = MXJOB(51) + 1
      IF (MXJOB(51) .GE. NMAX) GOTO 1550
      ISAFE = 0
      KC = K - 1
      DO 1140 J=1,KC
C
         IF (QPRTIM) CALL MONON (2)
         CALL MEXXEX (J, NP, NV, NL, NG, NU, FPROB,
     $        T, P, V, U, RLAM, SCAL, H,
     $        NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $        N2UDIM, N6UDIM, NVLDIM,
     $        LDA, LDG, LDF, MDA, MDG, MDF,
     $        NJ, TAB, KD, ITOL, RTOL, ATOL, ERR, HH, WK, AA,
     $        WF, WP, WV, WA, WRL, F0, WPDOT, PDOT0, G, GI,
     $        WU, WUDOT, UDOT0, WZ, FL0, FL1, AM, GP, FL,
     $        MXJOB, IFAIL, LIWKL, IWKL, LRWKL, RWKL,
     $        IDOFLG, ND, NDP, NDV, NDU, NDA, NDL,
     $        INDP, INDV, INDU, INDA, INDL, LSAFE, ISAFE, YSAFE)
         IF (QPRTIM) CALL MONOFF (2)
         IF (IFAIL .EQ. -9) GOTO 1610
C
         IF (QFAIL) THEN
            IF (MNINT .GE. 3) THEN
               WRITE (LUINT,*) ' - MEXXEX: Stability check fails',
     $              ' (MEXX21 DO 120)'
            ENDIF
            MXJOB(53) = MXJOB(53) + 1
            GOTO 1120
         ENDIF
 1140 CONTINUE
C ------------------------------------------------
C --- convergence monitor
C ------------------------------------------------
      IF (K .EQ. 2 .OR. QRJECT) GOTO 1150
      IF (ERR .LE. TOL) GOTO 1170
      IF ( ERR/TOL .GT. DBLE(NJ(K+1)*NJ(K))) GOTO 1460
C
 1150 CONTINUE
      J = K
      IF (QPRTIM) CALL MONON (2)
      CALL MEXXEX (J, NP, NV, NL, NG, NU, FPROB,
     $     T, P, V, U, RLAM, SCAL, H,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM,
     $     LDA, LDG, LDF, MDA, MDG, MDF,
     $     NJ, TAB, KD, ITOL, RTOL, ATOL, ERR, HH, WK, AA,
     $     WF, WP, WV, WA, WRL, F0, WPDOT, PDOT0 , G, GI,
     $     WU, WUDOT, UDOT0, WZ, FL0, FL1, AM, GP, FL,
     $     MXJOB, IFAIL, LIWKL, IWKL, LRWKL, RWKL,
     $     IDOFLG, ND, NDP, NDV, NDU, NDA, NDL,
     $     INDP, INDV, INDU, INDA, INDL, LSAFE, ISAFE, YSAFE)
      IF (QPRTIM) CALL MONOFF (2)
      IF (IFAIL .EQ. -9) GOTO 1610
C
      IF (QFAIL) THEN
         IF (MNINT .GE. 3) THEN
            WRITE (LUINT,*) ' - MEXXEX: Stability check fails (130)'
         ENDIF
         MXJOB(53) = MXJOB(53) + 1
C        nfail = nfail + 1
         GOTO 1120
      ENDIF
      KC = K
      IF (ERR .LE. TOL) GOTO 1170
C
C --- hope for convergence in line k+ 1 ------------------------
C
 1160 IF (ERR/TOL .GT. DBLE(NJ(K+1))) GOTO 1460
      KC = K + 1
C
      J = KC
      IF (QPRTIM) CALL MONON (2)
      CALL MEXXEX (J, NP, NV, NL, NG, NU, FPROB,
     $     T, P, V, U, RLAM, SCAL, H,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM,
     $     LDA, LDG, LDF, MDA, MDG, MDF,
     $     NJ, TAB, KD, ITOL, RTOL, ATOL, ERR, HH, WK, AA,
     $     WF, WP, WV, WA, WRL, F0, WPDOT, PDOT0, G, GI,
     $     WU, WUDOT, UDOT0, WZ, FL0, FL1, AM, GP, FL,
     $     MXJOB, IFAIL, LIWKL, IWKL, LRWKL, RWKL,
     $     IDOFLG, ND, NDP, NDV, NDU, NDA, NDL,
     $     INDP, INDV, INDU, INDA, INDL, LSAFE, ISAFE, YSAFE)
      IF (QPRTIM) CALL MONOFF (2)
      IF (IFAIL .EQ. -9) GOTO 1610
C
      IF (QFAIL) THEN
         IF (MNINT .GE. 3) THEN
            WRITE (LUINT,*) ' - MEXXEX: Stability check fails (150)'
         ENDIF
         MXJOB(53) = MXJOB(53) + 1
         GOTO 1120
      ENDIF
      IF (ERR .GT. TOL) GOTO 1460
C ------------------------------------------------
C ---  projection step
C ------------------------------------------------
 1170 CONTINUE
      T1 = T + H
      DO 1180 I=1,NP
         P1(I) = TAB(I,1)
 1180 CONTINUE
      NS = NP
      DO 1190 I=1,N
         V1(I) = TAB(NS+I,1)
 1190 CONTINUE
      NS = NS + NV
      DO 1200 I=1,NU
         U1(I) = TAB(NS+I, 1)
 1200 CONTINUE
      NS = NS + NU
      DO 1210 I=1,NV
         A1(I) = TAB(NS+I, 1)
 1210 CONTINUE
      NS = NS + NV
      DO 1220 I=1,NL
         RL1(I) = TAB(NS+I, 1)
 1220 CONTINUE
C
      IF (QPRTIM) CALL MONON (1)
      CALL MXPRJC(NP, NV, NL, NG, NU,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $     FPROB,
     $     T1, P1, V1, U1, A1, RL1, SCAL,
     $     RTOL, ATOL, ITOL,
     $     ITMAX, NITER, EPSREQ, EPSEST,
     $     WF, WP, WV, WZ, G, GI, WPDOT, WU, AM, GP, FL,
     $     MXJOB, RWKL, IWKL, LRWKL, LIWKL, IFAIL)
      IF (QPRTIM) CALL MONOFF (1)
      IF (IFAIL .EQ. -9) GOTO 1600
      IF (IFAIL .NE. 0) GOTO 1570
C
C ---    Test for changes in position    ----------------------
C
      CP = 0.D0
      DO 1230 I =1,NP
         IF (ITOL .EQ. 0) THEN
            WT = RTOL(1)*SCAL(I) + ATOL(1)
         ELSE
            WT = RTOL(I)*SCAL(I) + ATOL(I)
         ENDIF
         CP = CP + ((P1(I) - TAB(I,1))/WT)**2
 1230 CONTINUE
      CP = SQRT(CP/DBLE(NP))
C
      IF (CP .GT. TOL*FACP) THEN
         MXJOB(53) = MXJOB(53) + 1
         PFAC = MAX (SAFE4, FAC3*TOL/CP)
         H = H*PFAC
         QRJECT = .TRUE.
         QLAST = .FALSE.
         IF (MNINT .GE. 3) THEN
            WRITE (LUINT,*) ' - PROJEC: Changes in position too large'
         ENDIF
         MXJOB(53) = MXJOB(53) + 1
C        nfail = nfail + 1
         GOTO 1120
      ENDIF
C
C ---    Test for changes in velocity    ----------------------
C
      CV = 0.D0
      DO 1240 I=1,NV
         NI = NP + I
         IF (ITOL .EQ. 0) THEN
            WT = RTOL(1)*SCAL(NI) + ATOL(1)
         ELSE
            WT = RTOL(NI)*SCAL(NI) + ATOL(NI)
         ENDIF
         CV = CV + ((V1(I) - TAB(NI,1))/WT)**2
 1240 CONTINUE
      CV = SQRT(CV/DBLE(N))
      IF (CV .GT. TOL*FACV) THEN
         MXJOB(53) = MXJOB(53) + 1
         PFAC = MAX (SAFE4, FAC3*TOL/CV)
         H = H*PFAC
         QRJECT = .TRUE.
         QLAST = .FALSE.
         IF (MNINT .GE. 3) THEN
            WRITE (LUINT,*) ' - PROJEC: Changes in velocity too large'
         ENDIF
         MXJOB(53) = MXJOB(53) + 1
C        nfail = nfail + 1
         GOTO 1120
      ENDIF
C ------------------------------------------------
C ---  step is accepted
C ------------------------------------------------
C
      IF (MXJOB(52) .EQ. 0) THEN
C 
C ---  Initiations after first step
C
         IF (IALAM .EQ. 1) THEN
            NS = NP + 2*NV + NU + NL
            DO 1250 I=1,NV
               A(I) = TAB(NS+I, 1)
 1250       CONTINUE
            NS = NS + NV
            DO 1260 I=1,NL
               RLAM(I) = TAB(NS+I, 1)
 1260       CONTINUE
         ENDIF
C
C  Solution output
C
         IF (MOUT .GE. 1) THEN
C           Initialize solution output
            INFOS(41) = 0
            INFOS(42) = MXJOB(52)
            INFOS(44) = INFOS(44) + 1
            IRTRN = 0
            CALL SOLOUT(NP, NV, NU, NL,
     $           T, P, V, U, A, RLAM, INFOS, IRTRN)
            IF (IRTRN .LT. 0) GOTO 1580
         ENDIF
C
C  Dense output
C
         IF (MDOUT .GT. 0 .OR. MGLOB .GT. 0) THEN
C           Initialize dense output
            NS = 0
            NDPH = IWOPT(51)
            NDVH = IWOPT(52)
            NDUH = IWOPT(53)
            NDAH = IWOPT(54)
            NDLH = IWOPT(55)
            DO 1270 I=1,NDPH
               YIP(NS+I) = P(INDP(I))
 1270       CONTINUE
            NS = NDPH
            DO 1280 I=1,NDVH
               YIP(NS+I) = V(INDV(I))
 1280       CONTINUE
            NS = NS + NDVH
            DO 1290 I=1,NDUH
               YIP(NS+I) = U(INDU(I))
 1290       CONTINUE
            NS = NS + NDUH
            DO 1300 I=1,NDAH
               YIP(NS+I) = A(INDA(I))
 1300       CONTINUE
            NS = NS + NDAH
            DO 1310 I=1,NDLH
               YIP(NS+I) = RLAM(INDL(I))
 1310       CONTINUE
            NS = NS + NDLH
            NDH=NS
C
            INFOD(41) = 0
            INFOD(42) = MXJOB(52)
            IRTRN = 0 
C
            IF (MGLOB .GT. 0) THEN
               CALL MXDDU1(INFOD,NDH,NDPH,NDVH,NDUH,NDAH,NDLH,KDENS,
     $              ILAM, IRHO, T, T, H, DENS, YIP,
     $              T0, TFIN, DOUTS,
     $              NS, IWOPT(51), IWOPT(52), IWOPT(53), IWOPT(54),
     $              IWOPT(55),
     $              INDP, INDV, INDU, INDA, INDL,IRTRN)
            ENDIF
C
            IF (MDOUT .GT. 0) THEN
               CALL MXDOUT(INFOD, NDH, NDPH, NDVH, NDUH, NDAH, NDLH,
     $              KDENS, ILAM, IRHO, T, T, H, DENS, YIP,
     $              T0, TFIN, DOUTS, DENOUT,
     $              NS, IWOPT(51), IWOPT(52), IWOPT(53), IWOPT(54),
     $              IWOPT(55),
     $              INDP, INDV, INDU, INDA, INDL,QPRTIM,IRTRN)
               IF (IRTRN .LT. 0) GOTO 1580
            ENDIF
C
         ENDIF
C
C  root finding
C
         IF (NSWIF .GT. 0) THEN
            ICTYPE=0
            IRTRN=0
            CALL MXRTSC (MXJOB(52), ICTYPE, NP, NV, NL, NU, NSWIF,
     $                   NPDIM, NVDIM, NLDIM, NUDIM,
     $                   T, P, V, U, A, RLAM, GVAL0, GVAL1,
     $                   IRES0, IRES1, NSC, MXJOB(37), FSWIT,
     $                   MXJOB(17), MXJOB(18), IRTRN)
C
            CALL MXRTFD (MXJOB(52), ND, NDP, NDV, NDU, NDA, NDL,
     $                   KDENS, ILAM, IRHO, T, T, H, DENS, YIP,
     $                   NSWIF, MXJOB(37), GVALW1, GVAL0, GVAL1,
     $                   GVALW2, GVALW3, GVALW4,
     $                   IRES0, IRES1, IRESL, IRESR,
     $                   TTOL, TZERO, ISWIT, FSWIT, NSC, IRTRN,
     $                   NP, NV, NL, NG, NU,
     $                   NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $                   N2UDIM, N6UDIM, NVLDIM,
     $                   LDA, LDG, LDF, MDA, MDG, MDF,
     $                   FPROB,
     $                      P, V, U, A, RLAM, SCAL,
     $                   RTOL, ATOL, ITOL,
     $                   ITMAX, NITER, EPSREQ, EPSEST,
     $                   WF, WP, WV, WZ, G, GI, WPDOT, WU, AM, GP, FL,
     $                   MXJOB, RWKL, IWKL, LRWKL, LIWKL, IFAIL)
         ENDIF
C
C ---  End initiations after first step
C
      ENDIF
C
      MXJOB(52) = MXJOB(52) + 1
C
C     Function evaluations: f, pdot, udot, df/dl
C
      QFLAG(1) = .FALSE.
      QFLAG(2) = .FALSE.
      QFLAG(3) = .FALSE.
      QFLAG(4) = .TRUE.
      QFLAG(5) = .TRUE.
      QFLAG(6) = .TRUE.
      QFLAG(7) = .FALSE.
      QFLAG(8) = .FALSE.
      QFLAG(9) = (MDISC .EQ. 1) .AND. (.NOT. QLAST)
C
      JFAIL = 0
C
      IF (QFLAG(9)) THEN
         IMON = 10
      ELSE
         IMON=5
      ENDIF
      IF (QPRTIM) CALL MONON (IMON)
      CALL MXPROB (NP, NV, NL, NG, NU, T1, P1, V1, U1, RL1,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $     AM, GP, FL, F1, PDOT1, UDOT1, G, GI,
     $     QFLAG, FPROB, FPROB, FPROB,
     $     MXJOB, LIWKL, IWKL, JFAIL)
      IF (QPRTIM) CALL MONOFF (IMON)
C
      IF (JFAIL .NE. 0) GOTO 1620
C
      IF (MDISC .EQ. 1) CALL MMULTF (NV, NL, LDF, MDF,
     $                             FL, FL0, RLAM, MXJOB)
C
C --- root finding: check for sign change -----
C
      QSWIT = .FALSE.
      IF (NSWIF .GT. 0) THEN
         ICTYPE=1
         IRTRN=0
         CALL MXRTSC (MXJOB(52), ICTYPE, NP, NV, NL, NU, NSWIF,
     $                NPDIM, NVDIM, NLDIM, NUDIM,
     $                T1, P1, V1, U1, A1, RL1, GVAL0, GVAL1, 
     $                IRES0, IRES1, NSC, MXJOB(37), FSWIT,
     $                MXJOB(17), MXJOB(18), IRTRN)
      ENDIF
C
C --- dense output: Approximate derivatives, extrapolate, set
C                   coefficients
C                   (Called, if dense output option is on
C                               or
C                               root finding is on and sign changed
C
      IF (MDOUT .GT. 0 .OR. MGLOB .GT. 0 .OR.
     $     (NSWIF .GT. 0 .AND. NSC .GT. 0)) THEN
         IF (QPRTIM) CALL MONON(3)
C
         KEXMAX = KD
         KDOMAX = KD
         IDOCNT = IDOCNT + 1
C 
         IDOFLH = 1
         NDPH = IWOPT(51)
         NDVH = IWOPT(52)
         NDUH = IWOPT(53)
         NDAH = IWOPT(54)
         NDLH = IWOPT(55)
         NDH = NDPH+NDVH+NDUH+NDAH+NDLH
C
         IF (MDIND.EQ.0 .OR. NSC.GT.0) THEN
            IDOFLH= 2
            NDPH = NDP
            NDVH = NDV    
            NDUH = NDU    
            NDAH = NDA    
            NDLH = NDL     
            NDH = NDP+NDV+NDU+NDA+NDL
         ENDIF
C 
         IF (MDIND.EQ.1 .AND. NSC.EQ.0 .AND. NSWIF.GT.0) THEN
            NS = 0
            NSH = 0
            DO 1320 I=1,NDPH
            DO 1320 L=1,ISAFE
               YSAFE(NSH+I,L) = YSAFE(NS+INDP(I),L)
 1320       CONTINUE
            NS = NDP
            NSH = NDPH
            DO 1330 I=1,NDVH
            DO 1330 L=1,ISAFE
               YSAFE(NSH+I,L) = YSAFE(NS+INDV(I),L)
 1330       CONTINUE
            NS = NS + NDV
            NSH = NSH + NDVH
            DO 1340 I=1,NDUH
            DO 1340 L=1,ISAFE
               YSAFE(NSH+I,L) = YSAFE(NS+INDU(I),L)
 1340       CONTINUE
            NS = NS + NDU
            NSH = NSH + NDUH
            DO 1350 I=1,NDAH
            DO 1350 L=1,ISAFE
               YSAFE(NSH+I,L) = YSAFE(NS+INDA(I),L)
 1350       CONTINUE
            NS = NS + NDA
            NSH = NSH + NDAH
            DO 1360 I=1,NDLH
            DO 1360 L=1,ISAFE
               YSAFE(NSH+I,L) = YSAFE(NS+INDL(I),L)
 1360       CONTINUE
         ENDIF
         CALL MXDENS (NP, NV, NU, NL, T, P, V, U, A, RLAM,
     $                   PDOT0, A, UDOT0, PDOT1, A1, UDOT1,
     $                   H, KC, KDOMAX, KEXMAX, TAB, NJ,
     $                   P1, V1, U1, A1, RL1, ILAM, IRHO,
     $                   IDOFLH, N6UDIM, ISAFE, KDENS,
     $                   ND, NDH, NDPH, NDVH, NDUH, NDAH, NDLH,
     $                   ITOL, NTODIM, RTOL, ATOL, N2UDIM, SCAL, YIP,
     $                   INDP, INDV, INDU, INDA, INDL, YSAFE, DENS,
     $                   MXJOB, IDOCNT)
C
         IF (QPRTIM) CALL MONOFF(3)
C
      ENDIF
C 
C
C --- root finding: determine root with minimal distance to TOLD
C
      IF (NSWIF .GT. 0 .AND. NSC.GT.0) THEN  
         CALL MXRTFD (MXJOB(52), ND, NDP, NDV, NDU, NDA, NDL,
     $                KDENS, ILAM, IRHO, T1, T, H, DENS, YIP,
     $                NSWIF, MXJOB(37), GVALW1, GVAL0, GVAL1,
     $                GVALW2, GVALW3, GVALW4,
     $                IRES0,IRES1,IRESL,IRESR,
     $                TTOL, TZERO, ISWIT, FSWIT, NSC, IRTRN,
     $                NP, NV, NL, NG, NU,  
     $                NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $                N2UDIM, N6UDIM, NVLDIM,
     $                LDA, LDG, LDF, MDA, MDG, MDF,
     $                FPROB,
     $                   P1, V1, U1, A1, RL1, SCAL,
     $                RTOL, ATOL, ITOL,    
     $                ITMAX, NITER, EPSREQ, EPSEST,
     $                WF, WP, WV, WZ, G, GI, WPDOT, WU, AM, GP, FL,
     $                MXJOB, RWKL, IWKL, LRWKL, LIWKL, IFAIL)
         IF (IRTRN .LT. 0) GOTO 1590
         QSWIT = IRTRN .GT. 0
      ENDIF
C
      IF (MOUT .GE. 1) THEN
         INFOS(41) = 1
         IF(QSWIT .AND. MXJOB(35).EQ.2) INFOS(41) = 3
         IF(QSWIT .AND. MXJOB(35).EQ.1) INFOS(41) = 4
         INFOS(42) = MXJOB(52)
         IRTRN = 0
         IF(QSWIT) THEN
            INFOS(44) = INFOS(44) + 1
            CALL SOLOUT(NP, NV, NU, NL, TZERO, YIP, YIP(MIN(ND,NP+1)),
     $           YIP(MIN(ND,NP+NV+1)), YIP(MIN(ND,NP+NV+NU+1)),
     $           YIP(MIN(ND,NP+NV+NU+NV+1)), INFOS, IRTRN)
            IF (IRTRN .LT. 0) GOTO 1580
         ENDIF
         IF(INFOS(41).NE.4) THEN
         INFOS(41) = 1
         H1 = TFIN - T1
         INFOS(44) = INFOS(44) + 1
         IF (H1 .LE. UROUND) INFOS(41) = 2
         CALL SOLOUT(NP, NV, NU, NL, T1, P1, V1, U1, A1, RL1,
     $        INFOS, IRTRN)
         IF (IRTRN .LT. 0) GOTO 1580
         ENDIF
      ENDIF
C
      IF (MDOUT .GT. 0 .OR. MGLOB .GT. 0) THEN
         INFOD(41) = 1
         IF(QSWIT .AND. MXJOB(35).EQ.2) INFOD(41) = 3
         IF(QSWIT .AND. MXJOB(35).EQ.1) INFOD(41) = 4
         H1 = TFIN - T1
         IF (H1 .LE. UROUND) INFOD(41) = 2
         INFOD(42) = MXJOB(52)
         IF(QSWIT) INFOD(43) = 1
         IRTRN = 0
C
         IF (INFOD(41) .EQ. 4) THEN
            THH = TZERO
         ELSE
            THH = T1
         ENDIF
         IF (MGLOB .GT. 0) THEN
            CALL MXDDU1(INFOD, NDH, NDPH, NDVH, NDUH, NDAH, NDLH,
     $           KDENS, ILAM, IRHO, THH, T, H, DENS, YSAFE,
     $           T0, TFIN, DOUTS,
     $           NDH, IWOPT(51), IWOPT(52), IWOPT(53), IWOPT(54),
     $           IWOPT(55),
     $           INDP, INDV, INDU, INDA, INDL, IRTRN)
         ENDIF
         IF (MDOUT .GT. 0) THEN
            CALL MXDOUT(INFOD,NDH,NDPH,NDVH,NDUH,NDAH,NDLH,KDENS,
     $           ILAM, IRHO, THH, T, H, DENS, YSAFE,
     $           T0, TFIN, DOUTS, DENOUT,
     $         NDH,IWOPT(51),IWOPT(52),IWOPT(53),IWOPT(54),IWOPT(55),
     $           INDP, INDV, INDU, INDA, INDL, QPRTIM, IRTRN)
            IF (IRTRN .LT. 0) GOTO 1580
C
         ENDIF
C
      ENDIF
C
C ---   update t, p, v, a, u, rlam, f0, pdot0, udot0   -----------------
C
      T = T1
      DO 1370 I=1,NP
         P(I) = P1(I)
 1370 CONTINUE
      DO 1380 I=1,N
         V(I) = V1(I)
 1380 CONTINUE
      DO 1390 I=1,NU
         U(I) = U1(I)
 1390 CONTINUE
      DO 1400 I=1,N
         A(I) = A1(I)
 1400 CONTINUE
      DO 1410 I=1,NL
         RLAM(I) = RL1(I)
 1410 CONTINUE
C
C ---   exit if root was found and user wants mexx to stop  -----------
C
      IF(QSWIT .AND. MXJOB(35).EQ.1) GOTO 1480     
C
C  update  f0, pdot0, udot0
C
      DO 1420 I=1,NV
         F0(I) = F1(I)
 1420 CONTINUE
      DO 1430 I=1,NP
         PDOT0(I) = PDOT1(I)
 1430 CONTINUE
      DO 1440 I=1,NU
         UDOT0(I) = UDOT1(I)
 1440 CONTINUE
C
C --- compute optimal order   ---------------------------------
C
      IF (KC .EQ. 2) THEN
         KOPT = 3
         IF (QRJECT) KOPT = 2
         GOTO 1450
      ENDIF
      IF (KC .LE. K) THEN
         KOPT = KC
         IF (WK(KC-1) .LT. WK(KC)*FAC3) KOPT = KC - 1
         IF (WK(KC) .LT. WK(KC-1)*FAC5) KOPT = MIN (KC+1, KD-1)
      ELSE
         KOPT = KC - 1
         IF (KC .GT. 3 .AND. WK(KC-2) .LT. WK(KC-1)*FAC3) KOPT = KC -2
         IF (WK(KC) .LT. WK(KOPT)*FAC5) KOPT = MIN (KC, KD-1)
      ENDIF
C
C --- after a rejected step   ---------------------------------
C
 1450 IF (QRJECT) THEN
         K = MIN (KOPT, KC)
         H = MIN (H, HH(K))
         QRJECT = .FALSE.
         GOTO 1120
      ENDIF
C
C --- compute stepsize for next step  -------------------------
C
      IF (KOPT .LE. KC) THEN
         H = HH(KOPT)
      ELSE
         H = HH(KC) *AA(KC+1)/AA(KC)
      ENDIF
      K = KOPT
      GOTO 1120
C ------------------------------------------------
C --- step is rejected
C ------------------------------------------------
 1460 K = MIN (K, KC)
      IF (K .GT. 2 .AND. WK(K-1) .LT. WK(K)*FAC3) K = K - 1
      MXJOB(53) = MXJOB(53) + 1
      H = HH(K)
      QLAST = .FALSE.
      QRJECT = .TRUE.
      IF (MNINT .GE. 3) THEN
         WRITE (LUINT,*) ' - MEXX21: Discretization error too large'
      ENDIF
      GOTO 1120
C ------------------------------------------------
C --- solution exit
C ------------------------------------------------
C
 1470 CONTINUE
      GOTO 1630
C ------------------------------------------------
C --- special exit (root found)
C ------------------------------------------------
C
 1480 CONTINUE
      IERR=1
      T1=TZERO
      T=TZERO
      NS = 0
      DO 1490 I=1,NDP
         P(I) = YIP(NS+I)
 1490 CONTINUE
      NS = NDP
      DO 1500 I=1,NDV
         V(I) = YIP(NS+I)
 1500 CONTINUE
      NS = NS + NDV
      DO 1510 I=1,NDU
         U(I) = YIP(NS+I)
 1510 CONTINUE
      NS = NS + NDU
      DO 1520 I=1,NDA
         A(I) = YIP(NS+I)
 1520 CONTINUE
      NS = NS + NDA
      DO 1530 I=1,NDL
         RLAM(I) = YIP(NS+I)
 1530 CONTINUE

      GOTO 1630
C ------------------------------------------------
C --- fail exit
C ------------------------------------------------
C
 1540 CONTINUE
      IERR = -9
      GOTO 1560
 1550 CONTINUE
      IERR = -10
 1560 CONTINUE
      IF (MNGEN .GE. 2) THEN
         WRITE (LUGEN,*) ' - MEXX21: Integration fails:'
         WRITE (LUGEN,9010) T, H
 9010    FORMAT ('   at T = ', D15.7, '  for H = ', D15.7)
      ENDIF
      GOTO 1630
C
 1570 CONTINUE
      IERR = -11
      IF (MNGEN .GE. 2) THEN
         WRITE (LUGEN,*) ' - MEXX21: Projection step fails:'
         WRITE (LUGEN, 9020) '   with IFAIL = ', IFAIL
         WRITE (LUGEN, 9030) T
 9020    FORMAT (A, I5)
 9030    FORMAT ('   at T = ', D15.7)
         IF (IFAIL .EQ. -9) WRITE (LUGEN, 9020)
     $        '   PROJEC: error return from FPROB', MXJOB(73)
         IF (IFAIL .EQ. -15) THEN
            WRITE (LUGEN, 9020)
     $           '   PROJEC: error return from ADEC', MXJOB(73)
            WRITE (LUGEN, 9020)
     $           '   within routine                ', MXJOB(74)
         ENDIF
      ENDIF
      GOTO 1630
C
 1580 CONTINUE
      MXJOB(73) = IRTRN
      IERR = -19
      IF (MNGEN .GE. 2) THEN
         WRITE (LUGEN,*) ' - MEXX21: Integration interrupted:'
         WRITE (LUGEN,9030) T
      ENDIF
      GOTO 1630
C        
 1590 CONTINUE
      MXJOB(73) = IRTRN
      IERR = -20
      IF (MNGEN .GE. 2) THEN
         WRITE (LUGEN,*) ' - MEXX21: Root finder fail:'
         WRITE (LUGEN,*) '           return code:',IRTRN 
         WRITE (LUGEN,9030) T
      ENDIF
      GOTO 1630
C
 1600 CONTINUE
      IERR = -18
      IF (MNGEN .GE. 2) THEN
         WRITE (LUGEN,*) ' - MEXX21: Integration interrupted by user:'
         WRITE (LUGEN,*) ' - MEXX21: In FPROB called from PROJEC'
         WRITE (LUGEN,9030) T
      ENDIF
      GOTO 1630
C
 1610 CONTINUE
      IERR = -17
      IF (MNGEN .GE. 2) THEN
         WRITE (LUGEN,*) ' - MEXX21: Integration interrupted by user:'
         WRITE (LUGEN,*) ' - MEXX21: In FPROB called from MEXXEX'
         WRITE (LUGEN,9030) T
      ENDIF
      GOTO 1630
C
C
C ---    Error return from user function ------------------------------
C
 1620  MXJOB(73) = JFAIL
       IFAIL = -9
       RETURN
C
 1630 CONTINUE
C
      RETURN
C
C  end of subroutine MEXX21
C
      END
C
      SUBROUTINE MXPRJC (NP, NV, NL, NG, NU,
     $      NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $      N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $      FPROB,
     $      T, P, V, U, A, RLAM, SCAL,
     $      RTOL, ATOL, ITOL,
     $      ITMAX, NITER, EPSREQ, EPSEST,
     $      F, PP, VV, X, G, GI, PDOT, UDOT, AM, GP, FL,
     $      MXJOB, RWKL, IWKL, LRWKL, LIWKL, IFAIL)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER NP, NV, NL, NG, NU, NPDIM, NVDIM, NLDIM, NGDIM, NUDIM,
     $     NTODIM, N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF
C
      EXTERNAL FPROB
C
      DOUBLE PRECISION T, P(NPDIM), V(NVDIM), U(NUDIM), A(NVDIM),
     $     RLAM(NLDIM), SCAL(N2UDIM), RTOL(NTODIM), ATOL(NTODIM)
      INTEGER ITOL, ITMAX, NITER
      DOUBLE PRECISION EPSREQ, EPSEST
      DOUBLE PRECISION F(NVDIM), PP(NPDIM), VV(NVDIM), X(NVLDIM),
     $     G(NGDIM), GI(NLDIM), PDOT(NPDIM), UDOT(NUDIM), AM(LDA,MDA),
     $     GP(LDG,MDG), FL(LDF,MDF)
      INTEGER MXJOB(150)
      DOUBLE PRECISION RWKL(LRWKL)
      INTEGER IWKL(LIWKL), LRWKL, LIWKL, IFAIL
C
C ---------------------------------------------------------------------
C
C     This subroutine computes a projection of the vector (p,v)
C     onto the manifold defined by
C
C           g(t,p) = 0
C           G(t,p)*v+gI(t,p) = 0
C
C ---------------------------------------------------------------------
C
C* File              mxprjc.f
C* Version           1.1
C* Latest Change     96/03/06 (1.8)
C* Modification history
C    1.1           Correct error number 1.0-1 in initialization of PP
C                   (loop until NP instead of NV).  
C      
C ---------------------------------------------------------------------
C
C     Input: NP, NV, NL, NG, NU, NPDIM, NVDIM, NLDIM, NGDIM, NUDIM,
C            NTODIM, N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG,
C            MDF, 
C            FPROB, T, U, SCAL, RTOL, ATOL, ITOL,
C            ITMAX, EPSREQ, MXJOB,
C
C     InOut: P, V, A, RLAM
C
C     Out:   NITER, EPSEST, IFAIL
C
C     Work:  F, PP, VV, X, G, GI, PDOT, UDOT, AM, GP, FL, RWKL, IWKL,
C            LRWKL, LIWKL 
C
C ---------------------------------------------------------------------
C
C
      LOGICAL QFLAG(9)
      LOGICAL QPRTIM
C
      SAVE
C
      IUPDTE=0
      JFAIL = 0
      IFDEC = 0
      IFAIL = 0
      QPRTIM = MXJOB(19) .GE. 1
      MNINT = MXJOB(13)
      LUINT = MXJOB(14)
      EPS = EPSREQ
      N = NV
C
C ---    Set up the matrix and its decomposition   --------------------
C
C
C  function evaluation: am, dgdp*T, g
      QFLAG(1) = .TRUE.
      QFLAG(2) = .FALSE.
      QFLAG(3) = .TRUE.
      QFLAG(4) = .FALSE.
      QFLAG(5) = .FALSE.
      QFLAG(6) = .FALSE.
      QFLAG(7) = .TRUE.
      QFLAG(8) = .FALSE.
      QFLAG(9) = .FALSE.
C
      IF (QPRTIM) CALL MONON (7)
      CALL MXPROB (NP, NV, NL, NG, NU, T, P, V, U, RLAM,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM, N2UDIM,
     $     N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $     AM, GP, FL, F, PDOT, UDOT, G, GI, QFLAG,
     $     FPROB, FPROB, FPROB,
     $     MXJOB, LIWKL, IWKL, JFAIL)
      IF (QPRTIM) CALL MONOFF (7)
      IF (JFAIL .NE. 0) GOTO 1140
C
C  matrix decomposition
      IF (QPRTIM) CALL MONON (13)
      CALL ADEC (NV, NG, LDA, LDG, LDF, AM, GP, FL, MXJOB, NVLDIM, RWKL,
     $     IWKL, IFDEC, IUPDTE)
      IF (QPRTIM) CALL MONOFF (13)
      IF (IFDEC .NE. 0) GOTO 1150
C
C ---    Prepare for iteration    --------------------------------
C
      DEFOLD = 1.D20
      DO 1000 I=1,NP
         PP(I) = P(I)
 1000 CONTINUE
      DO 1010 I=1,N
         VV(I) = 0.D0
 1010 CONTINUE
C
      K = 1
C
C ---    Compute defect in  g    --------------------------------------
C
      DO 1020 I=1,NV
         X(I) = 0.D0
 1020 CONTINUE
      DO 1030 I=1,NG
         X(NV+I) = +G(I)
 1030 CONTINUE
C
 1040 CONTINUE
C
      IF (QPRTIM) CALL MONON (12)
      CALL ASOL (NV, NG, X, MXJOB, NVLDIM, RWKL, IWKL)
      IF (QPRTIM) CALL MONOFF (12)
C
      DO 1050 I=1,NV
         VV(I) = VV(I) - X(I)
 1050 CONTINUE
C
C     Evaluate pdot
C
      QFLAG(1) = .FALSE.
      QFLAG(2) = .FALSE.
      QFLAG(3) = .FALSE.
      QFLAG(4) = .FALSE.
      QFLAG(5) = .TRUE.
      QFLAG(6) = .FALSE.
      QFLAG(7) = .FALSE.
      QFLAG(8) = .FALSE.
      QFLAG(9) = .FALSE.
C
      IF (QPRTIM) CALL MONON (8)
      CALL MXPROB (NP, NV, NL, NG, NU, T, P, X, U, RLAM,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $     AM, GP, FL, F, PDOT, UDOT, G, GI, QFLAG,
     $     FPROB, FPROB, FPROB,
     $     MXJOB, LIWKL, IWKL, JFAIL)
      IF (QPRTIM) CALL MONOFF (8)
      IF (JFAIL .NE. 0) GOTO 1140
C
      DO 1060 I=1,NP
         PP(I) = PP(I) - PDOT(I)
 1060 CONTINUE
C
      DEF = 0.D0
      DO 1070 I=1,NP
          IF (ITOL .EQ. 0) THEN
             WT = RTOL(1)*SCAL(I) + ATOL(1)
          ELSE
             WT = RTOL(I)*SCAL(I) + ATOL(I)
          ENDIF
          DEF = DEF + (PDOT(I)/WT)**2
 1070 CONTINUE
      DEF = SQRT(DEF/DBLE(NP))
C
C ---    Test for convergence    --------------------------------------
C
      IF (MNINT .GE. 4) THEN
         ERRPRI = DEF/EPS
         WRITE (LUINT,*) ' - PROJEC: Est. normalized error:', K, ERRPRI
      ENDIF
      IF (DEF .LT. EPS) GOTO 1090
      IF (DEF .GE. DEFOLD*2.D0 .OR. K .GE. ITMAX) GOTO 1130
C
C ---    Next iteration    --------------------------------------------
C
      DEFOLD = DEF
      K = K + 1
C
C  right-hand side
C
      CALL MMULT (NV, NG, LDA, AM, VV, X, MXJOB)
C
C  function evaluation: g
C
      QFLAG(1) = .FALSE.
      QFLAG(2) = .FALSE.
      QFLAG(3) = .FALSE.
      QFLAG(4) = .FALSE.
      QFLAG(5) = .FALSE.
      QFLAG(6) = .FALSE.
      QFLAG(7) = .TRUE.
      QFLAG(8) = .FALSE.
      QFLAG(9) = .FALSE.
C
      IF (QPRTIM) CALL MONON (9)
      CALL MXPROB (NP, NV, NL, NG, NU, T, PP, V, U, RLAM,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $     AM, GP, FL, F, PDOT, UDOT, G, GI, QFLAG,
     $     FPROB, FPROB, FPROB,
     $     MXJOB, LIWKL, IWKL, JFAIL)
      IF (QPRTIM) CALL MONOFF (9)
      IF (JFAIL .NE. 0) GOTO 1140
C
      DO 1080 I=1,NG
         X(NV+I) = +G(I)
 1080 CONTINUE
C
      GOTO 1040
C
C ---    Projection of velocity     -----------------------------------
C
 1090 CONTINUE
      EPSEST = DEF
      NITER = K
      DO 1100 I=1,NP
         P(I) = PP(I)
 1100 CONTINUE
C
C  function evaluation: am, gp, gi
      QFLAG(1) = .TRUE.
      QFLAG(2) = .TRUE.
      QFLAG(3) = .FALSE.
      QFLAG(4) = .FALSE.
      QFLAG(5) = .FALSE.
      QFLAG(6) = .FALSE.
      QFLAG(7) = .FALSE.
      QFLAG(8) = .TRUE.
      QFLAG(9) = .FALSE.
C
      IF (QPRTIM) CALL MONON (6)
      CALL MXPROB (NP, NV, NL, NG, NU, T, P, V, U, RLAM,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $     AM, GP, FL, F, PDOT, UDOT, G, GI, QFLAG,
     $     FPROB, FPROB, FPROB,
     $     MXJOB, LIWKL, IWKL, JFAIL)
      IF (QPRTIM) CALL MONOFF (6)
      IF (JFAIL .NE. 0) GOTO 1140
C
      IF (QPRTIM) CALL MONON (13)
      CALL ADEC (NV, NL, LDA, LDG, LDF, AM, GP, FL, MXJOB, NVLDIM, RWKL,
     $     IWKL, IFDEC, IUPDTE)
      IF (QPRTIM) CALL MONOFF (13)
      IF (IFDEC .NE. 0) GOTO 1150
C
      CALL MMULT (NV, NL, LDA, AM, V, X, MXJOB)
C
      DO 1110 I=1,NL
         X(NV+I) = -GI(I)
 1110 CONTINUE
C
      IF (QPRTIM) CALL MONON (12)
      CALL ASOL (NV, NL, X, MXJOB, NVLDIM, RWKL, IWKL)
      IF (QPRTIM) CALL MONOFF (12)
C
      DO 1120 I=1, NV
         V(I) = X(I)
 1120 CONTINUE
C
      RETURN
C
C ---    No convergence in projection of p    -------------------------
C
 1130 CONTINUE
      EPSEST = DEF
      NITER = K
      IFAIL = -1
      RETURN
C
C ---    Error return from user function ------------------------------
C
 1140 MXJOB(73) = JFAIL
      IFAIL = -9
      RETURN
C
C ---    Error return from ADEC ---------------------------------------
C
 1150 MXJOB(73) = IFDEC
      IFAIL = -15
      RETURN
C
C  end of subroutine projec
C
      END
C
************************************************************************

      SUBROUTINE MXPROB (NP, NV, NL, NG, NU, T, P, V, U, RLAM,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $     AM, GP, FL, F, PDOT, UDOT, G, GI,
     $     QFLAG, FPROB1, FPROB2, FPROB3,
     $     MXJOB, LIWKL, IWKL, IFAIL)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C  input parameters:
      INTEGER NP, NV, NL, NG, NU
      DOUBLE PRECISION T, P(NPDIM), V(NVDIM), U(NUDIM), RLAM(NLDIM)
      INTEGER NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM, N2UDIM, N6UDIM,
     $     NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF, MXJOB(150)
      LOGICAL QFLAG(9)
C  output parameters:
      DOUBLE PRECISION AM(LDA,MDA), GP(LDG,MDG), FL(LDF,MDF), F(NVDIM),
     $     PDOT(NPDIM), UDOT(NUDIM),  G(NGDIM), GI(NLDIM)
      INTEGER IFAIL
C  worksapece parameters:
      INTEGER LIWKL, IWKL(LIWKL)
C
      EXTERNAL FPROB1, FPROB2, FPROB3
C
C
C* File              mxprob.f
C* Version           1.1
C* Latest Change     96/03/06 (1.6)
C* Modification history
C    1.1             Check for unnecessary flags.
C      
      SAVE
C
C  check for unnecessary flags
      IF (NG .EQ. 0) THEN
         QFLAG(3) = .FALSE.
         QFLAG(7) = .FALSE.
         IF (NL .EQ. 0) QFLAG(2) = .FALSE.
      ENDIF
      IF (NL .EQ. 0) THEN
         QFLAG(8) = .FALSE.
         QFLAG(9) = .FALSE.
      ENDIF
      IF (NU .EQ. 0) QFLAG(6) = .FALSE.
C
      JOB = MXJOB(2)
      MMODE = MXJOB(3)
C
      IF (JOB .LT. 10 .AND. MMODE .EQ. 0) THEN
         CALL FPROB1(NP, NV, NL, NG, NU, LDG,
     $               T, P, V, U, RLAM,
     $               AM, GP, F, PDOT, UDOT, G, GI, FL,
     $               QFLAG, IFAIL)
      ENDIF
C
      IF (JOB .LT. 10 .AND. MMODE .EQ. 1) THEN
         CALL FPROB2(NP, NV, NL, NG, NU, LDG,
     $               T, P, V, U, RLAM,
     $               AM, GP, F, PDOT, UDOT, G, GI, FL,
     $               QFLAG, MXJOB(4), MXJOB(5), IFAIL)
      ENDIF
C
      IF (JOB .EQ. 10 .AND. MMODE .EQ. 1) THEN
         MXJOB(102) = 0
         MXJOB(105) = 0
         IF (MXJOB(103) .EQ. -1) MXJOB(102) = 1
         IF (MXJOB(106) .EQ. -1) MXJOB(105) = 1
C
C        INTEGER work space for IROWG(NZGMAX), JCOLG(NZGMAX), 
C                               IROWF(NZFMAX), JCOLF(NZFMAX)
C        (compare linear algebra routine MXLIN)
C
         JIRG = 1
         JJCG = JIRG + MXJOB(7)
         JIRF = JJCG + MXJOB(7)
         JJCF = JIRF + MXJOB(8)
C
         CALL FPROB3(NP, NV, NL, NG, NU, LDG,
     $        T, P, V, U, RLAM, AM, GP, F,
     $        PDOT, UDOT, G, GI, FL, QFLAG,
     $        MXJOB(4), MXJOB(5), MXJOB(7), MXJOB(8),
     $        IWKL(JIRG), IWKL(JJCG), MXJOB(101), MXJOB(102),
     $        IWKL(JIRF), IWKL(JJCF), MXJOB(104), MXJOB(105),
     $        IFAIL)
         IF (MXJOB(103) .EQ. -1) MXJOB(102) = 1
         MXJOB(103) = 0
         IF (MXJOB(106) .EQ. -1) MXJOB(105) = 1
         MXJOB(106) = 0
      ENDIF
C
      MXJOB(54) = MXJOB(54) + 1
      IF (QFLAG(4)) MXJOB(55) = MXJOB(55) + 1
      IF (QFLAG(7)) MXJOB(56) = MXJOB(56) + 1
C
      RETURN
      END
C
      SUBROUTINE MEXXEX (J, NP, NV, NL, NG, NU, FPROB,
     $     T0, P0, V0, U0, RL0, SCAL, H,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM,
     $     LDA, LDG, LDF, MDA, MDG, MDF,
     $     NJ, TAB, KD, ITOL, RTOL, ATOL, ERR, HH, WK, AA,
     $     F, P, V, A, RLAM, F0, PDOT, PDOT0, G, GI,
     $     U, UDOT, UDOT0, Z, FL0, FL1, AM, GP, FL,
     $     MXJOB, IFAIL, LIWKL, IWKL, LRWKL, RWKL,
     $     IDOFLG, ND, NDP, NDV, NDU, NDA, NDL,
     $     INDP, INDV, INDU, INDA, INDL, LSAFE, ISAFE, YSAFE)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER J, NP, NV, NL, NG, NU
      EXTERNAL FPROB
      DOUBLE PRECISION T0, P0(NPDIM), V0(NVDIM), U0(NUDIM), RL0(NLDIM),
     $     SCAL(N2UDIM), H
      INTEGER NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM, N2UDIM, N6UDIM,
     $     NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF, NJ(KD)
      DOUBLE PRECISION TAB(N6UDIM,KD)
      INTEGER KD, ITOL
      DOUBLE PRECISION RTOL(NTODIM), ATOL(NTODIM), ERR, HH(KD), WK(KD),
     $     AA(KD), F(NVDIM), P(NPDIM), V(NVDIM), A(NVDIM), RLAM(NLDIM),
     $     F0(NVDIM), PDOT(NPDIM), PDOT0(NPDIM), G(NGDIM), GI(NLDIM),
     $     U(NUDIM), UDOT(NUDIM), UDOT0(NUDIM), Z(NVLDIM), FL0(LDF),
     $     FL1(LDF), AM(LDA,MDA), GP(LDG,MDG), FL(LDF,MDF)
      INTEGER MXJOB(150), IFAIL, LIWKL, IWKL(LIWKL), LRWKL
      DOUBLE PRECISION RWKL(LRWKL)
C -- dense output formula
      INTEGER IDOFLG, ND, NDP, NDV, NDU, NDA, NDL, INDP(NDP), INDV(NDV),
     $     INDU(NDU), INDA(NDA), INDL(NDL), LSAFE, ISAFE
      DOUBLE PRECISION YSAFE(ND,LSAFE)
C
C
C ----------------------------------------------------------------------
C
C       M E X X 2 1    (index 2,  h^1-extrapolation)
C
C     This subroutine computes the j-th line of the extrapolation
C     tableau and provides estimates of the error and the optimal step
C     size.
C
C ----------------------------------------------------------------------
C
C* File              mexxex.f
C* Version           1.0
C* Latest Change     96/03/04 (1.4)
C
C ----------------------------------------------------------------------
C
C
      LOGICAL QRJECT, QFAIL, QINIT, QFLAG(9)
      LOGICAL QPRTIM
C
      COMMON /ERRCOM/ TOLD4, UROUND, FAC, FAC1, FAC2,
     $      SAFE2, SAFE3, HMAX, ERROLD
      COMMON /LOGCOM/ QRJECT, QFAIL, QINIT
C
      SAVE
C
C
C ---    Preparations   -----------------------------------------------
C
      IUPDTE = MXJOB(1)
      IFAIL = 0
      MDISC = MXJOB(1)
      QPRTIM = MXJOB(19) .GE. 1
      MNINT = MXJOB(13)
      LUINT = MXJOB(14)
      IALAM = MXJOB(22)
      NN = NP + NV
      NX = NP + 2*NV + NU + NL
      NJJ = NJ(J)
      HJ = H/NJJ
C
      T = T0
      DO 1000 I=1,NP
         P(I) = P0(I)
 1000 CONTINUE
      DO 1010 I=1,NV
         V(I) = V0(I)
 1010 CONTINUE
      DO 1020 I=1,NL
         RLAM(I) = RL0(I)
 1020 CONTINUE
      DO 1030 I=1,NU
         U(I) = U0(I)
 1030 CONTINUE
      DO 1040 I=1,NV
         F(I) = F0(I)
 1040 CONTINUE
      DO 1050 I=1,NP
         PDOT(I) = PDOT0(I)
 1050 CONTINUE
      DO 1060 I=1,NU
         UDOT(I) = UDOT0(I)
 1060 CONTINUE
      IF (MDISC .EQ. 1) THEN
         DO 1070 I=1,NV
            FL1(I) = FL0(I)
 1070    CONTINUE
      ENDIF
C
C
C ---  General step  ---------------------------------------------------
C
      DO 1280 ISTEP=1,NJJ
C
         IF (ISTEP .GT. 1) THEN
C
C           Function evaluations: f, pdot, udot
C
            QFLAG(1) = .FALSE.
            QFLAG(2) = .FALSE.
            QFLAG(3) = .FALSE.
            QFLAG(4) = .TRUE.
            QFLAG(5) = .TRUE.
            QFLAG(6) = .TRUE.
            QFLAG(7) = .FALSE.
            QFLAG(8) = .FALSE.
            QFLAG(9) = .FALSE.
            JFAIL = 0
            IF (QPRTIM) CALL MONON (5)
            CALL MXPROB (NP, NV, NL, NG, NU, T, P, V, U, RLAM,
     $           NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $           N2UDIM, N6UDIM, NVLDIM,
     $           LDA, LDG, LDF, MDA, MDG, MDF,
     $           AM, GP, FL,
     $           F, PDOT, UDOT, G, GI,
     $           QFLAG, FPROB, FPROB, FPROB,
     $           MXJOB, LIWKL, IWKL, JFAIL)
            IF (QPRTIM) CALL MONOFF (5)
            IF (JFAIL .NE. 0) GOTO 1400
C
            IF (MDISC .EQ. 1) CALL MMULTF (NV, NL, LDF, MDF,
     $           FL, FL1, RLAM, MXJOB)
C
         ENDIF
C
C   position update
C
         T = T + HJ
         DO 1080 I=1,NP
            P(I) = P(I) + HJ*PDOT(I)
 1080    CONTINUE
C
C   additional dynamics
C
         DO 1090 I=1,NU
            U(I) = U(I) + HJ*UDOT(I)
 1090    CONTINUE
C
C  function evaluations: am, gp, gi
C
         QFLAG(1) = .TRUE.
         QFLAG(2) = .TRUE.
         QFLAG(3) = .FALSE.
         QFLAG(4) = .FALSE.
         QFLAG(5) = .FALSE.
         QFLAG(6) = .FALSE.
         QFLAG(7) = .FALSE.
         QFLAG(8) = .TRUE.
         QFLAG(9) = .FALSE.
C
         JFAIL = 0
C
         IF (QPRTIM) CALL MONON (6)
         CALL MXPROB (NP, NV, NL, NG, NU, T, P, V, U, RLAM,
     $        NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $        N2UDIM, N6UDIM, NVLDIM,
     $        LDA, LDG, LDF, MDA, MDG, MDF,
     $        AM, GP, FL,
     $        F, PDOT, UDOT, G, GI,
     $        QFLAG, FPROB, FPROB, FPROB,
     $        MXJOB, LIWKL, IWKL, JFAIL)
         IF (QPRTIM) CALL MONOFF (6)
         IF (JFAIL .NE. 0) GOTO 1400
C
C   right-hand side
C
         CALL MMULT (NV, NL, LDA, AM, V, Z, MXJOB)
C
         DO 1100 I=1,NV
            Z(I) = Z(I) + HJ*F(I)
 1100    CONTINUE
         IF (MDISC .EQ. 1) THEN
            DO 1110 I=1,NV
               Z(I) = Z(I) - HJ * FL1(I) 
 1110       CONTINUE
         ENDIF
         DO 1120 I=1,NL
            Z(NV+I) = -GI(I)
 1120    CONTINUE
C
C   linear system
C
         IFDEC = 0
         IF (QPRTIM) CALL MONON (13)
         CALL ADEC (NV, NL, LDA, LDG, LDF, AM, GP, FL, MXJOB, NVLDIM,
     $        RWKL, IWKL, IFDEC, IUPDTE)
         IF (QPRTIM) CALL MONOFF (13)
         IF (IFDEC .NE. 0) GOTO 1410
C
         IF (QPRTIM) CALL MONON (12)
         CALL ASOL (NV, NL, Z, MXJOB, NVLDIM, RWKL, IWKL)
         IF (QPRTIM) CALL MONOFF (12)
C
C   velocity, acceleration and Lagrange multiplier update
C
         DO 1130 I=1,NV
            A(I) = (Z(I) - V(I))/HJ
            V(I) = Z(I)
 1130    CONTINUE
C
         DO 1140 I=1,NL
            RLAM(I) = Z(NV+I)/HJ
 1140    CONTINUE
C
         IF (MXJOB(52) .EQ. 0 .AND. ISTEP .EQ. 1 .AND. IALAM .EQ. 1)
     $        THEN
            DO 1150 I=1,NV
               TAB(NX+I, J) = A(I)
 1150       CONTINUE
            DO 1160 I=1,NL
               TAB(NX+NV+I, J) = RLAM(I)
 1160       CONTINUE
         ENDIF
C
C        stability check
C
         IF (J .LE. 2 .AND. ISTEP .EQ. 2) THEN
            DEL1 = 0.D0
            DEL2 = 0.D0
            DO 1170 I=1,NV
               IF (ITOL .EQ. 0) THEN
                  WT = RTOL(1)*SCAL(NP+I) + ATOL(1)
               ELSE
                  WT = RTOL(NP+I)*SCAL(NP+I) + ATOL(NP+I)
               ENDIF
               DEL1 = DEL1 + (Z(I)/WT)**2
               DEL2 = DEL2 + (HJ*A(I)/WT)**2
 1170       CONTINUE
            IF (DEL2 .GT. DEL1*4.D0) THEN
               GOTO 1390
            ENDIF
         ENDIF
C     
         IF (IDOFLG .EQ. 1) THEN
            ISAFE = ISAFE + 1
            NS = 0
            DO 1180 I=1,NDP
               YSAFE(NS+I,ISAFE) = PDOT(INDP(I))
 1180       CONTINUE
            NS = NDP
            DO 1190 I=1,NDV
               YSAFE(NS+I,ISAFE) = A(INDV(I))
 1190       CONTINUE
            NS = NS + NDV
            DO 1200 I=1,NDU
               YSAFE(NS+I,ISAFE) = UDOT(INDU(I))
 1200       CONTINUE
            NS = NS + NDU
            DO 1210 I=1,NDA
               YSAFE(NS+I,ISAFE) = A(INDA(I))
 1210       CONTINUE
            NS = NS + NDA
            DO 1220 I=1,NDL
               YSAFE(NS+I,ISAFE) = RLAM(INDL(I))
 1220       CONTINUE
         ENDIF
C
         IF (IDOFLG .EQ. 2) THEN
            ISAFE = ISAFE + 1
            NS = 0
            DO 1230 I=1,NP
               YSAFE(NS+I,ISAFE) = PDOT(I) 
 1230       CONTINUE
            NS = NP
            DO 1240 I=1,NV
               YSAFE(NS+I,ISAFE) = A(I)
 1240       CONTINUE    
            NS = NS + NV
            DO 1250 I=1,NU
               YSAFE(NS+I,ISAFE) = UDOT(I)
 1250       CONTINUE    
            NS = NS + NU
            DO 1260 I=1,NV
               YSAFE(NS+I,ISAFE) = A(I)
 1260       CONTINUE    
            NS = NS + NV
            DO 1270 I=1,NDL
               YSAFE(NS+I,ISAFE) = RLAM(I)
 1270       CONTINUE
         ENDIF
C
 1280 CONTINUE
C
C
C ---    Compute j-th line of extrapolation tableau  -------------------
C
C Set element in the first column of the tableau
C
      DO 1290 I=1,NP
         TAB(I,J) = P(I)
 1290 CONTINUE
      NS=NP
      DO 1300 I=1,NV
         TAB(NS+I,J) = V(I)
 1300 CONTINUE
      NS=NS+NV
      DO 1310 I=1,NU
         TAB(NS+I,J) = U(I)
 1310 CONTINUE
      NS=NS+NU
      DO 1320 I=1,NV
         TAB(NS+I, J) = A(I)
 1320 CONTINUE
      NS=NS+NV
      DO 1330 I=1,NL
         TAB(NS+I, J) = RLAM(I)
 1330 CONTINUE
C
      IF (J .EQ.1) RETURN
      IF (MXJOB(52) .EQ. 0 .AND. IALAM .EQ. 1) NX = N6UDIM
C
C ---    h - extrapolation
C
      DO 1350 L=J,2,-1
         FAC = (DBLE(NJ(J))/DBLE(NJ(L-1))) - 1.D0
         DO 1340 I=1,NX
            TAB(I,L-1) = TAB(I,L) + (TAB(I,L)-TAB(I,L-1))/FAC
 1340    CONTINUE
 1350 CONTINUE
C
C ---    Error estimates      -----------------------------------------
C
      ERR = 0.D0
      DO 1360 I=1,NP
         SCAL(I) = MAX (ABS(P0(I)), ABS(TAB(I,1)))
         IF (ITOL .EQ. 0) THEN
            WT = RTOL(1)*SCAL(I) + ATOL(1)
         ELSE
            WT = RTOL(I)*SCAL(I) + ATOL(I)
         ENDIF
         ERR = ERR + MIN (ABS((TAB(I,1) - TAB(I,2)))/WT, 1.D20)**2
 1360 CONTINUE
      ERRP = ERR / DBLE(NP)
C
      ERR = 0.D0
      DO 1370 I=1,NV
         SCAL(NP+I) = MAX (ABS(V0(I)), ABS(TAB(NP+I,1)))
         IF (ITOL .EQ. 0) THEN
            WT = RTOL(1)*SCAL(NP+I) + ATOL(1)
         ELSE
            WT = RTOL(NP+I)*SCAL(NP+I) + ATOL(NP+I)
         ENDIF
         ERR = ERR + MIN (ABS((TAB(NP+I,1) -
     $        TAB(NP+I,2)))/WT, 1.D20)**2
 1370 CONTINUE
      ERRV = ERR / DBLE(NV)
C
      ERRU = 0.D0
      DO 1380 I=1,NU
         SCAL(NN+I) = MAX (ABS(U0(I)), ABS(TAB(NN+I,1)))
         IF (ITOL .EQ. 0) THEN
            WT = RTOL(1)*SCAL(NN+I) + ATOL(1)
         ELSE
            WT = RTOL(NN+I)*SCAL(NN+I) + ATOL(NN+I)
         ENDIF
         ERRU = ERRU + MIN (ABS((TAB(NN+I,1) - TAB(NN+I,2)))/WT,
     $        1.D20)**2
 1380 CONTINUE
      IF (NU .GT. 0) ERRU = ERRU / DBLE(NU)
C
      ERR = SQRT(ERRP + ERRV + ERRU)
      IF (MNINT .GE. 4) THEN
         ERRPRI = ERR
         WRITE (LUINT,*) ' - MEXXEX: Est. normalized error:', J, ERRPRI
      ENDIF
      IF (J .GT. 2 .AND. ERR .GE. ERROLD) GOTO 1390
      ERROLD = MAX (ERR, TOLD4)
C
C ---    Compute optimal step sizes    --------------------------------
C
      EXPO = 1.D0/DBLE(J-1)
      FACMIN = FAC1**EXPO
      FAC = MIN (FAC2/FACMIN, MAX (FACMIN, (ERR/TOLD4)**EXPO/SAFE2))
      FAC = 1.D0/FAC
      HH(J) = MIN (H*FAC, HMAX)
      WK(J) = AA(J)/HH(J)
      RETURN
C
C ---    Step rejected     --------------------------------------------
C
 1390 QFAIL = .TRUE.
      QRJECT = .TRUE.
      H = H*SAFE3
      RETURN
C
C ---    Error return from user function ------------------------------
C
 1400 MXJOB(73) = JFAIL
      IFAIL = -9
      RETURN
C
C ---    Error return from ADEC ---------------------------------------
C
 1410 MXJOB(73) = IFDEC
      IFAIL = -15
      RETURN
C
C
C  end of subroutine mexxex
C
      END
C
      SUBROUTINE MXDDU1 (INFO, ND, NDP, NDV, NDU, NDA, NDL,
     $     KDENS, ILAM, IRHO, T, TOLD, H, DENS, YIP,
     $     T0, TFIN, DOUTS,
     $     NDUMMY, NDPH, NDVH, NDUH, NDAH, NDLH,
     $     INDP, INDV, INDU, INDA, INDL, IRTRN)
C
C
C  This routine dumps the data which are usually passed to MXDOUT
C  to a file (the coefficients of hermite interpolation and the user 
C  specifications (selction of dense output components,NDOUT,...))
C  This allows dense output postprocessing with the special driver
C  routine MXPST1. MXPST1 read the data with MXDRE1 and calls MXDOUT
C  as if MEXX would do this, except for one case. 
C
C  ND,NDP,.. are the values used for the last call to MXDENS 
C  NDH,NDPH,.. are the values prescibed by the user
C  usually the values of ND,NDP,..  and NDH,NDPH,.. will coincide
C  but they will differ, if the user selected special components
C  and in the current step the root finder was invoked. In such
C  a situation the values must be arranged in this routine. In
C  all other cases, one can dump according to ND,NDP,..
C  and pass these values (and the solution) to the dump file.
C
C* File              mxddu1.f
C* Version           1.1
C* Latest Change     96/04/03 (1.6)
C* Modification history
C    1.1             Avoid free format because of compiler dependency.
C  
C
C  INPUT
C=======
C  
C   INFO(50)  Int   Info array
C                   INFO(1) .. INFO(40) see MXJOB(1) .. MXJOB(40) 
C                                           in MEXX
C                   INFO(41): Type of call: 
C                             = 0: the first call
C                             = 1: an intermediate call
C                             = 2: the last call at TFIN
C                             = 3: a root of the switching function was
C                                  found within the last integration
C                                  step and the integration continues
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
C  t     :  current integration point or root value (told+h .ne. t)
C
C  yip   :  for istep=0: the current "solution" of selected
C           components at t=t0
C           workspace otherwise
C
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER INFO(50), ND, NDP, NDV, NDU, NDA, NDL, KDENS,
     $        ILAM, IRHO, IRTRN,
     $        INDP(*), INDV(*), INDU(*), INDA(*), INDL(*)
      DOUBLE PRECISION T, TOLD, H, DENS(ND,KDENS), YIP(ND)
      DOUBLE PRECISION DOUTS(*)
C
      SAVE
C
      ITYPE = INFO(41)
      IRTRN = 0
      NDH = NDPH+NDVH+NDUH+NDAH+NDLH
C
C
C  Output at the first call (icall=0)
C========================================
C
      IF (ITYPE .EQ. 0) THEN
         LUNIT = INFO(16)
         WRITE (LUNIT,9010) ' '
         WRITE (LUNIT,9010) ' type_1 dumpfile from mexx21'
         WRITE (LUNIT,9010) ' '
         WRITE (LUNIT,9010) 'C  the options vector INFO:'
         WRITE (LUNIT,9020) (INFO(I),I=1,50)
         MDOUT = INFO(31)
         NDOUT = INFO(32)
         WRITE (LUNIT,9150)
     $        NDPH,NDVH,NDUH,NDAH,NDLH,' = ndp,ndv,ndu,nda,ndl'
C
         WRITE (LUNIT,9010) 'C  the selected components of P are:'
         WRITE (LUNIT,9020) (INDP(I),I=1,NDPH)
         WRITE (LUNIT,9010) 'C  the selected components of V are:'
         WRITE (LUNIT,9020) (INDV(I),I=1,NDVH)
         WRITE (LUNIT,9010) 'C  the selected components of U are:'
         WRITE (LUNIT,9020) (INDU(I),I=1,NDUH)
         WRITE (LUNIT,9010) 'C  the selected components of A are:'
         WRITE (LUNIT,9020) (INDA(I),I=1,NDVH)
         WRITE (LUNIT,9010) 'C  the selected components of LAMBDA are:'
         WRITE (LUNIT,9020) (INDL(I),I=1,NDLH)
C
         WRITE (LUNIT,9120) NDH,KDENS,' = nd,kdens'
         IF (MDOUT .EQ. 3) THEN 
            WRITE (LUNIT,9210) DOUTS(1),' = delta_t_max'
         ENDIF
         IF (MDOUT .EQ. 4) THEN 
            WRITE (LUNIT,9030) (DOUTS(L),L=1,NDOUT) 
         ENDIF
         WRITE (LUNIT,9220) T0,TFIN,' = t0,tfin'
         IF (NDH .EQ. ND) THEN
            WRITE (LUNIT,9030) (YIP(L),L=1,ND) 
         ELSE
            NS = 0
            NSH = 0
            DO 1000 I=1,NDPH
               YIP(NS+I) = YIP(NS+INDP(I))
 1000       CONTINUE
            NS = NDP
            NSH = NDPH
            DO 1010 I=1,NDVH
               YIP(NSH+I) = YIP(NS+INDV(I))
 1010       CONTINUE
            NS = NS + NDV
            NSH = NSH + NDVH
            DO 1020 I=1,NDUH
               YIP(NSH+I) = YIP(NS+INDU(I))
 1020       CONTINUE
            NS = NS + NDU
            NSH = NSH + NDUH
            DO 1030 I=1,NDAH
               YIP(NSH+I) = YIP(NS+INDA(I))
 1030       CONTINUE
            NS = NS + NDA
            NSH = NSH + NDAH
            DO 1040 I=1,NDLH
               YIP(NSH+I) = YIP(NS+INDL(I))
 1040       CONTINUE
            WRITE (LUNIT,9030) (YIP(L),L=1,NDH) 
         ENDIF
      ENDIF
C
C  OUTPUT at the other calls
C===========================
C
      IF (ITYPE .GT. 0) THEN
         KDDUMP = ILAM+IRHO+3
         WRITE (LUNIT,9110) INFO(41),' = itype'
         WRITE (LUNIT,9110) INFO(42),' = istep'
         WRITE (LUNIT,9230) TOLD,T,H,' = told,t,h'
         WRITE (LUNIT,9130) ILAM,IRHO,KDDUMP,' = ilam,irho,kddump'
         DO 1100 K=1,KDDUMP
            WRITE (LUNIT,9110) K,' = no of column of dens'
            IF (NDH .EQ. ND) THEN
               WRITE (LUNIT,9030) (DENS(L,K),L=1,NDH)
            ELSE
               NS = 0 
               NSH = 0 
               DO 1050 I=1,NDPH 
                  YIP(NS+I) = DENS(NS+INDP(I),K) 
 1050          CONTINUE 
               NS = NDP 
               NSH = NDPH 
               DO 1060 I=1,NDVH 
                  YIP(NSH+I) = DENS(NS+INDV(I),K) 
 1060          CONTINUE 
               NS = NS + NDV 
               NSH = NSH + NDVH 
               DO 1070 I=1,NDUH 
                  YIP(NSH+I) = DENS(NS+INDU(I),K) 
 1070          CONTINUE 
               NS = NS + NDU 
               NSH = NSH + NDUH 
               DO 1080 I=1,NDAH
                  YIP(NSH+I) = DENS(NS+INDA(I),K) 
 1080          CONTINUE 
               NS = NS + NDA 
               NSH = NSH + NDAH 
               DO 1090 I=1,NDLH 
                  YIP(NSH+I) = DENS(NS+INDL(I),K) 
 1090          CONTINUE
               WRITE (LUNIT,9030) (YIP(L),L=1,NDH) 
            ENDIF
 1100    CONTINUE
      ENDIF
C
      RETURN
C
 9010 FORMAT (A)
 9020 FORMAT ((10(2X,I6)))
 9030 FORMAT ((3(2X,1PD21.13)))
C
 9110 FORMAT (2X,I6, A)
 9120 FORMAT (2(2X,I6), A)
 9130 FORMAT (3(2X,I6), A)
 9150 FORMAT (5(2X,I6), A)
C
 9210 FORMAT (2X,1PD21.13, A)
 9220 FORMAT (2(2X,1PD21.13), A)
 9230 FORMAT (3(2X,1PD21.13), A)
C
      END
C
      SUBROUTINE MXDOUT(INFO, ND, NDP, NDV, NDU, NDA, NDL,
     $     KDENS, ILAM, IRHO, T, TOLD, H, DENS, YIP,
     $     T0, TFIN, DOUTS, DENOUT,
     $     NDH, NDPH, NDVH, NDUH, NDAH, NDLH,
     $     INDP, INDV, INDU, INDA, INDL, QPRTIM, IRTRN)
C
C  MXDOUT generates the dense output solution values according
C  to the user selected options and calls the dense output
C  routine DENOUT supplied by the user.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      PARAMETER (NPTSMX=1000)
C
      DOUBLE PRECISION THALT(NPTSMX)
C
      LOGICAL QPRTIM
      INTEGER INFO(50), ND, NDP, NDV, NDU, NDA, NDL, KDENS,
     $        ILAM, IRHO, IRTRN,
     $        INDP(*), INDV(*), INDU(*), INDA(*), INDL(*)
      DOUBLE PRECISION T, TOLD, H, DENS(ND,KDENS), YIP(ND+1)
      DOUBLE PRECISION DOUTS(*) 
C
      INTEGER IDEL, ILAMH, NDY, NDZ
      DOUBLE PRECISION X
C
      EXTERNAL DENOUT
C
      SAVE
C
C
C  nd,ndp,.. are the values used for the last call to mxdens
C  ndh,ndph,.. are the values prescibed by the user 
C  usually the values of nd,ndp,..  and ndh,ndph,.. will coincide
C  but they will differ, if the user selected special components
C  and in the current step the root finder was invoked. In such
C  a situation the values must be arranged in this routine. In
C  all other cases, one can interpolate according to nd,ndp,..
C  and pass these values (and the solution) to denout.
C  
C
C  INPUT
C=======
C  
C   INFO(50)  Int   Info array
C                   INFO(1) .. INFO(40) see MXJOB(1) .. MXJOB(40) 
C                                           in MEXX
C                   INFO(41): Type of call: 
C                             = 0: the first call
C                             = 1: an intermediate call
C                             = 2: the last call at TFIN
C                             = 3: a root of the switching function was
C                                  found within the last integration
C                                  step and the integration continues
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
C
C
C* File              mxdout.f
C* Version           1.1
C* Latest Change     96/03/14 (1.7)
C* Modification history
C    1.1             Correct setting of INFO(41)
C                    New description for INFO(3)
C                    Adapt to ZIB GUI.
C      
C
      IRTRN = 0
      MDOUT = INFO(31)
      NDOUT = INFO(32)
C
C
C  Initiation at the first call 
C==============================
C
      IF (INFO(41) .EQ. 0) THEN
         IDP=1
         IDV=IDP+NDPH
         IDU=IDV+NDVH
         IDA=IDU+NDUH
         IDL=IDA+NDAH
C
         IF(MDOUT.EQ.1) THEN
            TDEL=(TFIN-T0)/DBLE(NDOUT+1)
            IH1=2
            NPTS=1
            THALT(1)=T0
         ENDIF
C
         IF(MDOUT.EQ.2) THEN
            NPTS=1
            THALT(1)=T0
         ENDIF
C
         IF(MDOUT.EQ.3) THEN
            NPTS=1
            THALT(1)=T0
         ENDIF
C
         IF(MDOUT.EQ.4) THEN
            NPTS=0
            IH1=1
            IF(DOUTS(1).LE.T) THEN
               NPTS=1
               THALT(1)=T0
               IH1=2
            ENDIF
         ENDIF
C
      ENDIF
C
C  Subsequent calls
C==================
C
      IF (INFO(41) .GT. 0) THEN
C
         IF(MDOUT.EQ.1) THEN
            II=0 
            IH2=IH1
            DO 1000 I=IH1,NDOUT+2
               TI=T0+DBLE(I-1)*TDEL
               IF(TI.LE.T) THEN
                  II=II+1
                  THALT(II)=TI
               ELSE
                  IH2=I
                  GOTO 1010
               ENDIF
 1000       CONTINUE
 1010       CONTINUE
            IH1=IH2
            NPTS=II
         ENDIF
C
         IF(MDOUT.EQ.2) THEN
            TDEL=(T-TOLD)/DBLE(NDOUT+1)
            DO 1020 I=1,NDOUT+1
               THALT(I)=TOLD+DBLE(I)*TDEL
 1020       CONTINUE
            NPTS=NDOUT+1
         ENDIF
C
         IF(MDOUT.EQ.3) THEN
            XIP=(T-TOLD)/DOUTS(1)
            IF(XIP.LE.1.D0) THEN
               NPTS=1
               THALT(1)=T
            ELSE
               NPTS=(XIP+1.D0)
               TDEL=(T-TOLD)/DBLE(NPTS)    
               DO 1030 I=1,NPTS    
                  THALT(I)=TOLD+DBLE(I)*TDEL
 1030          CONTINUE
            ENDIF
         ENDIF
C
         IF(MDOUT.EQ.4) THEN
            II=0 
            IH2=IH1 
            DO 1040 I=IH1,NDOUT
               TI=DOUTS(I)
               IF(TI.LE.T) THEN 
                  II=II+1 
                  THALT(II)=TI 
               ELSE 
                  IH2=I 
                  GOTO 1050 
               ENDIF 
 1040       IH2=NDOUT+1
 1050       CONTINUE 
            IH1=IH2 
            NPTS=II 
         ENDIF
C
         IF (ILAM .GT. 0) THEN 
            IDEL = 0
            ILAMH = ILAM - 1
         ELSE
            IDEL = 1
            ILAMH = 0
         ENDIF
C
      ENDIF
C
C  For all calls
C===============
C
      NDY = NDP + NDV + NDU
      NDZ = NDA + NDL
      INFO1= INFO(41)
C
      DO 1110 IPTS=1,NPTS
C
         X = THALT(IPTS) 
C
         IF(INFO1.GT.0) THEN
         IF(QPRTIM) CALL MONON(4)
            IF(NDY.GT.0) THEN
              CALL MXIPOL (ND, DENS, NDY, X, TOLD, H, ILAM, IRHO, YIP)
            ENDIF
            IF(NDZ.GT.0) THEN
               CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, X, TOLD,
     $                      H, ILAMH, IRHO, YIP(NDY+1))
            ENDIF
         IF(QPRTIM) CALL MONOFF(4)
         ENDIF
C
         IF(ND.NE.NDH) THEN
            NS = 0
            NSH = 0
            DO 1060 I=1,NDPH
               YIP(NS+I) = YIP(NS+INDP(I))
 1060       CONTINUE
            NS = NDP
            NSH = NDPH
            DO 1070 I=1,NDVH
               YIP(NSH+I) = YIP(NS+INDV(I))
 1070       CONTINUE
            NS = NS + NDV
            NSH = NSH + NDVH
            DO 1080 I=1,NDUH
               YIP(NSH+I) = YIP(NS+INDU(I))
 1080       CONTINUE
            NS = NS + NDU
            NSH = NSH + NDUH
            DO 1090 I=1,NDAH
               YIP(NSH+I) = YIP(NS+INDA(I))
 1090       CONTINUE
            NS = NS + NDA
            NSH = NSH + NDAH
            DO 1100 I=1,NDLH
               YIP(NSH+I) = YIP(NS+INDL(I))
 1100       CONTINUE
         ENDIF
         IF (INFO1 .EQ. 0) THEN
            INFO(43) = 0
         ELSE
            INFO(43)=1
         ENDIF
         IF(INFO1.EQ.1 .OR. INFO1.EQ.2) THEN
            IF(ABS(T-X).LE.1.D-6*T) THEN      
               INFO(43)=0
               INFO(41)=INFO1
            ELSE
               INFO(41)=1    
            ENDIF
         ENDIF
         IF(INFO1.EQ.4) THEN
            IF(ABS(T-X).LE.1.D-6*T) THEN      
               INFO(41)=INFO1
            ELSE
               INFO(41)=1    
            ENDIF
         ENDIF
         INFO(44) = INFO(44) + 1
         IRTRN = 0
         CALL DENOUT (NDPH,NDVH,NDUH,NDAH,NDLH,
     &                X,YIP(IDP),YIP(IDV),YIP(IDU),YIP(IDA),YIP(IDL),
     &                INDP,INDV,INDU,INDA,INDL,INFO,IRTRN)
C
 1110 CONTINUE  
C
C
      RETURN
C
C  end of subroutine MXDOUT
C
      END
C
      SUBROUTINE MXLIN (N, M, LDA, MXJOB, IERR)
C
C ---------------------------------------------------------------------
C
C     This subroutine performs linear algebra for MEXX
C
C* File              mexxlin.f
C* Version           1.0
C* Latest Change     96/05/08 (1.6)
C
C ---------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER MXJOB(150), IWK(*)
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION AM(LDAM,*), GP(LDGP,*), FL(LDF,*), B(N+M), V(N),
     $     W(M), Z(N)
C
C
C
C  +++ Sparse LU-decomposition (MXJOB(2)=3) with MA28.. package +++
C
C  Problem depedent dimensions: 
C    
C    N      -  Dimension of square matrix AM
C    M      -  Row dimension of rectangular matrix G
C    NM     -  Dimension of square matrix A = |AM G^T-F0|  
C              ( = N + M )                    |G  0     |
C    NZGMAX -  Expected maximum number of nonzero elements of G
C    NZG    -  Actual maximum number of nonzero elements of G
C    NZFMAX -  Expected maximum number of nonzero elements of FL
C    NZF    -  Actual maximum number of nonzero elements of FL
C    NBLK   -  Number of blocks of AM
C    NMRC   -  Dimension of blocks of AM
C    NZA    -  Number of nonzero elements of matrix A
C              ( = 2*NZG + NBLK*NMRC**2)
C
C  Storage requirements
C
C    IROWG  -  NZGMAX  -  Row indices of G
C                         (input to MCOPY)
C    JCOLG  -  NZGMAX  -  Column indices of G
C                         (input to MCOPY)
C    IROWF  -  NZFMAX  -  Row indices of FL
C                         (input to MCOPY)
C    JCOLF  -  NZFMAX  -  Column indices of FL
C                         (input to MCOPY)
C    IROW   -  NZA     -  Row indices of A
C                         (output of MCOPY, input to MA28B)
C    JCOL   -  NZA     -  Column indices of A
C                         (output of MCOPY, input to MA28B)
C    IKEEP  -  MIKEEP  -  Intermediate storage
C                         (output of MA28A, input to MA28B/C)
C              MIKEEP  =  5*NM
C    IWMA   -  MIW     -  Work array 
C              MIW     =  8*NM
C    IRN    -  LIRN    -  For row indices of A/decomposed A
C                         (inout for MA28A, input to MA28B/C)
C              LIRN    =  NZA*FILLR, FILLR depends on problem and
C                                    algorithmic performance
C    ICN    -  LICN    -  For column indices of A/decomposed A
C                         (inout for MA28A, input to MA28B/C)
C              LICN    =  NZA*FILLC, FILLC depends on problem and
C                                    algorithmic performance
C                         
C    RWMA   -  NM      -  Work array
C
C    A      -  NZA/LICN-  On input the matrix A, on output the
C                         LU decomposition
C
C    FILLR, FILLC are set in parameter statement below
C 
C  Split of work space:
C
C  KIN=0, KIN=1:
C    IPVT(1, ..., NP+NG)  -  IWK(1, ..., NP+NG)
C
C  KIN=2:
C    JPVTc(1, ..., NP)    -  IWK(1, ..., NP)                 , if MODE=0
C                               in NBLK blocks of length NMRC, if MODE=1
C    JPVTq(1, ..., NG)    -  IWK(NP+1, ..., NP+NG)
C
C  KIN=3:
C    IPVTk(1, ..., NG)    -  IWK(1, ..., NG)
C    IPVTm(1, ..., NP-NG) -  IWK(NG+1, ..., NP)
C
C  KIN=10:
C    IROWG(1, ..., NZG)   -  IWK(IIROWG, ..., NZG)
C                            IIROWG = 1
C    JCOLG(1, ..., NZG)   -  IWK(IJCOLG, ...)
C                            IJCOLG = IIROWG + NZGMAX
C    IROWF(1, ..., NZF)   -  IWK(IIROWF, ..., NZF)
C                            IIROWF = 1
C    JCOLF(1, ..., NZF)   -  IWK(IJCOLF, ...)
C                            IJCOLF = IIROWF + NZFMAX
C    IROW(1, ..., NZA)    -  IWK(IIROW, ...)
C                            IIROW  = IJCOLG +NZGMAX
C    JCOL(1, ..., NZA)    -  IWK(IJCOL, ...)
C                            IJCOL  = IIROW + NZA
C    IKEEP(1, ..., MIKEEP)-  IWK(IIKEEP, ...)
C                            IIKEEP = IJCOL + NZA
C    IWMA(1, ..., MIW)    -  IWK(IIWMA, ...)
C                            IIWMA  = IIKEEP + 5*NM
C    IRN(1, ..., LIRN)    -  IWK(IIRN, ...)
C                            IIRN   = IIWMA + 8*NM
C    ICN(1, ..., LICN)    -  IWK(IICN, ...)
C                            IICN = 2*NZGMAX + 2*NZA + 13*NM + LIRN + 1
C
C    Minimum length of IWK:  LIWKMI = 2*NZGMAX + 2*NZFMAX + 2*NZA
C                                     + 13*NM + LIRN + LICN
C
C    RWMA(1, ..., NM)     -  A(1, ..., NM)
C    A(1, ..., NZA)       -  A(LP1, ...)
C                            LP1 = NM + 1
C
C    Minimum length of RWK:  LRWKMI = NM + LICN
C
C    For a description of the internal parameters of the package
C    see: Harwell Subroutine Library Specification
C                  
C
      INTEGER LP, MP, IRNCP, ICNCP, MINIRN, MINICN, IRANK, IDISP(2)
      LOGICAL QBLOCK, QGROW, QABRT1, QABRT2
      DOUBLE PRECISION EPSMA, RMIN, RESID
C
      COMMON /MA28ED/ LP, MP, QBLOCK, QGROW
      COMMON /MA28FD/ EPSMA, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
     $     IRANK, QABRT1, QABRT2
      COMMON /MA28GD/ IDISP
C
C                       Description see Harwell Subroutine
C                       Library Specification
C
      SAVE /MA28ED/, /MA28FD/, /MA28GD/
C
C  Essential parameters to control performance within MEXX:
C
C    QFACTH -  If .true. perform 'analyse factor' for each decomposition
C    FILLC  -  Fillin factor for column indices during decomposition
C              (ordinarily 2 to 4)
C    FILLR  -  Fillin factor for row indices during decomposition
C              (ordinarily not very much greater than 1)
C    THRSH1 -  Value between zero and one to control choice of pivots
C              in MA28AD
C    THRSH2 -  Smallest ratio of the pivot to the largest entry in the
C              row monitoring the stability of successive factorizations
C              in MA28BD
C
      LOGICAL QFACT, QFACTH
C
      PARAMETER ( QFACTH = .FALSE.       ,
     $            FILLC  = 2.0D0         ,
     $            FILLR  = 1.1D0         ,
     $            THRSH1 = 1.D-2         ,
     $            THRSH2 = 1.D-6         )
C
C
C* Local variables
C
      LOGICAL QREFCT
      LOGICAL QPRTIM
C
      SAVE
C
C
      IERR = 0
      NP = N
      NG = M
      MDISC = MXJOB(1)
      KIN = MXJOB(2)
      MODE = MXJOB(3)
      NBLK = MXJOB(4)
      NMRC = MXJOB(5)
      NZGMAX = MXJOB(7)
      NZFMAX = MXJOB(8)
      IMON = MXJOB(13)
      LOUT = MXJOB(14)
      QPRTIM = MXJOB(19) .GT. 0
      QREFCT = .FALSE.
      QFACT = .TRUE. 
C     
C  Check for sufficient storage
C
C    Workspace left for linear algebra
      LRWK = MXJOB(80)
      LIWK = MXJOB(81)
      NPGDIM = MXJOB(82)
C
C    Workspace required
      IF (KIN .LE. 1) THEN
         LIWKMI = N + M
         LRWKMI = (N + M)**2
      ELSE IF (KIN .EQ. 2) THEN
         LIWKMI = N + M
         LRWKMI = (N + M)**2 + 2*(N + M)
      ELSE IF (KIN .EQ. 3) THEN
         LIWKMI = N
         LRWKMI = (N + M)**2
      ELSE IF (KIN .EQ. 10) THEN
         NM = N + M
         NZA = 2*NZGMAX + 2*NZFMAX + NBLK*NMRC**2
         LIRN = DBLE(NZA)*FILLR
         LICN = DBLE(NZA)*FILLC
         LIWKMI = 2*NZGMAX + 2*NZFMAX + 2*NZA + 13*NM + LIRN + LICN
         LRWKMI = NM + LICN
      ENDIF
C
C    Minimum workspace required for decomposition
      MXJOB(86) = LRWKMI
      MXJOB(87) = LIWKMI
      MXJOB(88) = LRWKMI
      MXJOB(89) = LIWKMI
      MXJOB(90) = NZA    
      MXJOB(91) = LIRN   
      MXJOB(92) = LICN  
C
      IF (LIWKMI .GT. LIWK) THEN
         IERR = -1
         RETURN
      ENDIF
      IF (LRWKMI .GT. LRWK) THEN
         IERR = -2
         RETURN
      ENDIF
C
C  Workspace distribution
C------------------------
C
      IF (KIN .EQ. 0 .OR. KIN .EQ. 1) THEN
         IIPVT = 1
C     IFREE = IIPVT + NP + NG
C     
      ELSE IF (KIN .EQ. 2) THEN
         IJPVTC = 1
         IJPVTQ = IJPVTC + NP
C     IFREE  = IJPVTQ + NG
C     
      ELSE IF (KIN .EQ. 3) THEN
         IIPVTK = 1
         IIPVTM = 1 + NG
C        IFREE  = IIPVTM + NP - NG
C
      ELSE IF (KIN .EQ. 10) THEN
C
C  Adopt FILLR, FILLC to available workspace
C
         IFREEI = LIWK - LIWKMI
         FFACI  = DBLE(IFREEI)/DBLE(NZA)
         FFACI  = FFACI/2.D0
         IFREER = LRWK - LRWKMI
         FFACR  = DBLE(IFREER)/DBLE(NZA)
         FFAC   = DMIN1 (FFACR, FFACI)
         FILLRH = FILLR + FFAC
         FILLCH = FILLC + FFAC
         LIRN   = DBLE(NZA)*FILLRH
         LICN   = DBLE (NZA)*FILLCH
         LIWKMI = 2*NZGMAX + 2*NZFMAX + 2*NZA + 13*NM + LIRN + LICN
         LRWKMI = NM + LICN
C
         MXJOB (88) = LRWKMI
         MXJOB (89) = LIWKMI
         MXJOB (93) = LIRN   
         MXJOB (94) = LICN  
C
C  Set pointers
C
         IIROWG = 1
         IJCOLG = IIROWG + NZGMAX
         IIROWF = IJCOLG + NZGMAX
         IJCOLF = IIROWF + NZFMAX
         IIROW  = IJCOLF + NZFMAX
         IJCOL  = IIROW  + NZA
         IIKEEP = IJCOL  + NZA
         IIWMA  = IIKEEP + 5*NM
         IIRN   = IIWMA  + 8*NM
         IICN   = IIRN   + LIRN
C        IFREE  = IICN + LICN
         LP1 = NM + 1
      ENDIF
C
      RETURN
C     
      ENTRY MMULT (N, M, LDAM, AM, V, Z, MXJOB)
C
C ---------------------------------------------------------------------
C
C     This subroutine computes the matrix vector product
C
C     z = AM * v
C
C
C     Version of July 31, 1991.
C
C ---------------------------------------------------------------------
C
C
      MXJOB(57) = MXJOB(57) + 1
      GOTO (1010, 1010, 1010, 1040, 1000, 1000, 1000, 1000, 1000, 1000,
     $     1000, 1140, 1140, 1140, 1070, 1000, 1000, 1000, 1000, 1000,
     $     1000, 1140) MODE*11 + KIN + 1
C
C                     -------   MODE or KIN illegal   -------
 1000 CONTINUE
      WRITE (LOUT,*) ' Subroutine MMULT:',
     $     ' Illegal value of MODE and/or KIN ', MODE, KIN
      STOP
C                     -------   MODE = 0 / KIN = 0 or 1 or 2  -------
 1010 CONTINUE
      DO 1030 I=1,N
         Z(I) = 0.D0
         DO 1020 J=1,N
            Z(I) = Z(I) + AM(I,J)*V(J)
 1020    CONTINUE
 1030 CONTINUE
C
      RETURN
C
C                     -------   MODE = 0 / KIN = 3   -------
 1040 CONTINUE
      DO 1060 I=1,M
         Z(I) = 0.D0
         DO 1050 J=1,M
            Z(I) = Z(I) + AM(I,J)*V(J)
 1050    CONTINUE
 1060 CONTINUE
C
C
      RETURN
C
C                     -------   MODE = 1 / KIN = 3   -------
 1070 CONTINUE
      DO 1080 I=1,N
         Z(I) = 0.D0
 1080 CONTINUE
C
      DO 1130 IB=1,(M+1)/NMRC
         DO 1110 I=1,NMRC
            IF ((IB-1)*NMRC+I .GT. M) GOTO 1120
            SUM = 0.D0
            DO 1090 J=1,NMRC
               IF ((IB-1)*NMRC+J .GT. M) GOTO 1100
               SUM = SUM + AM((IB-1)*NMRC+I,J)*V((IB-1)*NMRC+J)
 1090       CONTINUE
 1100       CONTINUE
            Z((IB-1)*NMRC+I) = SUM
 1110    CONTINUE
 1120    CONTINUE
 1130 CONTINUE
C
      RETURN
C
C                -------   MODE = 1 / KIN = 0 or 1 or 2 or 10   -------
 1140 CONTINUE
C
      DO 1170 K=1,NBLK
         DO 1160 I=1,NMRC
            SUM = 0.D0
            DO 1150 J=1,NMRC
               SUM = SUM + AM((K-1)*NMRC+I,J)*V((K-1)*NMRC+J)
 1150       CONTINUE
            Z((K-1)*NMRC+I) = SUM
 1160    CONTINUE
 1170 CONTINUE
C
      RETURN
C
C
      ENTRY MMULTF (N, M, LDF, MDF, FL, V, W, MXJOB)
C
C ---------------------------------------------------------------------
C
C     This subroutine computes the matrix vector product
C
C     v = FL * w
C
C     in case  MODE = 0 and KIN = 4.
C
C     Version of April 29, 1992.
C
C ---------------------------------------------------------------------
C
C
      MXJOB(57) = MXJOB(57) + 1
C
      IF (MDISC .NE. 1 .OR. MODE .NE. 0) THEN
         WRITE (LOUT,*) ' Subroutine MMULTF:',
     $      ' Illegal value of MODE and/or MDISC ',MODE,MDISC
         STOP
      ENDIF
C
      DO 1180 I=1,N
         V(I) = 0.D0
 1180 CONTINUE
C
      DO 1200 J=1,M
         DO 1190 I=1,N
            V(I) = V(I) + FL(I,J)*W(J)
 1190    CONTINUE
 1200 CONTINUE
C
      RETURN
C
C* End of MMULTF
C
C
      ENTRY ADEC (N, M, LDAM, LDGP, LDF, AM, GP, FL, MXJOB, LDA, A, IWK,
     $     IFDEC, IUPDTE)
C
C ---------------------------------------------------------------------
C
C     This subroutine computes a decomposition of the matrix
C
C           ( AM   GP^T )
C           ( GP    0   )
C
C     Version of July 31, 1991.
C
C ---------------------------------------------------------------------
C
C
      IFDEC = 0
      EPSMA =DMIN1 (1.D0, 10.D0*THRSH2)
      MXJOB(58) = MXJOB(58) + 1
C
C
      IF (KIN .LE. 2) THEN
         IF (MODE .EQ. 1) THEN
            IF (MDISC .EQ. 1) GOTO 1280
C
            DO 1220 J=1,N
               DO 1210 I=1,N
                  A(I,J) = 0.D0
 1210          CONTINUE
 1220       CONTINUE
C
            DO 1250 K=1,NBLK
               DO 1240 J=1,NMRC 
                  DO 1230 I=1,NMRC
                     A((K-1)*NMRC+I,(K-1)*NMRC+J)=AM((K-1)*NMRC+I,J)
 1230             CONTINUE
 1240          CONTINUE
 1250       CONTINUE
C
         ELSE
C
            DO 1270 J=1,N
               IEND = J
               IF (KIN .EQ. 0) IEND = N
               DO 1260 I=1,IEND
                  A(I,J) = AM(I,J)
 1260          CONTINUE
 1270       CONTINUE
         ENDIF
      ENDIF
C
C
      IF (KIN .EQ. 10) GOTO 1800
      GOTO (1290, 1360, 1410, 1550) KIN+1
C                     -------   kin illegal   -------
      WRITE (LOUT,*) ' Subroutine ADEC aborted:',
     $     ' Illegal value of KIN ', KIN
      STOP
 1280 CONTINUE
      WRITE (LOUT,*) ' Subroutine ADEC aborted:',
     $     ' Illegal values MDISC, MODE ', MDISC, MODE
      STOP
C ---------------------------------------------------------------------
C                       ------   KIN = 0   ---------
 1290 CONTINUE
C
C  the lower left part
      NM = N + M
      DO 1310 I=1,M

         DO 1300 J=1,N
            A(N+I,J) = GP(I,J)
 1300    CONTINUE
 1310 CONTINUE
C
C  the upper right part
      NM = N + M
      DO 1330 I=1,N
         DO 1320 J=1,M
            A(I,N+J) = GP(J,I)
            IF (IUPDTE .EQ. 1) A(I,N+J) = A(I,N+J) - FL(I,J)
 1320    CONTINUE
 1330 CONTINUE
C
C  the lower right part
      DO 1350 J=1,M
         DO 1340 I=1,M
            A(N+I,N+J) =0.D0
 1340    CONTINUE
 1350 CONTINUE
      CALL DGEFA (A, NPGDIM, NM, IWK(IIPVT), IFDEC)
      IF (IFDEC .NE. 0) MXJOB(74) = 10
C
      RETURN
C
C ---------------------------------------------------------------------
C                       ------   KIN =  1  ---------
 1360 CONTINUE
      NM = N + M
      DO 1380 I=1,N
         DO 1370 J=1,M
            A(I,N+J) = GP(J,I)
 1370    CONTINUE
 1380 CONTINUE
      DO 1400 J=1,M
         DO 1390 I=1,J
            A(N+I,N+J) =0.D0
 1390    CONTINUE
 1400 CONTINUE
      CALL DSIFA (A, NPGDIM, NM, IWK(IIPVT), IFDEC)
      IF (IFDEC .NE. 0) MXJOB(74) = 1
      RETURN
C
C----------------------------------------------------------------------
C                       ------   KIN = 2   ---------
 1410 CONTINUE
      NM = N + M
C
      JOB = 0
      INFO = 0
C
      IF (MODE .EQ. 0) THEN
         IF (QPRTIM) CALL MONON(14)
         CALL DCHDC (A, NPGDIM, N, A(1,NM+1), IWK(IJPVTC), JOB, INFO)
         IF (QPRTIM) CALL MONOFF(14)
         IF (INFO .NE. N) THEN
            IFDEC = INFO
            MXJOB(74) = 2
            RETURN
         ENDIF
      ELSE IF (MODE .EQ. 1) THEN
         IF (NMRC .GT. 1) THEN
            DO 1420 KB=1,NBLK
               IF (QPRTIM) CALL MONON(14)
               CALL DCHDC (A((KB-1)*NMRC+1,(KB-1)*NMRC+1),
     $              NPGDIM, NMRC, A(1,NM+1), IWK((KB-1)*NMRC+1),
     $              JOB, INFO)
               IF (QPRTIM) CALL MONOFF(14)
               IF (INFO .NE. NMRC) THEN
                  IFDEC = INFO
                  MXJOB(74) = 3
                  RETURN
               ENDIF
 1420       CONTINUE
         ELSE IF (NMRC .EQ. 1) THEN
            DO 1430 I=1,N
               A(I,I) = SQRT (A(I,I))
 1430       CONTINUE
         ENDIF
      ENDIF
C
      DO 1450 I=1,N
         DO 1440 J=1,M
            A(I,N+J) = GP(J,I)
 1440    CONTINUE
 1450 CONTINUE
C
      IF (MODE .EQ. 0) THEN
         DO 1460 J=1,M
            IF (QPRTIM) CALL MONON(15)
            CALL DTRSL (A, NPGDIM, N, A(1,N+J), 11, IFDEC)
            IF (QPRTIM) CALL MONOFF(15)
            IF (IFDEC .NE. 0) THEN
               MXJOB(74) = 4
               RETURN
            ENDIF
 1460    CONTINUE
C
      ELSE IF (MODE .EQ. 1) THEN
         DO 1480 KB=1,NBLK
            DO 1470 J=1,M
               IF (QPRTIM) CALL MONON(15)
               CALL DTRSL (A((KB-1)*NMRC+1,(KB-1)*NMRC+1),
     $              NPGDIM, NMRC, A((KB-1)*NMRC+1,N+J), 11, IFDEC)
               IF (QPRTIM) CALL MONOFF(15)
               IF (IFDEC .NE. 0) THEN
                  MXJOB(74) = 5
                  RETURN
               ENDIF
 1470       CONTINUE
 1480    CONTINUE
      ENDIF
C
      DO 1500 J=1,M
         DO 1490 I=1,N
            A(N+J,I) = A(I,N+J)
 1490    CONTINUE
 1500 CONTINUE
C
      IF (QPRTIM) CALL MONON(16)
      CALL DQRDC (A(1,N+1), NPGDIM, N, M, A(1,NM+1), IWK(IJPVTQ),
     $     A(1,NM+2), JOB)
      IF (QPRTIM) CALL MONOFF(16)
C
      DO 1520 I=1,M
         DO 1510 J=I,M
            A(N+I,N+J) = A(I,N+J)
 1510    CONTINUE
 1520  CONTINUE
C
      DO 1540 I=1,N
         DO 1530 J=1,M
            A(I,N+J) = A(N+J,I)
 1530    CONTINUE
 1540 CONTINUE
C
      RETURN
C
C ---------------------------------------------------------------------
C                       ------   KIN = 3   ---------
 1550 CONTINUE
      NP1 = N + 1
      MP1 = M + 1
C
C  Store  K  and  J
C
      DO 1570 I=1,M
         DO 1560 J=1,N
            A(N+I,J) = GP(I,J)
 1560    CONTINUE
 1570 CONTINUE
C
C  Decompose  K
C
      CALL DGEFA (A(NP1,1), NPGDIM, M, IWK(IIPVTK), IFDEC)
      IF (IFDEC .NE. 0) THEN
         MXJOB(74) = 6
         RETURN
      ENDIF
C
C  Replace  J  by  K^(-1) * J
C
      DO 1580 J=MP1,N
         CALL DGESL (A(NP1,1), NPGDIM, M, IWK(IIPVTK), A(NP1,J), 0)
 1580 CONTINUE
C     
C  Compute  MM = J^T*AM*J
C
      IF (MODE .EQ. 0) THEN
C
         DO 1610 I=1,M
            DO 1600 J=MP1,N
               SUM = 0.D0
               DO 1590 K=1,M
                  SUM = SUM + AM(I,K)*A(N+K,J)
 1590          CONTINUE
               A(I,J) = SUM
 1600       CONTINUE
 1610    CONTINUE
C
      ELSE IF (MODE .EQ. 1) THEN
C
         DO 1670 IB=1,(M+1)/NMRC
            DO 1650 I=1,NMRC
               IF ((IB-1)*NMRC+I .GT. M) GOTO 1660
               DO 1640 J=MP1,N
                  SUM = 0.D0
                  DO 1620 K=1,NMRC
                     IF ((IB-1)*NMRC+K .GT. M) GOTO 1630
                     SUM = SUM + AM(I+(IB-1)*NMRC,K)*
     $                    A(N+(IB-1)*NMRC+K,J)
 1620             CONTINUE
 1630             CONTINUE
                  A(I+(IB-1)*NMRC,J) = SUM
 1640          CONTINUE
 1650       CONTINUE
 1660       CONTINUE
 1670    CONTINUE
C
      ENDIF
C
      DO 1700 I=MP1,N
         DO 1690 J=MP1,N
            SUM = 0.D0
            DO 1680 K=1,M
               SUM = SUM + A(N+K,I)*A(K,J)
 1680       CONTINUE
            A(I,J) = SUM
 1690    CONTINUE
 1700 CONTINUE
C
C  Decompose  MM
C
      NF = N - M
      CALL DGEFA (A(MP1,MP1), NPGDIM, NF, IWK(IIPVTM), IFDEC)
      IF (IFDEC .NE. 0) THEN
         MXJOB(74) = 7
         RETURN
      ENDIF
C
C  Store  AM
C
      IF (MODE .EQ. 0) THEN
         DO 1720 I=1,M
            DO 1710 J=1,M
               A(I,J) = AM(I,J)
 1710       CONTINUE
 1720    CONTINUE
C
      ELSE IF (MODE .EQ. 1) THEN
C
         DO 1740 I=1,M
            DO 1730 J=1,M
               A(I,J) = 0.D0
 1730       CONTINUE
 1740    CONTINUE
C
         DO 1790 IB=1,(M+1)/NMRC
            DO 1770 I=1,NMRC
               IF ((IB-1)*NMRC+I .GT. M) GOTO 1780
               DO 1750 J=1,NMRC
                  IF ((IB-1)*NMRC+J .GT. M) GOTO 1760
                  A(I+(IB-1)*NMRC,J+(IB-1)*NMRC) = AM(I+(IB-1)*NMRC,J)
 1750          CONTINUE
 1760          CONTINUE
 1770       CONTINUE
 1780       CONTINUE
 1790    CONTINUE
C
      ENDIF
C
      RETURN
C
C----------------------------------------------------------------------
C                       ------   KIN = 10   ---------
 1800 CONTINUE
C
C*    Start of decomposition
      QFACT = MXJOB(102) .EQ. 1 .OR. QFACTH
      NZG = MXJOB(101)
      NP = NMRC*NBLK
      NM = NP + M
      NZACT = NMRC*NMRC*NBLK + 2*NZG
C
 1810 CONTINUE
      CALL MCOPY (NMRC, NBLK, LDAM, NZG, NZACT, AM, GP, IWK(IIROWG),
     $     IWK(IJCOLG), A(LP1,1), IWK(IIROW), IWK(IJCOL))
      IF (QFACT .OR. QREFCT) THEN
C
         DO 1820 INZ=0,NZACT-1
            IWK(IIRN+INZ) = IWK(IIROW+INZ)
            IWK(IICN+INZ) = IWK(IJCOL+INZ)
 1820    CONTINUE
C
      ENDIF
C
      IF (QFACT .OR. QREFCT) THEN
         IF (IMON .GE. 4) THEN
            WRITE (LOUT,*) ' Call MA28A'
            WRITE (LOUT,*) ' ... NNZ(A) is:', NZACT
         ENDIF
         IF (QPRTIM) CALL MONON(14)
         CALL MA28AD (NM, NZACT, A(LP1,1), LICN, IWK(IIRN), LIRN,
     $        IWK(IICN), THRSH1, IWK(IIKEEP), IWK(IIWMA), A, IFDEC)
         IF (QPRTIM) CALL MONOFF(14)
         IF (IFDEC .LT. 0) THEN
            MXJOB(74) = 8
            CALL MADIAG (0, IFDEC, IMON .GE. 1, LOUT)
            IF (IFDEC .LT. 0) RETURN
         ELSE IF (IFDEC .GT. 0) THEN
            CALL MADIAG (0, IFDEC, IMON .GE. 2, LOUT)
         ENDIF
         MXJOB(95) = MAX (MINIRN, MXJOB(95))
         MXJOB(96) = MAX (MINICN, MXJOB(96))
         MXJOB(97) = MAX (IRNCP, MXJOB(97))
         MXJOB(98) = MAX (ICNCP, MXJOB(98))
         IF (IMON .GE. 4) THEN
            WRITE (LOUT,*) ' ... MINIRN/MINICN:', MINIRN, MINICN
            WRITE (LOUT,*) ' ... IRNCP/ICNCP  :', IRNCP, ICNCP  
            WRITE (LOUT,*) ' ... W(1) is:', A(1,1)
         ENDIF
         QREFCT = .FALSE.
C
      ELSE
C
         IF (IMON .GE. 5) WRITE (LOUT,*) ' Call MA28B'
         IF (QPRTIM) CALL MONON(15)
         CALL MA28BD (NM, NZACT, A(LP1,1), LICN, IWK(IIROW), IWK(IJCOL),
     $        IWK(IICN), IWK(IIKEEP), IWK(IIWMA), A, IFDEC)
         IF (QPRTIM) CALL MONOFF(15)
         IF (IFDEC .LT. 0) THEN
            CALL MADIAG (1, IFDEC, IMON .GE. 1, LOUT)
            IF (IFDEC .LT. 0) THEN
               QREFCT = .TRUE.
               IF (IMON .GE. 2) WRITE(LOUT, *)
     $              ' Analyze factor because of IFAIL =', IFDEC
               GOTO 1810
            ENDIF
         ELSE IF (IFDEC .GT. 0) THEN
            CALL MADIAG (1, IFDEC, IMON .GE. 2, LOUT)
         ENDIF
C
         IF (IMON .GE. 5) THEN
            WRITE (LOUT,*) ' ... RMIN is:', RMIN
            WRITE (LOUT,*) ' ... W(1) is:', A(1,1)
         ENDIF
C
         IF (RMIN .LT. THRSH2) THEN
            QREFCT = .TRUE.
            IF (IMON .GE. 4) WRITE (LOUT, *)
     $           ' Analyze factor because of RMIN =',
     $           RMIN, ' RW(1) =', A(1,1)
            GOTO 1810
         ENDIF
      ENDIF
C     
      RETURN
C
C* End of ADEC
C
C
C
      ENTRY ASOL (N, M, B, MXJOB, LDA, A, IWK)
C
C
C ----------------------------------------------------------------------
C
C     This subroutine computes the solution of a linear system
C     with matrix decomposed by subroutine adec .
C
C     Version of July 31, 1991.
C
C ----------------------------------------------------------------------
C
      MXJOB(59) = MXJOB(59) + 1
      IF (KIN .EQ. 10) GOTO 1980
      GOTO (1830, 1840, 1850, 1910) KIN+1
C                     -------   KIN illegal   -------
      WRITE (LOUT,*) ' Subroutine ASOL: Illegal value of KIN ', KIN
      STOP
C ---------------------------------------------------------------------
C                       ------   KIN = 0   ---------
 1830 CONTINUE
      NM = N + M

      CALL DGESL (A, NPGDIM, NM, IWK(IIPVT), B, 0)
      RETURN
C ---------------------------------------------------------------------
C                       ------   KIN = 1   ---------
 1840 CONTINUE
      NM = N + M
      CALL DSISL (A, NPGDIM, NM, IWK(IIPVT), B)
      RETURN
C ---------------------------------------------------------------------
C                       ------   KIN = 2   ---------
 1850 CONTINUE
      NM = N + M
      B(1) = B(1)/A(1,1)
      DO 1870 K=2,NM
         IF (MODE .EQ. 0 .OR. NMRC .GT. 1 .OR. K .GT. N) THEN
            DO 1860 I=1,K-1
               B(K) = B(K) - A(I,K)*B(I)
 1860       CONTINUE
         ENDIF
         B(K) = B(K)/A(K,K)
 1870 CONTINUE
C
      DO 1880 K=N+1,NM
         B(K) = -B(K)
 1880 CONTINUE
C
      B(NM) = B(NM)/A(NM,NM)
      DO 1900 K=NM-1,1,-1
         DO 1890 I=NM,K+1,-1
            B(K) = B(K) - A(K,I)*B(I)
 1890    CONTINUE
         B(K) = B(K)/A(K,K)
 1900 CONTINUE
C
      RETURN
C ---------------------------------------------------------------------
C                       ------   KIN = 3   ---------
 1910  CONTINUE
      NP1 = N + 1
      MP1 = M + 1
      CALL DGESL (A(NP1,1), NPGDIM, M, IWK(IIPVTK), B(NP1), 0)
      DO 1930 I=1,M
         SUM = -B(I)
         DO 1920 J=1,M
            SUM = SUM + A(I,J)*B(N+J)
 1920     CONTINUE
         B(I) = SUM
 1930  CONTINUE
C
      DO 1950 I=MP1,N
         SUM = 0.D0
         DO 1940 J=1,M
            SUM = SUM + A(N+J,I)*B(J)
 1940     CONTINUE
         B(I) = SUM
 1950  CONTINUE
      NF = N - M
C
      CALL DGESL (A(MP1,MP1), NPGDIM, NF, IWK(IIPVTM), B(MP1), 0)
C
      DO 1970 I=1,M
         SUM = B(N+I)
         DO 1960 J=MP1,N
            SUM = SUM - A(N+I,J)*B(J)
 1960     CONTINUE
         B(I) = SUM
 1970  CONTINUE
      RETURN
C
C ---------------------------------------------------------------------
C                       ------   KIN = 10   ---------
 1980 CONTINUE
      JOB = 1
      CALL MA28CD (NM, A(LP1,1), LICN, IWK(IICN), IWK(IIKEEP), B,
     $     A, JOB)
C
      RETURN
C
C* End of ASOL
C
      END
      SUBROUTINE MADIAG (ISUB, IFLAG, QPRINT, LUNIT)
C
C  Error diagnostics for MA28AD and MA28BD
C
C
C* Parameters:
C  -----------
C
C    ISUB      I   IN   Subroutine called before MADIAG
C                       = 0      MA28AD
C                       = 1      MA28BD
C  * IFLAG     I   IN   Return code from ISUB
C                  OUT  .LT. 0   Error exit
C                       .EQ. 0   Continue
C                       .GT. 0   Reduce stepsize
C    QPRINT    L   IN   .TRUE.   Print diagnostic messages
C                       .FALSE.  Don't print anything
C    LUNIT     I   IN   Logical unit for print output
C
C* File              mexxlin.f
C* Version           1.0
C* Latest Change     96/05/08 (1.6)
C
C* Type declaration:
C  -----------------
C
      IMPLICIT LOGICAL (A-P,R-Z), REAL (Q)
C
      INTEGER IFLAG, ICNCP, IRANK, IRNCP, LUNIT, MINICN, MINIRN
      INTEGER ISUB
C
      DOUBLE PRECISION EPS, RESID, RMIN
C
      LOGICAL QABRT1, QABRT2, QPRINT
C
      COMMON /MA28FD/ EPS, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
     $     IRANK, QABRT1, QABRT2
C
C     More than one non-zero in the same position is no error
      IF (IFLAG .EQ. -14) IFLAG = 0
C
      IF (IFLAG .LE. -3) THEN
C
C---     Analyze error situations
         IF (QPRINT) WRITE (LUNIT,'(A,I4)') ' - Error in MA28 - ', IFLAG
C
         IF (IFLAG .GE. -13 .AND. IFLAG .LE. -8) THEN
            IF (QPRINT) WRITE (LUNIT,'(23X,A)')
     $           ' Input matrix given is incorrect.'
C
         ELSE IF (IFLAG .EQ. -6) THEN
            IF (QPRINT) WRITE (LUNIT,'(23X,A,I8)')
     $           ' LICN is too small. Try at least ', MINICN
            IF (QPRINT) WRITE (LUNIT,'(23X,A,I8)')
     $           ' LIRN is too small. Try at least ', MINIRN
C
         ELSE IF (IFLAG .EQ. -5) THEN
            IF (QPRINT) WRITE (LUNIT,'(23X,A,I8)')
     $           ' LICN is too small. Try at least ', MINICN
C
         ELSE IF (IFLAG .EQ. -4) THEN
            IF (QPRINT) WRITE (LUNIT,'(23X,A)')
     $           ' LICN is far too small.'
C
         ELSE IF (IFLAG .EQ. -3) THEN
            IF (QPRINT) WRITE (LUNIT,'(23X,A,I8)')
     $           ' LIRN is too small. Try at least ', MINIRN
         ENDIF
C
      ELSE IF (IFLAG .EQ. -2 .OR. IFLAG .EQ. -1) THEN
         IF (QPRINT) WRITE (LUNIT,'(A,I4,/,23X,A)')
     $        '   - Error in MA28 -   ', IFLAG,
     $        ' Matrix is numerically or structurally singular.'
C
      ELSE IF (IFLAG .EQ. +1 .AND. ISUB .EQ. 0) THEN
         IF (QPRINT) WRITE (LUNIT,'(2A)') ' - Warning from MA28 - ',
     $        ' Matrix seems to be structurally singular.'
         IFLAG = 0
C
      ELSE IF (IFLAG .EQ. +2 .AND. ISUB .EQ. 0) THEN
         IF (QPRINT) WRITE (LUNIT,'(2A)') ' - Warning from MA28 - ',
     $        ' Matrix seems to be numerically singular.'
         IFLAG = 0
C
      ELSE IF (IFLAG .GT. 0 .AND. ISUB .EQ. 1) THEN
         IF (QPRINT) WRITE (LUNIT,'(2A,I8)')
     $        ' - Warning from MA28 - ',
     $        ' Very small pivot in row ', IFLAG
         IFLAG = 0
      ENDIF
      RETURN
C
C* End of MADIAG
C
      END
      SUBROUTINE MCOPY (NMRC, NBLK, LDAM, NZG, NZA,
     $     AM, G, IROWG, JCOLG, AMAT, IROW, JCOL)
C
C     Copy block-diagonal matrix AM to sparse matrix AMAT
C
      INTEGER NBLK, NZG, NZA, IROWG(NZG), JCOLG(NZG), IROW(NZA),
     $     JCOL(NZA)
      DOUBLE PRECISION AM(LDAM,*), G(NZG), AMAT(NZA)
C
C* File              mexxlin.f
C* Version           1.0
C* Latest Change     96/05/08 (1.6)
C
C
      SAVE
C
      L = 0
      DO 1020 K=1,NBLK
         DO 1010 J=1,NMRC
            DO 1000 I=1,NMRC
               L = L + 1
               AMAT(L) = AM(I+(K-1)*NMRC,J)
               IROW(L) = (K-1)*NMRC + I
               JCOL(L) = (K-1)*NMRC + J
 1000       CONTINUE
 1010    CONTINUE
 1020 CONTINUE
C
      DO 1030 LG=1,NZG
         L = L + 1
         AMAT(L) = G(LG)
         IROW(L) = NMRC*NBLK + IROWG(LG)
         JCOL(L) = JCOLG(LG)
C
         L = L + 1
         AMAT(L) = G(LG)
         IROW(L) = JCOLG(LG)
         JCOL(L) = NMRC*NBLK + IROWG(LG)
 1030 CONTINUE
C
C
      RETURN
C
C* End of MCOPY
C
      END
      SUBROUTINE RECDEC (NB,N,LDAM,MDAM,LDGP,MDGP,LDA,MDA,
     &     AM,GP,A,NZ,NF,IPT,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NZ(NB),NF(NB),IPT(NB+1)
      DOUBLE PRECISION A(LDA,MDA),AM(LDAM,MDAM),GP(LDGP,MDGP)
C
C  ---------------------------------------------------------------------
C
C     NB              NUMBER OF BODIES
C
C     N               DIMENSION OF POSITION VECTOR OF BODY "K" (SAME FOR
C                     ALL K)
C
C     NZ              NZ(NB)
C                     NZ(K) CONTAINS THE NUMBER OF CONSTRAINTS DUE TO
C                     IMPLICIT JOINT "K". ( 0 .LE. NZ(K) .LE. N )
C                     UNCHANGED ON EXIT
C
C     NF              NF(NB)
C                     DESCRIBES THE TOPOLOGY OF THE SYSTEM. 
C                     NF(K)=J MEANS J=FATHER(K). WE REQUIRE J.LE.K.
C                     UNCHANGED ON EXIT
C
C     AM              AM(LDAM,MDAM)
C                     HOLDS THE MASS MATRIX IN NB BLOCKS OF SIZE N*N.
C                     UNCHANGED ON EXIT
C
C     LDAM, MDAM      DIMENSIONS OF ARRAY AM
C                     LDAM .GE. NB*N
C                     MDAM .GE. N
C
C     GP              GP(LDGP,MDGP)
C                     HOLDS THE VELOCITY CONSTRAINT MATRIX GP. IT IS
C                     STORED  COLUMNWISE IN NB BLOCKS OF SIZE
C                     NZ(K) TIMES 2*N EACH (CORRESPONDING TO JOINT "K").
C                     THE K-TH BLOCK MUST HAVE THE FORM 
C                     GP_K = [ G_K, C_K ], WHERE G_K CORRESPONDS TO THE 
C                     POSITION VARIABLES P_K OF BODY "K" AND C_K TO P_J 
C                     OF BODY "J" WHERE J=FATHER(K).
C                     UNCHANGED ON EXIT
C
C     LDGP, MDGP      DIMENSIONS OF ARRAY GP
C                     LDGP .GE. NZT := NZ(1)+ ... + NZ(NB)
C                     MDGP .GE. 2*N
C
C     A               A(LDA,MDA)
C                     DECOMPOSED MATRIX 
C
C     LDA,  MDA       DIMENSIONS OF ARRAY AM
C                     LDA .GE. NZT+(NB*N)
C                     MDA .GE. N+NZMAX+1 WHERE
C                     NZMAX := MAX(NZ(1),...,NZ(NB)) IS THE MAXIMUM
C                     NUMBER OF CONSTRAINTS DUE TO A JOINT.
C
C     IPT             IPT(NB+1)
C                     INTEGER WORK ARRAY
C
C
C     THE (MAXIMUM) REQUIRED WORKSPACE FOR MATRICES AM, GP, AND A
C     MAY BE ESTIMATED AS FOLLOWS:
C
C     LDAM = NB*N,                      MDAM = N
C     LDGP = NZT       .LE.   NB*N      MDGP = 2*N
C     LDA  = NB*N+NZT  .LE. 2*NB*N      MDA  = N+NZMAX+1 .LE. 2*N+1
C
C     HENCE:   LRWK_MIN = 2*NB*N*N + 3*NZT*N + (NB*N+NZT)*(NZMAX+1)
C              LRWK_MAX = 7*NB*N*N+2*NB*N
C
C
C     FUNCTIONS AND SUBROUTINES USED: LINPACK DCHDC, DQRDC
C
C
C ----------------------------------------------------------------------
C     DATE:   3-11-1992 
C ----------------------------------------------------------------------
C
C --- Preparations
C
      IERR=0
      NZMAX=0
      IPT(1)=0
      DO 1000 I=1,NB
         IPT(I+1)=IPT(I)+NZ(I)
         NZMAX=MAX(NZ(I),NZMAX)
 1000 CONTINUE 
C
      NZT=IPT(NB+1)
      ND=N+NZMAX
C
C --- Copy M_k to upper left block of A_k
C
      DO 1030 K=1,NB
         J2=N*(K-1)
         J1=J2+IPT(K)
         DO 1020 J=1,N
            DO 1010 I=1,N
               A(J1+I,J)=AM(J2+I,J)
 1010       CONTINUE
 1020    CONTINUE
 1030 CONTINUE
C
      JOB=0
      INFO=0
C
      DO 1250 K=NB,1,-1
C
         NZK=NZ(K)
         J1=N*(K-1)+IPT(K)
         J2=J1+N
         JGP=IPT(K)
C
C     *  Choleski-decomposition of M_k
C
         CALL DCHDC (A(J1+1,1),LDA,N,A(J1+1,ND+1),IDUM,JOB,INFO)
         IF (INFO.NE.N) GOTO 1260
C
C     *  L_k^-1 * G_k^T
C
         DO 1050 J=1,NZK
            DO  1040 I=1,N
               A(J1+I,N+J)=GP(JGP+J,I)
 1040       CONTINUE
 1050    CONTINUE
C
         DO 1080 J=N+1,N+NZK
            DO 1070 I=J1+1,J1+N
               DO 1060 L=J1+1,I-1
                  A(I,J)=A(I,J)-A(L,I-J1)*A(L,J)
 1060          CONTINUE
               A(I,J)=A(I,J)/A(I,I-J1)
 1070       CONTINUE
 1080    CONTINUE
C
C     *  QR-decomposition
C
         DO 1100 J=1,NZK
            DO 1090 I=1,N
               A(J2+J,I)=A(J1+I,N+J)
 1090       CONTINUE
 1100    CONTINUE
C
         CALL DQRDC (A(J1+1,N+1),LDA,N,NZK,A(J1+1,ND+1),IDUM,RDUM,JOB)
C     
         DO 1120 J=N+1,N+NZK
            DO 1110 I=1,J-N
               A(J2+I,J)=A(J1+I,J)
 1110       CONTINUE
            IF (A(J1+J,J).EQ.0.D0) GOTO 1270
 1120    CONTINUE
C
         DO 1140 J=1,NZK
            DO 1130 I=1,N
               A(J1+I,N+J)=A(J2+J,I)
 1130       CONTINUE
 1140    CONTINUE
C
C      * R_k^-T * C_k
C
         KF=NF(K)
         IF (KF.NE.0) THEN
            DO 1160 J=1,N
               DO 1150 I=1,NZK
                  A(J2+I,J)=GP(JGP+I,N+J)
 1150          CONTINUE
 1160       CONTINUE
C
            DO 1190 J=1,N
               DO 1180 I=J2+1,J2+NZK
                  DO 1170 L=J2+1,I-1
                     A(I,J)=A(I,J)-A(L,J)*A(L,I-J1)
 1170             CONTINUE
                  A(I,J)=A(I,J)/A(I,I-J1)
 1180          CONTINUE
 1190       CONTINUE
C
C        *  Update M_kf
C
            JF=N*(KF-1)+IPT(KF)
            DO 1220 J=1,N
               DO 1210 I=1,J
                  SUM=0
                  DO 1200 L=J2+1,J2+NZK
                     SUM=SUM+A(L,I)*A(L,J)
 1200             CONTINUE
                  A(JF+I,J)=A(JF+I,J)+SUM
 1210          CONTINUE
 1220       CONTINUE
            DO 1240 J=2,N
               DO 1230 I=J,N
                  A(JF+I,J)=A(JF+J,I)
 1230          CONTINUE
 1240       CONTINUE
C
         ENDIF
C
 1250 CONTINUE
      RETURN
C
 1260 CONTINUE
C     WRITE (*,*) '# Error in RECDEC in equation of joint ',K
C     WRITE (*,*) '# Error return from DCHDC, return code INFO=',INFO
      IERR=INFO
      GOTO 1280
 1270 CONTINUE
C      WRITE (*,*) '# Error in RECDEC in equation of joint ',K
C      WRITE (*,*) '# Zero diagonal element in QR factor R_k: ',N-J
      IERR=N-J
      GOTO 1280
 1280 CONTINUE
      RETURN
      END
      SUBROUTINE RECSOL (NB,N,LDA,MDA,A,B,NZ,NF,IPT,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NZ(NB),NF(NB),IPT(NB+1)
      DOUBLE PRECISION A(LDA,MDA),B(LDA)
C
C  ---------------------------------------------------------------------
C
C     NB              NUMBER OF BODIES
C
C     N               DIMENSION OF POSITION VECTOR OF BODY "K"
C
C     NZ              NZ(NB)
C                     NZ(K) CONTAINS THE NUMBER OF CONSTRAINTS DUE TO
C                     IMPLICIT JOINT "K". ( 0 .LE. NZ(K) .LE. N )
C                     UNCHANGED ON EXIT
C
C     NF              NF(NB)
C                     NF(K)=J MEANS J=FATHER(K). WE REQUIRE J.LE.K
C                     UNCHANGED ON EXIT
C
C     A               A(LDA,MDA)
C                     HOLDS THE MATRIX DECOMPOSED BY RECDEC
C
C     LDA,  MDA       DIMENSIONS OF ARRAY A
C
C     B               B(LDA)
C                     RIGHT HAND SIDE VECTOR. HOLDS [F_K, G_K]
C                     OVERWRITTEN BY SOLUTION-VECTOR [W_K,LAM_K]
C
C     IPT             IPT(NB+1)
C                     INTEGER WORK ARRAY
C
C
C ----------------------------------------------------------------------
C     DATE:   3-12-1992 
C ----------------------------------------------------------------------
C
C --- Preparations
C
      NZMAX=0
      IPT(1)=0
      DO 1000 I=1,NB
         IPT(I+1)=IPT(I)+NZ(I)
         NZMAX=MAX(NZ(I),NZMAX)
 1000 CONTINUE
      NZT=IPT(NB+1)
      ND=N+NZMAX
      NNB=N*NB
      IF (NNB+NZT.GT.LDA .OR. ND.GT.MDA) GOTO 1190
C
C --- Forward substitution
C
      DO 1090 K=NB,1,-1
C
         NZK=NZ(K)
         JF=N*(K-1)
         JG=NNB+IPT(K)
         J1=JF+IPT(K)
         J2=J1+N
C
C      * f_k := L_k^-1 * f_k
C
          B(JF+1)=B(JF+1)/A(J1+1,1)
          DO 1020 J=2,N
             TEMP1=0
             DO 1010 I=1,J-1
                TEMP1=TEMP1+A(J1+I,J)*B(JF+I)
 1010        CONTINUE
             B(JF+J)=(B(JF+J)-TEMP1)/A(J1+J,J)
 1020     CONTINUE
C
C      * g_k := G_k * f_k - g_k 
C
         DO 1040 J=1,NZK
            TEMP2=0
            DO 1030 I=1,N
               TEMP2=TEMP2+A(J1+I,N+J)*B(JF+I)
 1030       CONTINUE
            B(JG+J)=TEMP2-B(JG+J)
 1040    CONTINUE
C
C      * g_k := R_k^-T * g_k
C
         B(JG+1)=B(JG+1)/A(J2+1,N+1)
         DO 1060 J=2,NZK
            TEMP3=0
            DO 1050 I=1,J-1
               TEMP3=TEMP3+A(J2+I,N+J)*B(JG+I)
 1050       CONTINUE
            B(JG+J)=(B(JG+J)-TEMP3)/A(J2+J,N+J)
 1060    CONTINUE
C
         KF=NF(K)
         IF (KF.NE.0) THEN
C
C         * Compute f_kf := f_kf - C_k^T g_k
C
            JK=N*(KF-1)
            DO 1080  J=1,N
               TEMP4=0
               DO 1070 I=1,NZK
                  TEMP4=TEMP4+A(J2+I,J)*B(JG+I)
 1070          CONTINUE
               B(JK+J)=B(JK+J)-TEMP4
 1080       CONTINUE
         ENDIF
C
 1090 CONTINUE
C
C --- Backward substitution
C
      DO 1180 K=1,NB
C
         NZK=NZ(K)
         JF=N*(K-1)
         JG=NNB+IPT(K)
         J1=JF+IPT(K)
         J2=J1+N
C
         KF=NF(K)
         IF (KF.NE.0) THEN
C
C        * g_k := g_k + C_k w_fk
C
            JK=N*(KF-1)
            DO 1110 I=J2+1,J2+NZK
               TEMP5=0
               DO 1100 J=1,N
                  TEMP5=TEMP5+A(I,J)*B(JK+J)
 1100          CONTINUE
               B(I-J2+JG)=B(I-J2+JG)+TEMP5
 1110       CONTINUE
         ENDIF
C
C        * g_k := R_k^-1 * g_k
C
         B(JG+NZK)=B(JG+NZK)/A(J2+NZK,N+NZK)
         DO 1130 J=NZK-1,1,-1
            TEMP6=-B(JG+J+1)
            DO 1120 I=1,J
               B(JG+I)=B(JG+I)+A(J2+I,N+J+1)*TEMP6
 1120       CONTINUE
            B(JG+J)=B(JG+J)/A(J2+J,N+J)
 1130    CONTINUE
C
C        * f_k := f_k - G_k^T * lambda_k
C
         DO 1150 J=1,NZK
            TEMP7=-B(JG+J)
            DO 1140 I=1,N
               B(JF+I)=B(JF+I)+TEMP7*A(J1+I,N+J)
 1140       CONTINUE 
 1150    CONTINUE
C
C        * f_k := L_k^-T * f_k
C
         DO 1170 J=N,1,-1
            TEMP8=0
            DO 1160 I=J+1,N
               TEMP8=TEMP8+A(J1+J,I)*B(JF+I)
 1160       CONTINUE
            B(JF+J)=(B(JF+J)-TEMP8)/A(J1+J,J)
 1170    CONTINUE
C
 1180 CONTINUE
C
      RETURN
C
C --- Error exists
C
 1190 WRITE (*,*) '# ERROR IN RECSOL - WRONG DIMENSION OF MATRIX A'
      WRITE (*,*) '# OR INCORRECT VALUES FOR N, NZ ARE GIVEN.'
      IERR = -1
      GOTO 1200
C
C --- Fail exit
C
 1200 CONTINUE
      RETURN
      END
C
      SUBROUTINE MXDENS (NP, NV, NU, NL, T, P0, V0, U0, A0, R0,
     $                  PDOT0, VDOT0, UDOT0, PDOT1, VDOT1, UDOT1,
     $                  H, KC, KDOMAX, KEXMAX, TAB, NJ,
     $                  P1, V1, U1, A1, R1, ILAM, IRHO,
     $                  IDOFLG,N6UDIM, LSAFE, LDENS,
     $                  LDYS, ND, NDP, NDV, NDU, NDA, NDR,
     $                  ITOL, NTODIM, RTOL, ATOL, N2UDIM, SCAL, YWRK, 
     $                  INDP, INDV, INDU, INDA, INDR, YSAFE, DENS,
     $                  MXJOB, IDOCNT)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER IDOCNT, IDOFLG
      INTEGER NP, NV, NU, NL
      INTEGER KC, KDOMAX, KEXMAX, ILAM, IRHO
      INTEGER N6UDIM, NTODIM, N2UDIM, LSAFE, LDENS
      INTEGER ND, NDP, NDV, NDU, NDA, NDR
      INTEGER ITOL, LUINT, MNINT
      INTEGER INDP, INDV, INDA, INDU, INDR, NJ
C
      DOUBLE PRECISION T, H
      DOUBLE PRECISION P0, V0, U0, A0, R0, PDOT0, VDOT0, UDOT0
      DOUBLE PRECISION P1, V1, U1, A1, R1, PDOT1, VDOT1, UDOT1
      DOUBLE PRECISION RTOL, ATOL, SCAL, TAB
      DOUBLE PRECISION DENS, YSAFE, YWRK
C
      DIMENSION INDP(NDP), INDV(NDV), INDU(NDU), INDA(NDA), INDR(NDR)
      DIMENSION NJ(KEXMAX), MXJOB(150)
C
      DIMENSION P0(NP), V0(NV), U0(NU), A0(NV), R0(NL)
      DIMENSION PDOT0(NP), VDOT0(NV), UDOT0(NU)
      DIMENSION P1(NP), V1(NV), U1(NU), A1(NV), R1(NL)
      DIMENSION PDOT1(NP), VDOT1(NV), UDOT1(NU)
      DIMENSION RTOL(NTODIM), ATOL(NTODIM), SCAL(N2UDIM)
      DIMENSION TAB(N6UDIM,KDOMAX)
      DIMENSION YSAFE(LDYS,LSAFE), DENS(ND,LDENS), YWRK(ND)
C 
C  Dense output 
C
C  KEXMAX : max order for extrapolation
C  KDOMAX : max order for dense output
C  KC : order of last accepted extrapolation step
C
C
C* File              mxdens.f
C* Version           1.1
C* Latest Change     96/03/14 (1.7)
C* Modification history
C    1.1             Prevent erroneous call to MXHERM if NDY = ND
C                    or NDY = 0
C      
C
      INTEGER KCMAX1
      PARAMETER (KCMAX1 = 21)
C
      INTEGER I, IL, IR, IDEL, NS, M, J, MAXLR, L0, L1, L, IEN, ILAMH,
     $     KMAX
      INTEGER NDY, NDZ, NN, NND
      INTEGER IPDIF
      DOUBLE PRECISION THERR, ERRFAC, ERR, ERRP, ERRV, ERRU, ERRPRI
      DOUBLE PRECISION WT, FAC
      DOUBLE PRECISION FACUL
      DIMENSION IPDIF(KCMAX1), FACUL(KCMAX1)
C
      EXTERNAL MXHERM
C
      SAVE FACUL, IPDIF
C
C     INITIAL PREPARATIONS
C
C
      MNINT = MXJOB(13)
      LUINT = MXJOB(14)
C
      MODDEN = MXJOB(25)
      MU     = MXJOB(27)
C
C
      IF (IDOCNT .EQ. 1) THEN
         IPDIF(1) = 0
         DO 1000 I=1,KEXMAX
            IPDIF(I+1) = IPDIF(I) + NJ(I)
 1000    CONTINUE
C
         FACUL(1) = 1.D0
         DO 1010 I=1,KEXMAX
            FACUL(I+1) = FACUL(I)*DBLE(I)
 1010    CONTINUE
      ENDIF
C
      NDY = NDP + NDV + NDU
      NDZ = NDA + NDR
C
      KMAX = MIN (KDOMAX, KC)
C
C  Number of left/right derivatives
C
C    ILAM ... Number of left derivatives
C    IRHO ... Number of right derivatives
C
C    Assume: ILAM .GE. 0, IRHO .GE. 1, IRHO .GE. ILAM
C
C  - UNPERTURBED CASE: (MODDEN .EQ. 0)
C
C      ILAM = (KMAX+MU-1)/2,
C      IRHO = (KMAX+MU)/2,   
C
C      MU = 0/1: DETERMINES THE DEGREE OF INTERPOLATION
C
C  - PERTURBED CASE (MODDEN .EQ. 1)
C
C      ILAM = 0
C      IRHO = KMAX-1
C
C      MU = 0/1: DETERMINES THE NUMBER OF EXTRAPOLATION STEPS
C                USED FOR COMPUTING THE DERIVATIVES
C
C*
      IDR = MODDEN
C*
      IF (IDR .EQ. 0) THEN
         ILAM = (KMAX + MU - 1)/2
         IRHO = (KMAX + MU)/2
      ELSE
         ILAM = 0
         IRHO = KMAX-1
      ENDIF
C
      IL = ILAM + 1
      IR = IL + 1
C
      IF (IDOFLG .EQ. 1) THEN
      NS = 0
      DO 1020 I=1,NDP
         DENS(NS+I,IL) = P0(INDP(I))
         DENS(NS+I,IL+1) = P1(INDP(I))
         DENS(NS+I,IL+2) = H*PDOT1(INDP(I))
 1020 CONTINUE
      NS = NS + NDP
      DO 1030 I=1,NDV
         DENS(NS+I,IL) = V0(INDV(I))
         DENS(NS+I,IL+1) = V1(INDV(I))
         DENS(NS+I,IL+2) = H*VDOT1(INDV(I))
 1030 CONTINUE
      NS = NS + NDV
      DO 1040 I=1,NDU
         DENS(NS+I,IL) = U0(INDU(I))
         DENS(NS+I,IL+1) = U1(INDU(I))
         DENS(NS+I,IL+2) = H*UDOT1(INDU(I))
 1040 CONTINUE
      NS = NS + NDU
      DO 1050 I=1,NDA
         DENS(NS+I,IL) = A0(INDA(I))
         DENS(NS+I,IL+1) = A1(INDA(I))
 1050 CONTINUE
      NS = NS + NDA
      DO 1060 I=1,NDR
         DENS(NS+I,IL) = R0(INDR(I))
         DENS(NS+I,IL+1) = R1(INDR(I))
 1060 CONTINUE
C     
      IF (ILAM .GE. 1) THEN
         NS = 0
         DO 1070 I=1,NDP
            DENS(NS+I,IL-1) = H*PDOT0(INDP(I))
 1070    CONTINUE
         NS = NS + NDP
         DO 1080 I=1,NDV
            DENS(NS+I,IL-1) = H*VDOT0(INDV(I))
 1080    CONTINUE
         NS = NS + NDV
         DO 1090 I=1,NDU
            DENS(NS+I,IL-1) = H*UDOT0(INDU(I))
 1090    CONTINUE
      ENDIF
      ENDIF
C
      IF (IDOFLG .EQ. 2) THEN
      NS = 0
      DO 1100 I=1,NDP
         DENS(NS+I,IL) = P0(I)
         DENS(NS+I,IL+1) = P1(I)
         DENS(NS+I,IL+2) = H*PDOT1(I)
 1100 CONTINUE
      NS = NS + NDP
      DO 1110 I=1,NDV
         DENS(NS+I,IL) = V0(I)
         DENS(NS+I,IL+1) = V1(I)
         DENS(NS+I,IL+2) = H*VDOT1(I)
 1110 CONTINUE
      NS = NS + NDV
      DO 1120 I=1,NDU
         DENS(NS+I,IL) = U0(I)
         DENS(NS+I,IL+1) = U1(I)
         DENS(NS+I,IL+2) = H*UDOT1(I)
 1120 CONTINUE
      NS = NS + NDU
      DO 1130 I=1,NDA
         DENS(NS+I,IL) = A0(I)
         DENS(NS+I,IL+1) = A1(I)
 1130 CONTINUE
      NS = NS + NDA
      DO 1140 I=1,NDR
         DENS(NS+I,IL) = R0(I)
         DENS(NS+I,IL+1) = R1(I)
 1140 CONTINUE
C     
      IF (ILAM .GE. 1) THEN
         NS = 0
         DO 1150 I=1,NDP
            DENS(NS+I,IL-1) = H*PDOT0(I)
 1150    CONTINUE
         NS = NS + NDP
         DO 1160 I=1,NDV
            DENS(NS+I,IL-1) = H*VDOT0(I)
 1160    CONTINUE
         NS = NS + NDV
         DO 1170 I=1,NDU
            DENS(NS+I,IL-1) = H*UDOT0(I)
 1170    CONTINUE
      ENDIF 
      ENDIF
C
C
C --- Unperturbed case ---
C
C
      IF (IDR .EQ. 0) THEN
C
         MAXLR = MAX (ILAM-1, IRHO)
         DO 1350 M=1,MAXLR
C
C           compute differences 
C
            DO 1200 J=M,KMAX
               L1 = IPDIF(J+1)
               L0 = IPDIF(J) + M + 1
               DO 1190 L=L1,L0,-1
                  DO 1180 I=1,ND
                     YSAFE(I,L) = YSAFE(I,L) - YSAFE(I,L-1)
 1180             CONTINUE
 1190          CONTINUE               
 1200       CONTINUE
C
            IF (M .LE. ILAM-1) THEN
C     
C              extrapolation of left differences
C
               DO 1220 J=M,KMAX
                  FAC = (DBLE(NJ(J))**M)/FACUL(M+1)
                  IEN = IPDIF(J) + M + 1
                  DO 1210 I=1,ND
                     TAB(I,J) = YSAFE(I,IEN)*FAC
 1210             CONTINUE
 1220          CONTINUE               
C
               DO 1250 J=M+1,KMAX
                  DO 1240 L=J,M+1,-1
                     FAC = DBLE(NJ(J))/DBLE(NJ(L-1)) - 1.D0
                     DO 1230 I=1,ND
                        TAB(I,L-1) = TAB(I,L) + (TAB(I,L) - TAB(I,L-1))/
     $                       FAC
 1230                CONTINUE
 1240             CONTINUE
 1250          CONTINUE
C
C              Store coefficients
C     
C 
               FAC = H/DBLE(M+1)
               DO 1260 I=1,NDY
                  DENS(I,IL-M-1) = TAB(I,M)*FAC
 1260          CONTINUE
C
               DO 1270 I=1,NDZ
                  DENS(NDY+I,IL-M) = TAB(NDY+I,M)
 1270          CONTINUE                              
            ENDIF
C
            IF (M .LE. IRHO) THEN
C     
C              Extrapolation of right differences
C
               DO 1290 J=M,KMAX
                  FAC = (DBLE(NJ(J))**M)/FACUL(M+1)
                  IEN = IPDIF(J+1)
                  DO 1280 I=1,ND
                     TAB(I,J) = YSAFE(I,IEN)*FAC
 1280             CONTINUE
 1290          CONTINUE               
C
               DO 1320 J=M+1,KMAX
                  DO 1310 L=J,M+1,-1
                     FAC = DBLE(NJ(J))/DBLE(NJ(L-1)) - 1.D0
                     DO 1300 I=1,ND
                        TAB(I,L-1) = TAB(I,L) + (TAB(I,L) - TAB(I,L-1))/
     $                       FAC
 1300                CONTINUE
 1310             CONTINUE
 1320          CONTINUE
C
C              Store coefficients
C
               FAC = H/DBLE(M+1)
               DO 1330 I=1,NDY
                  DENS(I,IR+M+1) = TAB(I,M)*FAC
 1330          CONTINUE
               DO 1340 I=1,NDZ
                  DENS(NDY+I,IR+M) = TAB(NDY+I,M)
 1340          CONTINUE                              
            ENDIF
C
 1350    CONTINUE
C
      ELSE
C
C
C --- Perturbed case ---
C
C
         DO 1470 M=1,IRHO
C
C           COMPUTE DIFFERENCES
C
            IF (M .GT. 1) THEN
               INC=0
               NX=ND
               I0 = 1
            ELSE
               INC = NDY
               NX = NDZ
               I0 = NDY + 1
            ENDIF
C
            DO 1380 J=M+MU,KMAX
               L1 = IPDIF(J+1)
               L0 = L1 + M + MU - J
               DO 1370 L=L1,L0,-1
                  DO 1360 I=I0,ND
                     YSAFE(I,L) = YSAFE(I,L) - YSAFE(I,L-1)
 1360             CONTINUE
 1370          CONTINUE
 1380       CONTINUE
C
C           Extrapolation of right differences
C
            DO 1410 J=M+MU,KMAX
               FAC = (DBLE(NJ(J))**(M-1))/FACUL(M+1)
               IEN = IPDIF(J+1)
               DO 1390 I=I0,NDY
                  TAB(I,J) = YSAFE(I,IEN)*FAC
 1390          CONTINUE
               FAC = FAC*DBLE(NJ(J))
               DO 1400 I=NDY+1,ND
                  TAB(I,J) = YSAFE(I,IEN)*FAC
 1400          CONTINUE
 1410       CONTINUE
C
            DO 1440 J=M+MU+1,KMAX
               DO 1430 L=J,M+MU+1,-1
                  FAC = DBLE(NJ(J))/DBLE(NJ(L-1)) - 1.D0
                  DO 1420 I=I0,ND
                     TAB(I,L-1) = TAB(I,L)+(TAB(I,L)-TAB(I,L-1))/FAC
 1420             CONTINUE
 1430          CONTINUE
 1440       CONTINUE
C
C           Store coefficients
C
            DO 1450 I=I0,NDY
               DENS(I,IR+M) = H*TAB(I,M+MU)
 1450       CONTINUE
            DO 1460 I=NDY+1,ND
               DENS(I,IR+M) = TAB(I,M+MU)
 1460       CONTINUE
C
 1470    CONTINUE
C
      ENDIF
C
C --- Hermite interpolation ---
C
C
      IF (ILAM .GT. 0) THEN 
         IDEL = 0
         ILAMH = ILAM - 1
      ELSE
         IDEL = 1
         ILAMH = 0
      ENDIF
      IF (IDR .EQ. 0) THEN
         IRHOH = IRHO+1
      ELSE
         IRHOH = IRHO
      ENDIF
C
C     Evaluate coefficients for differential components
C
      IF (NDY .GT. 0) CALL MXHERM (ND, DENS, NDY, ILAM, IRHO+1)
C
C     Evaluate coefficients for algebraic components
C
      IF (NDZ .GT. 0) CALL MXHERM (ND, DENS(NDY+1,2-IDEL), NDZ, ILAMH,
     $     IRHO)
C
C
C --- Error estimation for differential components ---
C
      THERR = DBLE(ILAM+1)/DBLE(ILAM+IRHOH+2)
      ERRFAC = THERR**(ILAM+1)*(1.D0-THERR)**(IRHO+1)
      DO 1480 I=1,NDY
         YWRK(I) = ABS(DENS(I,ILAM+IRHOH+2))*ERRFAC
 1480 CONTINUE 
C
      ERR = 0.D0
      DO 1490 J=1,NDP
         I = J
         IF (IDOFLG .EQ. 1) I = INDP(J)
         IF (ITOL .EQ. 0) THEN
            WT = RTOL(1)*SCAL(I) + ATOL(1)
         ELSE
            WT = RTOL(I)*SCAL(I) + ATOL(I)
         ENDIF 
         ERR = ERR + MIN (ABS(YWRK(J))/WT, 1.D20)**2
 1490 CONTINUE 
      ERRP = ERR
      IF (NDP .GT. 0) ERRP = ERRP / DBLE(NDP)
C
      ERR = 0.D0
      DO 1500 J=1,NDV
         I = J
         IF (IDOFLG .EQ. 1) I = INDV(J)
         IF (ITOL .EQ. 0) THEN
            WT = RTOL(1)*SCAL(NP+I) + ATOL(1)
         ELSE  
            WT = RTOL(NP+I)*SCAL(NP+I) + ATOL(NP+I)
         ENDIF 
         ERR = ERR + MIN (ABS(YWRK(NDP+J))/WT, 1.D20)**2
 1500 CONTINUE 
      ERRV = ERR 
      IF (NDV .GT. 0) ERRV = ERRV / DBLE(NDV)
C     
      NN = NP + NV
      NND = NDP + NDV
      ERR = 0.D0
      DO 1510 J=1,NDU
         I = J
         IF (IDOFLG .EQ. 1) I = INDU(J)
         IF (ITOL .EQ. 0) THEN
            WT = RTOL(1)*SCAL(NN+I) + ATOL(1)
         ELSE  
            WT = RTOL(NN+I)*SCAL(NN+I) + ATOL(NN+I)
         ENDIF
         ERR = ERR + MIN (ABS(YWRK(NND+J))/WT, 1.D20)**2
 1510 CONTINUE
      ERRU = ERR
      IF (NDU .GT. 0) ERRU = ERRU / DBLE(NDU)
C
      ERR = SQRT(ERRP + ERRV + ERRU)
C
      IF (MNINT .GE. 4) THEN
         ERRPRI = ERR
         WRITE (LUINT,*) ' - MEXDO: Est. normalized error:', ERRPRI
      ENDIF
      RETURN
      END
C
      SUBROUTINE MXDRE1 (FILNAM, INFO, ND, NDP, NDV, NDU, NDA, NDL,
     $     KDENS, ILAM, IRHO, T, TOLD, H, DENS, YIP,
     $     T0, TFIN, DOUTS, ICALL,
     $     INDP, INDV, INDU, INDA, INDL, IRTRN)
C
C  This routine reads the dump data produced by routine MXDDU1
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      CHARACTER FILNAM*(*)
      INTEGER ND, NDP, NDV, NDU, NDA, NDL, KDENS, ILAM, IRHO
      DOUBLE PRECISION T, TOLD, H, DENS(*), YIP(*)
      DOUBLE PRECISION DOUTS(*)
      INTEGER INFO(50)
      INTEGER INDP(*), INDV(*), INDU(*), INDA(*), INDL(*)
C
C* File              mxdre1.f
C* Version           1.1
C* Latest Change     96/04/03 (1.6)
C* Modification history
C    1.1             Avoid free format because of compiler dependency.
C
C
C  Machine constant (logical unit stdout)
      LUSOUT = 6
      LUNIT=10
      IRTRN=0
C
C
C  First call (read data at t0)
C==============================
C
      IF (ICALL .EQ. 0) THEN
         OPEN (UNIT=LUNIT,FILE=FILNAM)
         READ (LUNIT,*,END=1020,ERR=1020) 
         READ (LUNIT,*,END=1020,ERR=1020)
         READ (LUNIT,*,END=1020,ERR=1020)
         READ (LUNIT,*,END=1020,ERR=1020) 
         READ (LUNIT,9020,END=1020,ERR=1020) (INFO(I),I=1,50)
         READ (LUNIT,9150,END=1020,ERR=1020) NDP,NDV,NDU,NDA,NDL
C
         READ (LUNIT,*)
         READ (LUNIT,9020) (INDP(I),I=1,NDP)
         READ (LUNIT,*)
         READ (LUNIT,9020) (INDV(I),I=1,NDV)
         READ (LUNIT,*)
         READ (LUNIT,9020) (INDU(I),I=1,NDU)
         READ (LUNIT,*) 
         READ (LUNIT,9020) (INDA(I),I=1,NDV)
         READ (LUNIT,*)
         READ (LUNIT,9020) (INDL(I),I=1,NDL)
C
         READ (LUNIT,9120) ND,KDENS
         MDOUT = INFO(31)
         NDOUT = INFO(32)
C         
         IF (MDOUT .EQ. 3) THEN
            READ (LUNIT,9210) DOUTS(1)
         ENDIF
         IF (MDOUT .EQ. 4) THEN
            READ (LUNIT,9030) (DOUTS(L),L=1,NDOUT)
         ENDIF
         READ (LUNIT,9220) T0,TFIN
         READ (LUNIT,9030) (YIP(L),L=1,ND)
      ENDIF
C
C  Subsequent calls
C==============================
C
      IF (ICALL .GT. 0) THEN
         READ (LUNIT,9110,END=1010) INFO(41) 
         READ (LUNIT,*,END=1010) INFO(42)
         READ (LUNIT,9230) TOLD,T,H
         READ (LUNIT,9130) ILAM,IRHO,KDDUMP
         DO 1000 K=1,KDDUMP
            READ (LUNIT,9110) KK
            READ (LUNIT,9030) (DENS((K-1)*ND+L),L=1,ND)
 1000    CONTINUE
      ENDIF
C
      RETURN
C
C  End of file reached
 1010 CONTINUE
      IRTRN=1
      RETURN
C
C  Error/no data
 1020 CONTINUE 
      IRTRN=-1
      RETURN
C
 9010 FORMAT (A)
 9020 FORMAT ((10(2X,I6)))
 9030 FORMAT ((3(2X,1PD21.13)))
C
 9110 FORMAT (2X,I6)
 9120 FORMAT (2(2X,I6))
 9130 FORMAT (3(2X,I6))
 9150 FORMAT (5(2X,I6))
C
 9210 FORMAT (2X,1PD21.13)
 9220 FORMAT (2(2X,1PD21.13))
 9230 FORMAT (3(2X,1PD21.13))
C
      END
C
      SUBROUTINE MXHERM (NDCL, Y, N, ILAM, IRHO)
C     Computes the coefficients of the interpolation formula
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER NDCL, N
      DOUBLE PRECISION Y(NDCL,*)    
C
C
C* File              mxherm.f
C* Version           1.0
C* Latest Change     94/02/22 (1.2)
C      
C
      LOGICAL QARRET
      INTEGER I, J, IL, IR, L, ML, MR
C
      SAVE
C 
      IL = ILAM + 1
      IR = IRHO + 1
      ML = IL
      MR = IL
      QARRET = .TRUE.
      DO 1060 J=2,IL+IR
         DO 1010 I = MAX (1,IL-J+2), MIN (IL+IR-J+1, IL-(J+1)/2)
            DO 1000 L=1,N 
               Y(L,I) = Y(L,I+1) - Y(L,I)
 1000       CONTINUE
 1010    CONTINUE
         DO 1030 I=MIN (IL+IR, IL+J-1), MAX (J, IL+J/2+1), -1
            DO 1020 L=1,N 
               Y(L,I) = Y(L,I) - Y(L,I-1)
 1020       CONTINUE
 1030    CONTINUE
         IF (MOD(J,2) .EQ. 1 .AND. QARRET) THEN
            IF (ML .GT. 1 .AND. MR .LE. IL+IR) THEN 
               ML = ML - 1
               DO 1040 L=1,N 
                  Y(L,ML) = Y(L,MR) - Y(L,ML)
 1040          CONTINUE
            ELSE
               QARRET = .FALSE.
            ENDIF
         ENDIF
         IF (MOD(J,2) .EQ. 0 .AND. QARRET) THEN
            IF (ML .GE. 1 .AND. MR .LT. IL+IR) THEN 
               MR = MR + 1
               DO 1050 L=1,N 
                  Y(L,MR) = Y(L,MR) - Y(L,ML)
 1050          CONTINUE
            ELSE
               QARRET = .FALSE.
            ENDIF
         ENDIF
 1060 CONTINUE
      RETURN
      END
C
      SUBROUTINE MXIPOL (NDCL, DENS, N, T, TOLD, H, ILAM, IRHO,
     $                   CONTEX)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER NDCL, N, ILAM, IRHO
      DOUBLE PRECISION DENS, T, TOLD, H, CONTEX
      DIMENSION DENS(NDCL,*), CONTEX(N)
C ---------------------------------------------------------
C     This function can be used for continuous output in 
C     connection with the output subroutine for MEXX.
C 
C
C
C* File              mxipol.f
C* Version           1.0
C* Latest Change     94/02/22 (1.2)
C      
C
C ---------------------------------------------------------
C
      INTEGER J, L, IL, IR
      DOUBLE PRECISION THETA
C
      IL = ILAM + 1
      IR = IRHO + 1
      THETA = (T-TOLD)/H
      IF (ILAM .GT. IRHO) THEN
         DO 1000 L=1,N
            CONTEX(L) = DENS(L,1)
 1000    CONTINUE 
         DO 1020 J=2,IL-IR
            DO 1010 L=1,N
               CONTEX(L) = DENS(L,J) + CONTEX(L)*THETA
 1010       CONTINUE 
 1020    CONTINUE 
         DO 1040 J=1,IR
            DO 1030 L=1,N
               CONTEX(L) = DENS(L,IL+IR-J+1) + CONTEX(L)*(THETA-1.D0)
               CONTEX(L) = DENS(L,IL-IR+J) + CONTEX(L)*THETA
 1030       CONTINUE
 1040    CONTINUE 
      ELSE
         DO 1050 L=1,N
            CONTEX(L) = DENS(L,IL+IR)
 1050    CONTINUE 
         DO 1070 J=2,IR-IL+1
            DO 1060 L=1,N
               CONTEX(L) = DENS(L,IL+IR-J+1) + CONTEX(L)*(THETA-1.D0)
 1060       CONTINUE
 1070    CONTINUE 
         DO 1090 J=1,ILAM
            DO 1080 L=1,N
               CONTEX(L) = DENS(L,J) + CONTEX(L)*THETA
               CONTEX(L) = DENS(L,2*ILAM-J+2) + CONTEX(L)*(THETA-1.D0)
 1080       CONTINUE
 1090    CONTINUE 
         DO 1100 L=1,N
            CONTEX(L) = DENS(L,ILAM+1) + CONTEX(L)*THETA
 1100    CONTINUE 
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE MXRTSC (ISTEP, ICTYPE, NP, NV, NL, NU, NSWIF,
     $     NPDIM, NVDIM, NLDIM, NUDIM,
     $     T, P, V, U, A, RLAM, GVAL0, GVAL1,
     $     IRES0, IRES1, NSC, NSUB, GSWIT, MNSWI, LUSWI, IRTRN)
C
C  Initialize and check for sign change in switching function
C
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      DOUBLE PRECISION P(NPDIM), V(NVDIM), U(NUDIM), A(NVDIM), 
     $     RLAM(NLDIM)
      DOUBLE PRECISION GVAL0(*), GVAL1(*)
      INTEGER IRES0(*), IRES1(*)
      EXTERNAL GSWIT
C
C
C* File              mxroot.f
C* Version           1.2
C* Latest Change     96/07/17 (1.7)
C* Modification history
C    1.2             If more than one root is found report the solution
C                    values for the first root instead of the root
C                    detected last.
C
      PARAMETER ( GMIN = 1.D-6 )
C
      SAVE
C
      IRTRN=0
C
      IF (ICTYPE .EQ. 0) THEN
C
         IF (MNSWI .GE. 3) 
     $        WRITE (LUSWI,*)
     $        ' - MXRTSC:  L-evaluation of switching function'
         CALL MONON (11)
         CALL GSWIT (NP, NV, NL, NU, T, P, V, U, A, RLAM, NSWIF, GVAL1)
         CALL MONOFF(11)
         DO 1000 I=1,NSWIF
            IRES1(I)=0
            IF (ABS(GVAL1(I)) .LE. GMIN) IRES1(I)=1
 1000    CONTINUE
         RETURN
C
      ELSE
C
         DO 1010 I=1,NSWIF  
            GVAL0(I)=GVAL1(I)
            IRES0(I)=IRES1(I)
 1010    CONTINUE
C
         IF (MNSWI .GE. 3) 
     $        WRITE (LUSWI,*)
     $        ' - MXRTSC:  R-evaluation of switching function'
         CALL MONON(11)
         CALL GSWIT (NP, NV, NL, NU, T, P, V, U, A, RLAM, NSWIF,
     $        GVAL1)
         CALL MONOFF(11)
C
C  check for sign change
C
         DO 1020 I=1,NSWIF
            IRES1(I)=0
            IF (ABS(GVAL1(I)) .LE. GMIN) IRES1(I)=1
 1020    CONTINUE
C
         NSC=0
         DO 1030 I=1,NSWIF  
            IF (GVAL0(I)*GVAL1(I) .LT. 0.D0 .AND. 
     $           IRES0(I)+IRES1(I) .EQ. 0) NSC=NSC+1
 1030    CONTINUE
         IF (NSUB .GT. 1) NSC=NSC+1
         IF (MNSWI .GE. 1 .AND. NSC .EQ. 0) WRITE (LUSWI,*)
     $        ' - MXRTSC:  no sign change detected'
C
      ENDIF
C
      RETURN
C
C  end of subroutine MXRTSC
C
      END
      SUBROUTINE MXRTFD (ISTEP, ND, NDP, NDV, NDU, NDA, NDR,
     $     KDENS, ILAM, IRHO, T, TOLD, H, DENS, YIP,
     $     NSWIF, NSUB, GVANEW, GVAL0, GVAL1, GVALL, GVALR, GVALH,
     $     IRES0, IRES1, IRESL, IRESR,
     $     TOL, TZERO, IZERO, GSWIT, NSC, IRTRN,
     $     NP, NV, NL, NG, NU,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $     FPROB,
     $     P, V, U, A, RLAM, SCAL,
     $     RTOL, ATOL, ITOL,
     $     ITMAXP, NITER, EPSREQ, EPSEST,
     $     WF, WP, WV, WZ, WG, WGI, WPDOT, WUDOT, AM, GP, FL,
     $     MXJOB, RWK, IWK, LRWK, LIWK, IFAIL)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
      INTEGER ND, NDP, NDV, NDU, NDA, NDR, KDENS, 
     $     ILAM, IRHO
      DOUBLE PRECISION T, TOLD, H, DENS(ND,KDENS), YIP(ND)
      DOUBLE PRECISION GVAL0(*), GVAL1(*), GVANEW(*)
      DOUBLE PRECISION GVALL(*), GVALR(*), GVALH(*)
      DIMENSION IZERO(*)
      DIMENSION IRESL(*), IRESR(*), IRES0(*), IRES1(*)
      DOUBLE PRECISION TZI(100)
      DIMENSION INDTZ(100), IHELP(100)
      INTEGER IRTRN
C
      EXTERNAL GSWIT  
C
C 
      INTEGER NP, NV, NL, NG, NU, NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, 
     $     NTODIM, N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF 
C 
      EXTERNAL FPROB 
C 
      DOUBLE PRECISION    P(NPDIM), V(NVDIM), U(NUDIM), A(NVDIM), 
     $     RLAM(NLDIM), SCAL(N2UDIM), RTOL(NTODIM), ATOL(NTODIM) 
      INTEGER ITOL, ITMAXP, NITER 
      DOUBLE PRECISION EPSREQ, EPSEST 
      DOUBLE PRECISION WF(NVDIM), WP(NPDIM), WV(NVDIM), WZ(NVLDIM), 
     $     WG(NGDIM), WGI(NLDIM), WPDOT(NPDIM), WUDOT(NUDIM), AM(LDA,
     $     MDA), 
     $     GP(LDG,MDG), FL(LDF,MDF) 
      INTEGER MXJOB(150) 
      DOUBLE PRECISION RWK(LRWK) 
      INTEGER IWK(LIWK), LRWK, LIWK, IFAIL 
      LOGICAL QLEFT, QRIGHT
C
C
C     PARAMETER ( TOLREQ=1.D-12 )
      PARAMETER ( GMIN = 1.D-6 )
      PARAMETER ( ITMAXR = 10  )
      PARAMETER ( NJAC1 = 10 , NJAC2 = 10 )
C
C  This routine is called only, if either the sign change routine
C  (which works at the integration points only) has indicated an sign
C  change or the user has prescribed to check at further internal
C  points. 
C
C  internal parameters:
C
C    ITMAXR: the maximum number of newton iterations for root finding
C
C    NJAC1: the maximum number of jacobians for the newton scheme
C           which solves the unprojected problem
C           ( njac1 > 0 ,  njac1 = 1 : usual simplified Newton)
C
C    NJAC2: the maximum number of jacobians for the newton scheme
C           which solves the projected problem
C           ( njac2 > 0 ,  njac2 = 1 : usual simplified Newton)
C
C  INPUT
C=======
C
C  ISTEP :  no. of integration step (istep=0 indicates t=t0)
C
C  YIP   :  for istep=0: the current "solution" of selected
C           components at t=t0
C           workspace otherwise
C
C  OUTPUT
C========
C
C  IRTRN : =0 no root found
C          = 1 root found
C          < 0 fail of mxswit
C      
C  TZERO: the root (if irtrn=1: g(t,y(t))=0 for a component
C
C  IZERO: indicates the components of gswit for which g=0
C         ( izero(l)=1 if component l of g equals zero )
C
C
C* File              mxroot.f
C* Version           1.0
C* Latest Change     96/07/17 (1.7)
C      
      INTEGER NSWMAX
C
      SAVE
C
      DATA NSWMAX /0/
C
      MNSWI = MXJOB(17)
      IRTRN = 0
C
C
C  Initiation at the first call (icall=0)
C========================================
C
      IF (ISTEP .EQ. 0) THEN
         LUSWI = MXJOB(18)
         NSWMAX = NSWIF
         IF (MNSWI .GE. 3) WRITE (LUSWI,*) '   MXRTF:  Initiation'
         NDY = NP + NV + NU
         NDZ = NV + NL
         IPP = MIN(ND, 1)
         IPV = MIN(ND, IPP+NP)
         IPU = MIN(ND, IPV+NV)
         IPA = MIN(ND, IPU+NU)
         IPL = MIN(ND, IPA+NV)
         PERT=SQRT(TOL)
C
C  call GSWIT or copy values from MXRTSC
C
         DO 1000 I=1,NSWIF
            GVALR(I)=GVAL1(I)
            IRESR(I)=IRES1(I)
            IF (IRESR(I) .EQ. 1 .AND. MNSWI .GE. 1)
     $           WRITE (LUSWI,*)
     $           ' - MXRTFD:  small left residual detected',
     $           ' for component:', I
 1000    CONTINUE
C
         RETURN
      ENDIF
C
C  Subsequent calls
C==================
C
      TOLREQ=1.D-1*TOL*MAX(ABS(TOLD), ABS(TOLD+H))
      IF (MNSWI .GE. 3) WRITE (LUSWI,*) '   MXRTFD: '
      IF (MNSWI .GE. 3) WRITE (LUSWI,*) '   MXRTFD:  Standard call ',
     $     ' at t =', T, ' nsub =', NSUB
      IF (MNSWI .GE. 3) WRITE (LUSWI,*) 
     $     '   MXRTFD: abs. tol. req. =', TOLREQ
C
C  prepare dense output evaluation
C
      IF (ILAM .GT. 0) THEN
         IDEL = 0
         ILAMH = ILAM - 1
      ELSE
         IDEL = 1
         ILAMH = 0
      ENDIF
C
      TLEFT=TOLD
      TRIGHT=T
      DTSUB=(TRIGHT-TLEFT)/DBLE(NSUB)
      TRR=TLEFT
      NOZERO=0
C
      DO 1010 I=1,NSWIF   
         GVALR(I)=GVAL0(I)
         IRESR(I)=IRES0(I) 
 1010 CONTINUE
C
C
C  the outer loop over an eventual subdivision
C==========================================================
C
      DO 1130 ISUB=1,NSUB
C
C  the old right is now left
C   (used values from last call to MXRTFD or values from MXRTSC)
C
         DO 1020 I=1,NSWIF   
            GVALL(I)=GVALR(I)
            IRESL(I)=IRESR(I)
 1020    CONTINUE
C
         TLL=TRR
         TRR=TRR+DTSUB
         IF (MNSWI .GE. 2 .AND. NSUB .GT. 1)
     $        WRITE (LUSWI,*) ' - MXRTFD: in subinterval no.:', ISUB
         IF (MNSWI .GE. 3)
     $        WRITE (LUSWI,*) '   MXRTFD:  R-evaluation of switching',
     $        ' function at t =', TRR
C
         IF (ISUB .NE. NSUB) THEN
            CALL MONON(4)
            IF (NDY .GT. 0) THEN
               CALL MXIPOL (ND, DENS, NDY, TRR, TOLD, H, ILAM, IRHO,
     $              YIP)
            ENDIF
            IF (NDZ .GT. 0) THEN
               CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, TRR, TOLD,
     $              H, ILAMH, IRHO, YIP(NDY+1))
            ENDIF
            CALL MONOFF(4)
C
            CALL MONON(11)
            CALL GSWIT (NP, NV, NL, NU, TRR,
     $           YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $           NSWIF, GVALR)
            NSWIF = MIN (NSWIF, NSWMAX)
            CALL MONOFF(11)
C
C
         ELSE
C  copy values evaluated in MXRTSC
            DO 1030 I=1,NSWIF
               GVALR(I)=GVAL1(I)
               IRESR(I)=IRES1(I)
 1030       CONTINUE
         ENDIF
C
C
         NSR=0
         DO 1040 I=1,NSWIF
            IRESR(I)=0
            IF (ABS(GVALR(I)) .LE. GMIN) THEN
               NSR=NSR+1
               IRESR(I)=1
               IF (MNSWI .GE. 5)
     $              WRITE (LUSWI,*)
     $              '   MXRTFD:  small residual detected for',
     $              ' component:', I, GVALR(I)
            ENDIF
 1040    CONTINUE
C
C check for sign change
C
         NSC=0
         DO 1050 I=1,NSWIF
            IZERO(I)=0
            IF (GVALL(I)*GVALR(I) .LT. 0.D0 .AND. 
     $           IRESL(I)+IRESR(I) .EQ. 0) THEN
               NSC=NSC+1
               IZERO(I)=1
               IF (MNSWI .GE. 5)
     $              WRITE (LUSWI,*)
     $              '   MXRTFD:  sign change detected for',
     $              ' component:', I, GVALL(I), GVALR(I)
            ENDIF
 1050    CONTINUE
C
C
         IF (NSC .LE. 0) THEN
            IF (MNSWI .GE. 2)
     $           WRITE (LUSWI,*) ' - MXRTFD:  no sign change detected'
         ELSE
            IF (MNSWI .GE. 2) THEN 
               WRITE (LUSWI,*) 
     $              ' - MXRTFD:  Try to find root(s) within the',
     $              ' interval' 
               WRITE (LUSWI,*) ' - MXRTFD: Tleft =', TLL, ' Tright =',
     $              TRR
            ENDIF
         ENDIF
C
C
C  the second main loop: component for component
C=================================================
C
C
         MOZERO=NOZERO
         ICASE2=0
         ICASE3=0
         ICASE4=0
         DO 1090 ISWIF=1,NSWIF
C
            IF (IRESL(ISWIF) .EQ. 1 .OR. IRESR(ISWIF) .EQ. 1) THEN
               IF (IRESL(ISWIF) .EQ. 1 .AND. IRESR(ISWIF) .EQ. 0) THEN
                  IZERO(ISWIF)=2
                  ICASE2=ICASE2+1
                  IF (MNSWI .GE. 4) WRITE (LUSWI,*)
     $                 '   MXRTFD:  small left residual for component:',
     $                 ISWIF
                  GOTO 1090
               ENDIF
               IF (IRESL(ISWIF) .EQ. 0 .AND. IRESR(ISWIF) .EQ. 1) THEN
                  IZERO(ISWIF)=3
                  ICASE3=ICASE3+1
                  IF (MNSWI .GE. 4) WRITE (LUSWI,*)
     $                 '   MXRTFD:  small right residual for',
     $                 ' component:', ISWIF 
                  GOTO 1090
               ENDIF
               IF (IRESL(ISWIF) .EQ. 1 .AND. IRESR(ISWIF) .EQ. 1) THEN
                  IZERO(ISWIF)=4
                  ICASE4=ICASE4+1
                  IF (MNSWI .GE. 4) WRITE (LUSWI,*)
     $                 '   MXRTFD:  small left+right residual for',
     $                 ' component:', ISWIF
                  GOTO 1090
               ENDIF
            ENDIF
            IF (IZERO(ISWIF) .EQ. 0) THEN      
               IZERO(ISWIF)=1
               IF (MNSWI .GE. 4) WRITE (LUSWI,*)
     $              '   MXRTFD:  no sign change, no small residual',
     $              ' for component:', ISWIF
               GOTO 1090
            ENDIF
C
            IF (MNSWI .GE. 3)
     $           WRITE (LUSWI,*) '   MXRTFD:  search root for',
     $           ' component:', ISWIF 
C
C  prepare subdivision/Newton scheme
C
            TMAX=MAX(ABS(TLL), ABS(TRR))
            DTPERT=TMAX*PERT
C
            TL=TLL
            TR=TRR
            PHIL=GVALL(ISWIF)
            PHIR=GVALR(ISWIF)
            ISTART=0
C
C
 1060       CONTINUE
C
C  subdivide interval 
C
            ISTART=ISTART+1
            IF (ISTART .GT. 10) THEN
               IF (MNSWI .GE. 1)
     $              WRITE (LUSWI,*) ' - MXRTFD: too many initial',
     $              ' subdivisions: zero location fails' 
               IZERO(ISWIF)=7
               GOTO 1090
            ENDIF
            T0=TR-PHIR*(TR-TL)/(PHIR-PHIL)
C
C  evaluate switching function at t0
C
C
            CALL MONON(4)
            IF (NDY .GT. 0) THEN
               CALL MXIPOL (ND, DENS, NDY, T0, TOLD, H, ILAM, IRHO, YIP)
            ENDIF
            IF (NDZ .GT. 0) THEN
               CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, T0, TOLD,
     $              H, ILAMH, IRHO, YIP(NDY+1))
            ENDIF
            CALL MONOFF(4)
C
            CAll MONON(11)
            CALL GSWIT (NP, NV, NL, NU, T0,
     $           YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $           NSWIF, GVAL0)
            CALL MONOFF(11)
            NSWIF = MIN (NSWIF, NSWMAX)
C
            PHI0=GVAL0(ISWIF)
            IF (MNSWI .GE. 3) 
     $           WRITE (LUSWI,'(A,I4,3D12.3))')
     $           '   MXRTFD: Divide: iter,tl,tr,t0: ', ISTART, TL, TR,
     $           T0 
C
C  check for position of root
C
            QLEFT=PHIL*PHI0 .LT. 0.D0
            QRIGHT=PHIR*PHI0 .LE. 0.D0
C
C  more than one root ?
C
            IF (QLEFT .AND. QRIGHT) THEN
               IF (MNSWI .GE. 1) THEN
                  WRITE (LUSWI,*) ' - MXRTFD: more than one root',
     $                 ' expected'
                  WRITE (LUSWI,*)
     $                 '   MXRTFD: one should subdivide the initial',
     $                 ' interval'
                  WRITE (LUSWI,*) ' - MXRTFD: zero location fails '
               ENDIF
               IZERO(ISWIF)=6
               GOTO 1090
            ENDIF
C
C  no root ?
C
            IF ( .NOT. QLEFT .AND. .NOT. QRIGHT) THEN
               IF (MNSWI .GE. 1) THEN
                  WRITE (LUSWI,*) ' - MXRTFD: no root expected'
                  WRITE (LUSWI,*)
     $                 '   MXRTFD: one should check the initial',
     $                 ' interval' 
                  WRITE (LUSWI,*) ' - MXRTFD: zero location fails '
               ENDIF
               IZERO(ISWIF)=6
               GOTO 1090
            ENDIF
C
C  jacobian
C
            IF (QLEFT) THEN
               TH=T0-DTPERT
            ELSE
               TH=T0+DTPERT
            ENDIF
C
C  evaluate switching function at t0+/-dtpert
C
            CALL MONON(4)
            IF (NDY .GT. 0) THEN
               CALL MXIPOL (ND, DENS, NDY, TH, TOLD, H, ILAM, IRHO, YIP)
            ENDIF
            IF (NDZ .GT. 0) THEN
               CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, TH, TOLD,
     $              H, ILAMH, IRHO, YIP(NDY+1))
            ENDIF
            CALL MONOFF(4)
C
            CALL MONON(11)
            CALL GSWIT (NP, NV, NL, NU, TH,
     $           YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $           NSWIF, GVALH)
            CALL MONOFF(11)
            NSWIF = MIN (NSWIF, NSWMAX)
C
            PHIH=GVALH(ISWIF)
            PHIJAC=(PHIH-PHI0)/(TH-T0)
C
C  first newton correction
C
            DELTAT=-PHI0/PHIJAC
C
C  is convergence possible ?
C
            IF (QLEFT .AND. DELTAT .GT. 0) THEN
               TR=T0
               PHIR=PHI0
               GOTO 1060
            ENDIF
            IF (QRIGHT .AND. DELTAT .LT. 0) THEN
               TL=T0
               PHIL=PHI0
               GOTO 1060
            ENDIF
            NEWTIT=0
            TACC=T0
C
C  the newton iteration
C
 1070       CONTINUE
            NEWTIT=NEWTIT+1
            IF (MNSWI .GE. 3) WRITE (LUSWI,*)
     $           '   MXRTFD: Newton: iter, tacc, deltat ', NEWTIT, TACC,
     $           DELTAT
            IF (NEWTIT .GT. 10) THEN
               IF (MNSWI .GE. 1) WRITE (LUSWI,*)
     $              ' - MXRTFD: too many newton steps: zero location',
     $              ' fails'
               IZERO(ISWIF)=8
               GOTO 1090
            ENDIF
            IDAMP=0
            
 1080       TNEW=TACC+DELTAT
C
C  check for solution
C
            IF (IDAMP .EQ. 0 .AND. ABS(DELTAT) .LE. TOLREQ) THEN
               TZERO=TNEW
               IF (MNSWI .GE. 3) THEN 
                  WRITE (LUSWI,*) '   MXRTFD: root detected, tzero = ',
     $                 TZERO
                  WRITE (LUSWI,*) '   MXRTFD:  # subdivisions: ', ISTART
                  WRITE (LUSWI,*) '   MXRTFD:  # newton steps: ', NEWTIT
               ENDIF
               NOZERO=NOZERO+1
               
               IZERO(ISWIF)=5
               INDTZ(NOZERO)=ISWIF
               TZI(NOZERO)=TZERO
               GOTO 1090
            ENDIF
C
            IF (TNEW .LT. TL .OR. TNEW .GT. TR) THEN
               IDAMP=IDAMP+1
               IF (IDAMP .LT. 5) THEN
                  IF (MNSWI .GE. 3) 
     $                 WRITE (LUSWI,*) '   MXRTFD: emergency damping',
     $                 ' for step:', NEWTIT
                  DELTAT=DELTAT/5.D0
                  GOTO 1080
               ELSE
                  IF (MNSWI .GE. 1) THEN 
                     WRITE (LUSWI,*) ' - MXRTFD: damped newton',
     $                    ' iteration fails'
                     WRITE (LUSWI,*) ' - MXRTFD: emergency subdivision'
                  ENDIF
                  IF (QRIGHT) TL=T0
                  IF (QRIGHT) PHIL=PHI0
                  IF (QLEFT) TR=T0
                  IF (QLEFT) PHIR=PHI0
                  GOTO 1060
               ENDIF
            ENDIF
C
C prepare next newton step
C
            TACC=TNEW
C
C  evaluate switching function at tacc
C
C
            CALL MONON(4)
            IF (NDY .GT. 0) THEN
               CALL MXIPOL (ND, DENS, NDY, TACC, TOLD, H, ILAM, IRHO,
     $              YIP)
            ENDIF
            IF (NDZ .GT. 0) THEN
               CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, TACC, TOLD,
     $              H, ILAMH, IRHO, YIP(NDY+1))
            ENDIF
            CALL MONOFF(4)
C
            CALL MONON(11)
            CALL GSWIT (NP, NV, NL, NU, TACC,
     $           YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $           NSWIF, GVANEW)
            CALL MONOFF(11)
            NSWIF = MIN (NSWIF, NSWMAX)
C
            PHINEW=GVANEW(ISWIF)
            TH=TNEW+DTPERT
            IF (TH .GT. TR) TH=TNEW-DTPERT
C
C  evaluate switching function at tacc+/-dtpert
C
C
            CALL MONON(4)
            IF (NDY .GT. 0) THEN
               CALL MXIPOL (ND, DENS, NDY, TH, TOLD, H, ILAM, IRHO, YIP)
            ENDIF
            IF (NDZ .GT. 0) THEN
               CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, TH, TOLD,
     $              H, ILAMH, IRHO, YIP(NDY+1))
            ENDIF
            CALL MONOFF(4)
C
            CALL MONON(11)
            CALL GSWIT (NP, NV, NL, NU, TH,
     $           YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $           NSWIF, GVALH)
            CALL MONOFF(11)
            NSWIF = MIN (NSWIF, NSWMAX)
            PHIH=GVALH(ISWIF)
            PHIJAC=(PHIH-PHINEW)/(TH-TNEW)
            DELTAT=-PHINEW/PHIJAC
            GOTO 1070
C

 1090    CONTINUE 
C
C  report small residuals
C
         IF (ICASE3 .GT. 0) THEN
            IF (MNSWI .GE. 1) THEN
               WRITE (LUSWI,*) ' - MXRTFD: at t =', TRR
               WRITE (LUSWI,*) ' - MXRTFD: # of small residual',
     $              ' components:', ICASE3
            ENDIF 
            IF (MNSWI .GE. 2) THEN
               II=0
               DO 1100 I=1,NSWIF
                  II=II+1
                  IF (IZERO(I) .EQ. 3) IHELP(II)=I
 1100          CONTINUE
               WRITE (LUSWI,'(9X,10I5)') (IHELP(II), II=1,ICASE3)
            ENDIF 
         ENDIF 
C
         IF (ICASE4 .GT. 0) THEN 
            IF (MNSWI .GE. 1) THEN
               WRITE (LUSWI,*) ' - MXRTFD: at t =', TRR
               WRITE (LUSWI,*) ' - MXRTFD: # of small left+right',
     $              ' residual components:', ICASE4   
            ENDIF    
            IF (MNSWI .GE. 2) THEN
               II=0 
               DO 1110 I=1,NSWIF 
                  II=II+1 
                  IF (IZERO(I) .EQ. 4) IHELP(II)=I
 1110          CONTINUE 
               WRITE (LUSWI,'(9X,10I5)') (IHELP(II), II=1,ICASE4) 
            ENDIF
         ENDIF

C     
C report the zero values
C
         IF (MOZERO .LT. NOZERO) THEN
            IF (MNSWI .GE. 2) THEN
               III=NOZERO-MOZERO
               WRITE (LUSWI,*)
     $              ' - MXRTFD:   # of located (unprojected) zeros:',
     $              III 
               IF (MNSWI .GE. 2) THEN
                  DO 1120 II=MOZERO+1,NOZERO
                     WRITE (LUSWI,*)
     $                    '   MXRTFD:   component', INDTZ(II),
     $                    ' at t = ', TZI(II)
 1120             CONTINUE
               ENDIF 
            ENDIF 
         ENDIF
C
C report fail of zero location
C
 1130 CONTINUE 
C
C  return if no zero has been detected
C  
      IF (NOZERO .EQ. 0) THEN
         DO 1140 I=1,NSWMAX
            IZERO(I) = 0
 1140    CONTINUE
         IRTRN=0
         RETURN
      ENDIF
C  
C==================================================================
C  find all roots for the projected solution 
C==================================================================
C
C
      IF (MNSWI .GE. 3) WRITE (LUSWI,*)
     $     '   MXRTFD:  find all roots for the projected values'
C
      IRTRN=1  
C
      NZH=0
      DO 1170 II=1,NOZERO
C
         IZEROC=INDTZ(II)
         TITER=TZI(II) 
         ITER=0
C
C  Function evaluation for Newton scheme
C
         CALL MONON(4)
         IF (NDY .GT. 0) THEN
            CALL MXIPOL (ND, DENS, NDY, TITER, TOLD, H, ILAM, IRHO, YIP)
         ENDIF
         IF (NDZ .GT. 0) THEN
            CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, TITER, TOLD,
     $           H, ILAMH, IRHO, YIP(NDY+1))
         ENDIF
         CALL MONOFF(4)
C
         CALL MONON(1)
         CALL MXPRJC(NP, NV, NL, NG, NU,
     $        NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $        N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $        FPROB,
     $        TITER, YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $        SCAL,
     $        RTOL, ATOL, ITOL,
     $        ITMAXP, NITER, EPSREQ, EPSEST,
     $        WF, WP, WV, WZ, WG, WGI, WPDOT, WUDOT, AM, GP, FL,
     $        MXJOB, RWK, IWK, LRWK, LIWK, IFAIL)
         CALL MONOFF(1)
C
         CALL MONON(11)
         CALL GSWIT (NP, NV, NL, NU, TITER,
     $        YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $        NSWIF, GVAL0)
         CALL MONOFF(11)
            NSWIF = MIN (NSWIF, NSWMAX)
C
C  Jacobian at starting point of Newton scheme
C
         TH=TITER+DTPERT
C
         CALL MONON(4)
         IF (NDY .GT. 0) THEN
            CALL MXIPOL (ND, DENS, NDY, TH, TOLD, H, ILAM, IRHO, YIP)
         ENDIF
         IF (NDZ .GT. 0) THEN
            CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, TH, TOLD,
     $           H, ILAMH, IRHO, YIP(NDY+1))
         ENDIF
         CALL MONOFF(4)
C
         CALL MONON(1)
         CALL MXPRJC(NP, NV, NL, NG, NU,
     $        NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $        N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $        FPROB,
     $        TH, YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $        SCAL,
     $        RTOL, ATOL, ITOL,
     $        ITMAXP, NITER, EPSREQ, EPSEST,
     $        WF, WP, WV, WZ, WG, WGI, WPDOT, WUDOT, AM, GP, FL,
     $        MXJOB, RWK, IWK, LRWK, LIWK, IFAIL)
         CALL MONOFF(1)
C
         CALL MONON(11)
         CALL GSWIT (NP, NV, NL, NU, TH,
     $        YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $        NSWIF, GVALR)
         CALL MONOFF(11)
            NSWIF = MIN (NSWIF, NSWMAX)
         GJAC=(GVALR(IZEROC)-GVAL0(IZEROC))/(TH-TITER)
C
C
C-----------------------
C  Newton iteration loop
C-----------------------
C
 1150    ITER=ITER+1
         TDEL=-GVAL0(IZEROC)/GJAC
         TITER=TITER+TDEL
C
         IF (MNSWI .GE. 3)
     $        WRITE (LUSWI,*) '   MXRTFD:  Iter,Tdel,Titer:', ITER,
     $        TDEL, TITER
C
         CALL MONON(4)
         IF (NDY .GT. 0) THEN
            CALL MXIPOL (ND, DENS, NDY, TITER, TOLD, H, ILAM, IRHO, YIP)
         ENDIF
         IF (NDZ .GT. 0) THEN
            CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, TITER, TOLD,
     $           H, ILAMH, IRHO, YIP(NDY+1))
         ENDIF
         CALL MONOFF(4)
C
         CALL MONON(1)
         CALL MXPRJC(NP, NV, NL, NG, NU,
     $        NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $        N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $        FPROB,
     $        TITER, YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $        SCAL,
     $        RTOL, ATOL, ITOL,
     $        ITMAXP, NITER, EPSREQ, EPSEST,
     $        WF, WP, WV, WZ, WG, WGI, WPDOT, WUDOT, AM, GP, FL,
     $        MXJOB, RWK, IWK, LRWK, LIWK, IFAIL)
         CALL MONOFF(1)
C
         CALL MONON(11)
         CALL GSWIT (NP, NV, NL, NU, TITER,
     $        YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $        NSWIF, GVAL0)
         CALL MONOFF(11)
            NSWIF = MIN (NSWIF, NSWMAX)
C
C  convergence check
C
         IF (ABS(TDEL) .LE. 1.D1*TOLREQ) THEN
            IF (MNSWI .GE. 3) 
     $           WRITE (LUSWI,*)
     $           '   MXRTFD:  Root found, refined Tswit =', TITER
            TZERO=TITER
            GOTO 1160
         ENDIF
C
         IF (ITER .LT. NJAC2) THEN
            TH=TITER+DTPERT
            CALL MONON(4)
            IF (NDY .GT. 0) THEN
               CALL MXIPOL (ND, DENS, NDY, TH, TOLD, H, ILAM, IRHO, YIP)
            ENDIF
            IF (NDZ .GT. 0) THEN
               CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, TH, TOLD,
     $              H, ILAMH, IRHO, YIP(NDY+1))
            ENDIF
            CALL MONOFF(4)
C
            CALL MONON(1)
            CALL MXPRJC(NP, NV, NL, NG, NU,
     $           NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $           N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $           FPROB,
     $           TH, YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $           SCAL,
     $           RTOL, ATOL, ITOL,
     $           ITMAXP, NITER, EPSREQ, EPSEST,
     $           WF, WP, WV, WZ, WG, WGI, WPDOT, WUDOT, AM, GP, FL,
     $           MXJOB, RWK, IWK, LRWK, LIWK, IFAIL)
            CALL MONOFF(1)
C
            CALL MONON(11)
            CALL GSWIT (NP, NV, NL, NU, TH,
     $           YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $           NSWIF, GVALR)
            CALL MONOFF(11)
            NSWIF = MIN (NSWIF, NSWMAX)
            GJAC=(GVALR(IZEROC)-GVAL0(IZEROC))/(TH-TITER)
         ENDIF

C     
C---------------------------------------------------------
C  end of Newton iteration loop, check for iteration count
C---------------------------------------------------------
C
         IF (ITER .LT. ITMAXR) GOTO 1150
C
C
         IF (MNSWI .GE. 1) THEN
            WRITE (LUSWI,*) ' - MXRTFD:  Root finding fails'
            WRITE (LUSWI,*)
     $           ' - MXRTFD:  More than itmax =', ITMAXR,
     $           ' Newton iterations'
         ENDIF
         IRTRN=-3
         RETURN

C     
 1160    CONTINUE
         NZH=NZH+1
         TZI(NZH)=TZERO
         INDTZ(NZH)=IZEROC    
 1170 CONTINUE
C
C  get the leftmost root
C
 1180 CONTINUE
      ICHANG=0
      DO 1190 I=1,NZH-1
         IF (TZI(I) .GT. TZI(I+1)) THEN
            ICHANG=1
            TZERO=TZI(I)
            TZI(I)=TZI(I+1)
            TZI(I+1)=TZERO    
            IZEROC=INDTZ(I)
            INDTZ(I)=INDTZ(I+1)
            INDTZ(I+1)=IZEROC    
         ENDIF
 1190 CONTINUE
      IF (ICHANG .GT. 0) GOTO 1180
C
      IF (NZH .GT. 0) THEN
         IF (MNSWI .GE. 2) WRITE (LUSWI,*) 
     $        ' - MXRTFD: The roots with the projected variables:'
         IF (MNSWI .GE. 1) THEN
            III=NZH 
            DO 1200 II=1,NZH 
               WRITE (LUSWI,'(A,I5,A,1PD20.10)')
     $              ' - MXRTFD: got zero of switching function no.',
     $              INDTZ(II), ' at t = ', TZI(II)
 1200       CONTINUE
         ENDIF
      ENDIF
C
C  return the interpolated and projected solution at the minimal root
C
      TZERO=TZI(1)
      DO 1210 I=1,NSWIF
         IZERO(I)=0
 1210 CONTINUE
      IZERO(INDTZ(1))=1
C
      T = TZERO
      NS = 0
      DO 1220 I=1,NP
         P(I) = YIP(NS+I)
 1220 CONTINUE
      NS = NP
      DO 1230 I=1,NV
         V(I) = YIP(NS+I)
 1230 CONTINUE
      NS = NS + NV
      DO 1240 I=1,NU
         U(I) = YIP(NS+I)
 1240 CONTINUE
      NS = NS + NU
      DO 1250 I=1,NV
         A(I) = YIP(NS+I)
 1250 CONTINUE
      NS = NS + NV
      DO 1260 I=1,NL
         RLAM(I) = YIP(NS+I)
 1260 CONTINUE
C
      CALL MONON(4)
      IF (NDY .GT. 0) THEN
         CALL MXIPOL (ND, DENS, NDY, TZERO, TOLD, H, ILAM, IRHO, YIP)
      ENDIF
      IF (NDZ .GT. 0) THEN
         CALL MXIPOL (ND, DENS(NDY+1,2-IDEL), NDZ, TZERO, TOLD,
     $        H, ILAMH, IRHO, YIP(NDY+1))
      ENDIF
      CALL MONOFF(4)
C
      CALL MONON(1)
      CALL MXPRJC(NP, NV, NL, NG, NU,
     $     NPDIM, NVDIM, NLDIM, NGDIM, NUDIM, NTODIM,
     $     N2UDIM, N6UDIM, NVLDIM, LDA, LDG, LDF, MDA, MDG, MDF,
     $     FPROB,
     $     TZERO, YIP(IPP), YIP(IPV), YIP(IPU), YIP(IPA), YIP(IPL),
     $     SCAL,
     $     RTOL, ATOL, ITOL,
     $     ITMAXP, NITER, EPSREQ, EPSEST,
     $     WF, WP, WV, WZ, WG, WGI, WPDOT, WUDOT, AM, GP, FL,
     $     MXJOB, RWK, IWK, LRWK, LIWK, IFAIL)
      CALL MONOFF(1)
C
      CALL MONON(11)
      CALL GSWIT (NP, NV, NL, NU, TZERO,
     $     P, V, U, A, RLAM,
     $     NSWIF, GVAL0)
      CALL MONOFF(11)
      NSWIF = MIN (NSWIF, NSWMAX)
C
      RETURN
C
C  end of subroutine MXRTFD
C
      END
