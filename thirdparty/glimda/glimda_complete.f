C**************************************************************************
C
C     GLIMDA version 0.1
C     General LInear Methods for Differential Agebraic equations
C
C     Copyright (C) 2006,  Steffen Voigtmann
C     
C     This program is free software; you can redistribute it and/or
C     modify it under the terms of the GNU Lesser General Public
C     License as published by the Free Software Foundation; either
C     version 2.1 of the License, or (at your option) any later version.
C     
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C     Lesser General Public License for more details.
C
C     You should have received a copy of the GNU Lesser General Public
C     License along with this library; if not, write to the Free Software
C     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
C
C**************************************************************************

C**************************************************************************
C     ****f* src/glimda
C
C  NAME
C     GLIMDA -- General LInear Methods for Differential Algebraic equations
C
C  DESCRIPTION
C     This is the solver GLIMDA using general linear methods (GLMs) to solve
C     differential algebraic equations of the form
C       
C       f( q'(x,t), x, t ) = 0,       f : R^n x R^m x R --> R^m,
C
C     having a properly stated leading term. Currently methods of order
C     1<=p<=3 are used. The methods satisfy
C
C     * s = p, i.e. p internal stages are computed (this is similar to an
C       s-stage Runge-Kutta method)
C     * r = p, i.e. p external stages are computed (this is similar to an
C       r-step linear multistep scheme, the external stages are passed on 
C       from step to step using a Nordsieck vector)
C     * q = p, i.e. the stage order agrees with the order
C
C     Notice that all methods are stiffly accurate. The methods for p=1,2
C     are A- and L-stable while the methods with p=3 are not A-stable. 
C     More details are given in 
C
C     Steffen Voigtmann, General linear methods for integrated circuit design,
C                        PhD thesis, Humboldt University Berlin (2006)
C
C  ARGUMENTS
C     
C     integer  m         .. dimension of f : R^n x R^m x R --> R^m
C     integer  n         .. dimension of q : R^m x R       --> R^n
C     external fevl      .. evaluates F = f(y,x,t)
C     external qevl      .. evaluates Q = q(x,t)
C     external dfyevl    .. computes the matrix A = df/dy
C     external dfxevl    .. computes the matrix B = df/dx
C     external dqxevl    .. computes the matrix D = dq/dx
C     real*8   t0        .. start point of time integration
C     real*8   tend      .. end point of time integration
C     real*8   xin(m)    .. initial value xin = x_0 = x(t0)
C     real*8   qprime(n) .. derivative q'(x,t) 
C     real*8   h0        .. initial stepsize
C     real*8   atolX(m)  .. absolute tolerances for x
C     real*8   rtolX(m)  .. relative tolerances for x
C     logical  tolvec    .. true  => atolX, rtolX are vectors such that
C                                 different tolerances for each component
C                                 are possible
C                           false => only atolX(1), rtolX(1) are relevant,
C                                 i.e. the same tolerances is used for 
C                                 each component: atolX(i)=atolX(1),
C                                 rtolX(i)=rtolX(1)
C     logical  num_A     .. true  => evaluate A = df/dy numerically
C                                 (provide a dummy argument for dfyeval)
C                           false => A is calculated using dfyevl
C     logical  num_D     .. true  => evaluate D = dq/dx numerically
C                                 (provide a dummy argument for dqxeval)
C                           false => D is calculated using dqxevl  
C     logical  num_B     .. true  => evaluate B = df/dx numerically
C                                 (provide a dummy argument for dfxeval)
C                           false => B is calculated using dfxevl  
C     logical  ode       .. true  => problem is an ODE, false => DAE 
C     logical  ADcnst    .. true  => the leading term (A,D) is constant
C                                 (A,D evaluations can be saved)
C                           false => A/D depends on x or t 
C     integer  icount(6) .. statistical informtion about integration
C     integer  iopt(9)   .. integer options for the simulator
C     real*8   ropt(11)  .. real    options for the simulator 
C     integer  ipar(*)   .. additional integer arguments for f, q, ...
C     real*8   rpar(*)   .. additional real    arguments for f, q, ...
C     integer  ierr      .. return code
C
C  COPYRIGHT
C
C  SOURCE

      subroutine glimda(m     , n     , fevl  , qevl  , dfyevl, dfxevl,
     $                  dqxevl, t0    , tend  , xin   , qprime, h0    ,
     $                  atolX , rtolX , tolvec, num_A , num_D , num_B ,
     $                  ode   , ADcnst, icount, iopt  , ropt  , ipar  ,
     $                  rpar  , ierr  , solout)

C#define DEBUG
C
C     parameters
C
      implicit none
      logical   tolvec  , num_A   , num_D    , num_B    , ode    ,ADcnst
      integer   m       , n       , icount(6), iopt(9)  , ipar(*),ierr
      real*8    t0      , tend    , xin(m)   , qprime(n), h0     , 
     $          atolX(m), rtolX(m), ropt(11) , rpar(*)
      external  fevl    , qevl    , dfyevl   , dfxevl   , dqxevl ,solout
C      
C     options of the integrator (see iniopt for details)
C
      integer   MD      , PRTLEV  , NEWTIT, MAXORD, MINORD,
     $          P3METH  , ORDER   , KEEP_P, NSTMAX, NRETRY 
      parameter ( MD=3 )
      real*8    TOLF    , TOLX    , DELACL, CUTFAC, INCFAC, DECFAC,
     $          HMAX    , HMIN    , UPTH  , DOWNTH, RCOND ,
     $          atolQ(n), rtolQ(n)
C
C     statistics of the integration
C
      logical   last
      integer   nsteps, nrefus, nnocon, nfe, njac, nlu, nsol, orders(3)
      real*8    eps
C      
C     integration process
C
      integer   BAKIFO
      parameter ( BAKIFO=4 )
      logical   varord, goup    , godown     , change      , caljac  ,
     $          jacnew, tolfok  , tolxok     
      integer   retry , same_p  , conv       , ipiv(m)     , X_use   ,
     $          X_perm(BAKIFO)
      real*8    h , t , hnew    , mat_A(m,n) , mat_D(n,m)  , Jac(m,m),
     $          Q(n)  , Qp(n,MD), qn_in(n,MD), qn_out(n,MD), fzxt(m) ,
     $          res(m), X(m,MD) , scdrv0(n)  , scdrv1(n)   , scalQ(n),
     $          lam   , scalX(m), omega(n)   , nresx       , nresx1  ,  
     $          psi   , psi1    , theta      , theta1      , err     ,
     $          resold, nrmres  , errvec(n)   ,
     $          X_vals(m,BAKIFO), X_time(BAKIFO)
C
C     the general linear method
C
      integer   p
      real*8    A(MD,MD)  , U(MD,MD)     , B(MD,MD), V(MD,MD), c(MD),
     $          c1beta(MD), delta(MD+1,2)
C
C     keep track of last accepted values
C
      integer   p_acc
      real*8    sd0_ac(n), scaacQ(n), qn_acc(n,MD), h_acc, hnw_ac,   
     $          sd1_ac(n), scaacX(m), qp_acc(n)   
C
C     misc
C
      integer   k, l
      real*8    alpha, fac
C
C     functions
C
      logical      chktol, naninf
      real*8       dlamch, dnrm2
      character*2  ifmt

C#ifdef DEBUG
C     ---------------------------------------------------------------------
C      PRTLEV defined?
C      if (PRTLEV.ge.2) then
C         print*,' ',('-',k=1,69)
C         write(*,*) ' 0     initialise the simulator'
C      endif
C     ---------------------------------------------------------------------
C#endif

      call iniopt(t0    , tend  , iopt  , ropt  , MD    , atolX ,
     $            rtolX , tolvec, m     , PRTLEV, NEWTIT, MAXORD,
     $            MINORD, P3METH, ORDER , KEEP_P, NSTMAX, NRETRY, 
     $            TOLF  , TOLX  , DELACL, CUTFAC, INCFAC, DECFAC, 
     $            HMAX  , HMIN  , UPTH  , DOWNTH, RCOND )
C     
C     trick to treat h0 = -n as h0 = 10**(-n)
C
      if (h0 .le. 0d0) then
         h = 10d0**h0
      else
         h = h0
      endif

      t         = t0
      ierr      = 0
      nsteps    = 0
      nrefus    = 0
      nnocon    = 0
      nfe       = 0
      njac      = 0
      nlu       = 0
      nsol      = 0
      orders(1) = 0
      orders(2) = 0
      orders(3) = 0
      eps       = dlamch( 'Precision' )
      varord    = ( ORDER .eq.0   )
      ADcnst    = ( ADcnst.or.ode )
      LAST      = .false.
C     
C     start with order 1
C
      p = 1
      call getmet(p,MD,A,U,B,V,c,c1beta,delta,P3METH,ierr)
      if (ierr.ne.0) goto 9010
C
C     initialise Q = q(xin,t0)
C
      call qevl(m,n,xin,t,Q,rpar,ipar,ierr)
      if (ierr.ne.0) goto 9090
C
C     in case of a constant leading term, evaluate A=df/dy and 
C     D=df/dx only once
C
      if ( ADcnst ) then
         if ( ODE ) then
            call dcopy(m*n,  0d0, 0, mat_A,   1)
            call dcopy(n  , -1d0, 0, mat_A, m+1)
            call dcopy(n*m,  0d0, 0, mat_D,   1)
            call dcopy(m  ,  1d0, 0, mat_D, n+1)
         else
            if ( num_A ) then
               call fevl(m, n, qprime, xin, t, fzxt, rpar, ipar, ierr)
               nfe = nfe + 1
            endif
C
C     .     the user has to make sure that qprime is initialised with
C     .     sth. meaningful
C
            call get_A(m   , n   , qprime, xin, t   , num_A, dfyevl, 
     $                 fevl, fzxt, mat_A , eps, ipar, rpar , PRTLEV, 
     $                 ierr)
            if (ierr.ne.0) goto 9080
            call get_D(m    , n  , xin , t   , num_D , dqxevl, qevl, Q,
     $                 mat_D, eps, ipar, rpar, PRTLEV, ierr)
            if (ierr.ne.0) goto 9080
         endif
      endif
C
C     compute tolerances for q(x,t) (atolQ, rtolQ) from atolX, rtolX
C     note: for (ADcnst.eq..FALSE.) these tolerances should be adapted
C           during integration, this is not yet done.
C
      call tolX2Q(m,n,atolX,rtolX,t,atolQ,rtolQ,num_D,ADcnst,
     $            mat_D,dqxevl,qevl,Q,rpar,ipar,PRTLEV,eps,ierr)
      if (ierr.ne.0) goto 9070
C
C     compute initial Nordsieck vector. further initialisation
C
      call dcopy(n,Q     ,1,qn_in ,1)
      call dcopy(n,qprime,1,scdrv0,1)
      call dscal(n,h       ,scdrv0,1)
      call dcopy(n,0d0   ,0,scdrv1,1)
      goup   = .false.
      godown = .false.
      change = .false.
      caljac = .true.
      jacnew = .false.
      do k=1,n
         scalQ(k) = atolQ(k) + rtolQ(k)*abs(Q(k))
      enddo
      do k=1,m
         scalX(k) = atolX(k) + rtolX(k)*abs(xin(k))
      enddo
C
C     collect information about last accepted timepoint
C
      p_acc  = 1
      h_acc  = h
      hnw_ac = h
      retry  = 0
      same_p = 0
      call dcopy(n,qn_in ,1,qn_acc,1) 
      call dcopy(n,qprime,1,qp_acc,1)
      call dcopy(n,scdrv0,1,sd0_ac,1)
      call dcopy(n,scdrv1,1,sd1_ac,1)
      call dcopy(n,scalQ  ,1,scaacQ,1)
      call dcopy(m,scalX  ,1,scaacX,1)
C
C     prepare data structures predictor: interpolation is used such that we
C     need to keep track of the last BAKIFO stages X
C
      do k=1,BAKIFO
         X_perm(k) = k
      enddo
      X_time(1) = t
      call dcopy(m,xin,1,X_vals,1)
      X_use = 1

C#ifdef DEBUG
C     ---------------------------------------------------------------------
      if (PRTLEV.ge.2) then
         write(*,*) ' 1     start time integration'
      endif
C     ---------------------------------------------------------------------
C#endif

      do while (.true.)

 1000    continue

         nsteps = nsteps + 1
C        too many steps?
         if (nsteps .ge. NSTMAX) goto 9020
C        too many consecutive failures?
         if (retry  .gt. NRETRY) goto 9030 
C        stepsize too small?
         if ((h.lt.HMIN).or.(t+h.le.t).or.(h.lt.1d1*eps*t)) goto 9040
C        check for last timestep
         if (abs((t+h-tend)/(tend-t0)).le.eps) LAST = .true.
         
C#ifdef DEBUG
         if (PRTLEV.ge.2) then
            print*,' ',('-',k=1,69)
C     .  ------------------------------------------------------------------
            write(*,*) ' 2     take a step from t to t+h'
C     .  ------------------------------------------------------------------
            write(*,'(8x,''step no. '','//ifmt(nsteps)//
     $        ','',  order p = '',I1,'',  retry ='',I2'//
     $        ','',  last = '',L1)') nsteps, p, retry, last
            write(*,'(8x,1P,''t ='',E15.8,'',     t+h ='',E15.8)')
     $        t, t+h
            write(*,'(8x,1P,''h ='',E15.8,'',   to go ='',E15.8)')
     $        h, tend-(t+h)
         endif
C#endif


C#ifdef DEBUG
C     .  ------------------------------------------------------------------
         if (PRTLEV.ge.2) then
            write(*,*) ' 2.1   change the order (if appropriate)'
         endif
C     .  ------------------------------------------------------------------
C#endif
C
C     .  the order is allowed to change provided that
C     .    (i)   the last step was accepted (change.eq..true.)
C     .    (ii)  KEEP_P steps have been computed with the current order
C     .    (iii) the stepsize has settled
C
C     .  in case of constant order computations the order is built
C     .  up quicker by discarding (iii)
C
         if ( change.and.(same_p.ge.KEEP_P) ) then
            
            if ( (.not.varord).and.(p.lt.ORDER) ) then

               goup = .true.

            elseif ( varord.and.abs((h-h_acc)/h_acc).lt.0.2d0 ) then

               if ( (psi.lt.UPTH).and.(p.lt.MAXORD) ) then
                  goup = .true.
               elseif ( (psi.gt.DOWNTH).and.(p.gt.MINORD) ) then
                  godown = .true.
               endif

            endif
         endif
C
C     .  change the order (if required)
C
         if ( goup ) then

            call purify('p',qn_in,n,p,scdrv1,c1beta)
            call dcopy(n,scdrv0,1,qn_in(1,p+1),1)
            p = p + 1
            call getmet(p,MD,A,U,B,V,c,c1beta,delta,P3METH,ierr)
            if (ierr.ne.0) goto 9010
            call dcopy(n,scdrv1,1,scdrv0,1)
C           no approximation for the new scdrv1 available
            call dcopy(n,0d0   ,0,scdrv1,1)
C           modify not possible since we are lacking scdrv1
C           call purify('m',yn_in,m,p,scdrv1,c1beta)
            same_p  = 0
            caljac  = .true.

         elseif ( godown ) then

            call purify('p',qn_in,n,p,scdrv1,c1beta)
            call dcopy(n,scdrv0    ,1,scdrv1,1)
            call dcopy(n,qn_in(1,p),1,scdrv0,1)
            p = p - 1
            call getmet(p,MD,A,U,B,V,c,c1beta,delta,P3METH,ierr)
            if (ierr.ne.0) goto 9010
            call purify('m',qn_in,n,p,scdrv1,c1beta)
            same_p  = 0
            caljac  = .true.

         endif  

C#ifdef DEBUG
         if (PRTLEV.ge.2) then
            if (goup.or.godown) then
               write(*,*) '       ... changed the order: up =',
     $           goup,', down =',godown,', new order p = ',p 
            else
               write(*,*) '       ... no order change' 
            endif
         endif
C#endif
         goup   = .false.
         godown = .false.
         change = .false.


C#ifdef DEBUG
C        ------------------------------------------------------------------
         if (PRTLEV.ge.2) then
            write(*,*) ' 2.2   solve for the stages (sequentially)'
         endif
C        ------------------------------------------------------------------
C#endif
C
C     .  evaluate Jacobian at least once per step
C
         caljac = .true.
C
C     .  reset convergence rate
C
         psi    = 0d0
C
C     .  loop over all s = p stages
C
         do k=1,p

C#ifdef DEBUG
            if (PRTLEV.ge.2) then
 2006          format(1x,a,I1,a)
C     .     ---------------------------------------------------------------
               write(*,2006) ' 2.',k+2,'   solving for stage no. '//
     $           char(k+ICHAR('0'))
            endif
C     .     ---------------------------------------------------------------
C#endif
C     
C     .     save the diagonal element
C     
            lam = A(k,k)
            
C#ifdef DEBUG
C     .     ---------------------------------------------------------------
            if (PRTLEV.ge.2) then
               write(*,2006) ' 2.',k+2,'.1 computing omega'
            endif
C     .     ---------------------------------------------------------------
C#endif
C
C           this part of the residual is independent of the current stage X_k
C             omega = h*Qp(:,1:k-1)*m.A(k,1:k-1)' + qn_in*m.U(k,:)'
C
            call dgemv('n',n,p,1d0,qn_in,n,U(k,1),MD,0d0,omega,1)
            call dgemv('n',n,k-1,h,Qp,n,A(k,1),MD,1d0,omega,1)

C#ifdef DEBUG
C           ---------------------------------------------------------------
            if (PRTLEV.ge.2) then
               write(*,2006) ' 2.',k+2,'.2 predict stage value'
            endif
C           ---------------------------------------------------------------
C#endif
C
C           predict a value for X_k (interpolation of X_use+k-1 past stages:
C           X_use old (accepted) stages are saved in X_vals, k-1 more are 
C           already available from the current step
C
            call prdicX(m, k-1, t, h, c, X, X_vals, BAKIFO, X_time, 
     $                  X_use, X_perm,X(1,k))
C
C           compute the corresponding value for Q = q(X,t)
C
            call qevl(m, n, X(1,k), t+c(k)*h, Q, rpar, ipar, ierr)
C
C           predict a stage derivative Qp0 = Qp(:,k) = (1/h/lam)*(Q-omega)
C
            call dcopy(n, Q, 1, Qp(1,k), 1)
            call dscal(n,  1d0/h/lam          , Qp(1,k), 1)
            call daxpy(n, -1d0/h/lam, omega, 1, Qp(1,k), 1)
            
 3000       continue

C#ifdef DEBUG
C           ---------------------------------------------------------------
            if (PRTLEV.ge.2) then
               write(*,2006) ' 2.',k+2,'.3 Newton iteration'
            endif
C           ---------------------------------------------------------------
C#endif
            conv   = 0
            nresx  = 0d0
            theta  = 0d0
            resold = -1d0
            
            do l = 1,NEWTIT

               if (psi.gt.1d-3.and..not.jacnew) then
C     .           poor convergence ... update jacobian
                  psi    = 0d0
                  caljac = .true.
               endif
C
C     .        evaluate f in order to calculate the residual
C
               call fevl(m, n, Qp(1,k), X(1,k), t+c(k)*h, fzxt,
     $                   rpar, ipar, ierr)
               nfe = nfe + 1
C
C     .        check for errors during the evalutation of f
C
               if ( naninf(fzxt,m,1,m) ) then
                  if ( PRTLEV.ge.2 ) then
                     print*,'GLIMDA WARNING: f(z,x,t) yields NAN or INF'
                  endif
                  ierr = -10
               endif
               if ( ierr.ne.0 ) then
                  if (PRTLEV.ge.2) then
                     print*, 'GLIMDA WARNING: could not evaluate '
     $                    //'f(z,x,t) ... treating this as noconv'
                  endif
                  conv = 0
                  ierr = 0
                  jacnew = .true.
                  goto 4000
               end if  
C  
C     .        compute (scaled) residual res = -h*lam*fzxt
C
               call dcopy(m,fzxt,1,res,1)
               call dscal(m,-h*lam,res,1)
               nrmres = dnrm2(m,res,1)
C
C     .        check norm of the residual
C
               if ((resold.gt.0d0).and.(nrmres.ge.resold)) then
C#ifdef DEBUG
                  if (PRTLEV.ge.2) then
                     print *, '=> residual is growing '
     $                     //'(',nrmres,resold,')'
                  endif
C#endif
                  conv = 0
                  ierr = 0
                  jacnew = .true.
                  goto 4000
               endif

               tolfok =   chktol(res,m,scaacX,TOLF)

C#ifdef DEBUG
               if (PRTLEV.ge.2) then
                  write(*,'(8x,1P,''|res| ='',E15.8,'',   |resold|='//
     $              ',''E15.8)') nrmres, resold
               endif
C#endif

               resold = nrmres

C
C              evaluate and decompose the jacobian: Jac = AD+h*lam*B
C
               if (caljac) then
C#ifdef DEBUG
                  if (PRTLEV.ge.2) then
                     write(*,'(8x,a)') 'compute jacobian'
                  endif
C#endif
                  call get_j(m     , n     , Qp(1,k),X(1,k), t+c(k)*h,
     $                       h     , lam   , num_A , num_D , num_B   ,
     $                       dfyevl, dqxevl, dfxevl, fevl  , qevl    ,
     $                       ADcnst, mat_A , mat_D , fzxt  , eps     ,
     $                       Jac   , PRTLEV, ipar  , rpar  , ierr    )
                  njac = njac + 1
C
C     .           check for errors during computation of Jac
C
                  if (ierr.ne.0) then
                     if (PRTLEV.ge.2) then
                        print*, 'GLIMDA WARNING: could not evaluate '
     $                       //'f(z,x,t) ... treating this as noconv'
                     endif
                     conv = 0
                     ierr = 0
                     jacnew = .true.
                     goto 4000
                  endif       
C
C                 LU decomposition: Jac = P * L * U
C
                  call dgetrf(m, m, Jac, m, ipiv, ierr)
                  nlu = nlu + 1
                  if (ierr.ne.0) goto 9050
C
C                 estimate condition number (if required)
C
                  if (PRTLEV.ge.2) call chkrcd(Jac,m,RCOND)
                  jacnew = .true.
                  caljac = .false.

               endif

C
C     .        solve the linear system Jac * dX = res (Jac = P*L*U)
C
               call dgetrs('n',m,1,Jac,m,ipiv,res,m,ierr)
               nsol = nsol + 1
               if (ierr.ne.0) goto 9060
C
C     .        check accuracy of X update
C
               tolxok = chktol(res,m,scaacX,TOLX)
C
C     .        update stage: X = res + dX
C
               call daxpy(m,1d0,res,1,X(1,k),1)
C
C     .        compute the corresponding value for Q = q(X,t)
C
               call qevl(m, n, X(1,k), t+c(k)*h, Q, rpar, ipar, ierr)
C
C              compute the corresponding stage derivative 
C                Qp(:,k)=(1/h/lam)*(Q-omega)
C
               call dcopy(n, Q, 1, Qp(1,k), 1)
               call dscal(n,  1d0/h/lam          , Qp(1,k), 1)
               call daxpy(n, -1d0/h/lam, omega, 1, Qp(1,k), 1)
C
C              estimate the rate of convergence (for changing the order)
C
               call convrt(psi,l,res,m,scaacX,nresx,nresx1,theta,theta1)
               if (psi.ge.1d0) then
C#ifdef DEBUG
                  if (PRTLEV.ge.2) then
                     print*,'   psi = ',psi,' > 1'
                  endif
C#endif
                  conv = 0
                  ierr = 0
                  goto 4000
               endif 

C#ifdef DEBUG
 5             format(8x,a,i2,L2,L2,1x,a,E10.3,1x,a,E10.3,1x,
     $              a,E10.3)
               if (PRTLEV.ge.2) then
                  write(*,5) 'it=',l,tolfok,tolxok,
     $              'nresx=',nresx,'nresx1=',nresx1,'psi=',psi
               endif
C#endif
C     
C     .        check for convergence 
C     
               if ( tolfok .and. tolxok ) then
                  conv = l
                  psi1 = psi
                  goto 4000
               endif

            enddo
C     .     -----------------------
C     .     end of newton iteration
C     .     -----------------------


 4000       continue


            if (conv.eq.0) then
C#ifdef DEBUG
C     .        ------------------------------------------------------------
               if (PRTLEV.ge.2) then
                  write(*,2006) ' 2.',k+2,'.4 no convergence'
               endif
C     .        ------------------------------------------------------------
C#endif
               nnocon = nnocon + 1

               if (.not.jacnew.and.psi.le.1d0) then
C
C     .           no convergence for wrong jacobian ... update and retry
C
                  caljac = .true.
                  psi    = psi1
                  goto 3000

               endif
C
C     .        force a re-computation of the jacobian
C
               caljac = .true.
               LAST   = .false.
C
C              no convergence in spite of updated jacobian ... redo step 
C
               if ( (retry.eq.0).and.(p_acc.ne.p) ) then
C
C     .           no convergence after order change ... retry with old order
C
C#ifdef DEBUG
                  if (PRTLEV.ge.2) then
                     print*,'       ... use old order p =',p_acc
                  endif
C#endif
                  h = hnw_ac
                  p = p_acc
                  call getmet(p,MD,A,U,B,V,c,c1beta,delta,P3METH,
     $                        ierr)
                  if (ierr.ne.0) goto 9010
C     .           scale and modify
                  call dcopy(n*p,qn_acc,1,qn_in ,1)
                  call dcopy(n  ,sd0_ac,1,scdrv0,1)
                  call dcopy(n  ,sd1_ac,1,scdrv1,1)
                  call sclmod(qn_in,n,p,h,h_acc,scdrv0,scdrv1,c1beta)

               else
C
C                 same order or consecutive no convergence 
C                 ... reduce stepsize h and order p
C
                  hnew = CUTFAC * h
C#ifdef DEBUG
                  if (PRTLEV.ge.2) then
                     print*,'       ... reduce h: hnew =',hnew
                  endif
C#endif
                  if (varord.and.(p.gt.MINORD)) godown = .true.
C
C     .           scale and modify
C
                  call sclmod(qn_in,n,p,hnew,h,scdrv0,scdrv1,c1beta)
                  h = hnew

               endif

               retry = retry + 1
               goto 1000

            else

C#ifdef DEBUG
C     .        ------------------------------------------------------------
               if (PRTLEV.ge.2) then
                  write(*,2006) ' 2.',k+2,'.4 convergence'
               endif
C     .        ------------------------------------------------------------
C#endif
C
C     .        stage number k converged ... compute next stage
C
               jacnew = .false.

            endif
           
         enddo
C     .  -----------------------------
C     .  end of solving for the stages
C     .  -----------------------------


C#ifdef DEBUG
C     .  ------------------------------------------------------------------
         if (PRTLEV.ge.2) then
            write(*,*) ' 3.1   all stages converged ... compute qn_out'
         endif
C     .  ------------------------------------------------------------------
C#endif

C        compute output  qn_out = h*Qp*B' + qn_in * V'
         call dgemm('n','t',n,p,p,1d0,qn_in,n,V,MD,0d0,qn_out,n)
         call dgemm('n','t',n,p,p,h  ,Qp   ,n,B,MD,1d0,qn_out,n)      

C#ifdef DEBUG
C     .  ------------------------------------------------------------------
         if (PRTLEV.ge.2) then
            write(*,*) ' 3.2   compute scaled derivatives'
         endif
C     .  ------------------------------------------------------------------
C#endif
         call scdrvs(h,Qp,n,p,qn_in,delta,MD,qp_acc,scdrv0,scdrv1)

C#ifdef DEBUG
C     .  ------------------------------------------------------------------
         if (PRTLEV.ge.2) then
            write(*,*) ' 3.3   estimate local error'
         endif
C     .  ------------------------------------------------------------------
C#endif
C        error estimate
         do k=1,m
            scalX(k)  = atolX(k)+rtolX(k)*abs(X(k,p))
         enddo
         do k=1,n
            alpha     = atolQ(k)+rtolQ(k)*abs(qn_out(k,1))
            scalQ(k)  = alpha
            errvec(k) = scdrv1(k) / scalQ(k) 
         enddo
         err  = max(1d-2,abs(c1beta(1)) * dnrm2(n,errvec,1))
C         err  = abs(c1beta(1)) * dnrm2(n,errvec,1)

Cccccccc test a parameter for a modified error constant cccccccccc
Cc         err  = 1d0 * abs(c1beta(1)) * dnrm2(n,errvec,1)

C#ifdef DEBUG
C     .  ------------------------------------------------------------------
         if (PRTLEV.ge.2) then
            write(*,*) ' 3.4   predict new stepsize'
         endif
C     .  ------------------------------------------------------------------
C#endif
C        predict new stepsize
         fac  = 0.85d0 * err**(-1d0/real(p+1))
C        do not increase by more than INCFAC
         fac  = min(fac,INCFAC)
C        do not decrease by more than DECFAC
         fac  = max(fac,DECFAC)
C        do not use stepsize larger than HMAX
         hnew = min(fac*h,HMAX)

C#ifdef DEBUG
         if (PRTLEV.ge.2) then
            write(*,'(8x,''err ='',1P,E15.8,'',   hnew ='',E15.8)')
     $        err, hnew
         endif
C#endif
C
C     .  check for end of integration interval
C
         if ( (.not.LAST).and.(t+h+1.05*hnew.ge.tend-eps)) then

            hnew = tend - h - t
            LAST = .true.
            caljac = .true.
            
         endif

C#ifdef DEBUG
C     .  ------------------------------------------------------------------
         if (PRTLEV.ge.2) then
            write(*,*) ' 3.5   accept / refuse ?'
         endif
C     .  ------------------------------------------------------------------
C#endif

         if (err*DELACL.gt.1) then
            
C#ifdef DEBUG
C     .     ---------------------------------------------------------------
            if (PRTLEV.ge.2) then
               write(*,'(8x,a,$)') '... refuse'
            endif
C     .     ---------------------------------------------------------------
C#endif
            LAST = .false.
            nrefus = nrefus + 1

            if ( (retry.eq.0).and.(p.ne.p_acc) ) then
C
C     .        refuse after changing the order ... try again with old order
C     
C#ifdef DEBUG
               if (PRTLEV.ge.2) then
                  print*,'... switch back to old order'
               endif
C#endif
               h = hnw_ac
               p = p_acc
               call getmet(p,MD,A,U,B,V,c,c1beta,delta,P3METH,ierr)
               if (ierr.ne.0) goto 9010
C              scale and modify
               call dcopy(n*p,qn_acc,1,qn_in ,1)
               call dcopy(n  ,sd0_ac,1,scdrv0,1)
               call dcopy(n  ,sd1_ac,1,scdrv1,1)
               call sclmod(qn_in,n,p,h,h_acc,scdrv0,scdrv1,c1beta)

            else
C#ifdef DEBUG
               if (PRTLEV.ge.2) then
                  print*,'... carry on with same order and hnew = ',hnew
               endif
C#endif
C$$$  if(retry.eq.5) hnew = max(CUTFAC * hnew, HMIN)
               call dcopy(n,sd0_ac,1        ,scdrv0,1)
               call dcopy(n,sd1_ac,1        ,scdrv1,1)
               call dscal(n,(h/h_acc)**p    ,scdrv0,1)
               call dscal(n,(h/h_acc)**(p+1),scdrv1,1)
               call sclmod(qn_in,n,p,hnew,h,scdrv0,scdrv1,c1beta)
               h = hnew
            endif
            retry = retry + 1
            goto 1000

         else 

C#ifdef DEBUG
C     .     ---------------------------------------------------------------
            if (PRTLEV.ge.2) then
               write(*,'(8x,a)') '... accept'
            endif
C     .     ---------------------------------------------------------------
C#endif
C
C     .     update X_vals used for stage prediction (interpolation)
C
            call updxvl(m,p,t,h,c,BAKIFO,X_use,X_time,X_vals,X_perm,X)
C     
C     .     update status of the simulator
C
            t         = t + h
            orders(p) = orders(p) + 1
            same_p    = same_p + 1  
            change    = .true.
            retry     = 0
            p_acc     = p
            h_acc     = h    
            hnw_ac    = hnew
            call dcopy(n*p,qn_out ,1,qn_acc,1) 
            call dcopy(n  ,Qp(1,p),1,qp_acc,1)
            call dcopy(n  ,scdrv0 ,1,sd0_ac,1)
            call dcopy(n  ,scdrv1 ,1,sd1_ac,1)
            call dcopy(n  ,scalQ  ,1,scaacQ,1)
            call dcopy(m  ,scalX  ,1,scaacX,1)

C            if (PRTLEV.ge.3) then
C     .        write file output
            call solout(PRTLEV,t,h,p,X(1,p),m,Qp(1,p),n)
C            endif

C#ifdef DEBUG
C     .     ---------------------------------------------------------------
            if (PRTLEV.ge.2) then
               write(*,*) ' 3.6   all done?'
            endif
C     .     ---------------------------------------------------------------
C#endif
            if (t.ge.tend-eps) goto 9999

            call dcopy(n*p,qn_out ,1,qn_in,1) 
            call sclmod(qn_in,n,p,hnew,h,scdrv0,scdrv1,c1beta)
            h = hnew
            
         endif 
         
      enddo
C     -----------------------
C     end of time integration
C     -----------------------




C     ---------------------------------------------------------------------
C     error messages
C     ---------------------------------------------------------------------
 9001 format(1x,'GLIMDA ERROR: ',a,i2   ,' (ierr = ',i2,')')
 9002 format(1x,'GLIMDA ERROR: ',a,i6   ,' (ierr = ',i2,')')
 9003 format(1x,'GLIMDA ERROR: ',a,e11.4,' (ierr = ',i2,')')
 9004 format(1x,'GLIMDA ERROR: ',a,' (ierr = ',i4,' [',a,'])')
 9005 format(1x,'GLIMDA ERROR: ',a,' (ierr = ',i2,')')
 9010 if (PRTLEV.eq.1) goto 9999
      write(*,9001) 'could not get method of order ',p,ierr
      goto 9999
 9020 ierr = -2
      if (PRTLEV.eq.1) goto 9999
      write(*,9002) 'too many steps: ',NSTMAX,ierr
      goto 9999
 9030 ierr = -3
      if (PRTLEV.eq.1) goto 9999
      write(*,9001) 'too many consecutive failures: ',NRETRY,ierr
      goto 9999
 9040 ierr = -4
      if (PRTLEV.eq.1) goto 9999
      write(*,9003) 'stepsize too small: ',h,ierr
      goto 9999
 9050 if (PRTLEV.eq.1) goto 9999
      write(*,9004) 'could not compute LU factorisation',ierr,'DGETRF'
      goto 9999
 9060 if (PRTLEV.eq.1) goto 9999
      write(*,9004) 'could not solve the linear system',ierr,'DGETRS'
      goto 9999
 9070 ierr = -7
      if (PRTLEV.eq.1) goto 9999
      write(*,9005) 'could not compute the tolerances atolQ/rtolQ',ierr
      goto 9999
 9080 ierr = -8
      if (PRTLEV.eq.1) goto 9999
      write(*,9005) 'could not evaluate partial derivatives',ierr
      goto 9999
 9090 ierr = -9
      if (PRTLEV.eq.1) goto 9999
      write(*,9005) 'could not evaluate q(x,t)',ierr
      goto 9999

 9999 continue
      

C#ifdef DEBUG
      if (PRTLEV.ge.2) then
C     ---------------------------------------------------------------------
         write(*,*) ' 4. produce some statistical output'
C     ---------------------------------------------------------------------
         print*,' ',('-',k=1,69)
      endif
C#endif

 71   format(3x,a,I6,a,I6,a,I6,a,I6)
 73   format(7x,a,1P,E12.5,'  (min=',E12.5,', max=',E12.5,' )')
      if ( PRTLEV.ge.2 ) then 
         write(*,'(/,2X,''run characteristics'')')
         write(*,71) '   COMPUTED',nsteps,'    REFUSED',nrefus,
     $        '     NOCONV',nnocon
         write(*,71) '      FEVAL',nfe,   '    JACEVAL',njac,
     $        '       NSOL',nsol
         write(*,'(8x,a,3X,''1->'','//ifmt(orders(1))//
     $                ',3X,''2->'','//ifmt(orders(2))//
     $                ',3X,''3->'','//ifmt(orders(3))//')')
     $        'ORDERS',orders
         write(*,73) 'AVERAGE   h=',
     $        (tend-t0)/(nsteps-nrefus-nnocon),HMIN,HMAX
         print*
      endif
      icount(1) = nsteps
      icount(2) = nrefus
      icount(3) = nnocon
      icount(4) = nfe   
      icount(5) = njac  
      icount(6) = nlu   

C     copy output values
      t0 = t
      call dcopy(m,X(1,p),1,xin,1)
      call dcopy(n,Qp(1,p),1,qprime,1)
      call dscal(n,1d0/h,qprime,1)

      call sclmod(qn_in,n,p,h0,h,scdrv0,scdrv1,c1beta)

      end


C**************************************************************************
C     ****v* src/m
C  NAME
C     m
C     
C  DESCRIPTION
C     integer  m         .. dimension of f : R^n x R^m x R --> R^m
C                           note: x \in R^m, but q \in R^n
C
C**************************************************************************
C     ****v* src/n
C  NAME
C     n
C     
C  DESCRIPTION
C     integer  n         .. dimension of q : R^m x R       --> R^n
C                           note: x \in R^m, but q \in R^n
C
C**************************************************************************
C     ****v* src/t0
C  NAME
C     t0
C     
C  DESCRIPTION
C     real*8   t0        .. start point of time integration
C
C**************************************************************************
C     ****v* src/tend
C  NAME
C     tend
C     
C  DESCRIPTION
C     real*8   tend      .. end point of time integration
C
C**************************************************************************
C     ****v* src/xin
C  NAME
C     xin
C     
C  DESCRIPTION
C     real*8   xin(m)    .. initial value xin = x_0 = x(t0), during 
C                           integration xin represents the value of the
C                           dependent variable x at the beginning of a step
C
C**************************************************************************
C     ****v* src/qprime
C  NAME
C     qprime
C     
C  DESCRIPTION
C     real*8   qprime(n) .. derivative q'(x,t) at time t
C
C**************************************************************************
C     ****v* src/h0
C  NAME
C     h0
C     
C  DESCRIPTION
C     real*8   h0        .. initial stepsize
C
C**************************************************************************
C     ****v* src/atolX
C  NAME
C     atolX
C     
C  DESCRIPTION
C     real*8   atolX(m)  .. absolute tolerances for x
C
C**************************************************************************
C     ****v* src/rtolX
C  NAME
C     rtolX
C     
C  DESCRIPTION
C     real*8   rtolX(m)  .. relative tolerances for x
C
C**************************************************************************
C     ****v* src/tolvec
C  NAME
C     tolvec
C     
C  DESCRIPTION
C     logical  tolvec    .. true  => atolX, rtolX are vectors such that
C                                 different tolerances for each component
C                                 are possible
C                           false => only atolX(1), rtolX(1) are relevant,
C                                 i.e. the same tolerances is used for 
C                                 each component: atolX(i)=atolX(1),
C                                 rtolX(i)=rtolX(1)
C
C**************************************************************************
C     ****v* src/num_A
C  NAME
C     num_A
C     
C  DESCRIPTION
C     logical  num_A     .. true  => evaluate A = df/dy numerically
C                                 (provide a dummy argument for dfyeval)
C                           false => A is calculated using dfyevl
C
C**************************************************************************
C     ****v* src/num_D
C  NAME
C     num_D
C     
C  DESCRIPTION
C     logical  num_D     .. true  => evaluate D = dq/dx numerically
C                                 (provide a dummy argument for dqxeval)
C                           false => D is calculated using dqxevl  
C
C**************************************************************************
C     ****v* src/num_B
C  NAME
C     num_B
C     
C  DESCRIPTION
C     logical  num_B     .. true  => evaluate B = df/dx numerically
C                                 (provide a dummy argument for dfxeval)
C                           false => B is calculated using dfxevl  
C
C**************************************************************************
C     ****v* src/ode
C  NAME
C     ode
C     
C  DESCRIPTION
C     logical  ode       .. true  => problem is an ODE
C                           false => problem is a  DAE 
C
C**************************************************************************
C     ****v* src/ADcnst
C  NAME
C     ADcnst
C     
C  DESCRIPTION
C     logical  ADcnst    .. true  => the leading term (A,D) is constant
C                                 (A,D evaluations can be saved)
C                           false => A/D depends on x or t 
C
C**************************************************************************
C     ****v* src/icount
C  NAME
C     icount
C     
C  DESCRIPTION
C     integer  icount(6) .. statistical informtion about integration
C
C              icount(1) =  number of steps taken
C              icount(2) =  number of refused steps
C              icount(3) =  number of steps without convergence
C              icount(4) =  number of f-evaluations   
C              icount(5) =  number of jacobian evaluations  
C              icount(6) =  number of LU decompositions
C
C**************************************************************************
C     ****v* src/iopt
C  NAME
C     iopt
C     
C  DESCRIPTION
C     integer  iopt(9)   .. integer options for the simulator
C
C              iopt(1)   -> PRTLEV [-P], default: 2
C              iopt(2)   -> NEWTIT [-n], default: 5 
C              iopt(3)   -> MAXORD [-O], default: 3 
C              iopt(4)   -> MINORD [-o], default: 1 
C              iopt(5)   -> P3METH [-3], default: 1 
C              iopt(6)   -> ORDER  [-p], default: 0 
C              iopt(7)   -> KEEP_P [-k], default: 5 
C              iopt(8)   -> NSTMAX [-N], default: 500000 
C              iopt(9)   -> NRETRY [-r], default: 15 
C
C**************************************************************************
C     ****v* src/ropt
C  NAME
C     ropt
C     
C  DESCRIPTION
C     real*8   ropt(11)  .. floating point options for the simulator 
C
C              ropt( 1)  -> TOLF   [-F], default: 0.1
C              ropt( 2)  -> TOLX   [-X], default: 0.1
C              ropt( 3)  -> DELACL [-a], default: 0.1
C              ropt( 4)  -> CUTFAC [-c], default: 0.1
C              ropt( 5)  -> INCFAC [-i], default: 2.0
C              ropt( 6)  -> DECFAC [-d], default: 0.5
C              ropt( 7)  -> HMAX   [-M], default: abs(tend-t0)
C              ropt( 8)  -> HMIN   [-m], default: dlamch( 'Precision' )
C              ropt( 9)  -> UPTH   [-U], default: 1d-3
C              ropt(10)  -> DOWNTH [-D], default: 0.8
C              ropt(11)  -> RCOND  [-C], default: 0.0
C
C**************************************************************************
C     ****v* src/ipar
C  NAME
C     ipar
C     
C  DESCRIPTION
C     integer  ipar(*)   .. additional integer arguments for f, q, ...
C                           
C  SEE ALSO
C     fevl, qevl, dfyevl, dfxevl, dqxevl
C
C**************************************************************************
C     ****v* src/rpar
C  NAME
C     rpar
C     
C  DESCRIPTION
C     real*8   rpar(*)   .. additional floating point arguments for f, q, ...
C
C  SEE ALSO
C     fevl, qevl, dfyevl, dfxevl, dqxevl
C
C**************************************************************************
C     ****v* src/ierr
C  NAME
C     ierr
C     
C  DESCRIPTION
C     integer  ierr      .. return code
C
C              ierr      =   0  computation successful
C              ierr      =  -1  could not get method of order p
C              ierr      =  -2  too many steps, increase NSTMAX
C              ierr      =  -3  too many consecutive failures, increase NRETRY
C              ierr      =  -4  stepsize too small, decrease HMIN
C              ierr      =   ?  could not compute LU factorisation
C                               (DGETRF's error code is returned)
C              ierr      =   ?  could not solve the linear system
C                               (DGETRS's error code is returned)
C              ierr      =  -7  could not compute the tolerances atolQ/rtolQ
C              ierr      =  -8  could not evaluate parial derivatives
C              ierr      =  -9  could not evaluate q(x,t)
C
C**************************************************************************
C     ****v* src/TOLF
C  NAME
C     TOLF
C     
C  DESCRIPTION
C     real*8   TOLF   [-F] .. controls the error in the residuum for
C                             Newton's iteration: |res|_err < TOLF is 
C                             required for convergence
C                             |.|_err is a weighted norm incorporating
C                             both absolute and reative tolrance requirements
C
C                             default: 0.1
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/TOLX
C  NAME
C     TOLX
C     
C  DESCRIPTION
C     real*8   TOLX   [-X] .. controls the error in the Newton update for
C                             Newton's iteration: |dX|_err < TOLX is 
C                             required for convergence
C                             |.|_err is a weighted norm incorporating
C                             both absolute and reative tolrance requirements
C
C                             default: 0.1
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/DELACL
C  NAME
C     DELACL
C     
C  DESCRIPTION
C     real*8   DELACL [-a] .. controls whether a given step is accepted:
C                             |DELACL*err_const*scdrv1|_err < 1 is required
C                             for acceptance
C                             |.|_err is a weighted norm incorporating
C                             both absolute and reative tolrance requirements
C
C                             default: 0.1
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/CUTFAC
C  NAME
C     CUTFAC
C     
C  DESCRIPTION
C     real*8   CUTFAC [-c] .. in case of 'no convergence' use the stepsize
C                             CUTFAC*h for the next step           
C
C                             default: 0.1
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/INCFAC
C  NAME
C     INCFAC
C     
C  DESCRIPTION
C     real*8   INCFAC [-i] .. do not increase stepsize by no more than INCFAC
C
C                             default: 2.0
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/DECFAC
C  NAME
C     DECFAC
C     
C  DESCRIPTION
C     real*8   DECFAC [-d] .. do not decrease stepsize by no more than DECFAC
C
C                             default: 0.5
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/HMAX
C  NAME
C     HMAX
C     
C  DESCRIPTION
C     real*8   HMAX   [-M] .. maximal stepsize to use 
C
C                             default: abs(tend-t0)
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/HMIN
C  NAME
C     HMIN
C     
C  DESCRIPTION
C     real*8   HMIN   [-m] .. minimal stepsize to use
C
C                             default: dlamch( 'Precision' )
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/UPTH
C  NAME
C     UPTH
C     
C  DESCRIPTION
C     real*8   UPTH   [-U] .. threshold for the convergence rate when 
C                             changing the order upwards
C
C                             default: 1d-3
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/DOWNTH
C  NAME
C     DOWNTH
C     
C  DESCRIPTION
C     real*8   DOWNTH [-D] .. threshold for the convergence rate when 
C                             changing the order downwards
C
C                             default: 0.8
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/RCOND
C  NAME
C     RCOND
C     
C  DESCRIPTION
C     real*8   RCOND  [-C] .. the minimum (inverse of the) condition number
C                             that is allowed during computation 
C                             (RCOND=0 skips this test)
C
C                             default: 0.0
C  SEE ALSO
C     ropt
C
C**************************************************************************
C     ****v* src/PRTLEV
C  NAME
C     PRTLEV
C     
C  DESCRIPTION
C     integer  PRTLEV [-P] .. print level during integration
C
C              PRTLEV      =  0   no output during integration
C              PRTLEV      =  1   print error messages
C              PRTLEV      =  2   print warnings as well
C              PRTLEV      =  3   produce file output as well
C
C                             default: 2
C  SEE ALSO
C     iopt
C
C**************************************************************************
C     ****v* src/NEWTIT
C  NAME
C     NEWTIT
C     
C  DESCRIPTION
C     integer  NEWTIT [-n] .. maximum number of newton iterations
C
C                             default: 5 
C  SEE ALSO
C     iopt
C
C**************************************************************************
C     ****v* src/MAXORD
C  NAME
C     MAXORD
C     
C  DESCRIPTION
C     integer  MAXORD [-O] .. maximum order p to use
C
C                             default: 3 
C  SEE ALSO
C     iopt
C
C**************************************************************************
C     ****v* src/MINORD
C  NAME
C     MINORD
C     
C  DESCRIPTION
C     integer  MINORD [-o] .. minimum order p to use 
C
C                             default: 1 
C  SEE ALSO
C     iopt
C
C**************************************************************************
C     ****v* src/P3METH
C  NAME
C     P3METH
C     
C  DESCRIPTION
C     integer  P3METH [-3] .. the particular order 3 method to use
C
C              P3METH      =  1   use a method with s=r=p=q=3 and 
C                                 A(alpha) stability for alpha=74°
C                                 but being unconditionally stable
C                                 at zero and at infinity
C              P3METH      =  2   use a method with s=r=p=q=3 and 
C                                 A(alpha) stability for alpha=88.1°
C
C                             default: 1 
C  SEE ALSO
C     iopt
C
C**************************************************************************
C     ****v* src/ORDER
C  NAME
C     ORDER
C     
C  DESCRIPTION
C     integer  ORDER  [-p] .. ORDER \in {1,2,3} => use constant order
C                             ORDER = 0         => use variable order
C
C                             default: 0 
C  SEE ALSO
C     p, iopt
C
C**************************************************************************
C     ****v* src/KEEP_P
C  NAME
C     KEEP_P
C     
C  DESCRIPTION
C     integer  KEEP_P [-k] .. an order increase is possible only after 
C                             completing KEEP_P steps with the current
C                             order
C
C                             default: 5  
C  SEE ALSO
C     same_p, iopt
C
C**************************************************************************
C     ****v* src/NSTMAX
C  NAME
C     NSTMAX
C     
C  DESCRIPTION
C     integer  NSTMAX [-N]    maximum number of steps to calculate
C
C                             default: 500000 
C  SEE ALSO
C     nsteps, iopt
C
C**************************************************************************
C     ****v* src/NRETRY
C  NAME
C     NRETRY
C     
C  DESCRIPTION
C     integer  NRETRY [-r] .. maximum number of consecutive retries 
C
C                             default: 15 
C  SEE ALSO
C     iopt
C
C**************************************************************************
C     ****v* src/atolQ
C  NAME
C     atolQ
C     
C  DESCRIPTION
C     real*8   atolQ(n)  .. absolute tolerances for q=q(x,t)
C                           atolQ is computes from atolQ using tolX2Q  
C
C**************************************************************************
C     ****v* src/rtolQ
C  NAME
C     rtolQ
C     
C  DESCRIPTION
C     real*8   rtolQ(n)  .. relative tolerances for q=q(x,t)
C                           atolQ is computes from atolQ using tolX2Q  
C
C**************************************************************************
C     ****v* src/MD
C  NAME
C     MD
C     
C  DESCRIPTION
C     integer  MD        .. Maximal Dimension of the A,U,B,V matrices
C                           since methds with s=r=p=q are used, this
C                           number coincides with the maximum allowed order
C
C                           parameter ( MD=3 )
C
C**************************************************************************
C     ****v* src/last
C  NAME
C     last
C     
C  DESCRIPTION
C     logical  last      .. true , if the current step reaches tend
C                           false, otherwise
C
C**************************************************************************
C     ****v* src/nsteps
C  NAME
C     nsteps
C     
C  DESCRIPTION
C     integer  nsteps    .. the number of steps already taken
C
C  SEE ALSO
C     icount
C
C**************************************************************************
C     ****v* src/nrefus
C  NAME
C     nrefus
C     
C  DESCRIPTION
C     integer  nrefus    .. the number of refused steps
C
C  SEE ALSO
C     icount
C
C**************************************************************************
C     ****v* src/nnocon
C  NAME
C     nnocon
C     
C  DESCRIPTION
C     integer  nnocon    .. the number of steps without convergence
C
C  SEE ALSO
C     icount
C
C**************************************************************************
C     ****v* src/nfe
C  NAME
C     nfe
C     
C  DESCRIPTION
C     integer  nfe       .. the number of f-evaluations
C
C  SEE ALSO
C     icount
C
C**************************************************************************
C     ****v* src/njac
C  NAME
C     njac
C     
C  DESCRIPTION
C     integer  njac      .. the number of Jacobian-evaluations 
C
C  SEE ALSO
C     icount
C
C**************************************************************************
C     ****v* src/nlu
C  NAME
C     nlu
C     
C  DESCRIPTION
C     integer  nlu       .. the number of LU decompositions
C
C  SEE ALSO
C     icount
C
C**************************************************************************
C     ****v* src/nsol
C  NAME
C     nsol
C     
C  DESCRIPTION
C     integer  nsol      .. the number of linear systems that were solved
C
C  SEE ALSO
C     icount
C
C**************************************************************************
C     ****v* src/orders
C  NAME
C     orders
C     
C  DESCRIPTION
C     integer  orders(3) .. orders(p) holds the number of steps taken
C                           with order p
C
C**************************************************************************
C     ****v* src/eps
C  NAME
C     eps
C     
C  DESCRIPTION
C     real*8   eps       .. eps = dlamch( 'Precision' ) represents the
C                           machine precision
C
C**************************************************************************
C     ****v* src/BAKIFO
C  NAME
C     BAKIFO
C     
C  DESCRIPTION
C     integer   BAKIFO   .. use BAKIFO timepoints for predicting stage values
C
C                           parameter ( BAKIFO=4 )
C  SEE ALSO
C     polint, prdicX, updxvl
C
C**************************************************************************
C     ****v* src/varord
C  NAME
C     varord
C     
C  DESCRIPTION
C     logical  varord    .. true  => perform variable order computation
C                           false => perform constant order computation
C  SEE ALSO
C     ORDER
C
C**************************************************************************
C     ****v* src/goup
C  NAME
C     goup
C     
C  DESCRIPTION
C     logical  goup      .. true  => switch order upwards
C                           false => do not switch order upwards
C
C**************************************************************************
C     ****v* src/godown
C  NAME
C     godown
C     
C  DESCRIPTION
C     logical  godown    .. true  => switch order downwards
C                           false => do not switch order downwards
C
C**************************************************************************
C     ****v* src/change
C  NAME
C     change
C     
C  DESCRIPTION
C     logical  change    .. true  => an order change might be possible
C                                    (after accepting a step)
C                           false => changing the order is not allowed
C
C**************************************************************************
C     ****v* src/caljac
C  NAME
C     caljac
C     
C  DESCRIPTION
C     logical  caljac    .. true  => update the Jacobian as soon as possible
C                           false => do not update the Jacobian
C
C**************************************************************************
C     ****v* src/jacnew
C  NAME
C     jacnew
C     
C  DESCRIPTION
C     logical  jacnew    .. true  => the Jacobian was evaluated for the 
C                                    current timepoint
C                           false => the Jacobian belongs to some previous
C                                    timepoint
C
C**************************************************************************
C     ****v* src/tolfok
C  NAME
C     tolfok
C     
C  DESCRIPTION
C     logical  tolfok    .. true  => the residuum is sufficiently small in
C                                    Newton's iteration
C                           false => no convergence (yet) due to a large
C                                    residuum
C  SEE ALSO
C     TOLF, TOLX, tolxok
C
C**************************************************************************
C     ****v* src/tolxok
C  NAME
C     tolxok
C     
C  DESCRIPTION
C     logical  tolxok    .. true  => the update dX is sufficiently small
C                                    in Newton's iteration
C                           false => no convergence (yet) due to a large
C                                    Newton update dX
C  SEE ALSO
C     TOLX, TOLF, tolfok
C
C**************************************************************************
C     ****v* src/retry
C  NAME
C     retry
C     
C  DESCRIPTION
C     integer  retry     .. the number of times that the current step is 
C                           being re-tried (due to convergence failures or
C                           refused timepoints)
C                           Usualy retry should be 0.
C  SEE ALSO 
C     NRETRY, iopt
C
C**************************************************************************
C     ****v* src/same_p
C  NAME
C     same_p
C     
C  DESCRIPTION
C     integer  same_p    .. an order increase is possible only after 
C                           completing KEEP_P steps with the current
C                           order
C                           -> same_p keeps track of the number of steps
C                           that have already been completed with the 
C                           current order
C  SEE ALSO
C     KEEP_P, iopt     
C
C**************************************************************************
C     ****v* src/conv
C  NAME
C     conv
C     
C  DESCRIPTION
C     integer  conv      .. conv = 0 signals 'no convergence'
C                           conv > 0 signals 'convergence' after conv 
C                                    iterations      
C  SEE ALSO
C     NEWTIT, iopt
C
C**************************************************************************
C     ****v* src/ipiv
C  NAME
C     ipiv
C     
C  DESCRIPTION
C     integer  ipiv(m)   .. array used by DGETRF, DGETRS for solving linear
C                           systems (pivot indices)       
C
C**************************************************************************
C     ****v* src/X_use
C  NAME
C     X_use
C     
C  DESCRIPTION
C     integer  X_use     .. X_use past timepoints are available for 
C                           predicting stage values 
C
C  SEE ALSO
C     polint, prdicX, updxvl
C
C**************************************************************************
C     ****v* src/X_perm
C  NAME
C     X_perm
C     
C  DESCRIPTION
C     integer  X_perm(BAKIFO) .. permutation in order to keep track of
C                                the past timepoints using the array 
C                                X_vals 
C  SEE ALSO
C     polint, prdicX, updxvl
C
C**************************************************************************
C     ****v* src/h
C  NAME
C     h
C     
C  DESCRIPTION
C     real*8   h         .. current stepsize
C
C**************************************************************************
C     ****v* src/t
C  NAME
C     t
C     
C  DESCRIPTION
C     real*8   t         .. current timepoint, the step is taken from 
C                           t to t+h
C
C**************************************************************************
C     ****v* src/hnew
C  NAME
C     hnew
C     
C  DESCRIPTION
C     real*8   hnew      .. predictions for the new stepsize (next step) 
C
C**************************************************************************
C     ****v* src/mat_A
C  NAME
C     mat_A
C     
C  DESCRIPTION
C     real*8   mat_A(m,n) .. the matrix df/dy for computing the Jacobian
C
C**************************************************************************
C     ****v* src/mat_D
C  NAME
C     mat_D
C     
C  DESCRIPTION
C     real*8   mat_D(n,m) .. the matrix dq/dx for computing the Jacobian 
C
C**************************************************************************
C     ****v* src/Jac
C  NAME
C     Jac
C     
C  DESCRIPTION
C     real*8   Jac(m,m)  .. the Jacobian 
C                           Jac = AD+h*lam B = df/dy dq/dx + h*lam df/dx
C
C**************************************************************************
C     ****v* src/Q
C  NAME
C     Q
C     
C  DESCRIPTION
C     real*8   Q(n)      .. charges/fluxes Q=q(x,t) at the current timepoint
C
C  SEE ALSO
C     omega
C
C**************************************************************************
C     ****v* src/Qp
C  NAME
C     Qp
C     
C  DESCRIPTION
C     real*8   Qp(n,MD)  .. derivatives \dot Q=q'(x,t) of charges/fluxes at
C                           the current timepoint
C  SEE ALSO
C     omega
C
C**************************************************************************
C     ****v* src/qn_in
C  NAME
C     qn_in
C     
C  DESCRIPTION
C     real*8   qn_in(n,MD) .. the incoming Nordsieck vector
C                             qn_in(:,k+1) contains an approximation to the
C                             scaled derivative h^k d^k/dt^k ( q(x(t),t) ),
C                             k=0,1,..,p-1
C
C**************************************************************************
C     ****v* src/qn_out
C  NAME
C     qn_out
C     
C  DESCRIPTION
C     real*8   qn_out(n,MD).. the outgoing Nordsieck vector
C                             qn_out(:,k+1) contains an approximation to the
C                             scaled derivative h^k d^k/dt^k(q(x(t+h),t+h)),
C                             k=0,1,..,p-1
C
C**************************************************************************
C     ****v* src/fzxt
C  NAME
C     fzxt
C     
C  DESCRIPTION
C     real*8   fzxt(m)   .. the function value fzxt = f(z,x,t) 
C
C**************************************************************************
C     ****v* src/res
C  NAME
C     res
C     
C  DESCRIPTION
C     real*8   res(m)    .. the residuum -h*lam fzxt and (later) the 
C                           Newton update
C
C**************************************************************************
C     ****v* src/X
C  NAME
C     X
C     
C  DESCRIPTION
C     real*8   X(m,MD)   .. the stage values of the current step 
C
C**************************************************************************
C     ****v* src/scdrv0
C  NAME
C     scdrv0
C     
C  DESCRIPTION
C     real*8   scdrv0(n) .. an estimate for the scaled derivative   
C                           h^p d^p/dt^p ( q(x(t+h),t+h) )
C                           of order p at the end of the current step
C
C**************************************************************************
C     ****v* src/scdrv1
C  NAME
C     scdrv1
C     
C  DESCRIPTION
C     real*8   scdrv1(n) .. an estimate for the scaled derivative   
C                           h^(p+1) d^(p+1)/dt^(p+1) ( q(x(t+h),t+h) )
C                           of order p+1 at the end of the current step
C
C**************************************************************************
C     ****v* src/scalQ
C  NAME
C     scalQ
C     
C  DESCRIPTION
C     real*8   scalQ(n)  .. mixed tolerance for charges and fluxes
C                           scalQ(k) = atolQ(k) + rtolQ(k)*abs(Q(k))
C
C**************************************************************************
C     ****v* src/lam
C  NAME
C     lam
C     
C  DESCRIPTION
C     real*8   lam       .. diagonal element A(i,i) of the method currently
C                           used
C
C**************************************************************************
C     ****v* src/scalX
C  NAME
C     scalX
C     
C  DESCRIPTION
C     real*8   scalX(m)  .. mixed tolerance for the dependent variable
C                           scalX(k) = atolX(k) + rtolX(k)*abs(X(k,p))
C
C**************************************************************************
C     ****v* src/omega
C  NAME
C     omega
C     
C  DESCRIPTION
C     real*8   omega(n)  .. this part of the residual res is independent 
C                           of the current stage and is thus only computed
C                           once within the Newton iteration
C  NOTE
C
C     For f(q'(x,t), x, t) =0 define
C
C       Q(X) = h A Qp + U qn_in,    qn_out = h B Qp + V qn_in    
C 
C     with
C
C       f( Qp, X, t) = 0.
C
C     More precisely:
C   
C       Q(X_i) = h lam Qp_i + h \sum_{j=1}^{i-1} A(i,j) Qp_j + U qn_in 
C              = h lam Qp_i + omega_i
C     
C     and thus
C
C       F(X_i) = h * lam * f( (Q(X_i)-omega_i)/h/lam, X_i, t+c_i h) = 0
C
C     needs to be solved using Newton's method. The vector omega_i is saved 
C     in omega as it does not depend on X_i.
C 
C**************************************************************************
C     ****v* src/nresx
C  NAME
C     nresx
C     
C  DESCRIPTION
C     real*8   nresx     .. nresx = dnrm2(m,res,1) is the norm of the 
C                           residual for computing the convergence rate of
C                           Newton's method   
C  SEE ALSO
C     convrt
C
C**************************************************************************
C     ****v* src/nresx1
C  NAME
C     nresx1
C     
C  DESCRIPTION
C     real*8   nresx1    .. previous value of nresx
C
C  SEE ALSO
C     convrt
C
C**************************************************************************
C     ****v* src/psi
C  NAME
C     psi
C     
C  DESCRIPTION
C     real*8   psi       .. psi = max(psi,dsqrt( theta * theta1 )) is a 
C                           measure for the convergence rate of Newton's 
C                           method
C  SEE ALSO
C     convrt
C
C**************************************************************************
C     ****v* src/psi1
C  NAME
C     psi1
C     
C  DESCRIPTION
C     real*8   psi1      .. previous estimate for the convergence rate 
C                           Newton's method
C  SEE ALSO
C     convrt
C
C**************************************************************************
C     ****v* src/theta
C  NAME
C     theta
C     
C  DESCRIPTION
C     real*8   theta     .. theta  = nresx / nresx1 is used for computing
C                           the convergence rate of Newton's method 
C  SEE ALSO
C     convrt
C
C**************************************************************************
C     ****v* src/theta1
C  NAME
C     theta1
C     
C  DESCRIPTION
C     real*8   theta1    .. previous value of theta  
C
C  SEE ALSO
C     convrt
C
C**************************************************************************
C     ****v* src/err
C  NAME
C     err
C     
C  DESCRIPTION
C     real*8   err       .. measure of the local truncation error for q(x,t)
C                           in the current step (scaled by scalQ)
C
C**************************************************************************
C     ****v* src/resold
C  NAME
C     resold
C     
C  DESCRIPTION
C     real*8   resold    .. norm of the previous (old) residual
C 
C  SEE ALSO
C     nrmres
C
C**************************************************************************
C     ****v* src/nrmres
C  NAME
C     nrmres
C     
C  DESCRIPTION
C     real*8   nrmres    .. norm of the residual 
C
C**************************************************************************
C     ****v* src/errvec
C  NAME
C     errvec
C     
C  DESCRIPTION
C     real*8   errvec(n) .. vector of local truncation errors for the 
C                           different components of q(x,t) (scaled by scalQ)
C
C**************************************************************************
C     ****v* src/X_vals
C  NAME
C     X_vals
C     
C  DESCRIPTION
C     real*8   X_vals(m,BAKIFO) .. backward inforation used for predicting
C                                  stage values
C  SEE ALSO
C     polint, prdicX, updxvl
C
C**************************************************************************
C     ****v* src/X_time
C  NAME
C     X_time
C     
C  DESCRIPTION
C     real*8   X_time(BAKIFO) .. timepoints correxponding to X_vals
C
C  SEE ALSO
C     polint, prdicX, updxvl
C
C**************************************************************************
C     ****v* src/A
C  NAME
C     A
C     
C  DESCRIPTION
C     real*8   A(MD,MD)  .. matrix A of the general linear method
C
C  SEE ALSO
C     omega
C
C**************************************************************************
C     ****v* src/U
C  NAME
C     U
C     
C  DESCRIPTION
C     real*8   U(MD,MD)  .. matrix U of the general linear method    
C
C  SEE ALSO
C     omega
C
C**************************************************************************
C     ****v* src/B
C  NAME
C     B
C     
C  DESCRIPTION
C     real*8   B(MD,MD)  .. matrix B of the general linear method
C
C  SEE ALSO
C     omega
C
C**************************************************************************
C     ****v* src/V
C  NAME
C     V
C     
C  DESCRIPTION
C     real*8   V(MD,MD)  .. matrix V of the general linear method
C
C  SEE ALSO
C     omega
C
C**************************************************************************
C     ****v* src/c
C  NAME
C     c
C     
C  DESCRIPTION
C     real*8   c(MD)     .. vector c of the general linear method 
C                           determining intermediate timepoints
C  SEE ALSO
C     omega
C
C**************************************************************************
C     ****v* src/c1beta
C  NAME
C     c1beta             
C
C  DESCRIPTION
C     real*8   c1beta(MD) .. c1beta(1) is the methods error coefficients
C                            c1beta(2:p) are are higher order error 
C                            coefficients that need to be kept constant by 
C                            the scale and modify technique
C  SEE ALSO
C     sclmod, getcon
C
C**************************************************************************
C     ****v* src/delta
C  NAME
C     delta
C     
C  DESCRIPTION
C     real*8   delta(MD+1,2) .. coefficients used for estimating the p-th
C                               and (p+1)-st scaled derivatives
C                               delta depends on the method currently used
C SEE ALSO
C     getest, scdrvs
C
C**************************************************************************
C     ****v* src/sd0_ac
C  NAME
C     sd0_ac
C     
C  DESCRIPTION
C     real*8   sd0_ac(n) .. last accepted value of scdrv0
C
C**************************************************************************
C     ****v* src/scaacQ
C  NAME
C     scaacQ
C     
C  DESCRIPTION
C     real*8   scaacQ(n) .. last accepted value of scalQ
C
C**************************************************************************
C     ****v* src/qn_acc
C  NAME
C     qn_acc
C     
C  DESCRIPTION
C     real*8   qn_acc(n,MD) .. last accepted value of qn_out
C
C**************************************************************************
C     ****v* src/h_acc
C  NAME
C     h_acc
C     
C  DESCRIPTION
C     real*8   h_acc     .. stepsize used when computing the last accepted
C                           step
C
C**************************************************************************
C     ****v* src/hnw_ac
C  NAME
C     hnw_ac
C     
C  DESCRIPTION
C     real*8   hnw_ac    .. stepsize that was proposed for the next step 
C                           during computing the last accepted step
C
C**************************************************************************
C     ****v* src/sd1_ac
C  NAME
C     sd1_ac
C     
C  DESCRIPTION
C     real*8   sd1_ac(n) .. last accepted value of scdrv1
C
C**************************************************************************
C     ****v* src/scaacX
C  NAME
C     scaacX
C     
C  DESCRIPTION
C     real*8   scaacX(m) .. last accepted value of scalX
C
C**************************************************************************
C     ****v* src/qp_acc
C  NAME
C     qp_acc
C     
C  DESCRIPTION
C     real*8   qp_acc(n) .. last accepted value of qprime  
C
C**************************************************************************
C     ****v* src/p
C  NAME
C     p
C     
C  DESCRIPTION
C     integer  p         .. the integration order currently used
C
C**************************************************************************
C     ****v* src/p_acc
C  NAME
C     p_acc
C     
C  DESCRIPTION
C     integer  p_acc     .. the integration order used for the last 
C                           accepted step
C
C**************************************************************************

C**************************************************************************
C
C     GLIMDA version 0.1
C     General LInear Methods for Differential Agebraic equations
C
C     Copyright (C) 2006,  Steffen Voigtmann
C     
C     This program is free software; you can redistribute it and/or
C     modify it under the terms of the GNU Lesser General Public
C     License as published by the Free Software Foundation; either
C     version 2.1 of the License, or (at your option) any later version.
C     
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C     Lesser General Public License for more details.
C
C     You should have received a copy of the GNU Lesser General Public
C     License along with this library; if not, write to the Free Software
C     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
C
C**************************************************************************

C**************************************************************************
C     ****f* src/getmet
C
C  NAME
C     getmet -- GET a METhod
C
C  DESCRIPTION
C     This subroutine provide all data required for a methd of order p.
C
C  ARGUMENTS
C     integer  p            .. the integration order currently used
C     integer  MD           .. dimension of the matrices A,U,B,V
C     real*8   A(MD,MD)     .. matrix A of the general linear method
C     real*8   U(MD,MD)     .. matrix U of the general linear method
C     real*8   B(MD,MD)     .. matrix B of the general linear method
C     real*8   V(MD,MD)     .. matrix V of the general linear method
C     real*8   c(MD)        .. vector c of the general linear method 
C     real*8   c1beta(MD)   .. c1beta(1) is the methods error coefficients
C                              c1beta(2:p) are are higher order error 
C                              coefficients that need to be kept constant by 
C                              the scale and modify technique
C     real*8   delta(MD+1,2).. coefficients used for estimating the p-th
C                              and (p+1)-st scaled derivatives
C                              delta depends on the method currently used
C     integer  type         .. the particular order 3 method to use
C                              type = 1  use a method with s=r=p=q=3 and 
C                                        A(alpha) stability for alpha=74°
C                                        but being unconditionally stable
C                                        at zero and at infinity
C                              type = 2  use a method with s=r=p=q=3 and 
C                                        A(alpha) stability for alpha=88.1°
C     integer  ierr         .. return code (0 signals success)
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine getmet(p,MD,A,U,B,V,c,c1beta,delta,type,ierr)
      implicit none
      integer          p,MD,ierr,type
      real*8           A(MD,MD), U(MD,MD), B(MD,MD), V(MD,MD), c(MD),
     $                 c1beta(MD), delta(MD+1,2)

      if ( p .eq. 1 ) then      

         A(1,1)  = 1d0
         U(1,1)  = 1d0
         B(1,1)  = 1d0
         V(1,1)  = 1d0
         c(1)    = 1d0
         c1beta(1)  = -0.5d0
         call getest(MD,c,p,delta)

      elseif ( p .eq. 2 ) then

         call trbdf2(p,MD,A,U,B,V,c,c1beta,delta,ierr)

      elseif ( p .eq. 3 ) then

         if     (type .eq. 1) then
            call p3s3  (p,MD,A,U,B,V,c,c1beta,delta,ierr)
         elseif (type .eq. 2) then
            call p3s3_2(p,MD,A,U,B,V,c,c1beta,delta,ierr)
         else
            ierr = -1
            return
         end if

      else
         ierr = -1
         return
      endif

      ierr = 0

      end
C**************************************************************************
C     ****f* src/getcon
C
C  NAME
C     getcon -- GET error CONstants
C
C  DESCRIPTION
C     This subroutine computes the error constant c1beta(1) and the 
C     higher order error coefficients c1beta(2:p)
C 
C  COPYRIGHT
C
C  SOURCE
C
      subroutine getcon(MD,B,V,c,p,c1beta)
      implicit none
      integer          p, MD, i, ipiv(p), info
      real*8           B(MD,MD), V(MD,MD), c(MD), c1beta(MD),
     $                 M(MD,p), cp(p), facul
C
C     we need to solve the linear system
C         [ e1*e1^T + (I-V) ] c1beta = a_{p,p} - B c^p / p!
C      
C     --- get coefficient matrix ---
C     M = V
      call dcopy(MD*p,V,1,M,1)      
C     M = V-I
      call daxpy(p,-1d0,1d0,0,M,MD+1)
C     M = V - I - e1*e1^T
      M(1,1) = -1d0
C
C     --- get right hand side ---
C     c1beta = a_{p,p}
      do i=1,p
         cp(i)     = c(i)**p
         c1beta(i) = 1d0/facul(p+2-i) 
      end do
C
C     rhs = B c^p / p! - a_{p,p}
      call dgemv ('n',p,p,1d0/facul(p),B,MD,cp,1,
     $     -1d0,c1beta,1)
C
C     --- solve the linear system ---
      call dgetrf(p, p, M, MD, ipiv, info )
      call dgetrs('n', p, 1, M, MD, ipiv, c1beta, MD, info)

      end
C**************************************************************************
C     ****f* src/facul
C
C  NAME
C     facul -- compute n!
C 
C  COPYRIGHT
C
C  SOURCE
C
      real*8           function facul(n)
      implicit none
      integer i,n
      facul = 1d0
      do i=1,n
         facul = facul*real(i)
      end do
      end
C**************************************************************************
C     ****f* src/getest
C
C  NAME
C     getest -- GET error ESTimator
C
C  DESCRIPTION
C     This subroutine computes the coefficients used when approximating
C     scaled derivatives using linear combinations of stage derivatives. 
C 
C  COPYRIGHT
C
C  SOURCE
C
      subroutine getest(MD,c,p,delta)
      implicit none
      integer          p, MD, i, j, ipiv(p+1), info
      real*8           c(MD), delta(MD+1,2), W(p+1,p+1), facul
C
C     we need to solve the linear system
C                                 [ 0 0 ]
C          [ 1  1 ... 1 ]         [ 0 0 ]
C          [ 0          ] delta = [ : : ]
C          [ 0     W    ]         [ 0 0 ]
C                                 [ 1 0 ]
C                                 [ 1 1 ] 
C     where W = [ c , c^2/2! , ... , c^p/p! ]^T.
C
C     compute right hand side
C
      call dcopy(p,0d0,0,delta(1,1),1)
      call dcopy(p,0d0,0,delta(1,2),1)
      delta(p,1)   = 1d0
      delta(p+1,1) = 1d0
      delta(p+1,2) = 1d0
C
C     compute matrix W
C
      call dcopy(p+1,1d0,0,W(1,1),p+1)
      call dcopy(p  ,0d0,0,W(2,1),1  )
      do i=1,p
         do j=1,p
            W(j+1,i+1) = c(i)**j / facul(j)
         end do
      end do
C
C     solve the linear system
C
      call dgetrf(p+1, p+1, W, p+1, ipiv, info )
      call dgetrs('n', p+1, 2, W, p+1, ipiv, delta, MD+1, info)

      end
C**************************************************************************
C     ****f* src/trbdf2
C
C  NAME
C     trbdf2
C
C  DESCRIPTION
C     This is a stiffly accurate, A- and L-stable method with s=r=p=q=2
C     that can be interpreted as a concatenation of a trapezoidal step 
C     and a step with the BDF2 method.
C     The diagonal element lambda controls the damping behaviour.
C     For lambda->0 the trapezoidal rule is approximated.
C     For 1d0-dsqrt(2d0)/2d0 the two diagonal elements are equal. 
C 
C  COPYRIGHT
C
C  SOURCE
C
      subroutine trbdf2(p,MD,A,U,B,V,c,c1beta,delta,ierr)
      implicit none
      integer          p, MD, ierr
      real*8           lambda,A(MD,MD),U(MD,MD),B(MD,MD),V(MD,MD),c(MD),
     $                 c1beta(MD),delta(MD+1,2)
C
C     parameters
C
      lambda = 1d0-dsqrt(2d0)/2d0
C
C     number of internal and external stages is s=r=p
C
C     A \in \R^{2 x 2}
C
      A(1,1) = (2*lambda-1)/(lambda-1)/2
      A(1,2) = 0
      A(2,1) = 1.D0/2.D0-lambda/2
      A(2,2) = lambda
C
C     U \in \R^{2 x 2}
C
      U(1,1) = 1
      U(1,2) = (2*lambda-1)/(lambda-1)/2
      U(2,1) = 1
      U(2,2) = 1.D0/2.D0-lambda/2
C
C     B \in \R^{2 x 2}
C
      B(1,1) = 1.D0/2.D0-lambda/2
      B(1,2) = lambda
      B(2,1) = 0
      B(2,2) = 1
C
C     V \in \R^{2 x 2}
C
      V(1,1) = 1
      V(1,2) = 1.D0/2.D0-lambda/2
      V(2,1) = 0
      V(2,2) = 0
C
C     c \in \R^2
C
      c(1) = (2*lambda-1)/(lambda-1)
      c(2) = 1
C
C     compute remaining parts of the method 
C
      call getcon(MD,B,V,c,p,c1beta)
      call getest(MD,c,p,delta)

      end
C**************************************************************************
C     ****f* src/p3s3
C
C  NAME
C     p3s3
C
C  DESCRIPTION
C     This is a method with s=r=p=q=3 being stiffly accurate and 
C     unconditionally stable at zero and at infinity. The method is not
C     A-stable but A(alpha) stable and the only eigenvalue of M(\infty)
C     is zero.
C 
C  COPYRIGHT
C
C  SOURCE
C
      subroutine p3s3(p,MD,A,U,B,V,c,c1beta,delta,ierr)
      implicit none
      integer          p, MD, ierr
      real*8           lambda,A(MD,MD),U(MD,MD),B(MD,MD),V(MD,MD),c(MD),
     $                 c1beta(MD),delta(MD+1,2)

C     .. alpha = 61.7° ..
      lambda = 0.15d0
C     .. alpha = 74° ..     
      lambda = 4.d0/25.d0

      A(1,1) = lambda
      A(1,2) = 0
      A(1,3) = 0
      A(2,1) = -(6*lambda**2+1-6*lambda)**2*(-21*lambda**2-1+9*lambda+9*
     $lambda**3)/(-5*lambda+3*lambda**2+1)**3/lambda**2/27
      A(2,2) = lambda
      A(2,3) = 0
      A(3,1) = (-81*lambda**2-1+15*lambda+189*lambda**3-171*lambda**4+54
     $*lambda**5)/lambda**2/(-21*lambda**2-1+9*lambda+9*lambda**3)/27
      A(3,2) = -(-5*lambda+3*lambda**2+1)**3*lambda/(-21*lambda**2-1+9*l
     $ambda+9*lambda**3)/(6*lambda**2+1-6*lambda)
      A(3,3) = lambda
      U(1,1) = 1
      U(1,2) = 2*lambda
      U(1,3) = 3.D0/2.D0*lambda**2
      U(2,1) = 1
      U(2,2) = -(1+12798*lambda**7-306*lambda**3-21*lambda-5103*lambda**
     $8+7452*lambda**5+150*lambda**2+729*lambda**9-1224*lambda**4-14526*
     $lambda**6)/(-5*lambda+3*lambda**2+1)**3/lambda**2/27
      U(2,3) = -(6*lambda**2+1-6*lambda)*(1350*lambda**4+279*lambda**2-9
     $09*lambda**3-810*lambda**5-39*lambda+2+162*lambda**6)/(-5*lambda+3
     $*lambda**2+1)**3/lambda/18
      U(3,1) = 1
      U(3,2) = (81*lambda**6-378*lambda**5+459*lambda**4-144*lambda**3-2
     $1*lambda**2+12*lambda-1)/(6*lambda**2+1-6*lambda)/lambda**2/27
      U(3,3) = (18*lambda**3-48*lambda**2+21*lambda-2)/lambda/18
      B(1,1) = (-81*lambda**2-1+15*lambda+189*lambda**3-171*lambda**4+54
     $*lambda**5)/lambda**2/(-21*lambda**2-1+9*lambda+9*lambda**3)/27
      B(1,2) = -(-5*lambda+3*lambda**2+1)**3*lambda/(-21*lambda**2-1+9*l
     $ambda+9*lambda**3)/(6*lambda**2+1-6*lambda)
      B(1,3) = lambda
      B(2,1) = 0
      B(2,2) = 0
      B(2,3) = 1
      B(3,1) = (1-366*lambda**3-18*lambda-288*lambda**5+120*lambda**2+51
     $0*lambda**4+54*lambda**6)/lambda/(6*lambda-1)/(-1+3*lambda)/(3*lam
     $bda**2-6*lambda+1)/(3*lambda**3-9*lambda**2+6*lambda-1)/3
      B(3,2) = 3*(-5*lambda+3*lambda**2+1)**3*lambda**2/(6*lambda**2+1-6
     $*lambda)/(6*lambda-1)/(-1+3*lambda)/(3*lambda**2-6*lambda+1)/(3*la
     $mbda**3-9*lambda**2+6*lambda-1)
      B(3,3) = (18*lambda**4-66*lambda**3+66*lambda**2-21*lambda+2)/(3*l
     $ambda**3-9*lambda**2+6*lambda-1)/(6*lambda-1)
      V(1,1) = 1
      V(1,2) = (81*lambda**6-378*lambda**5+459*lambda**4-144*lambda**3-2
     $1*lambda**2+12*lambda-1)/(6*lambda**2+1-6*lambda)/lambda**2/27
      V(1,3) = (18*lambda**3-48*lambda**2+21*lambda-2)/lambda/18
      V(2,1) = 0
      V(2,2) = 0
      V(2,3) = 0
      V(3,1) = 0
      V(3,2) = -(324*lambda**7-1485*lambda**6+2394*lambda**5-1851*lambda
     $**4+771*lambda**3-177*lambda**2+21*lambda-1)/(6*lambda**2+1-6*lamb
     $da)/(3*lambda**3-9*lambda**2+6*lambda-1)/(6*lambda-1)/lambda/3
      V(3,3) = 0
      c(1) = 3*lambda
      c(2) = (6*lambda**2+1-6*lambda)/(-5*lambda+3*lambda**2+1)
      c(3) = 1

      call getcon(MD,B,V,c,p,c1beta)
      call getest(MD,c,p,delta)

      end
C**************************************************************************
C     ****f* src/p3s3_2
C
C  NAME
C     p3s3_2
C
C  DESCRIPTION
C     This is another stiffly accurate method with s=r=p=q=3. Again, the
C     method is not A-stable but A(alpha) stable and the only eigenvalue
C     of M(\infty) is zero. Compared to p3s3, the method is not 
C     unconditionally stable but the region of A(alpha) stability is 
C     considerably enlarged to alpha=88.1°.
C 
C  COPYRIGHT
C
C  SOURCE
C
      subroutine p3s3_2(p,MD,A,U,B,V,c,c1beta,delta,ierr)
      implicit none
      integer          p, MD, ierr
      real*8           lambda, c2, 
     $                 A(MD,MD), U(MD,MD),
     $                 B(MD,MD), V(MD,MD),
     $                 c(MD),    c1beta(MD),
     $                 delta(MD+1,2)
      real*8           s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,
     $     s11,s12,s13,s14,s15,s16,s17,s18,s19,s20

      A(1,1) = 7.D0/40.D0
      A(1,2) = 0
      A(1,3) = 0
      A(2,1) = 343.D0/1458.D0
      A(2,2) = 7.D0/40.D0
      A(2,3) = 0
      A(3,1) = 5225432851.D0/18627840000.D0+sqrt(11831751737089.D0)/1267
     $20000
      A(3,2) = 2346585921.D0/19317760000.D0-9.D0/2759680000.D0*sqrt(1183
     $1751737089.D0)
      A(3,3) = 7.D0/40.D0
      U(1,1) = 1
      U(1,2) = 7.D0/20.D0
      U(1,3) = 147.D0/3200.D0
      U(2,1) = 1
      U(2,2) = 11851.D0/29160.D0
      U(2,3) = 6517.D0/97200.D0
      U(3,1) = 1
      U(3,2) = 44126632861.D0/104315904000.D0-23.D0/4967424000.D0*sqrt(1
     $1831751737089.D0)
      U(3,3) = 159211901.D0/2027520000.D0-sqrt(11831751737089.D0)/675840
     $000
      B(1,1) = 5225432851.D0/18627840000.D0+sqrt(11831751737089.D0)/1267
     $20000
      B(1,2) = 2346585921.D0/19317760000.D0-9.D0/2759680000.D0*sqrt(1183
     $1751737089.D0)
      B(1,3) = 7.D0/40.D0
      B(2,1) = 0
      B(2,2) = 0
      B(2,3) = 1
      B(3,1) = 3022207.D0/670320.D0+sqrt(11831751737089.D0)/670320
      B(3,2) = -3.D0/1207360.D0*sqrt(11831751737089.D0)-17130621.D0/1207
     $360.D0
      B(3,3) = 8202367.D0/802560.D0+sqrt(11831751737089.D0)/802560
      V(1,1) = 1
      V(1,2) = 44126632861.D0/104315904000.D0-23.D0/4967424000.D0*sqrt(1
     $1831751737089.D0)
      V(1,3) = 159211901.D0/2027520000.D0-sqrt(11831751737089.D0)/675840
     $000
      V(2,1) = 0
      V(2,2) = 0
      V(2,3) = 0
      V(3,1) = 0
      V(3,2) = -2135167.D0/3951360.D0-sqrt(11831751737089.D0)/3951360
      V(3,3) = 0
      c(1) = 21.D0/40.D0
      c(2) = 49.D0/60.D0
      c(3) = 1

      call getcon(MD,B,V,c,p,c1beta)
      call getest(MD,c,p,delta)
      end
C**************************************************************************
C**************************************************************************
C
C     GLIMDA version 0.1
C     General LInear Methods for Differential Agebraic equations
C
C     Copyright (C) 2006,  Steffen Voigtmann
C     
C     This program is free software; you can redistribute it and/or
C     modify it under the terms of the GNU Lesser General Public
C     License as published by the Free Software Foundation; either
C     version 2.1 of the License, or (at your option) any later version.
C     
C     This program is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C     Lesser General Public License for more details.
C
C     You should have received a copy of the GNU Lesser General Public
C     License along with this library; if not, write to the Free Software
C     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
C
C**************************************************************************

C**************************************************************************
C     ****f* src/tolX2Q
C
C  NAME
C     tolX2Q -- convert TOLerances for X into those for Q
C
C  DESCRIPTION
C     This subroutine computes the absolute and relative tolerances
C     for the vector q(x,t) from those given for x
C     -> computes atolQ and rtolQ(m)
C 
C  COPYRIGHT
C
C  SOURCE
C
      subroutine tolX2Q(m,n,atolX,rtolX,t,atolQ,rtolQ,num_D,ADcnst,
     $                  mat_D,dqxevl,qevl,Q0,rpar,ipar,PRTLEV,eps,ierr)
      implicit none
      logical          num_D, ADcnst
      integer          m, n, ipar(*),PRTLEV,ierr
      real*8           atolX(m),rtolX(m),t,atolQ(n),rtolQ(n),mat_D(n,m),
     $                 rpar(*), eps,Q0(n)
      external         qevl, dqxevl

      logical          naninf
      integer          i,j
      real*8           absD(n,m),Q1(n)

      if (.not.ADcnst) then
         call get_D(m,n,atolX,t,num_D,dqxevl,qevl,Q0,mat_D,eps,
     $              ipar,rpar,PRTLEV,ierr)
         if (ierr.ne.0) return
      endif

      do i=1,n
         do j=1,m
            absD(i,j) = abs(mat_D(i,j))
         enddo
      enddo

      call dgemv ('n',n,m,1d0,absD,n,atolX,1,0d0,atolQ,1)

      if (ADcnst) then

         call dgemv ('n',n,m,1d0,absD,n,rtolX,1,0d0,rtolQ,1)

      else

         call get_D(m,n,rtolX,t,num_D,dqxevl,qevl,Q0,mat_D,eps,
     $              ipar,rpar,PRTLEV,ierr)

         if (ierr.ne.0) return

         do i=1,n
            do j=1,m
               absD(i,j) = abs(mat_D(i,j))
            enddo
         enddo

         call dgemv ('n',n,m,1d0,absD,n,rtolX,1,0d0,rtolQ,1)

      endif

      end
C**************************************************************************
C     ****f* src/get_D
C
C  NAME
C     get_D -- compute the matrix mat_D = D = dq/dx 
C              (either analytically or numerically)
C
C  DESCRIPTION
C     In case of (num_D.eq..true.) finite differences are used
C     to compute D, otherwise D is evaluated analytically using
C     the subroutine dqxevl.
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine get_D(m,n,x,t,num_D,dqxevl,qevl,Q0,mat_D,eps,
     $                 ipar,rpar,PRTLEV,ierr)
      implicit none
      logical          num_D
      integer          m, n, ipar(*), PRTLEV,ierr
      real*8           x(m), t, mat_D(n,m), Q0(n), eps, rpar(*)

      logical          naninf
      integer          i, j
      real*8           safe, delta, Q1(n)

      if( num_D ) then

         do i=1,m
            safe = x(i)
            delta = dsqrt( eps * max( 1.d-5, abs(safe) ) )
            delta = (x(i)+delta)-x(i)
            x(i)  = safe + delta
            call qevl(m,n,x,t,Q1,rpar,ipar,ierr)
            if (ierr.ne.0) return
            do j=1,n
               mat_D(j,i) = ( Q1(j) - Q0(j) ) / delta
            enddo
            x(i) = safe
         enddo

      else

         call dqxevl(m,n,x,t,mat_D,rpar,ipar,ierr)
         if ( ierr.ne.0 ) return

      endif

      if (naninf(mat_D,n,m,n)) then
         if (PRTLEV.ge.2) print*,'GLIMDA WARNING: '
     $        //'evaluation of dqx yields NAN or INF'
         ierr = -10
         return
      endif

      end
C**************************************************************************
C     ****f* src/get_A
C
C  NAME
C     get_A -- compute the matrix mat_A = A = df/dy 
C              (either analytically or numerically)
C
C  DESCRIPTION
C     In case of (num_A.eq..true.) finite differences are used
C     to compute A, otherwise A is evaluated analytically using
C     the subroutine dfyevl.

C
C  DESCRIPTION
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine get_A(m,n,z,x,t,num_A,dfyevl,fevl,F0,mat_A,eps,
     $                 ipar,rpar,PRTLEV,ierr)
      implicit none
      logical          num_A
      integer          m, n, ipar(*), PRTLEV,ierr
      real*8           z(n), x(m), t, mat_A(m,n), F0(m), eps, rpar(*)

      logical          naninf
      integer          i, j
      real*8           safe, delta, F1(m)

      if( num_A ) then

         do i=1,n
            safe = z(i)
            delta = dsqrt( eps * max( 1.d-5, abs(safe) ) )
            z(i)  = safe + delta
            call fevl(m,n,z,x,t,F1,rpar,ipar,ierr)
            if (ierr.ne.0) return
            do j=1,m
               mat_A(j,i) = ( F1(j) - F0(j) ) / delta
            enddo
            z(i) = safe
         enddo

      else

         call dfyevl(m,n,z,x,t,mat_A,rpar,ipar,ierr)
         if ( ierr.ne.0 ) return

      endif

      if (naninf(mat_A,m,n,m)) then
         if (PRTLEV.ge.2) print*,'GLIMDA WARNING: '
     $        //'evaluation of dfy yields NAN or INF'
         ierr = -10
         return
      endif

      end
C**************************************************************************
C     ****f* src/naninf
C
C  NAME
C     naninf -- check for NANs and INFs in a matrix
C
C  DESCRIPTION
C     This function attempts to determine whether there are any
C     NANs or INFs contained in the matrix Mat.
C     There must be better ways to do this, but this works good enough. 
C
C  ARGUMENTS
C     real*8  Mat(ld,n) .. matrix of dimension (m,n)
C     integer m         .. number of rows
C     integer n         .. number if column
C     integer ld        .. leading dimension of Mat
C
C  RESULT
C     .true.  -> NANs or INFs were found
C     .false. -> otherwise
C
C  COPYRIGHT
C
C  SOURCE
C
      logical function naninf(Mat,m,n,ld)
      implicit none
      integer          m, n, ld
      real*8           Mat(ld,n), dlange, work(4*m), nan
      nan = dlange('f',m,n,Mat,ld,work)
      if ( (nan.ne.nan) .or. (nan*1d1.eq.(nan-1d0)) ) then
         naninf = .true.
      else
         naninf = .false.
      end if      
      end
C**************************************************************************
C     ****f* src/updxvl
C
C  NAME
C     updxvl -- UPDate X VaLues used for stage prediction
C
C  DESCRIPTION
C     Extrapolation based on BAKIFO past timepoints is used to predict
C     stage values. This subroutine updates the backward information 
C     (stored in the matrix X_vals) after completing a step.
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine updxvl(m,p,t,h,c,BAKIFO,X_use,X_time,X_vals,X_perm,X)
      implicit none
      integer          m,p,BAKIFO,X_use,X_perm(BAKIFO) 
      real*8           t,h,c(p),X_time(BAKIFO),X_vals(m,BAKIFO),X(m,p)

      integer          newuse, first, last, k, pos

      newuse = min(X_use+p,BAKIFO)
      first  = mod(X_use+p-newuse+X_perm(1)-1,BAKIFO)+1
      last   = X_perm(X_use)
      X_use  = newuse

      do k = 1,max(p,BAKIFO)

         if (k.le.BAKIFO) X_perm(k) = mod(first+k-2,BAKIFO) + 1

         if (k.le.p) then
            pos = mod(last+k-1,BAKIFO)+1
            call dcopy(m,X(1,k),1,X_vals(1,pos),1)
            X_time(pos) = t+c(k)*h
         endif

      enddo

      end
C**************************************************************************
C     ****f* src/filout
C
C  NAME
C     filout -- write FILe OUTput
C
C  DESCRIPTION
C     After completing a step the current solution is stored in (t+h,Y).
C     This subroutine writes the values
C       * timepoint t+h
C       * stepsize  h
C       * order     p
C       * solution  Y(1:m)
C     into the file with ID f_sol.
C
C  ARGUMENTS
C
C     integer f_sol .. file ID
C     real*8  tph   .. timepoint t+h
C     real*8  h     .. stepsize
C     integer p     .. integration order
C     real*8  Y     .. solution Y(t+h)
C     integer m     .. dimension of Y
C
C  COPYRIGHT
C
C  SOURCE
C
C      subroutine filout(f_sol,tph,h,p,Y,m)
C      implicit none
C      integer          f_sol, p, m, i
C      real*8           tph, h, Y(m)
C 1    format(d14.8)
C 2    format(d14.8,'  ',$)
C      write(f_sol,2) tph
C      write(f_sol,2) h
C      write(f_sol,2) real(p)
C      do i=1,m-1
C         write(f_sol,2) Y(i)
C      end do
C      write(f_sol,1) Y(m)
C      end
C
C**************************************************************************
C     ****f* src/scdrvs
C
C  NAME
C     scdrvs -- compute SCaleD DeRiVativeS
C
C  DESCRIPTION
C     This subroutines computes approximations to the scaled derivatives 
C     h^k d^k/dt^k ( q( x(t),t ) ) of order k=p and k=p+1.
C     Due to the general linear method's high stage order linear combinations
C     of stage derivatives can be used.
C     Yp containes the stage derivatives of the current step, and 
C     yp_acc represents the derivative of the last accepted stage. This is 
C     similar (but a bit more general) to using methods having the FSAL
C     property (first same as last).
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine scdrvs(h,Yp,m,p,yn_in,delta,md,yp_acc,scdrv0,scdrv1)
      implicit none
      integer          m, p, md
      real*8           h, Yp(m,p), yn_in(m,p), delta(md+1,2), yp_acc(m),
     $                 scdrv0(m), scdrv1(m)
      if (p.eq.1) then
C        scdrv0 = h * Yp(:,s)
         call dcopy(m,Yp(1,p),1,scdrv0,1)
         call dscal(m,h,scdrv0,1)
C        scdrv1 = ( Yp(:,s) - acc_Yp(:,s) ) * h
         call dcopy(m,scdrv0,1,scdrv1,1)
         call daxpy(m,-h,yp_acc,1,scdrv1,1)
      else
C        scdrv0 = h * Yp * delta(2:end,1) + delta(1,1) * y_in(:,2)
         call dcopy(m,yn_in(1,2),1,scdrv0,1)
         call dgemv('n',m,p,h,Yp,m,delta(2,1),1,delta(1,1),scdrv0,1)
C        scdrv1 = h * Yp * delta(2:end,2) + delta(1,2) * y_in(:,2)
         call dcopy(m,yn_in(1,2),1,scdrv1,1)
         call dgemv('n',m,p,h,Yp,m,delta(2,2),1,delta(1,2),scdrv1,1)
      endif
      end

C**************************************************************************
C     ****f* src/prdicX
C
C  NAME
C     prdicX -- PReDICt a value for X
C
C  DESCRIPTION
C     This subroutine sets up the backward data and uses the routine
C     polint in order to predict a stage value at timepoint t+c(k+1)*h.
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine prdicX(m,k,t,h,c,X,X_vals,BAKIFO,
     $                  X_time,X_use,X_perm,X0)
      implicit none
      integer          m, k, BAKIFO, X_perm(BAKIFO), X_use
      real*8           t, h, c(k+1), X(m,k+1), X_vals(m,BAKIFO),
     $                 X_time(BAKIFO), X0(m)

      integer          i
      real*8           XX(m,X_use+k+1), tt(X_use+k+1), dX0(m)

C     compute XX = [ X_vals(in corrected order) X(:,1:k) ]
      do i=1,X_use
         call dcopy(m,X_vals(1,X_perm(i)),1,XX(1,i),1)
      enddo
      call dcopy(m*k,X,1,XX(1,X_use+1),1)
C     compute tt = [ X_time(in corrected order) t+c*h ]
      do i=1,X_use
         tt(i) = X_time(X_perm(i))
      enddo
      do i=1,k
         tt(X_use+i) = t+c(i)*h
      enddo

C     do the polynomial interpolation
      call polint(m,X_use+k,tt,XX,t+c(k+1)*h,X0,dX0)

      end

C**************************************************************************
C     ****f* src/polint
C
C  NAME
C     polint -- POLynomial INTerpolation
C
C  DESCRIPTION
C     This subroutine implements Neville's algorithm (adapted from 
C     'Numerical Recipes in Fortran 77') for polynomial interpolation.
C
C     Let P_k be the unique polynomial satisfying P_k(xa(i)) = ya(k,i)
C     for 1<=i<=s and 1<=k<=m. Then y(k) = P_k(x) is calculated.
C     dy(k) is the corresponding error estimate.
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine polint(m,s,xa,ya,x,y,dy)
      implicit none
      integer          m, s
      real*8           xa(s), ya(m,s), x, y(m), dy(m)

      integer          i, k, ns
      real*8           ho, hp, w(m), c(m,s), d(m,s) 

C     initialise the increments c and d
      call dcopy(m*s,ya,1,c,1)
      call dcopy(m*s,ya,1,d,1)

C     as we are interested in extrapolation, we start from the right
      call dcopy(m,ya(1,s),1,y,1)
      ns=s-1 
      do k=1,s-1
         do i=1,s-k
            ho=xa(i)-x 
            hp=xa(i+k)-x 
            call dcopy(m,c(1,i+1),1,w,1)
            call daxpy(m,-1d0,d(1,i),1,w,1)
            call dcopy(m,w,1,c(1,i),1)
            call dcopy(m,w,1,d(1,i),1)
            call dscal(m,ho/(ho-hp),c(1,i),1)
            call dscal(m,hp/(ho-hp),d(1,i),1)
         enddo 
         if (2*ns.lt.s-k)then
            call dcopy(m,c(1,ns+1),1,dy,1)
         else 
            call dcopy(m,d(1,ns)  ,1,dy,1)
            ns=ns-1
         endif
         call daxpy(m,1d0,dy,1,y,1)
      enddo 
      return
      end
C**************************************************************************
C     ****f* src/sclmod
C
C  NAME
C     sclmod -- SCaLe and MODify
C
C  DESCRIPTION
C     When using general linear methods in a variable stepsize fashon 
C     one has to make sure that the higher order error terms of the 
C     Nordsieck vector are maintained correctly.
C     When changing the stepsize, the Nordsieck vector has to be reSCALED
C     using powers of the stepsize ratio. The error terms are corrected 
C     by additing a suitable MODIFIcation.
C
C     On input yn_in, scdrv0 and scdrv1 are based on the stepsize h, i.e.
C       scdrv0       = h^p     y^(p)  (t) 
C       scdrv1       = h^(p+1) y^(p+1)(t)
C       yn_in(:,k+1) = h^k     y^(k)  (t) - beta_{k+1} h^(p+1) y^(p+1)(t)
C                    = h^k     y^(k)  (t) - beta_{k+1} scdrv1
C     On output the corresponding quantities for hnew are returned.
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine sclmod(yn_in,m,p,hnew,h,scdrv0,scdrv1,c1beta)
      implicit none
      integer          m, p, k
      real*8           yn_in(m,p), hnew, h, scdrv0(m), scdrv1(m),
     $                 c1beta(p) , sigma, alpha, modify(m,p)
C
C     (0) compute the stepsize ratio
C
      sigma = hnew / h
C
C     (1) perform 'scale' 
C         yn_in(:,k+1) = sigma^k * yn_in(:,k+1)
C                      = hnew^k y^(k)(t) - beta_{k+1} sigma^k scdrv1
      do k=1,p-1
         call dscal(m,sigma**k,yn_in(1,k+1),1)
      end do
C
C     (2) perform 'modify': since we want to calculate
C         yn_in(:,k+1) = hnew^k y^(k)(t) - beta_{k+1} hnew^(p+1) y^(p+1)(t)
C                      = hnew^k y^(k)(t) - beta_{k+1} sigma^(p+1) scdrv1
C     we need to modify the result computed in (1):
C      
C     (2.1) compute 'modify'
C           modify(:,k+1) = (sigma^k-sigma^(p+1)) * beta_{k+1} * scdrv1
C           (k=1..p-1 since the first component is assumed to be exact)
C
      call dcopy(m*p,0d0,0,modify,1)
      do k=1,p-1
         alpha = (sigma**k-sigma**(p+1)) * c1beta(k+1)
         call daxpy(m,alpha,scdrv1,1,modify(1,k+1),1)
      end do
C
C     (2.2) add modify in order to get the correct result 
C           yn_in(:,k+1) = yn_in(:,k+1) + modify(:,k+1)
C
      call daxpy(m*p,1d0,modify,1,yn_in,1)
C
C     (3) correct the scaled derivatives as well
C
      call dscal(m,sigma**p    ,scdrv0,1)
      call dscal(m,sigma**(p+1),scdrv1,1)

      end

C**************************************************************************
C     ****f* src/convrt
C
C  NAME
C     convrt -- estimate the CONVergence RaTe of Newton's method
C
C  DESCRIPTION
C     This subroutine used the techniques of Hairer/Wanner in order
C     to estimate the convergence rate used for an efficient order 
C     changing strategy.
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine convrt(psi,l,res,m,scal,nresx,nresx1,theta,theta1)
      implicit none
      integer          l, m, i
      real*8           psi, res(m), scal(m), nresx, nresx1, theta,
     $                 theta1, tmp(m), dnrm2
      nresx1 = nresx
      do i=1,m
         tmp(i) = res(i) / scal(i)
      end do
      nresx = dnrm2(m,res,1)
      if (l.ge.2) then
         theta1 = theta
         theta  = nresx / nresx1
         if (l.eq.2) then 
            psi = max(psi,theta)
         else
            psi = max(psi,dsqrt( theta * theta1 ))
         end if
      end if
      end

C**************************************************************************
C     ****f* src/chkrcd
C
C  NAME
C     chkrcd -- CHecK inveRse of the ConDition number
C
C  DESCRIPTION
C     In case of RCOND.gt.0.0 the condition number is estimated when 
C     solving  linear systems.
C     Warnings are issued when numerical problems are likely.
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine chkrcd(Mat,m,RCOND)
C     CHecK RConD: estimate (inverse of the) condition number 
      implicit none
      integer          m, iwork(m), info
      real*8           RCOND, rc, Mat(m,m), dlange, anorm, work(4*m)
      if ( RCOND.gt.0 ) then
         anorm = dlange('i',m,m,Mat,m,work)
         call dgecon('i',m,Mat,m,anorm,rc,work,iwork,info)
         if ( rc .lt. RCOND ) then
            print*, 'GLIMDA WARNING: Matrix singular or badly scaled '
     $           //'(RCOND=',rc,')'
         end if
      end if
      end

C**************************************************************************
C     ****f* src/get_j
C
C  NAME
C     get_j -- compute the Jacobian for Newtons method
C
C  DESCRIPTION
C     This subroutine assembles the Jacobian
C     Jac = A D + h lam B = df/dy dq/dx + h*lam df/dx
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine get_j(m    , n     , z     , x     , t     , h, lam,
     $                 num_A, num_D , num_B , dfyevl, dqxevl, dfxevl,
     $                 fevl , qevl  , ADcnst, mat_A , mat_D , F0    ,
     $                 eps  , Jac   , PRTLEV, ipar  , rpar  , ierr  )
      implicit none
      logical          num_A , num_D , num_B, ADcnst
      integer          m     , n     , PRTLEV,ipar(*), ierr
      real*8           z(n)  , x(m)  , t    , h      , lam  , 
     $                 mat_A(m,n), mat_D(n,m), F0(m) ,
     $                 eps   , Jac(m,m), rpar(*)
      external         dfyevl, dqxevl, dfxevl, fevl  , qevl
C     locals
      logical          naninf
      integer          i, j
      real*8           F1(m) , Q0(n) , Q1(n), safe, delta
C
C (1) compute mat_A
C
      if (.not.ADcnst) then
         call get_A(m,n,z,x,t,num_A,dfyevl,fevl,F0,mat_A,eps,
     $              ipar,rpar,PRTLEV,ierr)
         if (ierr.ne.0) return
      endif
C
C (2) compute mat_D
C
      if (.not.ADcnst) then
         call qevl(m,n,x,t,Q0,rpar,ipar,ierr)
         call get_D(m,n,x,t,num_D,dqxevl,qevl,Q0,mat_D,eps,
     $        ipar,rpar,PRTLEV,ierr)
         if (ierr.ne.0) return
      endif
C
C     (3) compute mat_B -> Jac       
C
      if( num_B ) then
         do i=1,m
            safe  = x(i)
            delta = dsqrt( eps * max( 1.d-5, abs(safe) ) )
            x(i)  = safe + delta
            call fevl(m,n,z,x,t,F1,rpar,ipar,ierr)
            if (ierr.ne.0) return
            do j=1,m
               jac(j,i) = ( F1(j) - F0(j) ) / delta
            end do
            x(i)  = safe
         enddo
      else
         call dfxevl(m,n,z,x,t,Jac,rpar,ipar,ierr)
         if ( ierr.ne.0 ) return
      endif

      if (naninf(Jac,m,m,m)) then
         if (PRTLEV.ge.2) print*,'GLIMDA WARNING: '
     $        //'evaluation of dfx yields NAN or INF'
         ierr = -10
         return
      endif
C
C (4) compute Jac = mat_A*mat_D + h*lam*mat_B
C
      call dgemm('n','n',m,m,n,1d0,mat_A,m,mat_D,n,h*lam,Jac,m)

      end

C**************************************************************************
C     ****f* src/chktol
C
C  NAME
C     chktol -- CHeK TOLerance
C
C  DESCRIPTION
C     Given a vector res(m) and a vector of weights scal(m)
C     this function checks whether each component satisfies
C
C        abs(res(i)) .le. tol * scal(i)
C
C  RESULT
C     If the above relation holds for i=1,..,m, .true. is returned.
C     Otherwise the function returns .false.
C
C  COPYRIGHT
C
C  SOURCE
C
      logical function chktol(res,m,scal,tol)
      implicit none
      integer          m, i
      real*8           res(m), scal(m), tol 
      do i=1,m
         if (abs(res(i)).gt.tol*scal(i)) then
            chktol = .false.
            return
         end if
      end do
      chktol = .true.
      return
      end

C**************************************************************************
C     ****f* src/purify
C
C  NAME
C     purify -- purify (or modify) a Nordsieck vector
C
C  DESCRIPTION
C     When changing the order, the Nordsieck vector has to be extended 
C     (order increase) or the last component has to be deleted (order
C     decrease). In any case: The higher order error terms used in the 
C     scale and modify technique will no longer be correct. 
C     Before modifying the Nordsieck vector, this subroutine attempts
C            to purify the Nordsieck vector by removing higher order 
C            error terms. 
C     After  modifying the Nordsieck vector, this subroutine attempts
C            to correct the Nordsieck vector by adding higher order 
C            error terms. 
C
C  ARGUMENTS
C
C     character c       .. c='p' -> purify, c='m' -> modify
C     real*8    qn(n,p) .. Nordsieck vector
C     integer   n       .. dimension of the Nordsieck vector's components
C     integer   p       .. integration order ( = length of Nordsieck vector)
C     real*8    scdrv1  .. estimate of the scaled derivative of order p+1
C     real*8    c1beta  .. error constant and vector of higher order error 
C                          terms ( coefficients of h^(p+1) q^(p+1) (x,t) )
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine purify(c,qn,n,p,scdrv1,c1beta)
      implicit none
      character        c
      integer          n, p, i
      real*8           qn(n,p), scdrv1(n), c1beta(p), modify(n,p)
C
C     compute 'modify' --> modify(:,i) = c1beta(i) * deriv_p1 (i=2..p)
C
      call dcopy(n*p,0d0,0,modify,1)
      do i=2,p
         call daxpy(n,c1beta(i),scdrv1,1,modify(1,i),1)
      end do
C
C     add 'modify' to the Nordsieck vector --> qn = qn +/- modify
C
      if (c.eq.'p') then
C        purify
         call daxpy(n*p, 1d0,modify,1,qn,1)
      elseif (c.eq.'m') then
C        modify
         call daxpy(n*p,-1d0,modify,1,qn,1)
      end if

      end
C**************************************************************************
C     ****f* src/iniopt
C
C  NAME
C     iniopt -- INItialise OPTions for the simulator
C
C  DESCRIPTION
C     The simulator is controlled using integer options iopt and 
C     floating point options ropt. This subroutine sets defaults and/or 
C     user defined options. The options are copied into more meaningful
C     variables.
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine iniopt(t0    , tend  , iopt  , ropt  , MD    , atolX , 
     $                  rtolX , tolvec,
     $                  m     , PRTLEV, NEWTIT, MAXORD, MINORD, P3METH,
     $                  ORDER , KEEP_P, NSTMAX, NRETRY, TOLF  , TOLX  ,
     $                  DELACL, CUTFAC, INCFAC, DECFAC, HMAX  , HMIN  ,
     $                  UPTH  , DOWNTH, RCOND )
      implicit none
      logical          tolvec
      integer          iopt(9) , MD     , m       , NEWTIT , MAXORD, 
     $                 MINORD  , P3METH , ORDER   , KEEP_P , NSTMAX,
     $                 NRETRY  , PRTLEV, i
      real*8           t0      , tend   , ropt(11), atolX(m), rtolX(m),
     $                 TOLF    , TOLX   ,
     $                 DELACL  , CUTFAC , INCFAC , DECFAC , HMAX  ,
     $                 HMIN    , UPTH   , DOWNTH , RCOND  , dlamch,
     $                 atmp, rtmp
C
C      PRTLEV:  print level
C      NEWTIT:  max number of newton iterations
C      MAXORD:  maximum order to use
C      MINORD:  minimum order to use 
C      P3METH:  the particular order 3 method to use
C      ORDER :  order \in {1,2,3} -> const. order; 0 means 'variable order'
C      KEEP_P:  order increase only after keep_p successful steps 
C      NSTMAX:  maximum number of steps to calculate
C      NRETRY:  maximum number of consecutive retries 
C
      if (iopt(1).eq.0) then
         PRTLEV = 2 
      else
         PRTLEV = iopt(1) 
      endif
      if (iopt(2).le.0) then
         NEWTIT = 5 
      else
         NEWTIT = iopt(2) 
      endif
      if (iopt(3).le.0) then
         MAXORD = md 
      else
         MAXORD = iopt(3) 
      endif
      if (iopt(4).le.0) then
         MINORD = 1 
      else
         MINORD = iopt(4) 
      endif
      if (iopt(5).le.0) then
         P3METH = 1 
      else
         P3METH = iopt(5) 
      endif
      if (iopt(6).le.0) then
         ORDER  = 0 
      else
         ORDER  = iopt(6) 
      endif
      if (iopt(7).le.0) then
         KEEP_P = 5 
      else
         KEEP_P = iopt(7) 
      endif
      if (iopt(8).le.0) then
         NSTMAX = 500000 
      else
         NSTMAX = iopt(8) 
      endif
      if (iopt(9).le.0) then
         NRETRY = 15 
      else
         NRETRY = iopt(9) 
      endif
C
C     TOLF  :  require |res|_err < tolf for convergence 
C     TOLX  :  require |dY |_err < tolx for convergence
C     DELACL:  accept if |delacl*err_const*scdrv1|_err < 1
C     CUTFAC:  use stepsize cutfac*h after 'no convergence' 
C     INCFAC:  increase stepsize by no more than incfac
C     DECFAC:  decrease stepsize by no mare than decfac
C     HMAX  :  maximal stepsize
C     HMIN  :  minimal stepsize
C     UPTH  :  threshold for changing the order upwards
C     DOWNTH:  threshold for changing the order downwards
C     RCOND :  minimum allowed (inverse of the) condition number
C             (0 skips this test)
C
      if (ropt( 1).le.0) then
         TOLF   = 1d-1
      else
         TOLF   = ropt( 1) 
      endif
      if (ropt( 2).le.0) then
         TOLX   = 1d-1
      else
         TOLX   = ropt( 2) 
      endif
      if (ropt( 3).le.0) then
         DELACL = 1d-1 
      else
         delacl = ropt( 3) 
      endif
      if (ropt( 4).le.0) then
         CUTFAC = 1d-1 
      else
         CUTFAC = ropt( 4) 
      endif
      if (ropt( 5).le.0) then
         INCFAC = 2d0  
      else
         INCFAC = ropt( 5) 
      endif
      if (ropt( 6).le.0) then
         DECFAC = 5d-1 
      else
         DECFAC = ropt( 6) 
      endif
      if (ropt( 7).le.0) then
         HMAX = abs(tend-t0)
      else
         HMAX   = ropt( 7)
      endif
      if (ropt( 8).le.0) then
         HMIN   = dlamch( 'Precision' )
      else
         HMIN   = ropt( 8) 
      endif
      if (ropt( 9).le.0) then
         UPTH   = 1d-3 
      else
         UPTH   = ropt( 9) 
      endif
      if (ropt(10).le.0) then
         DOWNTH = 0.8d0 
      else
         DOWNTH = ropt(10) 
      endif
      if (ropt(11).le.0) then
         RCOND  = 0d0  
      else
         RCOND  = ropt(11) 
      endif
C
C     handle tolerances and tolvec
C
      if ( tolvec ) then
         do i=1,m
            if (atolX(i).le.0d0) atolX(i) = 10**(atolX(i))
            if (rtolX(i).le.0d0) rtolX(i) = 10**(rtolX(i))
         enddo
      else
         if (atolX(1).le.0d0) then
            atmp = 10d0**(atolX(1))
         else
            atmp = atolX(1)
         endif
         if (rtolX(1).le.0d0) then
            rtmp = 10d0**(rtolX(1))
         else
            rtmp = rtolX(1)
         endif
         call dcopy(m,atmp,0,atolX(1),1)
         call dcopy(m,rtmp,0,rtolX(1),1)
      endif

      end
C**************************************************************************
C     ****f* src/ifmt
C
C  NAME
C     ifmt -- compute format string for printing integer numbers
C
C  DESCRIPTION
C     Given an integer number n, this function return the string
C     'Ik' where 1<=k<=9 is the number of digits for n.
C 
C  COPYRIGHT
C
C  SOURCE
C
      character*2 function ifmt(n)
      implicit none
      integer n
      ifmt = 'I'//char(int(log(1d0*max(1,n))/log(1d1)+1d0)+ichar('0'))
      end
C**************************************************************************
C     ****f* src/prt_mat
C
C  NAME
C     prt_mat -- PRinT a MATRIX
C 
C  COPYRIGHT
C
C  SOURCE
C
      subroutine prt_mat(m,n,A,neqn,str)
      implicit none
      integer m,n,neqn,i,j
      real*8           A(neqn,n)
      character*(*) str
      print*,str
      do i=1,m
         do j=1,n
C            write(*,'(E14.6,$)') A(i,j)
            write(*,'(1P,E8.1,$)') A(i,j)
         end do
         print*
      end do      
      end
C**************************************************************************
C     ****f* src/prt_vec
C
C  NAME
C     prt_vec -- PRinT a VECtor
C
C  COPYRIGHT
C
C  SOURCE
C
      subroutine prt_vec(m,vec,incr,str)
      implicit none
      integer m,i,incr
      real*8           vec(1+(m-1)*incr)
      character*(*) str
 1    format(5x,a,' =',1P,E10.2,$)
 2    format(1P,E10.2,$)
 3    format(1P,E10.2)
      write(*,1) str,vec(1)
      if (m.eq.1) then
         write(*,*) ''
      else
         do i=2,m-1
            write(*,2) vec(1+(i-1)*incr)
         enddo
         write(*,3) vec(1+(m-1)*incr)
      endif      
      end
C**************************************************************************
