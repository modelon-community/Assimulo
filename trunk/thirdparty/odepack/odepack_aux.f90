! Copyright (C) 2010 Modelon AB
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 3 of the License.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.


subroutine set_lsod_common(meth,nq,nqu,miter,maxord,meo,nqnyh,nst,nfe,nje,init,tn,conit,el, nge, hu, jstart)
! helper subroutine to fill the LSOD* common blocks DLS001, 
! This is needed to avoid the direct access of the common block in Python
!   set_1lsod_1common(**args)

implicit none

integer, intent(in), optional :: meth,nq,nqu,miter,maxord,meo,nqnyh,nst,nfe,nje,init, nge, jstart
double precision, intent(in):: tn,conit,el(13),  hu

! variables in the common blocks 
! those set from the parameter list got a "_1"-postfix

integer ::                  init_1,mxstep,mxhnil,nhnil,nslast,nyh_1,     &     
                            ialth, ipup, lmax, meo_1, nqnyh_1, nslp,     &    
                            icf, ierpj, iersl, jcur, jstart_1, kflag,  &
                            l,lyh, lewt, lacor, lsavf, lwm, liwm,    &
                            meth_1, miter_1,maxord_1, maxcor, msbp, mxncf, &
                            n, nq_1, nst_1, nfe_1, nje_1, nqu_1,   &
                            i_rootcommon, ngc, nge_1
double precision :: conit_1, crate, el_1, elco, hold,             &    ! rowns(209)
                            rmax, tesco,                       &
                            ccmax, el0, h, hmin,                     &
                            hmxi, hu_1, rc, tn_1, uround ,                &
                            r_rootcommon
common /dls001/ conit_1, crate, el_1(13), elco(13,12), hold,             &    ! rowns(209)
                            rmax, tesco(3,12),                       &
                            ccmax, el0, h, hmin,                     &
                            hmxi, hu_1, rc, tn_1, uround,                &
                            init_1,mxstep,mxhnil,nhnil,nslast,nyh_1,     &    !iownd(6), 
                            ialth, ipup, lmax, meo_1, nqnyh_1, nslp,     &    !iowns(6)
                            icf, ierpj, iersl, jcur, jstart_1, kflag,  &
                            l,lyh, lewt, lacor, lsavf, lwm, liwm,    &
                            meth_1,miter_1, maxord_1, maxcor, msbp, mxncf, &
                            n, nq_1, nst_1, nfe_1, nje_1, nqu_1
common /dlsr01/ r_rootcommon(5), i_rootcommon(7), ngc, nge_1  ! root function counters

meth_1  = meth
nq_1    = nq
if (present(nqu)) then
   nqu_1  = nqu
end if   
if (present(jstart)) then
   jstart_1=jstart
end if   
h=hu
hu_1=hu
hold=hu
jstart_1=jstart
miter_1 = miter
maxord_1 = maxord
meo_1    = meo
nqnyh_1 = nqnyh
init_1 = init
rc = rc*el(1)/el0
el0   = el(1)
el_1    = el
conit_1 = conit
tn_1=tn
nst_1=nst
nfe_1=nfe
nje_1=nje
nge_1=nge
return
end subroutine set_lsod_common
subroutine get_lsod_common(hu_,nqu_,nq_, nyh_, nqnyh_)
! helper subroutine to read the LSOD* common blocks DLS001, 
! This is needed to avoid the direct access of the common block in Python
!   get_lsod_common()

implicit none

double precision, intent(out):: hu_
integer, intent(out) :: nq_, nqu_, nyh_, nqnyh_

integer ::                  init,mxstep,mxhnil,nhnil,nslast,nyh,     &     
                            ialth, ipup, lmax, meo, nqnyh, nslp,     &    
                            icf, ierpj, iersl, jcur, jstart, kflag,  &
                            l,lyh, lewt, lacor, lsavf, lwm, liwm,    &
                            meth, miter,maxord, maxcor, msbp, mxncf, &
                            n, nq, nst, nfe, nje, nqu
double precision :: conit, crate, el, elco, hold,             &    ! rowns(209)
                            rmax, tesco,                       &
                            ccmax, el0, h, hmin,                     &
                            hmxi, hu, rc, tn, uround 
common /dls001/ conit, crate, el(13), elco(13,12), hold,             &    ! rowns(209)
                            rmax, tesco(3,12),                       &
                            ccmax, el0, h, hmin,                     &
                            hmxi, hu, rc, tn, uround,                &
                            init,mxstep,mxhnil,nhnil,nslast,nyh,     &    !iownd(6), 
                            ialth, ipup, lmax, meo, nqnyh, nslp,     &    !iowns(6)
                            icf, ierpj, iersl, jcur, jstart, kflag,  &
                            l,lyh, lewt, lacor, lsavf, lwm, liwm,    &
                            meth,miter, maxord, maxcor, msbp, mxncf, &
                            n, nq, nst, nfe, nje, nqu
hu_=hu
nq_=nq
nqu_=nqu
nyh_=nyh
nqnyh_=nqnyh
return 
end subroutine get_lsod_common
