#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from scipy import *
from assimulo.lib import mexax

class Mexaxpendulum(object):
    def __init__(self,index,gr=9.81):
        self.index=index
        # initial conditions
        self.p=array([1.,0.])
        self.np=2
        self.v=array([0.,0.])
        self.nv=2
        self.u=array([0.])   # to avoid zero dimensioning in Fortran
        self.nu=0
        self.gr=gr
        self.lam=array([0.])
        self.nl=1
        self.a=array([0.,-gr])
        
    def fprob(self,nl,ng,nu,t,p,v,u,lam,mass,gp,f,pdot,udot,g,gi,fl,qflag,np,nv,ldg):
        """
        problem function for index-2 single pendulum example
        """
        ifail=0
        gr=self.gr # gravitation constant
        exception_text='We should not need to evaluate this {}'

        if qflag[0]:  # set mass matrix
            mass=eye(2)
        if qflag[1]: # set velocity constraint matrix
            gp=array(p).reshape((1,-1))
        if qflag[2]: # raise error
            raise Exception(exception_text.format('C-matrix'))
        if qflag[3]:
            f=array([0.,-gr]).reshape([-1,1])
        if qflag[4]:
            pdot = v.copy()
        if qflag[5]:
            raise Exception(exception_text.format('udot'))
        if qflag[6]:
            raise Exception(exception_text.format('g-position residual'))
        if qflag[7]:
            gi=array([0.])
        if qflag[8]:
            raise Exception(exception_text.format('fl (df/dla)'))
        return mass,gp,f,pdot,udot,g,gi,ifail
    def solout(self,t,p,v,u,a,rlam,infos,irtrn):
        raise Exception('solout only in dummy mode')
    def denout(self,t,p,v,u,a,rlam,infos,irtrn):
        raise Exception('denout only in dummy mode')
    def fswit(self,t,p,v,u,a,rlam,g):
        raise Exception('fswit only in dummy mode')
    def simulate(self,t0,te):
            if self.index==2:
                ng=0
                nl=1
            else:
                raise Exception('Only index 2 implemented so far')
            # Work arrays and their dimensions
            np=self.np;nv=self.nv;nu=self.nu
            liwk=np + 4*nv + nl + nu + 60
            iwk=empty((liwk,),dtype=int)
            ngl=max(ng,nl)
            lrwk=(nv+nl)**2+np*(ngl+18)+nv*(nv+45)+28*nl+max(ng,1)+18*max(nu,1)+50
            rwk=empty((lrwk,),dtype=float)
            t=t0
            tfin=te
            itol=0
            rtol=atol=1.e-5
            h=1.e-4
            mxjob=zeros((150,))
            mxjob[11-1]=2
            print 'v',self.v
            [t, p, v, u, a, lam, h, mxjob, ierr, iwk, rwk]= \
                 mexax.mexx(nl,ng,nu,self.fprob,t,tfin,
                 self.p,self.v,self.u,self.a,self.lam,
                 itol,rtol,atol,h,mxjob,iwk,rwk,
                 self.solout,self.denout,self.fswit,lrwk=lrwk,liwk=liwk)
            print 'Simulated until {} with return code ierr {} and step size {}'.format(t,ierr,h)
            print 't,p,v,lam',t,p,v,lam
            print mxjob[50:74]
if __name__=='__main__':
    mexax_test=Mexaxpendulum(2,gr=13.750371636040738)
    mexax_test.simulate(0.,10.)



       
