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
from assimulo.thirdparty import mexax

class mexaxpendulum(object):
    def __init__(index,gr=9.81):
        self.index=index
        # initial conditions
        self.p=array([0.,1.])
        self.v=array([0.,0.])
        self.u=array([])
        self.gr=gr
        self.lam=array([0.])
        self.a=array([0.,-gr])
        
    def fprob(self,t,p,v,u,lam,mass,gp,f,pdot,udot,g,gi,fl,qflag,ifail):
        """
        problem function for index-2 single pendulum example
        """
        gr=self.gr # gravitation constant
        exception_text='We should not need to evaluate this {}'

        if qflag[0]:  # set mass matrix
            mass=eye(2)
        if qflag[1]: # set velocity constraint matrix
            gp=array(p).reshape((1,-1))
        if qflag[2]: # raise error
            raise Exception(exception_text.format('C-matrix'))
        if qflag[3]:
            f=array([0.,gr]).reshape([-1,1])
        if qflag[4]:
            pdot = v
        if qflag[5]:
            raise Exception(exception_text.format('udot'))
        if qflag[6]:
            raise Exception(exception_text.format('g-position residual'))
        if qflag[7]:
            gi=array([0.])
        if qflag[8]:
            raise Exception(exception_text.format('fl (df/dla)'))
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
            t=t0
            tfin=te
            itol=0
            rtol=atol=1.e-5
            h=1.e-4
            mxjob=zeros((150,))
            ierr = mexax.mexx(ng,self.fprob,t,tfin,
                 self.p,self.v,self.u,self.a,self.lam,
                 itol,rtol,atol,h,mxjob,ierr,iwk,rwk,
                 self.solout,self.denout,self.fswit)



       
