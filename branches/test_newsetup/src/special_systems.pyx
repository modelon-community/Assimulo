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

import assimulo.problem as ap
from assimulo.support import set_type_shape_array
import numpy as N
cimport numpy as N

cdef class cMechanical_System:
    u"""
        Special problem class for (constrained) mechanical systems:

        .. math::
            :nowrap:
            
            \\begin{eqnarray*}
                \dot{p} & = & v \\\\
                M(p) \dot{v} & = & f(t,p,v)-G(p)^\mathrm{T} \lambda  \\\\
                0 & = & g(p)
            \end{eqnarray*}
        
 
        Parameters ::
        
            n_p   number of position variables
            forces(t, p, v) applied forces with len(y)=2*np
            n_la  number of constraints
            GT(p)  n_p x n_la constraint matrix, p=y[:np]
                  (required if n_la > 0)
            pos0, vel0, lam0
                Defines the initial values for positions, velocities 
                and constraint forces
            posd0, veld0
                Defines the initial derivatives for positions and velocities 
            t0
                Defines the initial time
            mass_matrix(p)  n_p x n_p nonsingular mass matrix
                  (if not defined it is assumed to be the identity matrix)
            constr_3(t,y)  n_la index-3 constraints
            constr_2(t,y)  n_la index-2 constraints
                  (optional)
            constr_1(t,y) n_la index-1 constraints
            index  sets the type of equations to be solved
                 'ind3' Index-3 DAE (position constraints)
                 'ind2' Index-2 DAE (velocity constraints)
                 'ind1' Index-1 DAE (lambda constraints)
                 'ostab2' overdetermined stabilized index 2
                 'ostab1' overdetermined stabilized index 1
                 
    """         
    def __init__(self, int n_p,  object forces, int n_la, object pos0,
                 object vel0, object lam0, object posd0,
                 object veld0, object GT,  
                 double t0 = 0.0, object mass_matrix = None, 
                 object constr3 = None, object constr2 = None, 
                 object constr1 = None, p0=None, sw0=None):
                                        
        self.pos0 = set_type_shape_array(pos0)
        self.vel0 = set_type_shape_array(vel0)
        self.posd0= set_type_shape_array(posd0)
        self.veld0= set_type_shape_array(veld0)
        if [constr1,constr2,constr3]==3*[None]:
            self.constrained=False
        else:
            self.constrained=True    
        if self.constrained:
            if lam0 is None:
                raise ValueError('lam0 should not be None.')
            self.lam0 = set_type_shape_array(lam0)         
        self.n_p=n_p
        self.n_la=n_la
        self.forces = forces
        self.GT=GT
        self.t0=t0
        self.constr3=constr3
        self.constr2=constr2
        self.constr1=constr1
        self.mass_matrix=mass_matrix
        self.sw0=sw0
            
    def make_res(self,index):
        n_p,n_v,n_la=self.n_p, 2*self.n_p, self.n_la
        M=self.mass_matrix
        
        def set_constraints(func,index):
            if index=='ind1':
                constraint=self.constr1
            elif index=='ind2':
                constraint=self.constr2
            elif index=='ind3':
                constraint=self.constr3
            elif index in ('ovstab2', 'ggl2'):
                constraint=lambda t,y: N.hstack((self.constr2(t,y),
                                                 self.constr3(t,y)))
            elif index=='ovstab1':
                constraint=lambda t,y: N.hstack((self.constr1(t,y),
                                                self.constr3(t,y),
                                                self.constr2(t,y)))         
            else:
                raise Exception("index should be one of 'ind1', 'ind2', ind3'"+ 
                                "ovstab2", "ovstab1","ggl2")
            
            if not index =="ggl2":          
                def with_constraint(t,y,yd):
                    p,v,la=y[0:n_p], y[n_p:n_v], y[n_v:]
                    residual=func(t,y,yd)
                    return N.hstack((
                          residual[0:n_p],
                          residual[n_p:n_v]+N.dot(self.GT(p),la.reshape(-1,)), 
                          constraint(t,y)
                          ))  
            else:
                def with_constraint(t,y,yd):
                    p,v,la,mue=y[0:n_p], y[n_p:n_v], y[n_v:n_v+n_la],y[n_v+n_la:]
                    residual=func(t,y,yd)
                    return N.hstack((
                          residual[0:n_p]+N.dot(self.GT(p),mue.reshape(-1,)),
                          residual[n_p:n_v]+N.dot(self.GT(p),la.reshape(-1,)), 
                          constraint(t,y)
                          )) 
            
            return with_constraint

        def res(t,y,yd):
            p,pd=y[0:n_p], yd[0:n_p]
            v,vd=y[n_p:n_v], yd[n_p:n_v]
            Mvd = N.dot(M,vd) if M is not None else vd
            
            return N.hstack((pd - v, Mvd - self.forces(t,p,v)))
            
        if n_la==0:
            return res
        else:   
            return set_constraints(res,index)
            
    def generate_problem(self,index):
        # 0. Input check
        index_values=['ind0', 'ind1','ind2', 'ind3','ovstab1','ovstab2','ggl1','ggl2'] 
        index_error= 'index got not correct value.\n Should be one of {}'.format(index_values)         
        if not (index is None or index in index_values):
            raise ValueError(index_error)
        
        # 1. Treatment of initial conditions depending on the index
        y0=N.hstack((self.pos0,self.vel0))
        yd0=N.hstack((self.posd0,self.veld0))
        if self.constrained:
            if index is None or index=='ind0':
                raise ValueError(index_error)
            y0=N.hstack((y0,self.lam0))
            yd0=N.hstack((yd0,N.zeros(self.lam0.shape)))
        # 2. Indicating algebraic variables    
        if index == 'ind1':
            algvar = (self.pos0.size + self.vel0.size) * [1]\
                        + self.lam0.size*[1]
        elif index == 'ind2':
            algvar = (self.pos0.size + self.vel0.size) * [1] \
                        + self.lam0.size*[0]
        elif index == 'ind3':
            algvar = self.pos0.size * [1] \
                        +self.vel0.size* [0] + self.lam0.size*[0]   
        elif index == 'ovstab1':
            algvar = (self.pos0.size + self.vel0.size) * [1]\
                        + self.lam0.size*[1]
            neq=len(algvar)+2*self.lam0.size            
        elif index == 'ovstab2':
            algvar = (self.pos0.size + self.vel0.size) * [1] \
                        + self.lam0.size*[0] 
            neq=len(algvar)+self.lam0.size                          
        elif index == 'ggl1':
            mue = N.zeros(self.lam0.shape)
            gamma = N.zeros(self.lam0.shape)
            y0 = N.hstack((y0,mue,gamma))
            yd0 = N.hstack((yd0,mue,gamma))
            algvar = (self.pos0.size + self.vel0.size) * [1] \
                        + 3*self.lam0.size*[0]
        elif index == 'ggl2': 
            mue = N.zeros(self.lam0.shape)
            y0 = N.hstack((y0,mue))
            yd0 = N.hstack((yd0,N.zeros(mue.shape)))
            algvar = (self.pos0.size + self.vel0.size) * [1] \
                         + 2*self.lam0.size*[0]
        elif index is None:
            algvar = (self.pos0.size + self.vel0.size) * [1]
        if index in ('ovstab2','ovstab1'):
            problem=ap.Overdetermined_Problem(self.make_res(index), y0, yd0, self.t0, self.sw0)
            problem.neq=neq           
        else:
            problem=ap.Implicit_Problem(self.make_res(index), y0, yd0, self.t0, self.sw0)
        problem.algvar=algvar
        return problem          
                    
                       
class Mechanical_System(cMechanical_System):
    """
        Problem for our explicit integrators (ODEs). A problem
        consists of the right-hand-side and some initial conditions.
 
        Parameters::
            n_p   number of position variables
            n_la  number of constraints
            f(t,y) applied forces with len(y)=2*np
            G(p)  n_la x n_p constraint matrix, p=y[:np]
                  (required if n_la > 0)
            constr_3(t,p)  n_la index-3 constraints
            constr_2(t,p,v) n_la index-2 constraints
                  (optional)
            constr_1(t,p,v,la) n_la index-1 constraints
            M(p)  n_p x n_p nonsingular mass matrix
                  (if not defined it is assumed to be the identity matrix)
            x0
                Defines the starting values of [p0,v0,lambda0]
            t0
                Defines the starting time

    """
    pass
