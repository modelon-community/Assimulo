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

import numpy as N
import pylab as P
import nose
from assimulo.solvers import ImplicitEuler
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    eps = 5.e-3
    my = 1./eps
    
    #Define the rhs
    def f(t,y):
        yd_0 = y[1]
        yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
        
        return N.array([yd_0,yd_1])
    
    #Define the jacobian (If no jacobian, one is approximated)
    def jac(t,y):
        jd_00 = 0.0
        jd_01 = 1.0
        jd_10 = -1.0*my-2*y[0]*y[1]*my
        jd_11 = my*(1.-y[0]**2)
        
        return N.array([[jd_00,jd_01],[jd_10,jd_11]])
    
    y0 = [2.0,-0.6] #Initial conditions
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f,y0)
    exp_mod.name = 'Van der Pol (explicit)'
    exp_mod.jac = jac
    
    #Define an explicit solver
    exp_sim = ImplicitEuler(exp_mod) #Create a ImplicitEuler solver
    
    #Sets the parameters
    exp_sim.h = 1e-4 #Stepsize
    exp_sim.usejac = True #If the user defined jacobian should be used or not
    
    #Simulate
    t, y = exp_sim.simulate(2.0) #Simulate 2 seconds
    
    #Plot
    if with_plots:
        P.plot(t,y[:,0], marker='o')
        P.show()

    #Basic test
    x1 = y[:,0]
    assert N.abs(x1[-1]-1.8601438) < 1e-1 #For test purpose

if __name__=='__main__':
    run_example()

