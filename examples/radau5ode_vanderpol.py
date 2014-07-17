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
from assimulo.solvers import Radau5ODE
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    r"""
    Example for the use of the implicit Euler method to solve
    Van der Pol's equation
    
    .. math::
       
        \dot y_1 &= y_2 \\
        \dot y_2 &= \mu ((1.-y_1^2) y_2-y_1)

    with :math:`\mu= 10^6`.

    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance

    """    
    #Define the rhs
    def f(t,y):
        eps = 1.e-6
        my = 1./eps
        yd_0 = y[1]
        yd_1 = my*((1.-y[0]**2)*y[1]-y[0])
        
        return N.array([yd_0,yd_1])
    
    y0 = [2.0,-0.6] #Initial conditions
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f,y0)
    exp_mod.name = 'Van der Pol (explicit)'
    
    #Define an explicit solver
    exp_sim = Radau5ODE(exp_mod) #Create a Radau5 solver
    
    #Sets the parameters
    exp_sim.atol = 1e-4 #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    exp_sim.inith = 1.e-4 #Initial step-size
    
    #Simulate
    t, y = exp_sim.simulate(2.) #Simulate 2 seconds
    print(y)
    #Plot
    if with_plots:
        P.plot(t,y[:,0])#, marker='o')
        P.xlabel('Time')
        P.ylabel('State')
        P.title(exp_mod.name)
        P.show()

    #Basic test
    x1 = y[:,0]
    assert N.abs(x1[-1]-1.706168035) < 1e-3 #For test purpose
    return exp_mod, exp_sim

if __name__=='__main__':
    mod,sim = run_example()

