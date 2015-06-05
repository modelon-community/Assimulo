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
from assimulo.solvers import Dopri5
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    r"""
    Example to demonstrate the use of the Runge-Kutta solver DOPRI5
    for the linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
    
    """    
    #Defines the rhs
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])

    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f, 4.0,
              name = 'DOPRI5 Example: $\dot y = - y$')
    
    exp_sim = Dopri5(exp_mod) #Create a Dopri5 solver

    #Simulate
    t, y = exp_sim.simulate(5) #Simulate 5 seconds
    
    #Basic test
    nose.tools.assert_almost_equal(float(y[-1]),0.02695199,5)
    
    #Plot
    if with_plots:
        P.plot(t,y)
        P.title(exp_mod.name)
        P.xlabel('Time')
        P.ylabel('State')
        P.show()
        
    return exp_mod, exp_sim

if __name__=='__main__':
    mod,sim = run_example()
