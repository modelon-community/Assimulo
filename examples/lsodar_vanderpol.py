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

import numpy as np
from assimulo.solvers import LSODAR
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    r"""
    Example for the use of LSODAR method to solve
    Van der Pol's equation
    
    .. math::
       
        \dot y_1 &= y_2 \\
        \dot y_2 &= \mu ((1.-y_1^2) y_2-y_1)

    with :math:`\mu=\frac{1}{5} 10^3`.

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
        
        return np.array([yd_0,yd_1])
    
    y0 = [2.0,-0.6] #Initial conditions
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f,y0, name = "LSODAR: Van der Pol's equation")
    
    #Define an explicit solver
    exp_sim = LSODAR(exp_mod) #Create a LSODAR solver
    
    #Sets the parameters
    exp_sim.atol = 1e-4 #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
        
    #Simulate
    t, y = exp_sim.simulate(2.) #Simulate 2 seconds
    
     #Plot
    if with_plots:
        import pylab as pl
        pl.plot(t,y[:,0], marker='o')
        pl.title(exp_mod.name)
        pl.ylabel("State: $y_1$")
        pl.xlabel("Time")
        pl.show()
    
    #Basic test
    x1 = y[:,0]
    assert np.abs(x1[-1] - 1.706168035) < 1e-3
    
    return exp_mod, exp_sim

if __name__=='__main__':
    mod,sim = run_example()

