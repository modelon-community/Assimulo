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
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    r"""
    Demonstration of the use of CVode by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """

    
    #Define the rhs
    def f(t,y):
        ydot = -y[0]
        return np.array([ydot])
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f, y0=4, name = r'CVode Test Example: $\dot y = - y$')
    
    #Define an explicit solver
    exp_sim = CVode(exp_mod) #Create a CVode solver
    
    #Sets the parameters
    exp_sim.iter  = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = [1e-4] #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6

    #Simulate
    t1, y1 = exp_sim.simulate(5,100) #Simulate 5 seconds
    t2, y2 = exp_sim.simulate(7) #Simulate 2 seconds more
    
    #Plot
    if with_plots:
        import pylab as pl
        pl.plot(t1, y1, color="b")
        pl.plot(t2, y2, color="r")
        pl.title(exp_mod.name)
        pl.ylabel('y')
        pl.xlabel('Time')
        pl.show()
    
    #Basic test
    assert abs(y2[-1][0] - 0.00347746) < 1e-5
    assert abs(exp_sim.get_last_step() - 0.0222169642893) < 1e-3
    
    return exp_mod, exp_sim

if __name__=='__main__':
    mod,sim = run_example()
