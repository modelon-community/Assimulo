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
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    
    #Define the rhs
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f,t0=5, y0=0.02695)
    exp_mod.name = 'Simple CVode Example (backward)'
    
    #Define an explicit solver
    exp_sim = CVode(exp_mod) #Create a CVode solver
    
    #Sets the parameters
    exp_sim.iter  = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = [1e-8] #Default 1e-6
    exp_sim.rtol = 1e-8 #Default 1e-6
    exp_sim.backward = True

    #Simulate
    t, y = exp_sim.simulate(0) #Simulate 5 seconds (t0=5 -> tf=0)
    
    #Basic test
    print y
    nose.tools.assert_almost_equal(y[-1], 4.00000000, 3)
    
    #Plot
    if with_plots:
        P.plot(t, y, color="b")
        P.show()

if __name__=='__main__':
    run_example()
