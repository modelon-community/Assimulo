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
    """
    The same as example :doc:`EXAMPLE_cvode_basic`  but now integrated backwards in time.
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    #Define the rhs
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])
    
    #Define an Assimulo problem
    exp_mod = Explicit_Problem(f,t0=5, y0=0.02695, name = r'CVode Test Example (reverse time): $\dot y = - y$ ')
    
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

    #print 'y(5) = {}, y(0) ={}'.format(y[0][0],y[-1][0])
    
    #Basic test
    nose.tools.assert_almost_equal(float(y[-1]), 4.00000000, 3)
    
    #Plot
    if with_plots:
        P.plot(t, y, color="b")
        P.title(exp_mod.name)
        P.ylabel('y')
        P.xlabel('Time')
        P.show()
        
    return exp_mod, exp_sim

if __name__=='__main__':
    mod,sim = run_example()
