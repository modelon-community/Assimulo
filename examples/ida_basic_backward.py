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
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem

def run_example(with_plots=True):
    """
    The same as example :doc:`EXAMPLE_cvode_basic`  but now integrated backwards in time.
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    
    #Define the rhs
    def f(t,y,yd):
        res = yd[0] + y[0]
        return N.array([res])
    
    #Define an Assimulo problem
    imp_mod = Implicit_Problem(f,t0=5, y0=0.02695, yd0=-0.02695,
              name = 'IDA Example: $\dot y + y = 0$ (reverse time)')
    
    #Define an explicit solver
    imp_sim = IDA(imp_mod) #Create a IDA solver
    
    #Sets the parameters
    imp_sim.atol = [1e-8] #Default 1e-6
    imp_sim.rtol = 1e-8 #Default 1e-6
    imp_sim.backward = True

    #Simulate
    t, y, yd = imp_sim.simulate(0) #Simulate 5 seconds (t0=5 -> tf=0)
    
    #Basic test
    nose.tools.assert_almost_equal(float(y[-1]), 4.00000000, 3)
    
    #Plot
    if with_plots:
        P.plot(t, y, color="b")
        P.xlabel('Time')
        P.ylabel('State')
        P.title(imp_mod.name)
        P.show()
    
    return imp_mod, imp_sim

if __name__=='__main__':
    mod,sim = run_example()
