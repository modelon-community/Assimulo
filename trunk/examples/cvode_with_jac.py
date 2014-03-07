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
    r"""
    Example for demonstrating the use of a user supplied Jacobian
    
    ODE:
    
    .. math::
       
       \dot y_1 &= y_2 \\
       \dot y_2 &= -9.82
       
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    
    #Defines the rhs
    def f(t,y):
        yd_0 = y[1]
        yd_1 = -9.82
        #print y, yd_0, yd_1
        return N.array([yd_0,yd_1])
    
    #Defines the Jacobian
    def jac(t,y):
        j = N.array([[0,1.],[0,0]])
        return j
    
    #Defines an Assimulo explicit problem
    y0 = [1.0,0.0] #Initial conditions

    exp_mod = Explicit_Problem(f,y0, name = 'Example using analytic Jacobian')
    
    exp_mod.jac = jac #Sets the Jacobian
   
    
    exp_sim = CVode(exp_mod) #Create a CVode solver
    
    #Set the parameters
    exp_sim.iter = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = 1e-5 #Default 1e-6
    exp_sim.rtol = 1e-5 #Default 1e-6
    
    #Simulate
    t, y = exp_sim.simulate(5, 1000) #Simulate 5 seconds with 1000 communication points
    
    #Basic tests
    nose.tools.assert_almost_equal(y[-1][0],-121.75000000,4)
    nose.tools.assert_almost_equal(y[-1][1],-49.100000000)
        
    #Plot
    if with_plots:
        P.plot(t,y,linestyle="dashed",marker="o") #Plot the solution
        P.xlabel('Time')
        P.ylabel('State')
        P.title(exp_mod.name)
        P.show()
        
    return exp_mod, exp_sim


if __name__=='__main__':
    mod,sim = run_example()
