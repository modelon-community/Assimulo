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
import scipy.sparse as SP
import pylab as P
import nose
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem


def run_example(with_plots=True):
    r"""
    Example for demonstrating the use of a user supplied Jacobian (sparse).
    Note that this will only work if Assimulo has been configured with
    Sundials + SuperLU. Based on the SUNDIALS example cvRoberts_sps.c
    
    ODE:
    
    .. math::
       
       \dot y_1 &= -0.04y_1 + 1e4 y_2 y_3 \\
       \dot y_2 &= - \dot y_1 - \dot y_3 \\
       \dot y_3 &= 3e7 y_2^2
       
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    
    #Defines the rhs
    def f(t,y):
        yd_0 = -0.04*y[0] + 1e4*y[1]*y[2]
        yd_2 = 3e7*y[1]*y[1]
        yd_1 = -yd_0 - yd_2
        return N.array([yd_0,yd_1,yd_2])
    
    #Defines the Jacobian
    def jac(t,y):
        
        colptrs = [0,3,6,9]
        rowvals = [0, 1, 2, 0, 1, 2, 0, 1, 2]
        data = [-0.04, 0.04, 0.0, 1e4*y[2], -1e4*y[2]-6e7*y[1], 6e7*y[1], 1e4*y[1], -1e4*y[1], 0.0]

        J = SP.csc_matrix((data, rowvals, colptrs))
        return J
    
    #Defines an Assimulo explicit problem
    y0 = [1.0,0.0,0.0] #Initial conditions

    exp_mod = Explicit_Problem(f,y0, name = 'Example using analytic (sparse) Jacobian')
    
    exp_mod.jac = jac #Sets the Jacobian
    exp_mod.jac_nnz = 9
   
    
    exp_sim = CVode(exp_mod) #Create a CVode solver
    
    #Set the parameters
    exp_sim.iter = 'Newton' #Default 'FixedPoint'
    exp_sim.discr = 'BDF' #Default 'Adams'
    exp_sim.atol = [1e-8,1e-14,1e-6] #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    exp_sim.linear_solver = "sparse"
    
    #Simulate
    t, y = exp_sim.simulate(0.4) #Simulate 0.4 seconds
    
    #Basic tests
    nose.tools.assert_almost_equal(y[-1][0],0.9851,3)
        
    #Plot
    if with_plots:
        P.plot(t,y[:,1],linestyle="dashed",marker="o") #Plot the solution
        P.xlabel('Time')
        P.ylabel('State')
        P.title(exp_mod.name)
        P.show()
        
    return exp_mod, exp_sim


if __name__=='__main__':
    mod,sim = run_example()
