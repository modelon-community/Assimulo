#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2021 Modelon AB
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
import scipy.sparse as sps
from assimulo.solvers import Radau5ODE
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    r"""
    Example for demonstrating the use of a user supplied Jacobian (sparse).
    Based on the SUNDIALS example cvRoberts_sps.c
    
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
        return np.array([yd_0,yd_1,yd_2])
    
    #Defines the Jacobian
    def jac(t,y):
        
        colptrs = [0,2,5,7]
        rowvals = [0, 1, 0, 1, 2, 0, 1]
        data = [-0.04, 0.04, 1e4*y[2], -1e4*y[2]-6e7*y[1], 6e7*y[1], 1e4*y[1], -1e4*y[1]]

        J = sps.csc_matrix((data, rowvals, colptrs))
        return J
    
    #Defines an Assimulo explicit problem
    y0 = [1.0, 0.0, 0.0] #Initial conditions

    exp_mod = Explicit_Problem(f, y0, name = 'Example using analytic (sparse) Jacobian')
    
    exp_mod.jac = jac #Sets the Jacobian
    exp_mod.jac_nnz = 7
   
    
    exp_sim = Radau5ODE(exp_mod) #Create a Radau5 solver
    
    #Set the parameters
    exp_sim.atol = [1e-8,1e-14,1e-6] #Default 1e-6
    exp_sim.rtol = 1e-4 #Default 1e-6
    exp_sim.linear_solver = 'sparse'
    
    #Simulate
    t, y = exp_sim.simulate(0.4) #Simulate 0.4 seconds
    
    #Plot
    if with_plots:
        import pylab as pl
        pl.plot(t,y[:,1],linestyle="dashed",marker="o") #Plot the solution
        pl.xlabel('Time')
        pl.ylabel('State')
        pl.title(exp_mod.name)
        pl.show()
    
    #Basic tests
    assert abs(y[-1][0] - 0.9851) < 1e-3
    
    return exp_mod, exp_sim

if __name__=='__main__':
    mod,sim = run_example()
