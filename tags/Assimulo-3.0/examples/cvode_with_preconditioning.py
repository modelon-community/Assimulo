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

"""
Example by Johannes Herold
"""

import numpy as np
import nose
from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem


def run_example(with_plots=True): 
    r"""
    Example to demonstrate the use of a preconditioner
    
    .. math::
        
        \dot y_1 & = 2 t \sin y_1  + t \sin y_2 \\
        \dot y_2 & = 3 t \sin y_1  + 2 t \sin y_2
        
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """              

    #Define the rhs
    def rhs(t, y):
        A = np.array([[2.0, 1.0], [3.0, 2.0]])
        yd = np.dot(A * t, np.sin(y))
        return yd
        
    #Define the preconditioner setup function
    def prec_setup(t, y, fy, jok, gamma, data):
        A = np.array([[2.0, 1.0], [3.0, 2.0]])
        
        #If jok is false the jacobian data needs to be recomputed
        if jok == False:
                
            #Extract the diagonal of the jacobian to form a Jacobi preconditioner
            a0 = A[0, 0] * t * np.cos(y[0])
            a1 = A[1, 1] * t * np.cos(y[1])
            a = np.array([(1. - gamma * a0), (1. - gamma * a1)])
            
            #Return true (jacobian data was recomputed) and the new data
            return [True, a]
            
        #If jok is true the existing jacobian data can be reused 
        if jok == True:
            
            #Return false (jacobian data was reused) and the old data
            return [False, data]

    #Define the preconditioner solve function
    def prec_solve(t, y, fy, r, gamma, delta, data):
        
        #Solve the system Pz = r
        z0 = r[0]/data[0]
        z1 = r[1]/data[1]
           
        z = np.array([z0, z1])
        return z
        
    #Initial conditions
    y0 = [1.0, 2.0]

    #Define an Assimulo problem
    exp_mod = Explicit_Problem(rhs, y0, name = "Example of using a preconditioner in SUNDIALS")

    #Set the preconditioner setup and solve function for the problem
    exp_mod.prec_setup = prec_setup
    exp_mod.prec_solve = prec_solve

    #Create a CVode solver
    exp_sim = CVode(exp_mod)

    #Set the parameters for the solver
    exp_sim.iter = 'Newton'
    exp_sim.discr = 'BDF'
    exp_sim.atol = 1e-5
    exp_sim.rtol = 1e-5
    exp_sim.linear_solver = 'SPGMR'
    exp_sim.precond = "PREC_RIGHT" #Set the desired type of preconditioning

    #Simulate
    t, y = exp_sim.simulate(5)
    
    #Basic verification
    nose.tools.assert_almost_equal(y[-1,0],3.11178295,4)
    nose.tools.assert_almost_equal(y[-1,1],3.19318992,4)
    
    if with_plots:
        exp_sim.plot()

    return exp_mod, exp_sim
    
if __name__=='__main__':
    mod,sim = run_example()
