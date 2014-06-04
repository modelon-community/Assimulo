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
from assimulo.solvers import KINSOL
from assimulo.problem import Algebraic_Problem

def run_example(with_plots=True):
    r"""
    Example to demonstrate the use of the Sundials solver Kinsol with
    a user provided Jacobian.
    
    on return:
    
       - :dfn:`alg_mod`    problem instance
    
       - :dfn:`alg_solver`    solver instance
    
    """
    #Define the res
    def res(y):
        r1 = 2*y[0]+3*y[1]-6
        r2 = 4*y[0]+9*y[1]-15
        return N.array([r1,r2])
        
    def jac(y):
        return N.array([[2.,3.],[4.,9.]])
    
    #Define an Assimulo problem
    alg_mod = Algebraic_Problem(res, y0=[0,0], jac=jac, name='KINSOL example with Jac')
    
    #Define the KINSOL solver
    alg_solver = KINSOL(alg_mod)
    
    #Sets the parameters
    
    #Solve
    y = alg_solver.solve()
    
    #Basic test
    nose.tools.assert_almost_equal(y[0], 1.5, 5)
    nose.tools.assert_almost_equal(y[1], 1.0, 5)
    
    return alg_mod, alg_solver
    
if __name__=='__main__':
    mod, solv = run_example()

