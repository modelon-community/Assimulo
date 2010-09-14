#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from lib import sundials_kinsol_core
import numpy as N
import pylab as P
from assimulo.non_linear_problem import *

class KINSOL_Exception(Exception):
    pass
    

class KINSOL:
    
    def __init__(self):
        """
        Create the solver
        """
        
        self.solver = sundials_kinsol_core.KINSOL_wrap()
        
    def solve(self,problem,use_jac = True):
        """
        Function called when solving fuction rhs_fct
        
        Parameters:
            problem:
                instance of NL_problem found in non_linear_problem.py
        """
        # extract info from problem
        if hasattr(problem,'_x0'):
            try:
                x0 = problem.get_x0()
            except NL_Problem_Exception:
                raise KINSOL_Exception("Problem has not implemented method 'get_x0'")
        else:
            raise KINSOL_Exception("Problem has no instance '_x0'")
        
        # calculate dimension
        try:
            if isinstance(x0, int) or isinstance(x0, float):
                x0 = [x0]
            dim = len([N.array(x0, dtype=float)][0])
        except ValueError:
            raise KINSOL_Exception("Initial guess must be a Numpy.array with either ints or floats.")
        
        try:
            tmp = problem.f(x0)
            func = problem.f
        except NL_Problem_Exception:
            raise KINSOL_Exception("Problem has not implemented method 'f'")
        except IndexError:
            raise KINSOL_Exception("Problem has mismatching f and initial guess")
        
        if use_jac and hasattr(problem,'jac'):
            jac = problem.jac
        else:
            jac = None
        
        
        
        # Initialize solver and solve
        self.solver.KINSOL_init(func,x0,dim,jac)

        
        return self.solver.KINSOL_solve()
        