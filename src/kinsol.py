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
import operator as O
import re
from assimulo.problem_algebraic import *

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
            except ProblemAlg_Exception:
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
        except ProblemAlg_Exception:
            raise KINSOL_Exception("Problem has not implemented method 'f'")
        except IndexError:
            raise KINSOL_Exception("Problem has mismatching f and initial guess")
        
        if use_jac and hasattr(problem,'jac'):
            jac = problem.jac
        else:
            jac = None
            
        if hasattr(problem, 'get_constraints'):
            constraints = problem.get_constraints()
            if constraints != None:
                # test if constraints are of correct type
                if type(constraints).__name__ != 'ndarray':
                    raise KINSOL_Exception("Constraints must be of type numpy.ndarray")
                
                if len(constraints) != len(x0):
                    raise KINSOL_Exception("Constraints must have same length as x0")
                # Test if initial guess x0 is consistant with constraints
                for c,xi in zip(constraints,x0):
                    if re.search('float',type(c).__name__) == None:
                        print "Problem with: ", c, type(c).__name__
                        raise KINSOL_Exception("Constraints must contain floats.")
                    if abs(c) > 2:
                        raise KINSOL_Exception("Entries in constraint vector must be between -2 and 2, see documentation.")
                    if O.xor(c>=0,xi>=0):
                        raise KINSOL_Exception("Initial guess does not fulfill applied constraints.")
                
        else:
            constraints = None
        
        
        # Initialize solver and solve        
        self.solver.KINSOL_init(func,x0,dim,jac,constraints)

        return self.solver.KINSOL_solve()
        