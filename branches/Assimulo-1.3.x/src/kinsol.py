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

from lib.sundials_kinsol_core import KINSOL_wrap, KINError
from scipy.linalg import pinv2
from numpy.linalg import solve
import numpy as N
import pylab as P
import operator as O
import re
from assimulo.problem_algebraic import *

class KINSOL_Exception(Exception):
    pass
    

class KINSOL:
    
    def __init__(self,problem):
        """
        Create the solver
        
        Parameters::
            problem--
                instance of ProblemAlgebraic found in problem_algebraic.py
        """
        
        self.solver = KINSOL_wrap()
        
        # extract info from problem
        self.problem = problem
        if hasattr(self.problem,'_x0'):
            try:
                _x0 = self.problem.get_x0()
            except ProblemAlg_Exception:
                raise KINSOL_Exception("Problem has not implemented method 'get_x0'")
        else:
            raise KINSOL_Exception("Problem has no instance '_x0'")
        
        # calculate dimension
        try:
            if isinstance(_x0, int) or isinstance(_x0, float):
                self.x0 = [_x0]
            else:
                self.x0 = _x0
            self.dim = len([N.array(self.x0, dtype=float)][0])
        except ValueError:
            raise KINSOL_Exception("Initial guess must be a Numpy.array with either ints or floats.")
        
        # check for functions and test them
        try:
            tmp = self.problem.f(self.x0)
            self.func = self.problem.f
        except ProblemAlg_Exception:
            raise KINSOL_Exception("Problem has not implemented method 'f'")
        except IndexError:
            raise KINSOL_Exception("Problem has mismatching dimensions of f and initial guess")
        
        # check for constraints and test them
        if hasattr(self.problem, 'get_constraints'):
            self.constraints = self.problem.get_constraints()
            if self.constraints != None:
                # test if constraints are of correct type
                if type(self.constraints).__name__ != 'ndarray':
                    raise KINSOL_Exception("Constraints must be of type numpy.ndarray")
                
                if len(self.constraints) != len(self.x0):
                    raise KINSOL_Exception("Constraints must have same length as x0")
                # Test if initial guess x0 is consistant with constraints
                for c,xi in zip(self.constraints,self.x0):
                    if re.search('float',type(c).__name__) == None:
                        print "Problem with: ", c, type(c).__name__
                        raise KINSOL_Exception("Constraints must contain floats.")
                    if abs(c) > 2:
                        raise KINSOL_Exception("Entries in constraint vector must be between -2 and 2, see documentation.")
                    if O.xor(c>=0,xi>=0):
                        raise KINSOL_Exception("Initial guess does not fulfill applied constraints.")
        else:
            self.constraints = None
            
        self._use_jac = True
        self.pinv_count = 0
        self.lin_count = 0
                
    def set_jac_usage(self,use_jac):
        """
        Set whether to use the jacobian supplied by the model
        or if we are to calculate it numericaly
            
        Parameters::
        
            use_jac --
                Boolean set to True if the jacobian is to be 
                supplied by the model
                
        """
        if type(use_jac).__name__ == 'bool':
            self._use_jac = use_jac
        else:
            raise KINSOL_Exception("The variable sent to 'set_jac_usage' must be a boolean.")
    
    def solve(self):
        """
        Function called when solving function rhs_fct
        
        """
        # check for jacobian and set it if present and to be used
        if self._use_jac and hasattr(self.problem,'jac'):
            jac = self.problem.jac
        else:
            jac = None
            
        # Initialize solver and solve        
        
        solved = False
        res = N.zeros(self.x0.__len__())
        while not solved and self.pinv_count < 10:
            try:
                self.solver.KINSOL_init(self.func,self.x0,self.dim,jac,self.constraints)
                res = self.solver.KINSOL_solve()
                solved = True
            except KINError as error:
                if error.value == -11 :
                    # the problem is caused by a singular jacobian try a pinv step
                    self.pinv_count += 1
                    self._do_pinv_step()
                elif error.value == -6 or error.value == -7:
                    self._brute_force()    
                else:
                    # Other error, send onward as exception
                    raise KINSOL_Exception(error.msg[error.value])
        
        if not solved:
            raise KINSOL_Exception("Singular Jacobian. Tried using pseudo inverse but stopped after ten steps.")
           
        if self.pinv_count != 0:
            print self.pinv_count, " steps using the pseudo inverse performed."
            print self.lin_count, " steps using the numpy.linalg.solve function performed."
            
        return res
            
    def _do_pinv_step(self):
        """
        Method used to perform a step using the pseudo inverse
        """
        print "Trying to do a step with the pseudo inverse"
        if self._use_jac:
            if hasattr(self.problem,'jac'):
                # Extract data from problem
                x0 = self.x0
                fx = self.func(x0)
                J  = self.problem.jac(x0)
                
                # Calculate pseudo inverse and calculate new step
                Jpinv = pinv2(J)
                dx = N.dot(Jpinv,-fx)
                
                # Do step in problem
                self.x0 = x0 + dx
                self.problem._x0 = self.x0
                
            else:
                raise KINSOL_Exception("Singular jacobian. Trying to do a step using the pseudo inverse, but no jacobian supplied. Using the jacobian of KINSOL is not implemented yet.")
                
        else:
            raise KINSOL_Exception("Singular jacobian. Trying to do a step using the pseudo inverse but 'use_jac' is set to false.")
                
    def _brute_force(self):
        """
        Method used to use a bit of brute force.
        In the case of difficulties with reducing the residual, the method will newton iterate
        , without any linesearch, until the (norm of the) residual is once again declining.
        """
        print "Trying to solve the problem using a bit of brute force"
        if self._use_jac:
            if hasattr(self.problem,'jac'):
                # Extract data from problem
                x0 = self.x0
                fx = self.func(x0)
                                
                tol = N.linalg.norm(fx)
                for i in N.arange(10):
                    
                    fx = self.func(x0)
                    print "|fx|: ", abs(N.linalg.norm(fx))
                    if N.linalg.norm(fx) < tol :
                        self.lin_count += i
                        self.x0 = x0
                        self.problem._x0 = self.x0
                        break
                    
                    tol = N.linalg.norm(fx)
                    J  = self.problem.jac(x0)
                
                    dx = solve(J,-fx)
                
                    x0 = x0 + dx
                    self.problem._x0 = self.x0
                
            else:
                raise KINSOL_Exception("Singular jacobian. Trying to do a step using the pseudo inverse, but no jacobian supplied. Using the jacobian of KINSOL is not implemented yet.")
                
        else:
            raise KINSOL_Exception("Singular jacobian. Trying to do a step using the pseudo inverse but 'use_jac' is set to false.")
                
        