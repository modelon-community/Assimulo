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
from scipy.optimize import fminbound
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
            self.norm_of_res = N.linalg.norm(self.x0)
            self.func = self.problem.f
        except ProblemAlg_Exception:
            raise KINSOL_Exception("Problem has not implemented method 'f'")
        except IndexError:
            raise KINSOL_Exception("Problem has mismatching dimensions of f and initial guess")
        
        # check for constraints and test them
        broken_constraints = []
        if hasattr(self.problem, 'get_constraints'):
            self.constraints = self.problem.get_constraints()
            
            if self.constraints != None:
                # test if constraints are of correct type
                if type(self.constraints).__name__ != 'ndarray':
                    raise KINSOL_Exception("Constraints must be of type numpy.ndarray")
                
                if len(self.constraints) != len(self.x0):
                    raise KINSOL_Exception("Constraints must have same length as x0")
                # Test if initial guess x0 is consistant with constraints
                for c,xi,i in zip(self.constraints,self.x0,N.arange(0,self.x0.__len__())):
                    if re.search('float',type(c).__name__) == None:
                        print "Type problem with: ", c, type(c).__name__
                        raise KINSOL_Exception("Constraints must contain floats.")
                    if abs(c) > 2:
                        raise KINSOL_Exception("Entries in constraint vector must be between -2 and 2, see documentation.")
                    
                    if c != 0.0:
                        if c == 1.0:
                            if xi < 0.0:
                                broken_constraints.append(i)
                        elif c == 2.0:
                            if xi <= 0.0:
                                broken_constraints.append(i)
                        elif c == -1.0:
                            if xi > 0.0:
                                broken_constraints.append(i)
                        elif c == -2.0:
                            if xi >= 0.0:
                                broken_constraints.append(i)
                        else:
                            raise KINSOL_Exception("Constraint vector contains illegal elements.")
                
        else:
            self.constraints = None
            
        if broken_constraints != []:
            print "Variables breaking initial constraint: "
            for i in broken_constraints:
                self.problem.print_var_info(i)

            raise KINSOL_Exception("Initial guess does not fulfill applied constraints.")
        
        
        if hasattr(self.problem, 'check_constraints'):
            self.check_with_model = True
        else:
            self.check_with_model = False 
            
        self._use_jac = True
        self.reg_count = 0
        self.lin_count = 0
        self.verbosity = 0
        self.max_reg = 2.0
        
                
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
        
    def set_verbosity(self,verbosity):
        """
        Method used to set the verbosity of KINSOL
        
        Parameters::
            
            verbosity --
                Integer set to one of 0,1,2,3 where 0 is no info and 3 the
                most info. For more information please see the documentation for KINSOL.
        """
        type_name = type(verbosity).__name__
        if re.search('int',type_name) != None:
            
            # It is an integer, tes bounds
            if verbosity < 4 and verbosity > -1:
                self.verbosity = verbosity
            else:
                raise KINSOL_Exception("The variable sent to 'set_verbosity' must be either 0, 1, 2 or 3.")
        else:
            raise KINSOL_Exception("The variable sent to 'set_verbosity' must be an integer.")
            
        
    
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
        while not solved and self.reg_count < 2:
            try:
                self.solver.KINSOL_init(self.func,self.x0,self.dim,jac,self.constraints,self.verbosity,self.norm_of_res)
                res = self.solver.KINSOL_solve()
                solved = True
            except KINError as error:
                if error.value == 42:
                    # Try the heuristic
                    if hasattr(self.problem, 'get_heuristic_x0'):
                        print "----------------------------------------------------"
                        print "      Solver stuck with zero step-length."
                        print "----------------------------------------------------"
                        print "The following variables have start value zero"
                        print "and min set to zero causing the zero step-lenght."
                        print "These settings are either set by default or by user."
                        print ""

                        self.x0 = self.problem.get_heuristic_x0()
                        self.reg_count += 1
                        
                        print ""
                        print "This setting (start and min to zero) can often"
                        print "cause problem when initializing the system. "
                        print ""
                        print "To avoid this the above variables have"
                        print "their start attributes reset to one."
                        print ""
                        print "Trying to solve the system again..."
                    else:
                        raise KINSOL_Exception("Regularization failed due to constraints, tried getting heuristic initial guess but failed.")
                
                    """
                    Following functions commented out since they are moved to C code
                    elif error.value == -11 :
                        # the problem is caused by a singular jacobian try a regularized step
                        self.reg_count += 1
                        self._do_reg_step()
                    elif error.value == -6 or error.value == -8 :
                        self._brute_force()
                    
                    """ 
                else:
                    # Other error, send onward as exception
                    raise KINSOL_Exception(error.msg[error.value])
        
        if not solved:
            raise KINSOL_Exception("Singular Jacobian. Tried using Tikhonov regularization but stopped after ten steps.")
        """
        Functionality moved to C-code
        if self.reg_count != 0:
            print self.reg_count, " steps using the Tikhonov regularization performed."
        if self.lin_count != 0:
            print self.lin_count, " steps using the numpy.linalg.solve function performed."
        """
        if self.check_with_model:
            self.problem.check_constraints(res)
        return res
    

    
    def _do_reg_step(self):
        """
        Method used to perform a step using Tikhonov regularization
        """
        print "Trying to do a Tikhonov regularized step"
        if self._use_jac:
            if hasattr(self.problem,'jac'):
                # Extract data from problem
                x0 = self.x0
                fx = self.func(x0)
                J  = self.problem.jac(x0)
                try:
                    h = fminbound(self._norm_of_next_step,0.001,100.0)
                    print "Regularization parameter: ", h
                
                except:
                    h = 0.4
                    print "Could not find optimal regularization parameter. Using 0.4 instead."
                
                # Calculate regularisation step if the regularization prameter is not 'too big'
                if h < self.max_reg:
                    J_T = J.transpose()
                    rhs = N.dot(J_T,-fx)
                    A = N.dot(J_T,J)+(h**2)*N.eye(fx.__len__())
                    dx = solve(A,rhs)
                    
                    # Do step in problem
                    self.x0 = x0 + dx
                    self.problem._x0 = self.x0
                else:
                    print "Regularization parameter too big, trying a step using pseudo inverse."
                    self._do_pinv_step()
                
            else:
                raise KINSOL_Exception("Singular jacobian. Trying to do a step using Tikhonov regularization, but no jacobian supplied. Using the jacobian of KINSOL is not implemented yet.")
                
        else:
            raise KINSOL_Exception("Singular jacobian. Trying to do a step using Tikhonov regularization but 'use_jac' is set to false.")
    
    def _norm_of_next_step(self,h):
        """
        Function used to calculate the norm of the solution after the first step
        """
        x0 = self.x0
        J  = self.problem.jac(x0)
        fx = self.func(x0)
        
        # Calculate pseudo inverse and calculate new step
        J_T = J.transpose()
        rhs = N.dot(J_T,-fx)
        A = N.dot(J_T,J)+(h**2)*N.eye(fx.__len__())
        dx = solve(A,rhs)
        
        # Return norm of solution
        return N.linalg.norm(self.func(x0+dx))
    
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
                
        