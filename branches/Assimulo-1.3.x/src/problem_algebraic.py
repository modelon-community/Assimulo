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

class ProblemAlg_Exception(Exception):
    """
    A Problem Exception.
    """
    pass

class ProblemAlgebraic(object):
    """
    Base problem for KINSOL in assimulo.
    
    Must define the problem function f(self,x)
        
        Mandatory option::
        
            def f(self, input)
                Defines the right-hand-side of the problem.
                
                Returns:
                    A numpy array of size len(input).
                    
            def set_x0(self,x0)
                Sets an initial guess of the problem
                
                Returns:
                    -
                    
            def get_x0(self)
                Gets the initial guess of the problem
                
                Returns:
                    A numpy array of size len(self._x0)
        
        Available (optional) options::
        
            def jac(self,input)
                Defines the jacobian. J=df/dx.
                
                Returns:
                    A numpy matrix of size len(input)*len(input).
                    
            def get_constraints(self)
                Used to get the constraints of the problem
            
                Returns:
                    A numpy array, size len(_x0) containing
                    the constraints, if constraint[i] is:
                        0.0  - no constraint on x[i]
                        1.0  - x[i] greater or equal than 0.0
                        -1.0 - x[i] lesser or equal than 0.0
                        2.0  - x[i] greater  than 0.0
                        -2.0 - x[i] lesser than 0.0
                        
                    or
                    None if constraints should not be used
        
        Parameters ::
        
            _x0
                Defines the initial guess of x

    
    """

    
    def f(self, input):
        """
        The residual (right-hand-side) for a non linear problem.
        """
        raise ProblemAlg_Exception('The residual is not specified.')
    
    def set_x0(self,x0):
        """
        Sets the initial guess of the problem, sets it to x0
        """
        
        raise ProblemAlg_Exception('The routine used to set x0 is not specified.')

    def get_x0(self):
        """
        Gets the initial guess of the problem
        """
        
        raise ProblemAlg_Exception('The routine used to get x0 is not specified.')

    