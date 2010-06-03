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

class Problem_Exception(Exception):
    """
    A Problem Exception.
    """
    pass

class Problem(object):
    """
    Base problem for Assimulo.
    
    Available (optional) options::
    
        def reset(self)
            Overrides the reset function to a user defined. The restriction is that
            t0, and y0 (and yd0) will be the new starting values. If they are not changed
            in the reset() they are left to be the the starting values defined in
            the beginning.
            
        def handle_event(self, solver, event_info)
            Defines how to handle a discontinuity. This functions gets called when
            a discontinuity has been found in the supplied event functions.
            
        def initiate(self)
            Used to supply a specialized initialization construct. Gets called when
            (solver).initiate() is called.
        
        def completed_step(self, solver)
            Method for specifying certain options or handling after a step have been
            successful. Return -1 to break and reinitiate the solver. Return 0 if not.
            
    Parameters (optional)::
    
        problem_name
            The name of the problem.
    """
    problem_name = '---'
    
    def initiate(self, solver):
        """
        Method for specializing initiation.
        """
        if solver.verbosity == solver.NORMAL:
            print 'No initiation defined for the problem.'
    
    def reset(self):
        """
        Resets a problem to its default values.
        """
        pass
        
    
    def handle_event(self, solver, event_info):
        """
        Method that is called when an event has triggered.
        """
        if solver.verbosity == solver.NORMAL:
            print 'No event handling defined.'
    
    
class Implicit_Problem(Problem):
    """
    Problem for our implicit integrators.
    
        Must define the problem function f(self, t, y, yd, sw=None)
        
        Mandatory option::
        
            def f(self, t, y, yd) or f(self, t, y, yd, sw)
                Defines the residual of the problem.
                
                Returns:
                    A numpy array of size len(y).
        
        Available (optional) options::
        
            def event_fcn(self ,t ,y ,yd, sw)
                Defines the event (root) functions.
                
                Returns:
                    A numpy array.
                
            def jac(self, c, t, y, yd, sw)
                Defines the Jacobian, which should be of the form
                J = dF/dx + c*dF/dx'.
                
                Returns:
                    A numpy array of size len(y)*len(y).
          
        Parameters (optional)::
          
            t0
                Defines the starting time.
            y0
                Defines the starting values of y0.
            yd0
                Defines the starting values of yd0.
            switches0
                Defines the starting values of the switches.
                Should be a list of booleans.
            algvar
                Defines the differential and algebraic components of the problem.
                Should be a list of integers. For more information, see the
                property algvar in Implicit_ODE.IDA
    """
    
    def f(self, t, y, yd, sw=None):
        """
        The residual function for a DAE problem.
        """
        raise Problem_Exception('The residual function is not specified.')
        
class Explicit_Problem(Problem):
    """
    Problem for our explicit integrators.
 
        Must define the problem function f(self, t, y, sw=None)
        
        Mandatory option::
        
            def f(self, t, y) or f(self, t, y, sw)
                Defines the right-hand-side of the problem.
                
                Returns:
                    A numpy array of size len(y).
        
        Available (optional) options::
        
            def event_fcn(self ,t ,y, sw)
                Defines the event (root) functions.
                
                Returns:
                    A numpy array.
                
            def jac(self, t, y, sw=None)
                Defines the jacobian. J=df/dx.
                
                Returns:
                    A numpy matrix of size len(y)*len(y).
        
        Parameters (optional)::
        
            t0
                Defines the starting time
            y0
                Defines the starting values of y0
            switches0
                Defines the starting values of the switches. 
                Should be a list of booleans.
    """
    
    def f(self, t, y, sw=None):
        """
        The rhs (right-hand-side) for a ODE problem.
        """
        raise Problem_Exception('The right-hand-side is not specified.')

