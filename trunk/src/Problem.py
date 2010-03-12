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
    
    Avaliable (optional) options:
        def reset(self)
            Overrides the reset function to a user defined. The restriction is that
            t0, and y0 (and yd0) will be the new starting values. If they are not change
            in the reset() they are left to be the the starting values defined in
            the beginning.
        def handle_event(self, solver, event_info)
            Defines how to handle a discontinuity. This functions gets called when
            a discontinuity has been found in the supplied event functions.
        def init_mode(self, solver)
            Override this function to create your own initialization of the new mode
            of operation. Default assumes that switching does not break consistency.
        def initiate(self)
            Used to supply a specialized initialization construct. Gets called when
            (solver).initiate() is called.
    """
    Problem_Name = '---'
    
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
        
        Avaliable (optional) options:
            def event_fcn(self ,t ,y ,yd, sw)
                Defines the event (root) functions.
            def jac(self, c, t, y, yd, sw)
                Defines the jacobian, which should be of the form
                J = dF/dx + c*dF/dx'
          
          Parameters (optional):
            t0
                Defines the starting time
            y0
                Defines the starting values of y0
            yd0
                Defines the starting values of yd0
            switches0
                Defines the starting values of the switches
            algvar
                Defines the differential and algebraic components of the problem (Ones and Zeros)
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
        
        Avaliable (optional) options:
            def event_fcn(self ,t ,y, sw)
                Defines the event (root) functions.
        
        Parameters (optional):
            t0
                Defines the starting time
            y0
                Defines the starting values of y0
            switches0
                Defines the starting values of the switches
    """
    
    def f(self, t, y, sw=None):
        """
        The rhs (right-hand-side) for a ODE problem.
        """
        raise Problem_Exception('The right-hand-side is not specified.')

