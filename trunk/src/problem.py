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
            a discontinuity has been found in the supplied event functions. The solver
            is the solver attribute while the event_info is a list of length 2 where
            the first element is a list containing information about state events and
            the second element is a boolean for indicating if there have been an time
            event. If there have not been a state event the first element is an empty
            list. The state event list contains a set of integers of values (-1,0,1),
            the values indicates which state event have triggered (determined from 
            state_event(...) ) and the value indicates to where the state event is 'headed'. 
            
                Parameters::
                
                    solver
                        The solver object.
                    
                    event_info
                        List containing event information.
                            event_info[0] = Information about state events. (List)
                            event_info[1] = Information about time events. (Boolean)
                            
                                Example:
                                    Occured time event:
                                        event_info[0] = []
                                        event_info[1] = True
                                    Occured state event:
                                        event_info[0] = [0,1,0]  #Equal the length of state_events(..)
                                        event_info[1] = False
                                    Occured state event and time event:
                                        event_info[0] = [0,0,-1] #Equal the length of state_events(..)
                                        event_info[1] = True
            
        def initiate(self)
            Used to supply a specialized initialization construct. Gets called when
            (solver).initiate() is called.
        
        def finalize(self, solver)
            Used for specifying finalization. Gets called in the simulate method,
            after the simulation have been preformed.
        
        def completed_step(self, solver)
            This method is intended to be used for events which need not to be hit with high accuracy.
            It is called after each successful step taken by the solver. An example of use is
            that of simulating a systems which need to change the coordinate set to avoid singularities.
            The state, where this coordinate change has to be preformed is normally not very critical. 
            
            This is different from the event functions (state and time) in the sense that
            this method is completely related to the numerics of the solver, while
            the event functions, couple the numerics of the solver with physics
            of the problem.
            
                Parameters::
                
                    solver
                        The solver object.
                
                Return::
                
                    True
                        Indicating that the solver is to be reinitiated.
                    False
                        Indicating that the solver is NOT to be reinitiated.
            
    Parameters (optional)::
    
        problem_name
            The name of the problem.
    """
    problem_name = '---'
    
    def initiate(self, solver):
        """
        Method for specializing initiation.
        """
        if solver.verbosity >= solver.NORMAL:
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
        if solver.verbosity >= solver.NORMAL:
            print 'No event handling defined.'
    
    def finalize(self,solver):
        """
        Method for specifying the finalization options when the simulation have
        finished.
        """
        pass

    
class Implicit_Problem(Problem):
    """
    Problem for our implicit integrators.
    
        Must define the problem function f(self, t, y, yd, sw=None, p=None)
        
        Mandatory option::
        
            def f(self, t, y, yd)        or 
            def f(self, t, y, yd, sw)    or   (sw = Switches in case of discontinuities)
            def f(self, t, y, yd, p)     or   (p = Parameters for which sensitivites are to be calculated)
            def f(self, t, y, yd, sw, p)      (Both)
                Defines the residual of the problem.
                
                Returns:
                    A numpy array of size len(y).
        
        Available (optional) options::
        
            def state_events(self ,t ,y ,yd, sw)
                Defines the event (root) functions.
                
                Returns:
                    A numpy array.
                
            def time_events(self, t, y, yd, sw)
                Defines the time events. This function should return
                the next time-point for a time event. At a time-event
                the usual method handle_event is called for the specific
                handling. If there are no more time events. This function
                should return None.
                
                Returns:
                    Float
                        The time-point for the next time-event.
                    None
                        No time-event.
                
            def jac(self, c, t, y, yd, sw)
                Defines the Jacobian, which should be of the form
                J = dF/dx + c*dF/dx'.
                
                Returns:
                    A numpy array of size len(y)*len(y).
                    
            def handle_result(self, solver, t, y, yd)
                See docstring of handle_result.
          
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
            p0
                Defines the starting value of the parameters.
    """
    
    def f(self, t, y, yd, sw=None):
        """
        The residual function for a DAE problem.
        """
        raise Problem_Exception('The residual function is not specified.')
        
    def handle_result(self, solver, t, y, yd):
        """
        Method for specifying how the result is to be handled. As default the
        data is stored in three vectors, solver.(t/y/yd). If this is changed the
        plotting functionality will not work. The method works differently depending
        of certain options.
        
        Cases:
        
            Radau5 and IDA (no events)
            
                    Gets called whenever the simulation has finished with a call for
                    each time point specified with the number of communication points or
                    for each internal time point (if number of communication points == 0)
        
            IDA (with events)
            
                    Gets called as above with the extension that it is called also at each
                    event.
            
            IDA (with completed_step method or store_cont property)
            
                    Gets called at each internal time point.
        
            Parameters::
            
                solver
                    The solver object.
                t
                    Time.
                y
                    State.
                yd
                    State derivatives.
        """
        solver.t.extend([t])
        solver.y.extend([y])
        solver.yd.extend([yd])
        
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
        
            def state_events(self ,t ,y, sw)
                Defines the event (root) functions.
                
                Returns:
                    A numpy array.
                    
            def time_events(self, t, y, sw)
                Defines the time events. This function should return
                the next time-point for a time event. At a time-event
                the usual method handle_event is called for the specific
                handling. If there are no more time events. This function
                should return None.
                
                Returns:
                    Float
                        The time-point for the next time-event.
                    None
                        No time-event.
                
            def jac(self, t, y, sw=None)
                Defines the jacobian. J=df/dx.
                
                Returns:
                    A numpy matrix of size len(y)*len(y).
                    
            def jacv(self, t, y, fy, v)
                Defines a Jacobian Vector product. df/dx*v.
                
                Returns:
                    A numpy vector of size len(y).
            
            def handle_result(self, solver, t, y)
                See docstring of handle_result.
        
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
        
    def handle_result(self, solver, t, y):
        """
        Method for specifying how the result is to be handled. As default the
        data is stored in three vectors, solver.(t/y). If this is changed the
        plotting functionality will not work. The method works differently depending
        of certain options.
        
        Cases:
        
            Radau5, RungeKutta4, RungeKutta34, Explicit_Euler, CVode (no events)
            
                    Gets called whenever the simulation has finished with a call for
                    each time point specified with the number of communication points or
                    for each internal time point (if number of communication points == 0)
        
            CVode (with events)
            
                    Gets called as above with the extension that it is called also at each
                    event.
            
            CVode (with method completed_step or property store_cont)
            
                    Gets called at each internal time point.
        
            Parameters::
            
                solver
                    The solver object.
                t
                    Time.
                y
                    State.
                yd
                    State derivatives.
        """
        solver.t.extend([t])
        solver.y.extend([y])
    
