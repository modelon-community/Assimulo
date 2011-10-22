#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2011 Modelon AB
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

from ode import *
from problem import Explicit_Problem

from exception import *

realtype = N.float

class Explicit_ODE(ODE):
    """
    Baseclass for our explicit ODE integrators.
    """
    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' or 'cExplicit_Problem'
                              class.
        """
        ODE.__init__(self) #Sets general attributes
        
        if problem is None:
            raise Explicit_ODE_Exception('The problem needs to be a subclass of a Explicit_Problem')
        
        if isinstance(problem, Explicit_Problem):
            self.problem = problem
        else:
            raise Explicit_ODE_Exception('The problem needs to be a subclass of a Explicit_Problem.')
        
        if hasattr(problem, 'y0'):
            problem.y0 = N.array(problem.y0,dtype=realtype) if len(N.array(problem.y0,dtype=realtype).shape)>0 else N.array([problem.y0],dtype=realtype)
        else:
            raise Explicit_ODE_Exception('y0 must be specified. Either in the problem or in the initialization')
        
        if hasattr(problem, "p0"):
            problem.p0 = N.array(problem.p0,dtype=realtype) if len(N.array(problem.p0,dtype=realtype).shape)>0 else N.array([problem.p0],dtype=realtype)
        
        if hasattr(problem, "sw0"):
            problem.sw0 = N.array(problem.sw0,dtype=N.bool) if len(N.array(problem.sw0,dtype=N.bool).shape)>0 else N.array([problem.sw0],dtype=N.bool)
        
        if hasattr(problem, 't0'):
            problem.t0 = float(problem.t0)
        else:
            problem.t0 = 0.0
        
        #Check Problem for event functions
        if hasattr(self.problem, 'time_events'):
            self.internal_flags["time_events"] = True
        
        if hasattr(self.problem, 'state_events'):
            self.internal_flags["state_events"] = True
        
        if hasattr(self.problem, 'step_events'):
            self.internal_flags["step_events"] = True
        
        self.t_cur = problem.t0
        self.y_cur = problem.y0.copy()
        
        self.t = []
        self.y = []
            
    def reset(self):
        """
        
        Resets the problem. If the problem is defined with a reset method, its called
        and then the method re_init. The re_init method is called with the initial
        values set in the problem, problem.t0 and problem.y0.
        
        """
        self.problem.reset()
        
        #Resets the results
        self.t = []
        self.y = []
        
        self.re_init(self.problem.t0, self.problem.y0)
        
    def re_init(self,t0, y0):
        """
        Reinitiates the solver.
        
            Parameters::
                
                t0  
                    - The initial time.
                y0  
                    - The initial values for the states
                
        See information in the __init__ method.
        """
        if len(self.y_cur) != len(y0):
            raise Explicit_ODE_Exception('y0 must be of the same length as the original problem.')
        
        #Set the new values as the current values
        self.t_cur = float(t0)
        self.y_cur = N.array(y0) if len(N.array(y0).shape)>0 else N.array([y0])

    def __call__(self, t0, tfinal, output_list, ONE_STEP, INTERPOLATE_OUTPUT,
                 TIME_EVENT, STEP_EVENT):
        """
        INTERNAL FUNCTION, FOR SIMULATION USE METHOD SIMULATE.
        
        Calls the integrator to perform the simulation over the given time-interval.
        If a second call to simulate is performed, the simulation starts from the last
        given final time.
        
            Parameters::
            
                tfinal  
                        - Final time for the simulation
                
                        - Should be a float or integer greater than the initial time.
                        
                ncp     
                        - Default '0'. Number of communication points where the 
                          solution is returned. If '0', the integrator will return 
                          at its internal steps.
                          
                        - Should be an integer.
                          
                    Example:
                    
                        __call__(10.0, 100), 10.0 is the final time and 100 is the number
                                             communication points.
        """
        y0 = self.y_cur
        t_logg = t0
        
        #Logg the first point
        self.problem.handle_result(self,t0,y0)
        
        #Reinitiate the solver
        flag_initialize = True

        #Start flag
        flag = ID_OK

        while flag != ID_COMPLETE or tevent != tfinal:

            #Time event function is specified.
            if  TIME_EVENT:
                tevent = self.problem.time_events(self.t_cur, self.y_cur, self.switches)
                tevent = tfinal if tevent is None else (tevent if tevent < tfinal else tfinal)
            else:
                tevent = tfinal
            
            
            if ONE_STEP:
                #Run in One step mode
                [flag, t, y]         = self.one_step_mode(self.t_cur, self.y_cur, tevent, flag_initialize)
                self.t_cur, self.y_cur = t, y.copy()
                
                #Store data depending on situation
                if INTERPOLATE_OUTPUT:
                    while output_list[output_index] <= t:
                        self.problem.handle_result(self, output_list[output_index], self.interpolate(output_list[output_index]))
                        
                        #Last logging point
                        t_logg = output_list[output_index]
                        
                        output_index = output_index+1 
                else:
                    self.problem.handle_result(self,t,y)
                    
                    #Last logging point
                    t_logg = self.t_cur
                    
                if STEP_EVENT: #If the option completed step is set.
                    flag_initialize = self.problem.completed_step(self)
                else:
                    flag_initialize = False
            else:
                #Run in Normal mode
                [flag, tlist, ylist] = zip(*list(self.integrator(self.t_cur, self.y_cur, tevent, flag_initialize, output_list)))
                flag, self.t_cur, self.y_cur = flag[-1], tlist[-1], ylist[-1].copy()

                #Store data
                map(self.problem.handle_result,itertools.repeat(self,len(tlist)), tlist, ylist)
                
                #Last logging point
                t_logg = self.t_cur
            
            #Event handling
            if flag == ID_EVENT or (flag == ID_COMPLETE and tevent != tfinal): #Event have been detected
                
                #Get and store event information
                event_info = [[],flag == ID_COMPLETE]
                if flag == ID_EVENT:
                    event_info[0] = self.state_event_info()
                
                #Log the information
                self.logg_event(self.t_cur, event_info, NORMAL)
                self.logg_message("A discontinuity occured at t = %e."%self.t_cur,NORMAL)
                self.logg_message("Current Switches: " + str(self.switches), LOUD)
                self.logg_message('Event info: ' + str(event_info), LOUD) 
                
                #Print statistics
                self.print_statistics(LOUD)

                try:
                    self.problem.handle_event(self, event_info) #self corresponds to the solver
                except TerminateSimulation: #Terminating the simulation after indication from handle event
                    self.logg_message("Terminating simulation at t = %f after signal from handle_event."%self.t_cur, NORMAL)
                    break
                    
                flag_initialize = True
            
            #Logg after the event handling if there was a communication point there.
            if flag_initialize and t_logg == self.t_cur: 
                self.problem.handle_result(self, self.t_cur, self.y_cur)
        
    
    def plot(self, mask=None, **kwargs):
        """
        Plot the computed solution.
        
            Parameters::
            
                mask    
                        - Default 'None'. Used to determine which variables that is to be plotted.
                          Used as a list of integers, ones represents the variable that is to be
                          plotted and zeros that is not. 
                        
                        - Should be a list of integers.
                        
                            Example:
                                mask = [1,0] , plots the first variable.
                
                **kwargs
                        - See http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
                          for information about the available options for **kwargs.
        """
        P.xlabel('time')
        P.ylabel('state')
        P.title(self.problem.name)
        
        if not mask:
            P.plot(self.t, self.y, **kwargs)
        else:
            if not isinstance(mask, list):
                raise Explicit_ODE_Exception('Mask must be a list of integers')
            if not len(mask)==len(self.y[-1]):
                raise Explicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                             'the number of variables.')
            for i in range(len(mask)):
                if mask[i]:
                    P.plot(self.t, N.array(self.y)[:,i],**kwargs)
        
        
        P.show()
