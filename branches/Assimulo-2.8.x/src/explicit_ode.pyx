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

from ode cimport ODE     
from problem import Explicit_Problem, Delay_Explicit_Problem, SingPerturbed_Problem

import pylab as P
import itertools
import sys
import numpy as N
cimport numpy as N

from exception import *
from time import clock, time

include "constants.pxi" #Includes the constants (textual include)

realtype = N.float

cdef class Explicit_ODE(ODE):
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
        ODE.__init__(self, problem) #Sets general attributes
        
        if isinstance(problem, Explicit_Problem) or isinstance(problem, Delay_Explicit_Problem) or isinstance(problem, SingPerturbed_Problem):
            self.problem = problem
        else:
            raise Explicit_ODE_Exception('The problem needs to be a subclass of a Explicit_Problem.')
        
        #Check the dimension of the state event function
        if self.problem_info["state_events"]:
            self.problem_info["dimRoot"] = len(problem.state_events(self.t0,self.y0, self.sw0))
        
        self.t = self.t0
        self.y = self.y0.copy()
            
    def reset(self):
        """
        Resets the problem. If the problem is defined with a reset method, its called
        and then the method re_init. The re_init method is called with the initial
        values provided to the solver (solver.t0 and solver.y0).
        
        """
        self.problem.reset()

        self.re_init(self.t0, self.y0, self.sw0 if self.problem_info["switches"] else None)
        
    def re_init(self,t0, y0, sw0=None):
        """
        Reinitiates the solver.
        
            Parameters::
                
                t0  
                    - The initial time.
                y0  
                    - The initial values for the states
                
        See information in the __init__ method.
        """
        y0 = N.array(y0) if len(N.array(y0).shape)>0 else N.array([y0])
        
        if len(self.y) != len(y0):
            raise Explicit_ODE_Exception('y0 must be of the same length as the original problem.')
        
        #Set the new values as the current values
        self.t = float(t0)
        self.y = y0
        
        if sw0 is not None:
            self.sw = (N.array(sw0,dtype=N.bool) if len(N.array(sw0,dtype=N.bool).shape)>0 else N.array([sw0],dtype=N.bool)).tolist()
            
        #Clear logs
        self.clear_logs()

    cpdef _simulate(self, double t0, double tfinal, N.ndarray output_list, int REPORT_CONTINUOUSLY, int INTERPOLATE_OUTPUT,
                 int TIME_EVENT):
        """
        INTERNAL FUNCTION, FOR SIMULATION USE METHOD SIMULATE.
        
        Calls the integrator to perform the simulation over the given time-interval.
        If a second call to simulate is performed, the simulation starts from the last
        given final time.
        
            Parameters::
            
            
                t0
                        - Initial time 
                        
                        - Should be float or integer less than final time.
            
                tfinal  
                        - Final time for the simulation
                
                        - Should be a float or integer greater than the initial time.
                        
                output_list
                          
                        - List of communication points
                        
                        - Should be an array.
                        
                REPORT_CONTINUOUSLY
                
                        - integer flag: 1 indicates that results are reported at every internal time step
                
                INTERPOLATE_OUTPUT
                
                        - integer flag: 1 indicates that results at output points should be obtained by interpolation
                                        0 indicates that results at output points should be computed exactly.
                                          This might slow down integration. 
                
                
                TIME_EVENT
                          
                        - integer flag: 1 indicates that time events were defined in problem else 0
                    
        """
        cdef double tevent
        cdef int flag, output_index
        cdef dict opts
        cdef double eps = N.finfo(float).eps*100 #Machine Epsilon
        cdef backward = 1 if self.backward else 0
        
        y0 = self.y

        #Log the first point
        self.problem.handle_result(self,t0,y0)

        #Reinitiate the solver
        flag_initialize = True

        #Start flag
        flag = ID_OK
        tevent = tfinal
        
        #Internal solver options
        opts = {}
        opts["initialize"] = flag_initialize
        opts["output_list"] = output_list
        opts["output_index"] = 0
        opts["report_continuously"] = 1 if REPORT_CONTINUOUSLY else 0
        output_index = 0
        
        self.time_limit_activated = 1 if self.time_limit > 0 else 0
        self.time_integration_start = time()
        
        while (flag == ID_COMPLETE and tevent == tfinal) is False and (self.t-eps > tfinal) if backward else (self.t+eps < tfinal):

            #Time event function is specified
            if TIME_EVENT == 1:
                tret = self.problem.time_events(self.t, self.y, self.sw)
                tevent = tfinal if tret is None else (tret if tret < tfinal else tfinal)
            else:
                tevent = tfinal
            
            #Initialize the clock, enabling storing elapsed time for each step
            if REPORT_CONTINUOUSLY and self.options["clock_step"]:
                self.clock_start = clock()
            
            flag, tlist, ylist = self.integrate(self.t, self.y, tevent, opts)
            
            #Store data if not done after each step
            if REPORT_CONTINUOUSLY is False and len(tlist) > 0:
                self.t, self.y = tlist[-1], ylist[-1].copy()
                list(map(self.problem.handle_result,itertools.repeat(self,len(tlist)), tlist, ylist))
            
            #Initialize flag to false
            flag_initialize = False
            
            #Event handling
            if flag == ID_EVENT or (flag == ID_COMPLETE and tevent != tfinal) or (flag == ID_COMPLETE and TIME_EVENT and tret==tevent): #Event has been detected
                
                if self.store_event_points and output_list is not None and abs(output_list[opts["output_index"]-1]-self.t) > eps:
                    self.problem.handle_result(self, self.t, self.y.copy())
                
                #Get and store event information
                event_info = [[],flag == ID_COMPLETE]
                if flag == ID_COMPLETE:
                    self.statistics["ntimeevents"] += 1#Time event detected
                if flag == ID_EVENT:
                    event_info[0] = self.state_event_info()
                    if REPORT_CONTINUOUSLY:
                        self._chattering_check(event_info)
                
                #Log the information
                if LOUD >= self.options["verbosity"]:
                    self.log_event(self.t, event_info, LOUD)
                    if SCREAM >= self.options["verbosity"]:
                        self.log_message("A discontinuity occured at t = %e."%self.t,SCREAM)
                        self.log_message("Current switches: " + str(self.sw), SCREAM)
                        self.log_message('Event info: ' + str(event_info), SCREAM) 
                
                    #Print statistics
                    self.print_statistics(LOUD)
                
                try:
                    self.problem.handle_event(self, event_info) #self corresponds to the solver
                except TerminateSimulation: #Terminating the simulation after indication from handle event
                    self.log_message("Terminating simulation at t = %f after signal from handle_event."%self.t, NORMAL)
                    break
                
                flag_initialize = True

            #Update options
            opts["initialize"] = flag_initialize
            
            #Logg after the event handling if there was a communication point there.
            if flag_initialize and (output_list is None or self.store_event_points):#output_list[opts["output_index"]] == self.t):
                self.problem.handle_result(self, self.t, self.y.copy())
                
            if self.t == tfinal: #Finished simulation (might occur due to event at the final time)
                break
        
    def report_solution(self, t, y, opts):
        '''Is called after each successful step in case the complete step
        option is active. Here possible interpolation is done and the result 
        handeled. Furthermore possible step events are checked.
        '''
        self.t, self.y = t, y.copy()
        
        #Store the elapsed time for a single step
        if self.options["clock_step"]:
            self.elapsed_step_time = clock() - self.clock_start
            self.clock_start = clock()
        
        #Check elapsed timed
        if self.time_limit_activated:
            if self.time_limit-(time()-self.time_integration_start) < 0.0:
                raise TimeLimitExceeded("The time limit was exceeded at integration time %.8E."%self.t)
                
        self.chattering_clear_counter += 1
        if self.chattering_clear_counter > 3:
            self.chattering_check = None
            self.chattering_ok_print = 1
                
        if self.display_progress:
            if (time() - self.time_integration_start) > self.display_counter*10:
                self.display_counter += 1
                
                sys.stdout.write(" Integrator time: %e" % self.t)
                sys.stdout.write('\r')
                sys.stdout.flush()
        
        #Store data depending on situation
        if opts["output_list"] is not None:
            output_list = opts["output_list"]
            output_index = opts["output_index"]
            try:
                while output_list[output_index] <= t:
                    self.problem.handle_result(self, output_list[output_index], 
                                    self.interpolate(output_list[output_index]))
                    output_index = output_index + 1
            except IndexError:
                pass
            opts["output_index"] = output_index
        else:
            self.problem.handle_result(self,t,y.copy())
        
        #Callback to the problem
        if self.problem_info["step_events"]:
            flag_initialize = self.problem.step_events(self) #completed step returned to the problem
            if flag_initialize:
                self.statistics["nstepevents"] += 1
        else:
            flag_initialize = False
            
        return flag_initialize
        
    def event_locator(self, t_low, t_high, y_high):
        '''Checks if an event occurs in [t_low, t_high], if that is the case event 
        localization is started. Event localization finds the earliest small interval 
        that contains a change in domain. The right endpoint of this interval is then 
        returned as the time to restart the integration at.
        '''
        
        g_high = self.event_func(t_high, y_high)
        g_low = self.g_old
        self.statistics["nstatefcns"] += 1
        n_g = self.problem_info["dimRoot"]
        TOL = max(abs(t_low), abs(t_high)) * 1e-13
        #Check for events in [t_low, t_high].
        for i in xrange(n_g):
            if (g_low[i] > 0) != (g_high[i] > 0):
                break
        else:
            self.g_old = g_high
            return (ID_PY_OK, t_high, y_high)
        
        side = 0
        sideprev = -1
        
        while abs(t_high - t_low) > TOL:
            #Adjust alpha if the same side is choosen more than once in a row.
            if (sideprev == side):
                if side == 2:
                    alpha = alpha * 2.0
                else:
                    alpha = alpha / 2.0
            #Otherwise alpha = 1 and the secant rule is used.
            else:
                alpha = 1
            
            #Decide which event function to iterate with.
            maxfrac = 0
            imax = 0 #Avoid compilation problem
            for i in xrange(n_g):
                if ((g_low[i] > 0) != (g_high[i] > 0)):
                    gfrac = abs(g_high[i]/(g_low[i] - g_high[i]))
                    if gfrac >= maxfrac:
                        maxfrac = gfrac
                        imax = i
            
            #Hack for solving the slow converging case when g is zero for a large part of [t_low, t_high].
            if g_high[imax] == 0 or g_low[imax] == 0:
                t_mid = (t_low + t_high)/2
            else:
                t_mid = t_high - (t_high - t_low)*g_high[imax]/ \
                                 (g_high[imax] - alpha*g_low[imax])
        
            #Check if t_mid is to close to current brackets and adjust inwards if so is the case.
            if (abs(t_mid - t_low) < TOL/2):
                fracint = abs(t_low - t_high)/TOL
                if fracint > 5:
                    delta = (t_high - t_low) / 10.0
                else:
                    delta = (t_high - t_low) / (2.0 * fracint)
                t_mid = t_low + delta
        
            if (abs(t_mid - t_high) < TOL/2):
                fracint = abs(t_low - t_high)/TOL
                if fracint > 5:
                    delta = (t_high - t_low) / 10.0
                else:
                    delta = (t_high - t_low) / (2.0 * fracint)
                t_mid = t_high - delta
            
            #Calculate g at t_mid and check for events in [t_low, t_mid].
            g_mid = self.event_func(t_mid, self.interpolate(t_mid))
            self.statistics["nstatefcns"] += 1
            sideprev = side
            for i in xrange(n_g):
                if (g_low[i] > 0) != (g_mid[i] > 0):
                    (t_high, g_high) = (t_mid, g_mid[0:n_g])
                    side = 1
                    break
            #If there are no events in [t_low, t_mid] there must be some event in [t_mid, t_high].
            else:
                (t_low, g_low) = (t_mid, g_mid[0:n_g])
                side = 2
        
        event_info = N.array([0] * n_g)
        for i in xrange(n_g):
            if (g_low[i] > 0) != (g_high[i] > 0):
                event_info[i] = 1 if g_high[i] > 0 else -1
                
        self.set_event_info(event_info)
        self.statistics["nstateevents"] += 1
        self.g_old = g_high
        return (ID_PY_EVENT, t_high, self.interpolate(t_high))
    
    def plot(self, mask=None, **kwargs):
        """
        Plot the computed solution.
        
            Parameters::
            
                mask    
                        - Default 'None'. Used to determine which solution components should be plotted.
                          Used as a list of integers, ones represent the components to be
                          plotted and zeros that is not. 
                        
                        - Should be a list of integers.
                        
                            Example:
                                mask = [1,0] , plots the first variable.
                
                **kwargs
                        - See http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
                          for information about the available options for **kwargs.
        """
        if len(self.t_sol) > 0:
            P.xlabel('time')
            P.ylabel('state')
            P.title(self.problem.name)
            
            if not mask:
                P.plot(self.t_sol, self.y_sol, **kwargs)
            else:
                if not isinstance(mask, list):
                    raise Explicit_ODE_Exception('Mask must be a list of integers')
                if not len(mask)==len(self.y_sol[-1]):
                    raise Explicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                                 'the number of variables.')
                for i in range(len(mask)):
                    if mask[i]:
                        P.plot(self.t_sol, N.array(self.y_sol)[:,i],**kwargs)
            
            
            P.show()
        else:
            self.log_message("No result for plotting found.",NORMAL)
