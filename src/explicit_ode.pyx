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

# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

cimport cython

import itertools
import sys
import numpy as np
cimport numpy as np
from timeit import default_timer as timer

from assimulo.ode cimport ODE
from assimulo.explicit_ode cimport Explicit_ODE, f_event_locator

from assimulo.problem import Explicit_Problem, Delay_Explicit_Problem, SingPerturbed_Problem, cExplicit_Problem
from assimulo.exception import Explicit_ODE_Exception, TimeLimitExceeded, TerminateSimulation

include "constants.pxi" #Includes the constants (textual include)

realtype = float

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c_d(double* dest, object source, int dim):
    """Copy 1D numpy (double) array to (double *) C vector."""
    if not (isinstance(source, np.ndarray) and source.flags.contiguous and source.dtype == np.float64):
        source = np.ascontiguousarray(source, dtype=np.float64)
    assert source.size >= dim, "The dimension of the vector is {} and not equal to the problem dimension {}. Please verify the output vectors from the min/max/nominal/evalute methods in the Problem class.".format(source.size, dim)
    memcpy(dest, <double*>np.PyArray_DATA(source), dim*sizeof(double))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c2py_d(np.ndarray[double, ndim=1, mode='c'] dest, double* source, int dim):
    """Copy (double *) C vector to 1D numpy array."""
    memcpy(np.PyArray_DATA(dest), source, dim*sizeof(double))

cdef int callback_event(int n_y, int n_g, double t, double* y_in, double* g_out, void* f_event_EXT) noexcept:
    """Event indicator callback function to event_locator.c"""
    cdef np.ndarray[double, ndim=1, mode="c"]y_py = np.empty(n_y, dtype = np.double)
    c2py_d(y_py, y_in, n_y)
    ret, g_high = (<object>f_event_EXT)(t, y_py)
    if ret < 0: ## immediate return, no not try to copy g_high
        return ret
    py2c_d(g_out, g_high, n_g)
    return ret

cdef int callback_interp(int n, double t, double* y_out, void* f_interp_EXT) noexcept:
    """Interpolation callback function to event_locator.c"""
    y_interp = (<object>f_interp_EXT)(t)
    py2c_d(y_out, y_interp, n)
    return 0

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
        
        if isinstance(problem, (cExplicit_Problem, Delay_Explicit_Problem, SingPerturbed_Problem)):
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
        y0 = np.array(y0) if len(np.array(y0).shape)>0 else np.array([y0])
        
        if len(self.y) != len(y0):
            raise Explicit_ODE_Exception('y0 must be of the same length as the original problem.')
        
        #Set the new values as the current values
        self.t = float(t0)
        self.y = y0
        
        if sw0 is not None:
            self.sw = (np.array(sw0,dtype=bool) if len(np.array(sw0,dtype=bool).shape)>0 else np.array([sw0],dtype=bool)).tolist()
            
        #Clear logs
        self.clear_logs()

    cpdef _simulate(self, double t0, double tfinal, np.ndarray output_list, int REPORT_CONTINUOUSLY, int INTERPOLATE_OUTPUT,
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
        cdef double eps = np.finfo(float).eps*100 #Machine Epsilon
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
        
        self.display_progress_activated = 1 if self.display_progress else 0
        self.time_integration_start = timer()
        
        while (flag == ID_COMPLETE and tevent == tfinal) is False and (self.t-eps > tfinal) if backward else (self.t+eps < tfinal):

            #Time event function is specified
            if TIME_EVENT == 1:
                tret = self.problem.time_events(self.t, self.y, self.sw)
                tevent = tfinal if tret is None else (tret if tret < tfinal else tfinal)
            else:
                tevent = tfinal
            
            #Initialize the clock, enabling storing elapsed time for each step
            if REPORT_CONTINUOUSLY and self.options["clock_step"]:
                self.clock_start = timer()
            
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
        
    cpdef report_solution(self, double t, np.ndarray y, opts):
        '''Is called after each successful step in case the complete step
        option is active. Here possible interpolation is done and the result 
        handeled. Furthermore possible step events are checked.
        '''
        self.t = t
        self.y = y
        
        #Store the elapsed time for a single step
        if self.options["clock_step"]:
            self.elapsed_step_time = timer() - self.clock_start
            self.clock_start = timer()
        
        #Check elapsed timed
        if self.time_limit_activated:
            if self.time_limit-(timer()-self.time_integration_start) < 0.0:
                raise TimeLimitExceeded("The time limit was exceeded at integration time %.8E."%self.t)
                
        self.chattering_clear_counter += 1
        if self.chattering_clear_counter > 3:
            self.chattering_check = None
            self.chattering_ok_print = 1
                
        if self.display_progress_activated:
            if (timer() - self.time_integration_start) > self.display_counter*10:
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
        
    cpdef event_locator(self, double t_low, double t_high, np.ndarray y_high):
        '''Checks if an event occurs in [t_low, t_high], if that is the case event 
        localization is started. Event localization finds the earliest small interval 
        that contains a change in domain. The right endpoint of this interval is then 
        returned as the time to restart the integration at.
        '''
        cdef int n_g = self.problem_info["dimRoot"]
        cdef np.ndarray[double, mode="c", ndim=1] g_low_c  = np.array(self.g_old)
        cdef np.ndarray[double, mode="c", ndim=1] g_mid_c  = np.empty(n_g, dtype = np.double)
        cdef np.ndarray[double, mode="c", ndim=1] g_high_c = np.empty(n_g, dtype = np.double)
        cdef np.ndarray[double, mode="c", ndim=1] y_high_c = np.array(y_high)
        cdef int nstatefcns = 0
        cdef int ret = f_event_locator(len(y_high), n_g, 1.e-13, t_low, &t_high,
                                       &y_high_c[0], &g_low_c[0], &g_mid_c[0], &g_high_c[0],
                                       callback_event, <void*>self.event_func,
                                       callback_interp, <void*>self.interpolate,
                                       &nstatefcns)
        self.statistics["nstatefcns"] += nstatefcns

        if ret == ID_PY_EVENT:
            event_info = np.zeros(n_g, dtype = int)
            for i in range(n_g):
                if (g_low_c[i] > 0) != (g_high_c[i] > 0):
                    event_info[i] = 1 if g_high_c[i] > 0 else -1
            self.set_event_info(event_info)

        self.g_old = g_high_c[:]
        y_high = y_high_c[:]
        return ret, t_high, y_high
    
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
        import pylab as pl
        
        if len(self.t_sol) > 0:
            pl.xlabel('time')
            pl.ylabel('state')
            pl.title(self.problem.name)
            
            if not mask:
                pl.plot(self.t_sol, self.y_sol, **kwargs)
            else:
                if not isinstance(mask, list):
                    raise Explicit_ODE_Exception('Mask must be a list of integers')
                if not len(mask)==len(self.y_sol[-1]):
                    raise Explicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                                 'the number of variables.')
                for i in range(len(mask)):
                    if mask[i]:
                        pl.plot(self.t_sol, np.array(self.y_sol)[:,i],**kwargs)
            
            pl.show()
        else:
            self.log_message("No result for plotting found.",NORMAL)
