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

import numpy as N
cimport numpy as N
import time
import pylab as P

import itertools
import multiprocessing

from exception import *
from problem import Explicit_Problem, Delay_Explicit_Problem, Implicit_Problem, SingPerturbed_Problem
from support import Statistics

include "constants.pxi" #Includes the constants (textual include)

realtype = N.float 

cdef class ODE:
    """
    Base class for all our integrators.
    """
    
    def __init__(self, problem):
        """
        Defines general starting attributes for a simulation
        problem.
        """
        self.statistics = Statistics() #Initialize the statistics dictionary
        self.options = {"report_continuously":False,
                        "display_progress":True,
                        "verbosity":NORMAL,
                        "backward":False, 
                        "store_event_points":True, 
                        "time_limit":0, 
                        "clock_step":False, 
                        "num_threads":1} #multiprocessing.cpu_count()
        #self.internal_flags = {"state_events":False,"step_events":False,"time_events":False} #Flags for checking the problem (Does the problem have state events?)
        self.supports = {"state_events":False,"interpolated_output":False,"report_continuously":False,"sensitivity_calculations":False,"interpolated_sensitivity_output":False} #Flags for determining what the solver supports
        self.problem_info = {"dim":0,"dimRoot":0,"dimSens":0,"state_events":False,"step_events":False,"time_events":False
                             ,"jac_fcn":False, "sens_fcn":False, "jacv_fcn":False,"switches":False,"type":0,"jaclag_fcn":False,'prec_solve':False,'prec_setup':False
                             ,"jac_fcn_nnz": -1}
        #Type of the problem
        #0 = Explicit
        #1 = Implicit
        
        #Data object for storing the event data
        self.event_data = []
        self._event_info = N.array([])
        
        if problem is None:
            raise ODE_Exception('The problem needs to be a subclass of a Problem.')
        
        #Check Problem for event functions
        if hasattr(problem, 'time_events'):
            self.problem_info["time_events"] = True
        
        if hasattr(problem, 'state_events'):
            self.problem_info["state_events"] = True
        
        if hasattr(problem, 'step_events'):
            self.problem_info["step_events"] = True
        
        if hasattr(problem, 'y0'):
            self.y0 = N.array(problem.y0,dtype=realtype) if len(N.array(problem.y0,dtype=realtype).shape)>0 else N.array([problem.y0],dtype=realtype)
            self.problem_info["dim"] = len(self.y0)
        else:
            raise ODE_Exception('y0 must be specified in the problem.')
            
        if hasattr(problem, 'neq'):  # relevant for overdetermined problems: neq=number of equations >= dim
            self.problem_info["neq"] = problem.neq
        else:
            self.problem_info["neq"] = self.problem_info["dim"]
        
        if hasattr(problem, "p0"):
            self.p0 = N.array(problem.p0,dtype=realtype) if len(N.array(problem.p0,dtype=realtype).shape)>0 else N.array([problem.p0],dtype=realtype)
            self.problem_info["dimSens"] = len(self.p0)
            self.p = self.p0.copy()
        
        if hasattr(problem, "sw0"):
            self.sw0 = N.array(problem.sw0,dtype=N.bool) if len(N.array(problem.sw0,dtype=N.bool).shape)>0 else N.array([problem.sw0],dtype=N.bool)
            self.problem_info["switches"] = True
            self.sw = self.sw0.tolist()
        
        if hasattr(problem, 't0'):
            self.t0 = float(problem.t0)
        else:
            self.t0 = 0.0
            
        if hasattr(problem, "jac"):
            self.problem_info["jac_fcn"] = True
        if hasattr(problem, "jac_nnz"):
            self.problem_info["jac_fcn_nnz"] = problem.jac_nnz
        if hasattr(problem, "jacv"):
            self.problem_info["jacv_fcn"] = True
        if hasattr(problem, "jaclag"):
            self.problem_info["jaclag_fcn"] = True
        if hasattr(problem, "prec_solve"):
            self.problem_info["prec_solve"] = True
        if hasattr(problem, "prec_setup"):
            self.problem_info["prec_setup"] = True
            
        #Reset solution variables
        self._reset_solution_variables()
        
        #Specify storing of sensitivity to 0
        problem._sensitivity_result = 0
        
        #Initialize timer
        self.elapsed_step_time = -1.0
        self.clock_start = -1.0
        self.display_counter = 1
        self.chattering_clear_counter = 0
        self.chattering_ok_print = 1
        
        #Add common statistics
        self.statistics.add_key("nsteps", "Number of steps")
        self.statistics.add_key("nfcns", "Number of function evaluations")
        self.statistics.add_key("njacs", "Number of Jacobian evaluations")
        self.statistics.add_key("njacvecs", "Number of Jacobian*vector evaluations")
        self.statistics.add_key("nfcnjacs", "Number of function eval. due to Jacobian eval.")
        self.statistics.add_key("nerrfails", "Number of error test failures")
        self.statistics.add_key("nlus", "Number of LU decompositions")
        self.statistics.add_key("nniters", "Number of nonlinear iterations")
        self.statistics.add_key("nnfails", "Number of nonlinear convergence failures")
        self.statistics.add_key("nstatefcns", "Number of state function evaluations")
        self.statistics.add_key("nstateevents", "Number of state events")
        self.statistics.add_key("ntimeevents", "Number of time events")
        self.statistics.add_key("nstepevents", "Number of step events")
        self.statistics.add_key("nprecs", "Number of pre-conditioner solves")
        self.statistics.add_key("nprecsetups", "Number of pre-conditioner setups")
        self.statistics.add_key("nsensfcns", "Number of sensitivity evaluations")
        self.statistics.add_key("nsensfcnfcns", "Number of function eval. due to sensitivity eval.")
        self.statistics.add_key("nsensniters", "Number of sensitivity nonlinear iterations")
        self.statistics.add_key("nsensnfails", "Number of sensitivity nonlinear convergence failures")
        self.statistics.add_key("nsenserrfails", "Number of sensitivity error test failures")

        
    def __call__(self, double tfinal, int ncp=0, list cpts=None):
        return self.simulate(tfinal, ncp, cpts)
        
    cdef _reset_solution_variables(self):
        """
        Resets solution variables.
        """
        self.t_sol = []
        self.y_sol = []
        self.yd_sol = []
        self.p_sol = [[] for i in range(self.problem_info["dimSens"])]
        
        
    cpdef simulate(self, double tfinal, int ncp=0, object ncp_list=None):
        """
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
                        
                ncp_list
                        - Default None. A list of time points where the solution
                          should be returned. Note, requires that ncp == 0.
                          
                    Example:
                    
                        simulate(10.0, 100), 10.0 is the final time and 100 is the number
                                             communication points.
                 
        """
        t0 = self.t
        
        #Reset solution variables
        self._reset_solution_variables()
        
        #Error checking
        try:
            tfinal = float(tfinal)
        except ValueError:
            raise AssimuloException('Final time must be an integer or float.')
            
        if (self.t > tfinal) and not self.options["backward"]:
            raise AssimuloException('Final time {} must be greater than start time {}.\n Perhaps you should consider to reset the integration.'.format(tfinal,self.t))
        
        if not isinstance(ncp, int):
            raise AssimuloException('Number of communication points must be an integer')
        
        if ncp < 0:
            ncp = 0
            self.log_message('Number of communication points must be a positive integer, setting ncp = 0.',WARNING)
        
        #Check solver support against current problem
        if self.problem_info["state_events"] and self.supports["state_events"] is False:
            self.log_message("The current solver does not support state events (root functions). Disabling and continues.", WHISPER)
            self.problem_info["state_events"] = False
        
        if self.problem_info["step_events"] and self.supports["report_continuously"] is False:
            self.log_message("The current solver does not support step events (report continuously). Disabling step events and continues.", WHISPER)
            self.problem_info["step_events"] = False
        
        if self.supports["report_continuously"] is False and self.options["report_continuously"]:
            self.log_message("The current solver does not support to report continuously. Setting report_continuously to False and continues.", WHISPER)
            self.options["report_continuously"] = False
        
        if (ncp != 0 or ncp_list is not None) and (self.options["report_continuously"] or self.problem_info["step_events"]) and self.supports["interpolated_output"] is False:
            self.log_message("The current solver does not support interpolated output. Setting ncp to 0 and ncp_list to None and continues.", WHISPER)
            ncp = 0
            ncp_list = None
            
        if (ncp != 0 or ncp_list is not None) and self.problem_info["state_events"] and self.supports["report_continuously"] is False:
            self.log_message("The current solver does not support interpolated output together with state events. Setting ncp to 0 and ncp_list to None and continues.", WHISPER)
            ncp = 0
            ncp_list = None
        elif (ncp != 0 or ncp_list is not None) and self.problem_info["step_events"] and self.supports["report_continuously"]:
            if not self.report_continuously:
                 self.log_message("The problem contains step events: report_continuously is set to True", WHISPER)
            self.report_continuously = True
            
        #Determine the output list
        if ncp != 0:
            output_list = N.linspace(t0,tfinal,ncp+1)[1:]
            output_index = 0
        elif ncp_list is not None:
            output_list = N.array(ncp_list, dtype=realtype, ndmin=1)[N.logical_and(N.array(ncp_list, dtype=realtype, ndmin=1)>t0,N.array(ncp_list, dtype=realtype, ndmin=1)<=tfinal)]
            if output_list[-1] < tfinal: #Add the last point if necessary!
                output_list = N.append(output_list, tfinal)
            output_index = 0
        else:
            output_list = None
            output_index = 0
        
        #Determine if we are using one step mode or normal mode
        if self.problem_info['step_events'] or self.options['report_continuously']:
            REPORT_CONTINUOUSLY = 1
        else:
            REPORT_CONTINUOUSLY = 0
        
        #Determine if the output should be interpolated or not
        if output_list is None:
            INTERPOLATE_OUTPUT = 0
        else:
            INTERPOLATE_OUTPUT = 1

        #Time and Step events
        TIME_EVENT = 1 if self.problem_info['time_events'] is True else 0

        #Simulation starting, call initialize
        self.problem.initialize(self)
        self.initialize()
        
        #Start of simulation, start the clock
        time_start = time.clock()
        
        #Start the simulation
        self._simulate(t0, tfinal, output_list, REPORT_CONTINUOUSLY, INTERPOLATE_OUTPUT, TIME_EVENT)
        
        #End of simulation, stop the clock
        time_stop = time.clock()
        
        #Simulation complete, call finalize
        self.finalize()
        self.problem.finalize(self)
        
        #Print the simulation statistics
        self.print_statistics(NORMAL)
        
        #Log elapsed time
        self.log_message('Simulation interval    : ' + str(t0) + ' - ' + str(self.t) + ' seconds.', NORMAL)
        self.log_message('Elapsed simulation time: ' + str(time_stop-time_start) + ' seconds.', NORMAL)
        
        #Return the results
        if isinstance(self.problem, Explicit_Problem) or isinstance(self.problem, Delay_Explicit_Problem) or isinstance(self.problem, SingPerturbed_Problem):
            return self.t_sol, N.array(self.y_sol)
        else:
            return self.t_sol, N.array(self.y_sol), N.array(self.yd_sol)
        
    def _simulate(self,t0, tfinal, output_list, REPORT_CONTINUOUSLY, INTERPOLATE_OUTPUT, TIME_EVENT):
         pass
         
    cpdef initialize(self):
        pass
    
    cpdef finalize(self):
        pass
    
    def _set_verbosity(self, verb):
        try:
            self.options["verbosity"] = int(verb)
        except:
            raise AssimuloException("Verbosity must be an integer.")
    
    def _get_verbosity(self):
        """
        This determines the level of the output. A smaller value
        means more output. The following values can be set:
            
            QUIET   = 50
            WHISPER = 40
            NORMAL  = 30
            LOUD    = 20
            SCREAM  = 10
        
            Parameters::
            
                verb  
                        - Default 30 (NORMAL)
                    
                        - Should be a integer.

        """
        return self.options["verbosity"]
    
    verbosity = property(_get_verbosity,_set_verbosity)
    
    def _set_time_limit(self, time_limit):
        if time_limit < 0:
            raise AssimuloException("The time limit must be positive or zero.")
        self.options["time_limit"] = time_limit
        
    def _get_time_limit(self):
        """
        This option can be used to limit the time of an integration. I.e
        to set an upper bound on the time allowed for the integration
        to be completed. The time limit is specified in seconds. For the
        limit to be checked, the option report_continuously must be True.
        
            Parameters::
            
                time_limit
                            - Default 0, i.e. NO limit.
        """
        return self.options["time_limit"]
    
    time_limit = property(_get_time_limit, _set_time_limit)
    
    def _set_display_progress(self, display_progress):
        self.options["display_progress"] = bool(display_progress)
    
    def _get_display_progress(self):
        """
        This option actives output during the integration in terms of
        that the current integration is periodically printed to the
        stdout. Note though that report_continuously needs to be 
        activated.
        
            Parameters::
            
                display_progress
                            - Default True
        """
        return self.options["display_progress"]
    
    display_progress = property(_get_display_progress, _set_display_progress)
    
    def _set_report_continuously(self, report_continuously):
        self.options["report_continuously"] = bool(report_continuously)
    
    def _get_report_continuously(self):
        """
        This options specifies if the solver should report the solution 
        continuously after steps.
        
            Parameters::
            
                report_continuously
                  
                        - Default False
                    
                        - Should be a boolean.

        """
        return self.options["report_continuously"]
    
    report_continuously = property(_get_report_continuously,_set_report_continuously)
    
    def _set_number_threads(self, num_threads):
        self.options["num_threads"] = int(num_threads)
    
    def _get_number_threads(self):
        """
        This options specifies the number of threads to be used for those
        solvers that supports it.
        
            Parameters::
            
                num_threads
                  
                        - Default is the number of cores
                    
                        - Should be a integer.

        """
        return self.options["num_threads"]
    
    num_threads = property(_get_number_threads,_set_number_threads)
    
    
    def _set_store_event_points(self, store_event_points):
        self.options["store_event_points"] = bool(store_event_points)
    
    def _get_store_event_points(self):
        """
        This options specifies if the solver should save additional points
        at the events, :math:`t_e^-, t_e^+`.
        
            Parameters::
            
                store_event_points
                  
                        - Default True
                    
                        - Should be a Boolean.

        """
        return self.options["store_event_points"]
    
    store_event_points = property(_get_store_event_points,_set_store_event_points)
    
    def _set_clock_step(self, clock_step):
        self.options["clock_step"] = clock_step
    
    def _get_clock_step(self):
        """
        Specifies if the elapsed time of an integrator step should be
        timed or not. Not that this is only possible if running in 
        report continuously mode.
        """
        return self.options["clock_step"]
        
    clock_step = property(_get_clock_step, _set_clock_step)
    
    def _set_backward(self, backward):
        self.options["backward"] = bool(backward)
        
    def _get_backward(self):
        """
        Specifies if the simulation is done in reverse time. (NOTE:
        experimental!)
        
            Parameters::
            
                backward
                
                    - Default False
                    - Boolean
        """
        return self.options["backward"]
        
    backward = property(_get_backward,_set_backward)
    
    cpdef log_message(self, message,int level):
        if level >= self.options["verbosity"]:
            print(message)
            
    cpdef log_event(self,double time,object event_info, int level):
        if level >= self.options["verbosity"]:
            self.event_data.append([time,event_info])
            
    cpdef clear_logs(self):
        """
        Clears the currently stored log messages.
        """
        self.event_data = []
            
    cpdef get_options(self):
        """
        Returns the current solver options.
        """
        return self.options
        
    cpdef get_supports(self):
        """
        Returns the functionality which the solver supports.
        """
        return self.supports
        
    cpdef get_statistics(self):
        """
        Returns the run-time statistics (if any).
        """
        return self.statistics
        
    cpdef get_event_data(self):
        """
        Returns the event information (if any). If there exists information
        about events, it will be returned as a list of tuples where the
        first value is the time of the occured event and the second is the
        event information.
        """
        return self.event_data
        
    cpdef print_event_data(self):
        """
        Prints the event information (if any).
        """
        cdef i = 0
        for i in self.event_data:
            print 'Time, t = %e'%i[0]
            print '  Event info, ', i[1]
        print 'Number of events: ', len(self.event_data)
        
    def print_statistics(self, verbose=NORMAL):
        """
        This method should print the statistics of the solver.
        """
        if verbose >= self.options["verbosity"]:
            self.log_message('Final Run Statistics: %s ' % self.problem.name,        verbose)
            self.statistics.print_stats()
    
    cpdef get_elapsed_step_time(self):
        """
        Returns the elapsed time of a step. I.e. how long a step took.
        Note that this is only possible if running in report_continuously
        mode and should only be used to get a general idea of how long
        a step really took. 
        
        Returns::
        
            Elapsed time (note -1.0 indicates that it was not used)
        """
        return self.elapsed_step_time
        
    def _compact_atol(self):
        """
        Reduces atol to a scalar if it is an ndarray and  all entries are the same.
        Used for print solver options in a more compact way
        """
        if isinstance(self.atol,N.ndarray) and (self.atol==self.atol[0]).all():
                return self.atol[0]
        else:
                return self.atol
    
    cpdef _chattering_check(self, object event_info):
        self.chattering_clear_counter = 0
        if event_info[0] is not None and event_info[0] != []:
            if self.chattering_check is None:
                self.chattering_check  = abs(N.array(event_info[0]))
            else:
                self.chattering_check += abs(N.array(event_info[0]))

                if max(self.chattering_check) > 5 and self.chattering_ok_print:
                    self.chattering_ok_print = 0
                    self.log_message("Warning: Possible chattering detected at t = %e in state event(s): "%self.t + 
                                     str(N.where(self.chattering_check == max(self.chattering_check))[0]), NORMAL)
