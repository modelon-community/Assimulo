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

import numpy as N
cimport numpy as N
import time
import pylab as P

import itertools

from exception import *

include "constants.pxi" #Includes the constants (textual include)
 
cdef class ODE:
    """
    Base class for all our integrators.
    """
    
    def __init__(self):
        """
        Defines general starting attributes for a simulation
        problem.
        """
        self.options = {"solver": "","type": "", "verbosity": NORMAL} #Options dict
        self.solver_options = {"continuous_output":False}
        self.internal_flags = {"state_events":False,"step_events":False,"time_events":False} #Flags for checking the problem (Does the problem have state events?)
        self.solver_support = {"state_events":False,"interpolated_output":False,"one_step_mode":False} #Flags for determining what the solver supports
        
        #Data object for storing the event data
        self.event_data = []
        
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
                          
                    Example:
                    
                        __call__(10.0, 100), 10.0 is the final time and 100 is the number
                                             communication points.
                 
        """
        t0 = self.t_cur
        
        #Error checking
        try:
            tfinal = float(tfinal)
        except ValueError:
            raise AssimuloException('Final time must be an integer or float.')
            
        if self.t_cur > tfinal:
            raise AssimuloException('Final time must be greater than start time.')
        
        if not isinstance(ncp, int):
            raise AssimuloException('Number of communication points must be an integer')
        
        if ncp < 0:
            ncp = 0
            self.logg_message('Number of communication points must be a positive integer, setting ncp = 0.',WARNING)
        
        #Check solver support against current problem
        if self.internal_flags["step_events"] and self.solver_support["one_step_mode"] is False:
            self.logg_message("The current solver does not support step events (completed steps). Disabling step events and continues.", WHISPER)
            self.internal_flags["step_events"] = False
        
        if self.solver_support["one_step_mode"] is False and self.solver_options["continuous_output"]:
            self.logg_message("The current solver does not support continuous output. Setting continuous_output to False and continues.", WHISPER)
            self.solver_options["continuous_output"] = False
        
        if (ncp != 0 or ncp_list != None) and (self.solver_options["continuous_output"] or self.internal_flags["step_events"]) and self.solver_support["interpolated_output"] is False:
            self.logg_message("The current solver does not support interpolated output. Setting ncp to 0 and ncp_list to None and continues.", WHISPER)
            ncp = 0
            ncp_list = None
            
        #Determine the output list
        if ncp != 0:
            output_list = N.linspace(t0,tfinal,ncp+1)
            output_index = 0
        elif ncp_list != None:
            output_list = N.array(ncp_list, dtype=realtype, ndmin=1)
            output_index = 0
        else:
            output_list = None
            output_index = 0
        
        #Determine if we are using one step mode or normal mode
        if self.internal_flags['step_events'] or self.solver_options['continuous_output']:
            ONE_STEP = True
        else:
            ONE_STEP = False
        
        #Determine if the output should be interpolated or not
        if output_list == None:
            INTERPOLATE_OUTPUT = False
        else:
            INTERPOLATE_OUTPUT = True

        #Time and Step events
        TIME_EVENT = self.internal_flags['time_events']
        STEP_EVENT = self.internal_flags["step_events"]

        #Simulation starting, call initialize
        self.problem.initialize(self)
        
        #Start of simulation, start the clock
        time_start = time.clock()
        
        #Start the simulation
        self.__call__(t0, tfinal, output_list, ONE_STEP, INTERPOLATE_OUTPUT, TIME_EVENT, STEP_EVENT)
        
        #End of simulation, stop the clock
        time_stop = time.clock()
        
        #Simulation complete, call finalize
        self.problem.finalize(self)
        
        #Print the simulation statistics
        self.print_statistics(NORMAL)
        
        #Log elapsed time
        self.logg_message('Simulation interval    : ' + str(t0) + ' - ' + str(self.t_cur) + ' seconds.', NORMAL)
        self.logg_message('Elapsed simulation time: ' + str(time_stop-time_start) + ' seconds.', NORMAL)

    cpdef logg_message(self, message,int level):
        if level >= self.options["verbosity"]:
            print message
            
    cpdef logg_event(self,double time,object event_info, int level):
        if level >= self.options["verbosity"]:
            self.event_data.append([time,event_info])
            
    cpdef get_options(self):
        """
        Returns the solver options.
        """
        return self.solver_options
