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

from exception import *
from problem import Algebraic_Problem

include "constants.pxi" #Includes the constants (textual include)

realtype = N.float

cdef class Algebraic:
    
    def __init__(self, problem):
        """
        Defines general starting attributes for a simulation
        problem.
        """
        self.statistics = {} #Initialize the statistics dictionary
        self.options = {"verbosity":NORMAL, "y_nominal":None, "y_min":None, "y_max":None}
        self.problem_info = {"dim":0,"jac_fcn":False,"jacv_fcn":False,'prec_solve':False,'prec_setup':False}
        
        if problem is None:
            raise Algebraic_Exception('The problem needs to be a subclass of a Problem.')
        self.problem = problem
        
        if hasattr(problem, 'y0'):
            self.y0 = N.array(problem.y0,dtype=realtype) if len(N.array(problem.y0,dtype=realtype).shape)>0 else N.array([problem.y0],dtype=realtype)
            self.y = N.array(problem.y0,dtype=realtype) if len(N.array(problem.y0,dtype=realtype).shape)>0 else N.array([problem.y0],dtype=realtype)
            self.problem_info["dim"] = len(self.y0)
        else:
            raise Algebraic_Exception('y0 must be specified in the problem.')
            
        if hasattr(problem, "jac"):
            self.problem_info["jac_fcn"] = True
        if hasattr(problem, "jacv"):
            self.problem_info["jacv_fcn"] = True
        if hasattr(problem, "prec_solve"):
            self.problem_info["prec_solve"] = True
        if hasattr(problem, "prec_setup"):
            self.problem_info["prec_setup"] = True
        
        if hasattr(problem, 'y0_nominal'):
            self.options["y_nominal"] = problem.y0_nominal
        if hasattr(problem, 'y0_min'):
            self.options["y_min"] = problem.y0_min
        if hasattr(problem, 'y0_max'):
            self.options["y_max"] = problem.y0_max
        
            
        #Reset solution variables
        self._reset_solution_variables()
        
        #Initialize timer
        self.elapsed_step_time = -1.0
        self.clock_start = -1.0

    cdef _reset_solution_variables(self):
        pass
    
    cpdef initialize(self):
        pass
    
    cpdef finalize(self):
        pass
    
    cpdef solve(self, y=None):
        """
        Calls the solver to solve the problem at hand.
        
            Parameters::
            
                y  
                        - The initial guess. If none, the stored value
                          of y will be used as initial guess.
        """
        #Reset solution variables
        self._reset_solution_variables()
         
        #Simulation starting, call initialize
        self.problem.initialize(self)
        self.initialize()
        
        #Start of simulation, start the clock
        time_start = time.clock()
        
        #Start the simulation
        self.y = self._solve(y)
        
        #End of simulation, stop the clock
        time_stop = time.clock()
        
        #Simulation complete, call finalize
        self.finalize()
        self.problem.finalize(self)
        
        #Print the simulation statistics
        self.print_statistics(NORMAL)
        
        #Log elapsed time
        self.log_message('Elapsed simulation time: ' + str(time_stop-time_start) + ' seconds.', NORMAL)
        
        #Return the results
        return self.y
    
    cpdef log_message(self, message,int level):
        if level >= self.options["verbosity"]:
            print(message)
    
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
    
