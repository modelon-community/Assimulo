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

import numpy as N
import pylab as P

class ODE_Exception(Exception):
    """ An integrator exception. """
    pass

class ODE(object):
    """
    Base class for all our integrators.
    """
    
    def __init__(self):
        """
        Defines general starting attributes for a simulation
        problem.
        """
        
        self.verbosity = self.NORMAL #Output level
        self.maxsteps = 10000 #Max number of steps
        self.atol = 1.0e-6 #Absolute tolerance
        self.rtol = 1.0e-6 #Relative tolerance
        self.h = 0.01 #Stepsize used for methods with fixed stepsize
        #self.max_eIter = 50 #Default number of allowed event iterations
        
        #Internal values
        self._log_event_info = []
    
    def _set_max_steps(self, maxsteps):
        """
        Sets the maximum number of steps the integrator is allowed
        to take to finnish the simulation.
        """
        if not isinstance(maxsteps,int):
            raise ODE_Exception('The maximum number of steps must be an integer.')
        if maxsteps < 1:
            raise ODE_Exception('The maximum number of steps must be a positive integer.')
        self.__max_steps = maxsteps
    
    def _get_max_steps(self):
        """
        Returns the maximum number of steps used.
        """
        return self.__max_steps
    maxstepsdocstring='Value to restrict the number of steps the integrator is allowed to take'
    maxsteps = property(_get_max_steps, _set_max_steps, doc=maxstepsdocstring)
    
    #def _get_eps(self):
    #    """
    #    Returns the epsilon used in the event indicators.
    #    """
    #    return self.atol*0.001
    #epsdocstring='Value used for adjusting the event indicators'
    #eps = property(_get_eps, doc=epsdocstring)
    
    def _set_problem_name(self, name):
        """
        Sets the name describing the problem
        """
        self.__name = name
        
    def _get_problem_name(self):
        """
        Returns the problem name.
        """
        return self.__name
    problemnamedocstring = 'String to describe the problem'
    problemname = property(_get_problem_name, _set_problem_name, doc=problemnamedocstring)
    
    #def _set_max_eIteration(self, max_eIter):
    #    """
    #    Sets the maximum number of iterations allowed in the event iteration.
    #    """
    #    if not isinstance(max_eIter, int) or max_eIter < 0:
    #        raise ODE_Exception('max_eIter must be a positive integer.')
    #    self.__max_eIter = max_eIter
    #    
    #def _get_max_eIteration(self):
    #    """
    #    Returns max_eIter.
    #    """
    #    return self.__max_eIter
    #    
    #max_eIterdocstring='Maximum number of event iterations allowed.'
    #max_eIter = property(_get_max_eIteration, _set_max_eIteration, doc=max_eIterdocstring)

    @property
    def is_disc(self):
        """Method to test if we are at an event."""
        return False
        
    #def check_eIter(self,before,after):
    #    """
    #    Helper function for event_switch to determine if we have event
    #    iteration.
    #    
    #        Input: Values of the event indicator functions (event_fcn)
    #        before and after we have changed mode of operations.
    #    """
    #    
    #    eIter = [False]*len(before)
    #    
    #    for i in range(len(before)):
    #        if (before[i] < 0.0 and after[i] > 0.0) or (before[i] > 0.0 and after[i] < 0.0):
    #            eIter[i] = True
    #            
    #    return eIter
    
    #Verbosity levels
    QUIET = 0
    WHISPER = 1
    NORMAL = 2
    LOUD = 3
    SCREAM = 4
    VERBOSE_VALUES = [QUIET, WHISPER, NORMAL, LOUD, SCREAM]
    def _get_verbosity(self):
        """Return the verbosity of the integrator."""
        return self.__verbosity
        
    def _set_verbosity(self, verbosity):
        """
        Sets the verbosity of the integrator. The verbosity levels are used to
        determine the amount of output from the integrater the user wishes to receive.
        
        These are the options:
            QUIET = 0
            WHISPER = 1
            NORMAL = 2
            LOUD = 3
            SCREAM = 4
        """
        if verbosity not in self.VERBOSE_VALUES:
            raise ODE_Exception('Verbosity values must be within %s - %s'%(self.QUIET, self.SCREAM))
        self.__verbosity = verbosity
        
    verbositydocstring = 'Determine the output level from the integrator.'
    verbosity = property(_get_verbosity, _set_verbosity,doc=verbositydocstring)
    
    def initiate(self):
        """
        Initiates the problem (if defined in the problem class).
        """
        self._problem.initiate(self)

    def print_event_info(self):
        """
        Prints the event information.
        """
        for i in range(len(self._log_event_info)):
            print 'Time, t = %e'%self._log_event_info[i][0]
            print '  Event info, ', self._log_event_info[i][1]
        print 'Number of events: ', len(self._log_event_info)
