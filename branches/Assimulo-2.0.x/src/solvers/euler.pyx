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

cimport numpy as N
import numpy as N

#from assimulo.ode import *
from assimulo.explicit_ode cimport Explicit_ODE

include "constants.pxi" #Includes the constants (textual include)

cdef class ExplicitEuler(Explicit_ODE):
    """
    Explicit Euler.
    """
    cdef N.ndarray yd1
    cdef object f
    
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["h"] = 0.01
        
        #Internal temporary result vector
        self.yd1 = N.array([0.0]*len(problem.y0))
        
        #RHS-Function
        self.f = problem.rhs_internal
        
        #Solver support
        self.supports["one_step_mode"] = True
    
    cpdef step(self,double t,N.ndarray y,double tf,dict opts):
        cdef double h
        h = self.options["h"]
        
        if t+h < tf:
            t, y = self._step(t,y,h)
            return ID_OK, t, y
        else:
            h = min(h, abs(tf-t))
            t, y = self._step(t,y,h)
            return ID_COMPLETE, t, y
    
    cpdef integrate(self, double t,N.ndarray y,double tf, dict opts):
        cdef double h
        cdef list tr,yr
        
        h = self.options["h"]
        h = min(h, abs(tf-t))
        
        tr = []
        yr = []
        
        while t+h < tf:
            t, y = self._step(t,y,h)
            tr.append(t)
            yr.append(y)
            h=min(h, abs(tf-t))
        else:
            t, y = self._step(t, y, h)
            tr.append(t)
            yr.append(y)
        
        return ID_COMPLETE, tr, yr
    
    cdef tuple _step(self,double t,N.ndarray y,double h):
        """
        This calculates the next step in the integration.
        """
        self.f(self.yd1,t,y) #The output is stored in yd
        return t + h, y + h*self.yd1
        
    def print_statistics(self, verbose=NORMAL):
        """
        Should print the statistics.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        self.log_message(' Step-length          : %s '%(self.options["h"]), verbose)
        self.log_message('\nSolver options:\n',                                    verbose)
        self.log_message(' Solver            : ExplicitEuler',                     verbose)
        self.log_message(' Solver type       : Fixed step\n',                      verbose)
