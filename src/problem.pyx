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

include "constants.pxi" #Includes the constants (textual include)

realtype = N.float

cdef class cProblem:
    
    name = '---'
    
    cpdef initialize(self, solver):
        """
        Method for specializing initiation.
        """
        solver.log_message("No initialization defined for the problem.", LOUD)
    
    cpdef reset(self):
        """
        Resets a problem to its default values.
        """
        pass
        
    cpdef handle_event(self, object solver, event_info):
        """
        Method that is called when an event has triggered.
        """
        solver.log_message("No event handling defined.", NORMAL)
    
    cpdef finalize(self,object solver):
        """
        Method for specifying the finalization options when the simulation have
        finished.
        """
        solver.log_message("No finalization defined for the problem.", LOUD)

cdef class cImplicit_Problem(cProblem):
    
    def __init__(self, object f=None, y0=None, yd0=None,double t0=0.0):
        
        self.f   = f
        self.y0  = None if y0 is None else (N.array(y0,dtype=realtype) if len(N.array(y0,dtype=realtype).shape)>0 else N.array([y0],dtype=realtype))
        self.yd0 = None if yd0 is None else (N.array(yd0,dtype=realtype) if len(N.array(yd0,dtype=realtype).shape)>0 else N.array([yd0],dtype=realtype))
        self.t0  = t0
        
        #Switches for discontinuities 
        pass
    
    cpdef handle_result(self, solver, t, N.ndarray[double, ndim=1] y, N.ndarray[double, ndim=1] yd):
        
        solver.t.extend([t])
        solver.y.extend([y])
        solver.yd.extend([yd])
        
    cpdef res_internal(self, N.ndarray[double, ndim=1] res, double t, N.ndarray[double, ndim=1] y, N.ndarray[double, ndim=1] yd):
        try:
            res[:] = self.f(t,y,yd)
        except:
            return ID_FAIL
        return ID_OK
    
cdef class cExplicit_Problem(cProblem):

    def __init__(self, object f=None, y0=None,double t0=0.0):
        
        self.f   = f
        self.y0  = None if y0 is None else (N.array(y0,dtype=realtype) if len(N.array(y0,dtype=realtype).shape)>0 else N.array([y0],dtype=realtype))
        self.t0  = t0
    
    cpdef handle_result(self, solver, double t, N.ndarray[double, ndim=1] y):
        solver.t.extend([t])
        solver.y.extend([y])
        
    cpdef int rhs_internal(self, N.ndarray[double, ndim=1] yd, double t, N.ndarray[double, ndim=1] y):
        try:
            yd[:] = self.f(t,y)
        except:
            return ID_FAIL
        return ID_OK

class Implicit_Problem(cImplicit_Problem):
    pass
    
class Explicit_Problem(cExplicit_Problem):
    pass
