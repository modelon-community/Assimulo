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

import numpy as np
cimport numpy as np

cdef class cProblem:
    cdef public int _sensitivity_result
    
    cpdef initialize(self, solver)
    cpdef reset(self)
    cpdef handle_event(self, object solver, event_info)
    cpdef finalize(self,object solver)

cdef class cImplicit_Problem(cProblem):
    cpdef res_internal(self, np.ndarray[double, ndim=1] res, double t, np.ndarray[double, ndim=1] y, np.ndarray[double, ndim=1] yd)
    
cdef class cOverdetermined_Problem(cProblem):
    cpdef res_internal(self, np.ndarray[double, ndim=1] res, double t, np.ndarray[double, ndim=1] y, np.ndarray[double, ndim=1] yd)
    
cdef class cExplicit_Problem(cProblem):
    cpdef int rhs_internal(self, np.ndarray[double, ndim=1] yd, double t, np.ndarray[double, ndim=1] y)
    cpdef np.ndarray res(self, t, y, yd, sw=*)
        
cdef class cDelay_Explicit_Problem(cExplicit_Problem):
    pass
    
cdef class cSingPerturbed_Problem(cExplicit_Problem):
    pass

cdef class cAlgebraic_Problem:
    cpdef initialize(self, solver)
    cpdef finalize(self,object solver)
