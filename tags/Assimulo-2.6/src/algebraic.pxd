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

cdef class Algebraic:
    cdef public object problem
    cdef public dict options, solver_options, problem_info
    cdef public dict statistics
    
    cdef public N.ndarray y
    cdef public N.ndarray y0

    cdef _reset_solution_variables(self)
    
    cdef double elapsed_step_time
    cdef double clock_start
    
    cpdef log_message(self, message, int level)
    cpdef solve(self, object y=*)
    
    cpdef finalize(self)
    cpdef initialize(self)
