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

cdef class ODE:
    cdef public dict options, solver_options, problem_info
    cdef public dict supports, 
    cdef public object statistics
    
    cdef public list event_data
    
    cdef public object problem
    cdef public object chattering_check
    
    cdef public double t, t0
    cdef public int display_counter
    cdef public int chattering_clear_counter
    cdef public int chattering_ok_print
    cdef public N.ndarray y,yd, p
    cdef public N.ndarray y0, yd0, p0, sw0
    cdef double elapsed_step_time, time_integration_start
    cdef int time_limit_activated
    cdef double clock_start
    cdef public object _event_info
    
    #cdef public list t,y,yd,p,sw_cur
    cdef public list t_sol, y_sol, yd_sol, p_sol, sw
        
    cpdef log_message(self, message, int level)
    cpdef log_event(self, double time, object event_info, int level)
    cpdef clear_logs(self)
    cpdef simulate(self, double tfinal, int ncp=*, object ncp_list=*)
    cpdef get_options(self)
    cpdef get_supports(self)
    cpdef get_statistics(self)
    cpdef get_event_data(self)
    cpdef print_event_data(self)
    cpdef finalize(self)
    cpdef initialize(self)
    cdef _reset_solution_variables(self)
    cpdef get_elapsed_step_time(self)
    cpdef _chattering_check(self, object event_info)
