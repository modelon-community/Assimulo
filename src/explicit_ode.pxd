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

from assimulo.ode cimport ODE

cdef class Explicit_ODE(ODE):

    cpdef _simulate(self, double t0, double tfinal,N.ndarray output_list,int COMPLETE_STEP, int INTERPOLATE_OUTPUT,int TIME_EVENT)
    cpdef report_solution(self, double t, N.ndarray y, opts)
    cpdef event_locator(self, double t_low, double t_high, N.ndarray y_high)

cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)

cdef extern from "ode_event_locator.h":
    ctypedef int (*FP_event_func)(int, int, double, double*, double*, void*);
    ctypedef int (*FP_interpolation)(int, double, double*, void*);

    int f_event_locator(int n_y, int n_g, double TOL, double t_low, double *t_high,
                        double *y_high, double *g_old, double *g_mid, double *g_high, 
                        FP_event_func f_event, void *f_event_EXT,
                        FP_interpolation f_interp, void *f_interp_EXT,
                        int *f_event_cnt);
