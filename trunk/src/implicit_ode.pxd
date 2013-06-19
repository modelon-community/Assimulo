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


from ode cimport ODE
import numpy as N
cimport numpy as N


cdef class Implicit_ODE(ODE):
    cpdef _simulate(self, double t0, double tfinal,N.ndarray output_list,int COMPLETE_STEP, int INTERPOLATE_OUTPUT,int TIME_EVENT)
    
