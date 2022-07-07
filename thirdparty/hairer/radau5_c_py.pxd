#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2021 Modelon AB
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

from numpy cimport int32_t

cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)

cdef extern from "radau_decsol_c.h":
    ctypedef int32_t integer
    ctypedef double doublereal
    
    ## FunctionPointer_CallBack
    ctypedef int (*FP_CB_f)(integer, doublereal*, doublereal*, doublereal*,
                            doublereal*, integer*, void*)
    ctypedef int (*FP_CB_jac)(integer, doublereal*, doublereal*, doublereal*,
                            integer*, doublereal*, integer*, void*)
    ctypedef int (*FP_CB_mas)(integer, doublereal*, integer*, doublereal*,
                            integer*, void*)
    ctypedef int (*FP_CB_solout)(integer*, doublereal*, doublereal*, doublereal*,
                                doublereal*, doublereal*, integer*, integer*,
                                doublereal*, integer*, integer*, void*)
    ctypedef int (*FP_CB_jac_sparse)(int, double*, double*, int*, double*, int*, int*, doublereal*, integer*, void*)

    int radau5_c(integer, FP_CB_f, void*, doublereal*, doublereal*,
                 doublereal*, doublereal*, doublereal*, doublereal*,
                 integer*, FP_CB_jac, FP_CB_jac_sparse, void*, integer*, integer*, integer*,
                 integer*, integer*, integer*, FP_CB_solout,
                 void*, integer*, doublereal*, integer*, integer*, integer*,
                 doublereal*, integer*, integer*, integer)

    doublereal contr5_c(integer*, doublereal*, doublereal*, integer*)
