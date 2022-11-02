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

cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)

cdef extern from "radau5_c.h":
    ctypedef struct radau_mem_t

    int RADAU_CALLBACK_OK
    int RADAU_CALLBACK_ERROR_RECOVERABLE
    int RADAU_CALLBACK_ERROR_NONRECOVERABLE
    int RADAU_CALLBACK_ERROR_INVALID_JAC_FORMAT
    int RADAU_CALLBACK_ERROR_INVALID_NNZ

    ## FunctionPointer_CallBack
    ctypedef int (*FP_CB_f)(int, double*, double*, double*, void*)
    ctypedef int (*FP_CB_jac)(int, double*, double*, double*, void*)
    ctypedef int (*FP_CB_solout)(int*, double*, double*, double*,
                                double*, double*, int*, int*,
                                int*, void*)
    ctypedef int (*FP_CB_jac_sparse)(int, double*, double*, int*, double*, int*, int*, void*)

    int setup_radau_mem(int n, int sparseLU, int nprocs, int nnz, void **mem_out)

    int reset_radau_stats(void* mem)
    int radau_get_stats(void *radau_mem, int *nfcn, int *njac, int *nsteps, int *naccpt, int *nreject, int * ludecomps, int *lusolves)

    int radau_set_para_nmax(void *radau_mem, int val)
    int radau_set_para_nmax_newton(void *radau_mem, int val)

    int radau5_c(void*, FP_CB_f, void*,
			 double*, double*, double*, double*,
			 double*, double*,
			 FP_CB_jac, FP_CB_jac_sparse, void*, int*,
			 FP_CB_solout, void*, int*,
			 double*, int*)

    double contr5_c(int*, double*, double*, int*)

    void free_radau_mem(void **radau_mem)
