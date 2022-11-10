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

cdef extern from "radau5_impl.h":
    ## FunctionPointer_CallBack
    ctypedef int (*FP_CB_f)(int, double, double*, double*, void*)
    ctypedef int (*FP_CB_jac)(int, double, double*, double*, void*)
    ctypedef int (*FP_CB_solout)(int, double, double, double*, double*, int, void*)
    ctypedef int (*FP_CB_jac_sparse)(int, double, double*, int*, double*, int*, int*, void*)
    
    int RADAU_OK
    int RADAU_ERROR_CALLBACK_RECOVERABLE
    int RADAU_ERROR_CALLBACK_UNRECOVERABLE
    int RADAU_ERROR_CALLBACK_JAC_FORMAT
    int RADAU_ERROR_CALLBACK_INVALID_NNZ

cdef extern from "radau5_io.h":
    int radau_setup_mem(int n, int sparseLU, int nprocs, int nnz, void **mem_out)
    int radau_reinit(void *radau_mem)
    int radau_get_stats(void *radau_mem, int *nfcn, int *njac, int *nsteps, int *naccpt, int *nreject, int *ludecomps, int *lusolves)
    char *radau_get_err_msg(void *radau_mem)

    int radau_set_nmax              (void *radau_mem, int val)
    int radau_set_nmax_newton       (void *radau_mem, int val)

    int radau_set_step_size_safety  (void *radau_mem, double val)
    int radau_set_theta_jac_recomp  (void *radau_mem, double val)
    int radau_set_fnewt             (void *radau_mem, double val)
    int radau_set_quot1             (void *radau_mem, double val)
    int radau_set_quot2             (void *radau_mem, double val)
    int radau_set_hmax              (void *radau_mem, double val)
    int radau_set_fac_lower         (void *radau_mem, double val)
    int radau_set_fac_upper         (void *radau_mem, double val)

    void radau_free_mem(void **radau_mem)

cdef extern from "radau5_c.h":
    int radau5_solve(void*, FP_CB_f, void*,
			 double*, double*, double*, double*,
			 double*, double*,
			 FP_CB_jac, FP_CB_jac_sparse, void*, int,
			 FP_CB_solout, void*, int, int*)

    int radau_get_cont_output(void *radau_mem, double x, double *out)
