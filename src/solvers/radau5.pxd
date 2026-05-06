#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2010-2024 Modelon AB
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

from assimulo.explicit_ode cimport Explicit_ODE

# ---------------------------------------------------------------------------
# C declarations
# ---------------------------------------------------------------------------

cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)

cdef extern from "radau5_impl.h":
    ctypedef int (*FP_CB_f)(int, double, double*, double*, void*) except? -1
    ctypedef int (*FP_CB_jac)(int, double, double*, double*, void*) except? -1
    ctypedef int (*FP_CB_solout)(int, double, double*, double*, double*, int, void*) except? -1
    ctypedef int (*FP_CB_jac_sparse)(int, double, double*, int*, double*, int*, int*, void*) except? -1

    int RADAU_OK
    int RADAU_ERROR_CALLBACK_RECOVERABLE
    int RADAU_ERROR_CALLBACK_UNRECOVERABLE
    int RADAU_ERROR_CALLBACK_JAC_FORMAT
    int RADAU_ERROR_CALLBACK_INVALID_NNZ

cdef extern from "radau5_io.h":
    int radau_setup_mem(int n, int sparseLU, int nprocs, int nnz, void **mem_out)
    int radau_reinit(void *radau_mem)
    int radau_get_stats(void *radau_mem, int *nfcn, int *njac, int *nsteps,
                        int *naccpt, int *nreject, int *ludecomps, int *lusolves)
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

cdef extern from "radau5.h":
    int radau5_solve(void*, FP_CB_f, void*,
                     double*, double*, double*, double*,
                     double*, double*,
                     FP_CB_jac, FP_CB_jac_sparse, void*, int,
                     FP_CB_solout, void*, int, int*)
    int radau_get_cont_output(void *radau_mem, double x, double *out)

# ---------------------------------------------------------------------------
# Class declarations
# ---------------------------------------------------------------------------

cdef class RadauMemory:
    cdef void* rmem
    cdef int n

    cpdef int initialize(self, int n, int superLU, int nprocs, int nnz)
    cpdef int set_nmax(self, int val)
    cpdef int set_nmax_newton(self, int val)
    cpdef int set_step_size_safety(self, double val)
    cpdef int set_theta_jac(self, double val)
    cpdef int set_fnewt(self, double val)
    cpdef int set_quot1(self, double val)
    cpdef int set_quot2(self, double val)
    cpdef int set_hmax(self, double val)
    cpdef int set_fac_lower(self, double val)
    cpdef int set_fac_upper(self, double val)
    cpdef str get_err_msg(self)
    cpdef int reinit(self)
    cpdef int interpolate(self, double t, np.ndarray output_array)
    cpdef list get_stats(self)
    cpdef void finalize(self)


cdef class Radau5ODE(Explicit_ODE):
    # __dict__ needed because Radau_Common (plain Python base) has instance attributes
    cdef dict __dict__
    # Solver state - typed for O(1) C-slot access
    cdef RadauMemory      rad_memory
    cdef public int       _leny
    cdef public str       _type
    cdef public np.ndarray _werr      # pre-allocated weighted-error vector
    cdef np.ndarray        _y_work    # pre-allocated callback work array
    cdef list              _tlist
    cdef list              _ylist
    cdef dict              _opts
    # Python callable references
    cdef object            pt_fcn
    cdef object            pt_jac
    cdef object            pt_root
    cdef public object     event_func

    cpdef initialize(self)
    cpdef interpolate(self, double time)
    cdef int _solout(self, int nrsol, double told, double t,
                     double* y_ptr, double* werr_ptr, int n) except ? -1
    cpdef integrate(self, double t, np.ndarray[ndim=1, dtype=np.float64_t] y,
                    double tf, dict opts)
    cpdef finalize(self)
