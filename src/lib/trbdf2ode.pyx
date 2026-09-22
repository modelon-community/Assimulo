#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2026 Modelon AB
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

"""Cython wrapper of the TR-BDF2 stepper in trbdf2.c: rhs, Jacobian and solution-output
callbacks go from C through function pointers into Python callables."""

cimport cython
import numpy as np
cimport numpy as np

cimport trbdf2ode  # .pxd

from numpy cimport PyArray_DATA


cdef class _Callbacks:
    """The Python callables and the problem size, passed to C as the user pointer."""
    cdef object rhs, jac, solout
    cdef int n


cdef int cb_rhs(int n, double t, const double *y, double *ydot, void *user) except? -1:
    cdef _Callbacks cb = <_Callbacks>user
    cdef np.ndarray[double, ndim=1, mode="c"] y_py = np.empty(n, dtype=np.double)
    memcpy(PyArray_DATA(y_py), y, n * sizeof(double))
    rhs, ret = cb.rhs(t, y_py)
    if ret[0] == 0:
        rhs = np.ascontiguousarray(rhs, dtype=np.double)
        memcpy(ydot, <double*>PyArray_DATA(rhs), n * sizeof(double))
    return ret[0]


cdef int cb_jac(int n, double t, const double *y, double *jac, void *user) except? -1:
    cdef _Callbacks cb = <_Callbacks>user
    cdef np.ndarray[double, ndim=1, mode="c"] y_py = np.empty(n, dtype=np.double)
    cdef np.ndarray[double, ndim=2, mode="fortran"] J_f
    memcpy(PyArray_DATA(y_py), y, n * sizeof(double))
    J, ret = cb.jac(t, y_py)
    if ret[0] != 0:
        return ret[0]
    J_f = np.asfortranarray(J, dtype=np.double)          # column-major, as the C code expects
    memcpy(jac, <double*>PyArray_DATA(J_f), n * n * sizeof(double))
    return 0


cdef int cb_solout(int naccpt, double told, double t, const double *y, const double *werr, void *user) except? -1:
    cdef _Callbacks cb = <_Callbacks>user
    cdef np.ndarray[double, ndim=1, mode="c"] y_py = np.empty(cb.n, dtype=np.double)
    cdef np.ndarray[double, ndim=1, mode="c"] werr_py = np.empty(cb.n, dtype=np.double)
    memcpy(PyArray_DATA(y_py), y, cb.n * sizeof(double))
    memcpy(PyArray_DATA(werr_py), werr, cb.n * sizeof(double))
    return cb.solout(naccpt, told, t, y_py, werr_py)


cdef class TRBDF2Memory:
    """The C stepper's memory; persists over integrate calls."""
    cdef trbdf2_mem *mem
    cdef int n

    def __cinit__(self, int n):
        self.n = n
        self.mem = NULL
        if trbdf2_create(n, &self.mem) != 0:
            raise MemoryError("TR-BDF2: could not allocate the solver memory")

    def __dealloc__(self):
        if self.mem != NULL:
            trbdf2_free(&self.mem)

    cpdef int set_tolerances(self, np.ndarray rtol, np.ndarray atol):
        cdef np.ndarray[double, ndim=1, mode="c"] r = np.ascontiguousarray(rtol, dtype=np.double)
        cdef np.ndarray[double, ndim=1, mode="c"] a = np.ascontiguousarray(atol, dtype=np.double)
        return trbdf2_set_tolerances(self.mem, &r[0], &a[0])

    cpdef int set_hmax(self, double hmax):
        return trbdf2_set_hmax(self.mem, hmax)

    cpdef int set_max_steps(self, long nmax):
        return trbdf2_set_max_steps(self.mem, nmax)

    cpdef int set_newton(self, int max_iter, double tol):
        return trbdf2_set_newton(self.mem, max_iter, tol)

    cpdef int set_step_control(self, double safety, double fac_min, double fac_max, double lu_lo, double lu_hi):
        return trbdf2_set_step_control(self.mem, safety, fac_min, fac_max, lu_lo, lu_hi)

    cpdef int set_jac_policy(self, int max_steps_between_jac, int use_user_jac):
        return trbdf2_set_jac_policy(self.mem, max_steps_between_jac, use_user_jac)

    cpdef int set_failure_policy(self, double step_factor, int max_consecutive):
        return trbdf2_set_failure_policy(self.mem, step_factor, max_consecutive)

    cpdef int reinit(self):
        return trbdf2_reinit(self.mem)

    cpdef int interpolate(self, double t, np.ndarray out):
        cdef np.ndarray[double, ndim=1, mode="c"] o = out
        return trbdf2_interpolate(self.mem, t, &o[0])

    cpdef double get_h(self):
        return trbdf2_get_h(self.mem)

    cpdef dict get_stats(self):
        cdef trbdf2_stats s
        trbdf2_get_stats(self.mem, &s)
        return {"nsteps": s.nsteps, "naccpt": s.naccpt, "nreject": s.nreject, "nfcn": s.nfcn, "nfcnjac": s.nfcnjac,
                "njac": s.njac, "nlu": s.nlu, "nsolve": s.nsolve, "nnewton": s.nnewton, "nnfail": s.nnfail, "nrhsfail": s.nrhsfail,
                "nnfail_div": s.nnfail_div, "nnfail_iter": s.nnfail_iter}

    cpdef int reset_stats(self):
        return trbdf2_reset_stats(self.mem)

    cpdef str get_err_msg(self):
        return trbdf2_get_err_msg(self.mem).decode("utf-8")


cpdef tuple trbdf2_py_solve(TRBDF2Memory memory, rhs_py, jac_py, solout_py, double t, np.ndarray y, double tend, double h0):
    """
    Integrate from (t, y) to tend. Returns (flag, t, y): flag 0 = reached tend, 1 = solout asked
    to stop, < 0 = error (message from memory.get_err_msg()).

    rhs_py(t, y) -> (ydot, [ret]); jac_py(t, y) -> (J, [ret]) with J (n, n) dense, or None for
    internal finite differences; solout_py(naccpt, told, t, y, werr) -> int (0 continue, 1 stop,
    < 0 unrecoverable). ret: 0 ok, > 0 recoverable (the trial point is refused), < 0 unrecoverable.
    """
    cdef _Callbacks cb = _Callbacks()
    cdef np.ndarray[double, ndim=1, mode="c"] y_c = np.ascontiguousarray(y, dtype=np.double).copy()
    cdef double t_c = t
    cdef int flag
    cdef trbdf2_jac_fn jac_c = NULL
    cdef trbdf2_solout_fn solout_c = NULL
    cb.rhs = rhs_py
    cb.jac = jac_py
    cb.solout = solout_py
    cb.n = memory.n
    if jac_py is not None:
        jac_c = cb_jac
    if solout_py is not None:
        solout_c = cb_solout
    flag = trbdf2_solve(memory.mem, cb_rhs, jac_c, solout_c, <void*>cb, &t_c, &y_c[0], tend, h0)
    return flag, t_c, y_c
