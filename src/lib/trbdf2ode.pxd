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

cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, size_t n)

cdef extern from "trbdf2.h":
    ctypedef int (*trbdf2_rhs_fn)(int n, double t, const double *y, double *ydot, void *user) except? -1
    ctypedef int (*trbdf2_jac_fn)(int n, double t, const double *y, double *jac, void *user) except? -1
    ctypedef int (*trbdf2_solout_fn)(int naccpt, double told, double t, const double *y, const double *werr, void *user) except? -1

    ctypedef struct trbdf2_stats:
        long nsteps, naccpt, nreject, nfcn, nfcnjac, njac, nlu, nsolve, nnewton, nnfail, nrhsfail
        long nnfail_div, nnfail_iter
    ctypedef struct trbdf2_mem:
        pass

    int TRBDF2_OK
    int TRBDF2_STOP
    int TRBDF2_ERROR_CALLBACK

    int    trbdf2_create(int n, trbdf2_mem **mem_out)
    void   trbdf2_free(trbdf2_mem **mem)
    int    trbdf2_set_tolerances(trbdf2_mem *mem, const double *rtol, const double *atol)
    int    trbdf2_set_hmax(trbdf2_mem *mem, double hmax)
    int    trbdf2_set_max_steps(trbdf2_mem *mem, long nmax)
    int    trbdf2_set_newton(trbdf2_mem *mem, int max_iter, double tol)
    int    trbdf2_set_step_control(trbdf2_mem *mem, double safety, double fac_min, double fac_max, double keep_lo, double keep_hi)
    int    trbdf2_set_jac_policy(trbdf2_mem *mem, int max_steps_between_jac, int use_user_jac)
    int    trbdf2_set_failure_policy(trbdf2_mem *mem, double step_factor, int max_consecutive)
    int    trbdf2_reinit(trbdf2_mem *mem)
    int    trbdf2_reset_stats(trbdf2_mem *mem)
    int    trbdf2_solve(trbdf2_mem *mem, trbdf2_rhs_fn rhs, trbdf2_jac_fn jac, trbdf2_solout_fn solout, void *user,
                        double *t, double *y, double tend, double h0)
    int    trbdf2_interpolate(trbdf2_mem *mem, double t, double *y_out)
    double trbdf2_get_h(trbdf2_mem *mem)
    void   trbdf2_get_stats(trbdf2_mem *mem, trbdf2_stats *stats)
    const char *trbdf2_get_err_msg(trbdf2_mem *mem)
