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

# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

cimport cython
import numpy as np
cimport numpy as np
import scipy.sparse as sps

from numpy cimport PyArray_DATA

from assimulo.explicit_ode cimport Explicit_ODE
from assimulo.lib.radau_core import Radau_Common, Radau_Exception
from assimulo.exception import (
    AssimuloRecoverableError,
    TimeLimitExceeded,
)

include "constants.pxi"

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c_d(double* dest, object source, int dim) noexcept:
    """Copy 1D numpy (double) array to (double *) C vector."""
    if not (isinstance(source, np.ndarray) and source.flags.contiguous and source.dtype == np.float64):
        source = np.ascontiguousarray(source, dtype=np.float64)
    assert source.size >= dim
    memcpy(dest, <double*>PyArray_DATA(source), dim*sizeof(double))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c_d_matrix_flat_F(double* dest, object source, int nrow, int ncol) noexcept:
    """Copy 2D numpy array (C order) to (double *) C matrix (Fortran column-major)."""
    cdef np.ndarray[double, ndim=2] source_np = np.array(source, copy=False, dtype=np.float64)
    cdef int i, j
    for i in range(ncol):
        for j in range(nrow):
            dest[j + i*nrow] = source_np[j][i]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c2py_d(np.ndarray[double, ndim=1, mode='c'] dest, double* source, int dim) noexcept:
    """Copy (double *) C vector to 1D numpy array."""
    memcpy(PyArray_DATA(dest), source, dim*sizeof(double))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c_i(int* dest, object source, int dim) noexcept:
    """Copy 1D numpy (int) array to (int *) C vector."""
    if not (isinstance(source, np.ndarray) and source.flags.contiguous and source.dtype == np.intc):
        source = np.ascontiguousarray(source, dtype=np.intc)
    assert source.size >= dim
    memcpy(dest, <int*>PyArray_DATA(source), dim*sizeof(int))


cdef class RadauMemory:
    """Auxiliary data structure to persist C solver state across integrate calls."""

    cpdef int initialize(self, int n, int superLU, int nprocs, int nnz):
        self.n = n
        return radau_setup_mem(n, superLU, nprocs, nnz, &self.rmem)

    cpdef int set_nmax(self, int val):
        return radau_set_nmax(self.rmem, val)

    cpdef int set_nmax_newton(self, int val):
        return radau_set_nmax_newton(self.rmem, val)

    cpdef int set_step_size_safety(self, double val):
        return radau_set_step_size_safety(self.rmem, val)

    cpdef int set_theta_jac(self, double val):
        return radau_set_theta_jac_recomp(self.rmem, val)

    cpdef int set_fnewt(self, double val):
        return radau_set_fnewt(self.rmem, val)

    cpdef int set_quot1(self, double val):
        return radau_set_quot1(self.rmem, val)

    cpdef int set_quot2(self, double val):
        return radau_set_quot2(self.rmem, val)

    cpdef int set_hmax(self, double val):
        return radau_set_hmax(self.rmem, val)

    cpdef int set_fac_lower(self, double val):
        return radau_set_fac_lower(self.rmem, val)

    cpdef int set_fac_upper(self, double val):
        return radau_set_fac_upper(self.rmem, val)

    cpdef str get_err_msg(self):
        cdef char* ret = radau_get_err_msg(self.rmem)
        return ret.decode('UTF-8')

    cpdef int reinit(self):
        return radau_reinit(self.rmem)

    cpdef int interpolate(self, double t, np.ndarray output_array):
        cdef np.ndarray[double, ndim=1, mode="c"] output_array_c = output_array
        return radau_get_cont_output(self.rmem, t, &output_array_c[0])

    cpdef list get_stats(self):
        cdef int nfcn = 0, njac = 0, nsteps = 0, naccpt = 0, nreject = 0, ludecomps = 0, lusolves = 0
        radau_get_stats(self.rmem, &nfcn, &njac, &nsteps, &naccpt, &nreject, &ludecomps, &lusolves)
        return [nfcn, njac, nsteps, naccpt, nreject, ludecomps, lusolves]

    cpdef void finalize(self):
        radau_free_mem(&self.rmem)

from assimulo.solvers._radau5_py import Radau5Error, Radau5DAE, _Radau5ODE, _Radau5DAE

cdef class Radau5ODE(Explicit_ODE, Radau_Common):
    """
    Radau IIA fifth-order three-stages with step-size control and
    continuous output. Based on the FORTRAN code RADAU5 by E.Hairer and
    G.Wanner, which can be found here:
    http://www.unige.ch/~hairer/software.html

    Details about the implementation (FORTRAN) can be found in the book::

        Solving Ordinary Differential Equations II,
        Stiff and Differential-Algebraic Problems

        Authors: E. Hairer and G. Wanner
        Springer-Verlag, ISBN: 3-540-60452-9
    """
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem)

        self.options["inith"]         = 0.01
        self.options["newt"]          = 7
        self.options["thet"]          = 1.e-3
        self.options["fnewt"]         = None
        self.options["quot1"]         = 1.0
        self.options["quot2"]         = 1.2
        self.options["fac1"]          = 0.2
        self.options["fac2"]          = 8.0
        self.options["maxh"]          = None
        self.options["safe"]          = 0.9
        self.options["atol"]          = 1.0e-6 * np.ones(self.problem_info["dim"])
        self.options["rtol"]          = 1.0e-6
        self.options["usejac"]        = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"]      = 100000
        self.options["linear_solver"] = "DENSE"

        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = True
        self.supports["state_events"]        = True

        self._leny   = len(self.y)
        self._type   = '(explicit)'
        self._werr   = np.zeros(self._leny)
        self._y_work = np.empty(self._leny, dtype=np.double)

    def _get_linear_solver(self):
        return self.options["linear_solver"]

    def _set_linear_solver(self, linear_solver):
        """
        Which type of linear solver to use, "DENSE" or "SPARSE".
        """
        try:
            linear_solver_upper = linear_solver.upper()
        except Exception:
            raise Radau_Exception(
                "'linear_solver' parameter needs to be the STRING 'DENSE' or 'SPARSE'. "
                "Set value: {}, type: {}".format(linear_solver, type(linear_solver))
            ) from None
        if linear_solver_upper not in ["DENSE", "SPARSE"]:
            raise Radau_Exception(
                "'linear_solver' parameter needs to be either 'DENSE' or 'SPARSE'. "
                "Set value: {}".format(linear_solver)
            ) from None
        self.options["linear_solver"] = linear_solver.upper()

    linear_solver = property(_get_linear_solver, _set_linear_solver)

    cpdef initialize(self):
        cdef int ret, sparseLU
        self.statistics.reset()

        if self.usejac and not hasattr(self.problem, "jac"):
            raise Radau_Exception(
                "Use of an analytical Jacobian is enabled, but problem does not contain a 'jac' function.")

        if self.options["linear_solver"] == "SPARSE":
            if not self.usejac:
                self.log_message(
                    "Switching to 'DENSE' linear solver since a Jacobian method has not been provided.", LOUD)
                self.linear_solver = "DENSE"

        if self.options["linear_solver"] == "SPARSE":
            if not isinstance(self.problem_info["jac_fcn_nnz"], int):
                raise Radau_Exception(
                    "Number of non-zero elements of sparse Jacobian must be an integer, received: {}.".format(
                        self.problem_info["jac_fcn_nnz"]))
            if self.problem_info["jac_fcn_nnz"] < 0:
                if self.problem_info["jac_fcn_nnz"] == -1:
                    raise Radau_Exception(
                        "Number of non-zero elements of sparse Jacobian must be non-negative. "
                        "Detected default value of '-1', has 'problem.jac_fcn_nnz' been set?")
                raise Radau_Exception(
                    "Number of non-zero elements of sparse Jacobian must be non-negative, given value = {}.".format(
                        self.problem_info["jac_fcn_nnz"]))
            if self.problem_info["jac_fcn_nnz"] > self.problem_info["dim"]**2 + self.problem_info["dim"]:
                raise Radau_Exception(
                    "Number of non-zero elements of sparse Jacobian infeasible, "
                    "must be smaller than the problem dimension squared.")

        self.rad_memory = RadauMemory()
        sparseLU = int(self.options["linear_solver"] == "SPARSE")
        ret = self.rad_memory.initialize(
            self.problem_info["dim"], sparseLU,
            self.options["num_threads"], self.problem_info["jac_fcn_nnz"])
        if ret == -3:
            self.finalize()
            raise Radau5Error(value=ret, err_msg="Radau5 solver has not been compiled with superLU enabled.")
        if ret < 0:
            self.finalize()
            raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())

        ret = self.rad_memory.set_nmax(self.maxsteps)
        if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())
        ret = self.rad_memory.set_nmax_newton(self.newt)
        if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())
        ret = self.rad_memory.set_step_size_safety(self.safe)
        if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())
        ret = self.rad_memory.set_theta_jac(self.thet)
        if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())

        if self.options["fnewt"] is not None:
            ret = self.rad_memory.set_fnewt(self.fnewt)
            if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())

        ret = self.rad_memory.set_quot1(self.quot1)
        if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())
        ret = self.rad_memory.set_quot2(self.quot2)
        if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())

        if self.options["maxh"] is not None:
            ret = self.rad_memory.set_hmax(self.maxh)
            if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())

        ret = self.rad_memory.set_fac_lower(self.fac1)
        if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())
        ret = self.rad_memory.set_fac_upper(self.fac2)
        if ret < 0: self.finalize(); raise Radau5Error(value=ret, err_msg=self.rad_memory.get_err_msg())

    def set_problem_data(self):
        if self.problem_info["state_events"]:
            def event_func(t, y):
                try:
                    res = self.problem.state_events(t, y, self.sw)
                except BaseException as E:
                    self._py_err = E
                    return -1, None
                return 0, res

            def f(t, y):
                ret = 0
                try:
                    rhs = self.problem.rhs(t, y, self.sw)
                    return rhs, [ret]
                except BaseException as E:
                    rhs = y.copy()
                    if isinstance(E, (np.linalg.LinAlgError, ZeroDivisionError, AssimuloRecoverableError)):
                        ret = 1
                    else:
                        self._py_err = E
                        ret = -1
                return rhs, [ret]

            self.pt_fcn   = f
            self.pt_root  = event_func
            self.event_func = event_func
            self._event_info = [0] * self.problem_info["dimRoot"]
            ret, self.g_old = event_func(self.t, self.y)
            self.g_old = np.array(self.g_old)
            if ret < 0:
                raise self._py_err
            self.statistics["nstatefcns"] += 1
        else:
            def f(t, y):
                ret = 0
                try:
                    rhs = self.problem.rhs(t, y)
                    return rhs, [ret]
                except BaseException as E:
                    rhs = y.copy()
                    if isinstance(E, (np.linalg.LinAlgError, ZeroDivisionError, AssimuloRecoverableError)):
                        ret = 1
                    else:
                        self._py_err = E
                        ret = -1
                return rhs, [ret]
            self.pt_fcn = f

        if self.usejac:
            self.pt_jac = self._jacobian

    cpdef interpolate(self, double time):
        y = np.empty(self._leny)
        self.rad_memory.interpolate(time, y)
        return y

    def get_weighted_local_errors(self):
        """Returns the vector of weighted estimated local errors at the current step."""
        return np.abs(self._werr)

    # --- _solout as cdef: called directly from callback_solout without Python overhead ---

    cdef int _solout(self, int nrsol, double told, double t,
                     double* y_ptr, double* werr_ptr, int n) except ? -1:
        cdef int ret = 0, flag = 0
        cdef int output_index
        cdef np.ndarray y
        try:
            c2py_d(self._y_work, y_ptr, n)
            c2py_d(self._werr, werr_ptr, n)
            y = self._y_work

            if self.problem_info["state_events"]:
                flag, t, y = self.event_locator(told, t, y)
                if flag == ID_PY_EVENT:
                    ret = 1
                if flag < 0:
                    ret = -1

            if self._opts["report_continuously"]:
                try:
                    initialize_flag = self.report_solution(t, y.copy(), self._opts)
                    if initialize_flag:
                        ret = 1
                except TimeLimitExceeded as e:
                    self._py_err = e
                    ret = -2
            else:
                if self._opts["output_list"] is None:
                    self._tlist.append(t)
                    self._ylist.append(y.copy())
                else:
                    output_list  = self._opts["output_list"]
                    output_index = self._opts["output_index"]
                    try:
                        while output_list[output_index] <= t:
                            self._tlist.append(output_list[output_index])
                            self._ylist.append(self.interpolate(output_list[output_index]))
                            output_index += 1
                    except IndexError:
                        pass
                    self._opts["output_index"] = output_index

                    if (self.problem_info["state_events"] and flag == ID_PY_EVENT
                            and len(self._tlist) > 0 and self._tlist[-1] != t):
                        self._tlist.append(t)
                        self._ylist.append(y)
        except TimeLimitExceeded as e:
            self._py_err = e
            ret = -2
        except BaseException as E:
            self._py_err = E
            ret = -1
        return ret

    def _jacobian(self, t, y):
        """Compute the Jacobian (analytical or return error flags)."""
        ret = 0
        try:
            jac = self.problem.jac(t, y)
            if isinstance(jac, sps.csc_matrix) and (self.options["linear_solver"] == "DENSE"):
                jac = jac.toarray()
        except BaseException as E:
            jac = np.eye(len(y))
            if isinstance(E, (np.linalg.LinAlgError, ZeroDivisionError, AssimuloRecoverableError)):
                ret = 1
            else:
                self._py_err = E
                ret = -1
        return jac, [ret]

    cpdef integrate(self, double t, np.ndarray[ndim=1, dtype=np.float64_t] y,
                    double tf, dict opts):
        cdef int IJAC       = 1 if self.usejac else 0
        cdef int IOUT       = 1
        cdef int solout_ret = 1   # receives last return value from callback_solout
        cdef int flag                 # C solver exit code: 0=complete, 1=event, <0=error
        cdef double h  = self.inith
        cdef double x  = t, xend = tf

        if self.usejac and not hasattr(self.problem, "jac"):
            raise Radau_Exception(
                "Use of an analytical Jacobian is enabled, but problem does not contain a 'jac' function.")

        if opts["initialize"]:
            self.set_problem_data()
            self._tlist = []
            self._ylist = []

        self._py_err = None
        self._opts   = opts
        self.rad_memory.reinit()

        cdef np.ndarray[double, mode="c", ndim=1] y_vec    = y.copy()
        cdef np.ndarray[double, mode="c", ndim=1] rtol_vec = self.rtol * np.ones(self.problem_info["dim"])
        cdef np.ndarray[double, mode="c", ndim=1] atol_vec = self.atol

        # flag = C solver exit code (0=complete, 1=interrupted by solout, <0=error)
        # solout_ret = last return value of callback_solout (stored in *solout_ret by C)
        flag = radau5_solve(
            self.rad_memory.rmem,
            callback_fcn,      <void*>self,
            &x, &y_vec[0], &xend, &h,
            &rtol_vec[0], &atol_vec[0],
            callback_jac, callback_jac_sparse, <void*>self,
            IJAC,
            callback_solout, <void*>self,
            IOUT, &solout_ret)

        t = x
        y[:] = y_vec

        # Retrieve statistics via typed cpdef calls (no Python dispatch)
        cdef list stats = self.rad_memory.get_stats()
        cdef int nfcns    = stats[0]
        cdef int njacs    = stats[1]
        cdef int nsteps   = stats[3]
        cdef int nerrfails = stats[4]
        cdef int nLU      = stats[5]
        self.statistics["nsteps"]    += nsteps
        self.statistics["nfcns"]     += nfcns
        self.statistics["njacs"]     += njacs
        self.statistics["nfcnjacs"]  += (njacs * self.problem_info["dim"] if not self.usejac else 0)
        self.statistics["nerrfails"] += nerrfails
        self.statistics["nlus"]      += nLU

        if (flag >= 0 and opts["output_list"] is not None
                and (not self._tlist or self._tlist[-1] != t)):
            self._tlist.append(t)
            self._ylist.append(y)

        if flag == 0:
            return ID_PY_COMPLETE, self._tlist, self._ylist
        elif flag == 1:
            return ID_PY_EVENT, self._tlist, self._ylist
        else:
            msg = self.rad_memory.get_err_msg()
            self.finalize()
            if isinstance(self._py_err, BaseException):
                raise self._py_err from None
            raise Radau5Error(value=flag, t=t, err_msg=msg) from None

    def state_event_info(self):
        return self._event_info

    def set_event_info(self, event_info):
        self._event_info = event_info

    def print_statistics(self, verbose=NORMAL):
        """Prints the run-time statistics for the problem."""
        Explicit_ODE.print_statistics(self, verbose)

        log_message_verbose = lambda msg: self.log_message(msg, verbose)
        log_message_verbose('\nSolver options:\n')
        log_message_verbose(' Solver                  : Radau5' + self._type)
        log_message_verbose(' Linear solver           : ' + str(self.options["linear_solver"]))
        log_message_verbose(' Tolerances (absolute)   : ' + str(self._compact_tol(self.options["atol"])))
        log_message_verbose(' Tolerances (relative)   : ' + str(self.options["rtol"]))
        log_message_verbose('')

    cpdef finalize(self):
        """De-allocate internal C solver memory."""
        self.rad_memory.finalize()

# ---------------------------------------------------------------------------
# Module-level C callbacks - defined AFTER Radau5ODE so the type is known.
# Each receives <void*>self (a Radau5ODE instance) and casts it back.
# This eliminates the Python object-call overhead on every integration step.
# ---------------------------------------------------------------------------

cdef int callback_fcn(int n, double x, double* y_in, double* y_out,
                      void* solver_ptr) except? -1:
    """RHS callback: convert arrays, call solver.pt_fcn (Python user code)."""
    cdef Radau5ODE solver = <Radau5ODE>solver_ptr
    cdef np.ndarray[double, ndim=1, mode="c"] y_py = np.empty(n, dtype=np.double)
    c2py_d(y_py, y_in, n)
    rhs, ret = solver.pt_fcn(x, y_py)
    py2c_d(y_out, rhs, len(rhs))
    return ret[0]

cdef int callback_jac(int n, double x, double* y, double* fjac,
                      void* solver_ptr) except? -1:
    """Dense Jacobian callback."""
    cdef Radau5ODE solver = <Radau5ODE>solver_ptr
    cdef np.ndarray[double, ndim=1, mode="c"] y_py = np.empty(n, dtype=np.double)
    c2py_d(y_py, y, n)
    J, ret = solver.pt_jac(x, y_py)
    if ret[0]:
        return ret[0]
    py2c_d_matrix_flat_F(fjac, J, J.shape[0], J.shape[1])
    return RADAU_OK

cdef int callback_jac_sparse(int n, double x, double* y, int* nnz,
                              double* data, int* indices, int* indptr,
                              void* solver_ptr) except? -1:
    """Sparse (CSC) Jacobian callback."""
    cdef Radau5ODE solver = <Radau5ODE>solver_ptr
    cdef np.ndarray[double, ndim=1, mode="c"] y_py = np.empty(n, dtype=np.double)
    c2py_d(y_py, y, n)
    J, ret = solver.pt_jac(x, y_py)
    if ret[0]:
        return ret[0]
    if not isinstance(J, sps.csc_matrix):
        return RADAU_ERROR_CALLBACK_JAC_FORMAT
    if J.nnz > nnz[0]:
        return RADAU_ERROR_CALLBACK_INVALID_NNZ - J.nnz
    cdef np.ndarray[double, mode="c", ndim=1] jac_data    = J.data.astype(np.double)
    cdef np.ndarray[int,    mode="c", ndim=1] jac_indices = J.indices.astype(np.intc)
    cdef np.ndarray[int,    mode="c", ndim=1] jac_indptr  = J.indptr.astype(np.intc)
    nnz[0] = J.nnz
    py2c_d(data,    jac_data,    nnz[0])
    py2c_i(indices, jac_indices, nnz[0])
    py2c_i(indptr,  jac_indptr,  n + 1)
    return RADAU_OK

cdef int callback_solout(int nrsol, double xosol, double* xsol,
                         double* y, double* werr, int n,
                         void* solver_ptr) except? -1:
    """Solution-output callback: direct cdef call into solver._solout - no Python overhead."""
    cdef Radau5ODE solver = <Radau5ODE>solver_ptr
    return solver._solout(nrsol, xosol, xsol[0], y, werr, n)
