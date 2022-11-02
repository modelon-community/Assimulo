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

cimport radau5ode # .pxd
cimport cython

import numpy as np
cimport numpy as np
import scipy.sparse as sps

from numpy cimport PyArray_DATA

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c_d(double* dest, object source, int dim):
    """
    Copy 1D numpy (double) array to (double *) C vector
    """
    if not (isinstance(source, np.ndarray) and source.flags.contiguous and source.dtype == np.float):
        source = np.ascontiguousarray(source, dtype=np.float)
    assert source.size >= dim, "The dimension of the vector is {} and not equal to the problem dimension {}. Please verify the output vectors from the min/max/nominal/evalute methods in the Problem class.".format(source.size, dim)
    memcpy(dest, <double*>PyArray_DATA(source), dim*sizeof(double))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c_d_matrix_flat_F(double* dest, object source, int nrow, int ncol):
    """
    Copy (square) 2D numpy array (order = c) to (double *) C matrix (with Fortran-style column major ordering)
    """
    cdef np.ndarray[double, ndim=2] source_np = np.array(source, copy=False, dtype = np.float)
    for i in range(ncol):
        for j in range(nrow):
            dest[j + i*nrow] = source_np[j][i]
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c2py_d(np.ndarray[double, ndim=1, mode='c'] dest, double* source, int dim):
    """
    Copy (double *) C vector to 1D numpy array
    """
    memcpy(PyArray_DATA(dest), source, dim*sizeof(double))
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c_i(int* dest, object source, int dim):
    """
    Copy 1D numpy (int) array to (int *) C vector
    """
    if not (isinstance(source, np.ndarray) and source.flags.contiguous and source.dtype == np.intc):
        source = np.ascontiguousarray(source, dtype=np.intc)
    assert source.size >= dim, "The dimension of the vector is {} and not equal to the problem dimension {}. Please verify the output vectors from the min/max/nominal/evalute methods in the Problem class.".format(source.size, dim)
    memcpy(dest, <int*>PyArray_DATA(source), dim*sizeof(int))

cdef int callback_fcn(int n, double* x, double* y_in, double* y_out, void* fcn_PY):
    """
    Internal callback function to enable call to Python based rhs function from C
    """
    cdef np.ndarray[double, ndim=1, mode="c"]y_py_in = np.empty(n, dtype = np.double)
    c2py_d(y_py_in, y_in, n)
    rhs, ret = (<object>fcn_PY)(x[0], y_py_in)

    py2c_d(y_out, rhs, len(rhs))

    # RADAU_CALLBACK_OK || RADAU_CALLBACK_ERROR_RECOVERABLE || RADAU_CALLBACK_ERROR_NONRECOVERABLE
    return ret[0] 

cdef int callback_jac(int n, double* x, double* y, double* fjac, void* jac_PY):
    """
    Internal callback function to enable call to Python based Jacobian function from C
    """
    cdef np.ndarray[double, ndim=1, mode="c"]y_py_in = np.empty(n, dtype = np.double)
    c2py_d(y_py_in, y, n)
    J, ret = (<object>jac_PY)(x[0], y_py_in)

    if ret[0]:
        return ret[0] # RADAU_CALLBACK_ERROR_RECOVERABLE || RADAU_CALLBACK_ERROR_NONRECOVERABLE

    py2c_d_matrix_flat_F(fjac, J, J.shape[0], J.shape[1])
    return RADAU_CALLBACK_OK

cdef int callback_solout(int* nrsol, double* xosol, double* xsol, double* y,
                         double* cont, double* werr, int* lrc, int* nsolu,
                         int* irtrn, void* solout_PY):
    """
    Internal callback function to enable call to Python based solution output function from C
    """
    cdef np.ndarray[double, ndim=1, mode="c"]y_py = np.empty(nsolu[0], dtype = np.double)
    cdef np.ndarray[double, ndim=1, mode="c"]cont_py = np.empty(4*nsolu[0], dtype = np.double)
    cdef np.ndarray[double, ndim=1, mode="c"]werr_py = np.empty(nsolu[0], dtype = np.double)
    c2py_d(y_py, y, nsolu[0])
    c2py_d(cont_py, cont, 4*nsolu[0])
    c2py_d(werr_py, werr, nsolu[0])

    irtrn[0] = (<object>solout_PY)(nrsol[0], xosol[0], xsol[0],
                                   y_py, cont_py, werr_py,
                                   lrc[0], irtrn[0])

    return irtrn[0]

cdef class RadauMemory:
    """Auxiliary data structure required to have C structs persists over multiple integrate calls."""
    cdef void* rmem

    cpdef int initialize(self, int n, int superLU, int nprocs, int nnz):
        """
        n = problem size
        superLU = 0 || 1, flag if using superLU
        nprocs = number of processors/threads in superLU
        nnz = number of non-zero elements with sparse LU
        """
        return radau5ode.setup_radau_mem(n, superLU, nprocs, nnz, &self.rmem)

    cpdef int set_nmax(self, int val):
        return radau5ode.radau_set_para_nmax(self.rmem, val)

    cpdef int set_nmax_newton(self, int val):
        return radau5ode.radau_set_para_nmax_newton(self.rmem, val)

    cpdef int set_step_size_safety(self, double val):
        return radau5ode.radau_set_para_step_size_safety(self.rmem, val)

    cpdef int reset_internal(self):
        """
        Reset internal memory and stats in Radau5
        """
        return radau5ode.reset_radau_internal_mem(self.rmem)

    cpdef list get_stats(self):
        """
        Return runtime stats logged in Radau5.
        """
        cdef int nfcn, njac, nsteps, naccpt, nreject, ludecomps, lusolves
        radau5ode.radau_get_stats(self.rmem, &nfcn, &njac, &nsteps, &naccpt, &nreject, &ludecomps, &lusolves)
        return [nfcn, njac, nsteps, naccpt, nreject, ludecomps, lusolves]

    cpdef void finalize(self):
        free_radau_mem(&self.rmem)

cdef int callback_jac_sparse(int n, double *x, double *y, int *nnz,
                             double *data, int *indices, int *indptr,
                             void* jac_PY):
    """Internal callback function to enable call to Python based evaluation of sparse (csc) jacobians."""
    cdef np.ndarray[double, ndim=1, mode="c"]y_py = np.empty(n, dtype = np.double)
    c2py_d(y_py, y, n)

    J, ret = (<object>jac_PY)(x[0], y_py)

    if ret[0]:
        return ret[0] # RADAU_CALLBACK_ERROR_RECOVERABLE || RADAU_CALLBACK_ERROR_NONRECOVERABLE

    if not isinstance(J, sps.csc.csc_matrix):
        return RADAU_CALLBACK_ERROR_INVALID_JAC_FORMAT

    if J.nnz > nnz[0]:
        return RADAU_CALLBACK_ERROR_INVALID_NNZ - J.nnz

    cdef np.ndarray[double, mode="c", ndim=1] jac_data_py = J.data.astype(np.double)
    cdef np.ndarray[int, mode="c", ndim=1] jac_indices_py = J.indices.astype(np.intc)
    cdef np.ndarray[int, mode="c", ndim=1] jac_indptr_py = J.indptr.astype(np.intc)

    ## copy data to output structures
    nnz[0] = J.nnz
    py2c_d(data, jac_data_py, nnz[0])
    py2c_i(indices, jac_indices_py, nnz[0])
    py2c_i(indptr, jac_indptr_py, n + 1)
    return RADAU_CALLBACK_OK

cpdef radau5(fcn_PY, double x, np.ndarray y,
             double xend, double h__, np.ndarray rtol, np.ndarray atol,
             jac_PY, int ijac, solout_PY,
             int iout, np.ndarray work,
             RadauMemory rad_memory):
    """
    Python interface for calling the C based Radau solver

        Parameters::

            fcn_PY
                        - Right-hand side function f(x, y), where 'x' is time, returning the evaluated value, 
            x
                        - Start time
            y
                        - Array, initial value
            xend
                        - End time
            h__
                        - Initial step-size guess
            rtol
                        - Array (len == len(y)) or scalar, Relative error tolerance in step-size control
            atol
                        - Array (len == len(y)) or scalar, Absolute error tolerance in step-size control
            jac_PY
                        - Jacobian function jac(x, y), where 'x' is time
            ijac
                        - Switch for Jacobian computation:
                          ijac == 0: C based finite differences
                          ijac == 1: Calls supplied 'jac_PY' function 
            solout_PY
                        - Callback function for logging solution of time-integration:
                          solout_PY(nrsol, told, t, y, cont, werr, lrc, irtrn)
                            - nrsol: number of solution point
                            - told:  Previous time-point
                            - t:     Current time-point
                            - y:     Solution at current time-point
                            - cont:  Output to be used to obtain high-order dense output, via the 'contr5' function
                            - werr:  Local error estimate
                            - lrc:   Unsused optional parameter
                            - irtrn: Optional parameter for interrupting time-integation if irtrn < 0
            iout
                        - Switch for using solout_PY:
                          iout == 0: solout_PY is never called
                          iout == 1: solout_PY is called after each successful time-integration step
            work
                        - Advanced tuning parameters of Radau solver, see radau5_c.c for details
            rad_memory
                        - instance of RadauMemory, needs to be initialized via RadauMemory.initialize
        Returns::
            
            x
                        - Final time for which a solution has been computed, x == xend, if succesful
            y
                        - Final solution
            idid
                        - Return flag, see radau_decsol.c for details (1 == Successful computation)
    """
    cdef int ret
    cdef int idid = 1
    
    cdef np.ndarray[double, mode="c", ndim=1] y_vec = y
    cdef np.ndarray[double, mode="c", ndim=1] rtol_vec = rtol
    cdef np.ndarray[double, mode="c", ndim=1] atol_vec = atol
    cdef np.ndarray[double, mode="c", ndim=1] work_vec = work

    ret = radau5ode.radau5_c(rad_memory.rmem, callback_fcn, <void*>fcn_PY, &x, &y_vec[0], &xend,
                        &h__, &rtol_vec[0], &atol_vec[0], callback_jac, callback_jac_sparse, <void*> jac_PY,
                        &ijac, callback_solout, <void*>solout_PY, &iout, &work_vec[0], &idid)

    return x, y, ret

cpdef contr5(int i, double x, np.ndarray cont):
    """
        Python interface for calling the C based interpolation function using dense output
        Returns 'i'-th component (index 0 based) at time 'x'.

            Parameters::

                i
                        - Which component to compute the solution for.
                x
                        - time-point at which the solution is requested. Needs to be within the time-interval defined by the last successful step.
                cont
                        - 'cont' input parameter to 'solout_PY' callback in 'radau5' function
            
            Returns::
                        - See function description

    """
    cdef np.ndarray[double, mode="c", ndim=1] cont_vec = cont
    return radau5ode.contr5_c(i, x, &cont_vec[0])
