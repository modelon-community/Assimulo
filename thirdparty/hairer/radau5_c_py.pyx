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

cimport radau5_c_py
cimport cython

import numpy as np
cimport numpy as np
import scipy.sparse as sps

from numpy cimport PyArray_DATA

from assimulo.exception import AssimuloException

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
cdef void c2py_d_mat_F(np.ndarray[double, ndim=2, mode='fortran'] dest, double* source, int dim):
    """
    Copy (double *) C matrix (Fotran-style column major ordering) to 2D numpy array
    """
    memcpy(PyArray_DATA(dest), source, dim*sizeof(double))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c_z(doublecomplex* dest, object source, int dim):
    """
    Copy 1D numpy (double complex) array to (doublecomplex *) C vector
    """
    if not (isinstance(source, np.ndarray) and source.flags.contiguous and source.dtype == np.complex128):
        source = np.ascontiguousarray(source, dtype=np.float)
    assert source.size >= dim, "The dimension of the vector is {} and not equal to the problem dimension {}. Please verify the output vectors from the min/max/nominal/evalute methods in the Problem class.".format(source.size, dim)
    memcpy(dest, <doublecomplex*>PyArray_DATA(source), dim*sizeof(doublecomplex))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c2py_i(np.ndarray[int, ndim=1, mode='c'] dest, int* source, int dim):
    """
    Copy (int *) C vector to 1D numpy array
    """
    memcpy(PyArray_DATA(dest), source, dim*sizeof(int))

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

cdef int callback_fcn(integer n, doublereal* x, doublereal* y_in, doublereal* y_out,
                      doublereal* rpar, integer* ipar, void* fcn_PY):
    """
    Internal callback function to enable call to Python based rhs function from C
    """
    cdef np.ndarray[double, ndim=1, mode="c"]y_py_in = np.empty(n, dtype = np.double)
    c2py_d(y_py_in, y_in, n)
    res = (<object>fcn_PY)(x[0], y_py_in)
    py2c_d(y_out, res[0], len(res[0]))
    ipar[0] = res[1][0]
    return 0

cdef int callback_jac(integer n, doublereal* x, doublereal* y, doublereal* fjac,
                      integer* ldjac, doublereal* rpar, integer* ipar, void* jac_PY):
    """
    Internal callback function to enable call to Python based Jacobian function from C
    """
    cdef np.ndarray[double, ndim=1, mode="c"]y_py_in = np.empty(n, dtype = np.double)
    c2py_d(y_py_in, y, n)
    res = (<object>jac_PY)(x[0], y_py_in)
    if type(res) != np.ndarray:
        raise AssimuloException("Error: Encountered a sparse Jacobian while expecting a dense Jacobian, type of Jacobian: {}".format(type(res)))
    py2c_d_matrix_flat_F(fjac, res, res.shape[0], res.shape[1])
    return 0

cdef int callback_mas(integer n, doublereal* am, integer* lmas, doublereal* rpar,
                      integer* ipar, void* mas_PY):
    """
    Internal callback function to enable call to Python based mass matrix function from C
    """
    cdef np.ndarray[double, mode="fortran", ndim=2]am_py = np.empty((lmas[0], n), order = 'F', dtype = np.double)
    c2py_d_mat_F(am_py, am, n*lmas[0])
    res = (<object>mas_PY)(am_py)
    py2c_d_matrix_flat_F(am, res, res.shape[0], res.shape[1])
    return 0

cdef int callback_solout(integer* nrsol, doublereal* xosol, doublereal* xsol, doublereal* y,
                         doublereal* cont, doublereal* werr, integer* lrc, integer* nsolu,
                         doublereal* rpar, integer* ipar, integer* irtrn, void* solout_PY):
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

cdef int callback_jac_sparse(int n, double *x, double *y, int *nnz,
                             double * data, int *indices, int *indptr,
                             doublereal* rpar, integer* ipar,
                             void* jac_PY):
    ## TODO: Add docstring
    cdef np.ndarray[double, ndim=1, mode="c"]y_py = np.empty(n, dtype = np.double)
    c2py_d(y_py, y, n)

    J = (<object>jac_PY)(x[0], y_py)
    if not isinstance(J, sps.csc.csc_matrix):
        raise AssimuloException("The Jacobian must be provided as scipy.sparse.csc_matrix. Given type: {}".format(type(J)))

    if J.nnz > nnz[0]:
        raise AssimuloException("Mismatch of nnz in the jacobian, specified in problem: {}, jacobian evaluation: {}".format(nnz[0], J.nnz))

    cdef np.ndarray[double, mode="c", ndim=1] data_py = J.data
    cdef np.ndarray[int, mode="c", ndim=1] indices_py = J.indices.astype(np.intc)
    cdef np.ndarray[int, mode="c", ndim=1] indptr_py = J.indptr.astype(np.intc)

    py2c_d(data, data_py, len(J.data))
    py2c_i(indices, indices_py, len(J.indices))
    py2c_i(indptr, indptr_py, len(J.indptr))

    nnz[0] = J.nnz
    return 0

cdef int assemble_sparse_system_d(int n, double fac, int *nnz,
                                  double *data_in, int *indices_in, int *indptr_in,
                                  double *data_out, int *indices_out, int *indptr_out, 
                                  int flag_mass, double *mass_diag):
    ## TODO: Add docstring
    cdef np.ndarray[double, mode="c", ndim=1] data_J_py = np.empty(nnz[0], dtype = np.double)
    cdef np.ndarray[int, mode="c", ndim=1] indices_J_py = np.empty(nnz[0], dtype = np.intc)
    cdef np.ndarray[int, mode="c", ndim=1] indptr_J_py = np.empty(n + 1, dtype = np.intc)

    c2py_d(data_J_py, data_in, nnz[0])
    c2py_i(indices_J_py, indices_in, nnz[0])
    c2py_i(indptr_J_py, indptr_in, n + 1)

    # Reconstruct Jacobian
    J = sps.csc_matrix((data_J_py, indices_J_py, indptr_J_py), shape = (n, n))

    cdef np.ndarray[double, mode="c", ndim=1] M_diag = np.empty(n, dtype = np.double)
    if flag_mass:
        c2py_d(M_diag, mass_diag, n)
        M = sps.diags(M_diag, offsets = 0, shape = J.shape, format = 'csc')
    else:
        M = sps.eye(*J.shape, k = 0, dtype = np.double, format = 'csc')

    A = fac * M - J

    cdef np.ndarray[double, mode="c", ndim=1] data_A_py = A.data
    cdef np.ndarray[int, mode="c", ndim=1] indices_A_py = A.indices.astype(np.intc)
    cdef np.ndarray[int, mode="c", ndim=1] indptr_A_py = A.indptr.astype(np.intc)

    py2c_d(data_out, data_A_py, len(A.data))
    py2c_i(indices_out, indices_A_py, len(A.indices))
    py2c_i(indptr_out, indptr_A_py, len(A.indptr))
    nnz[0] = A.nnz

    return 0

cdef int assemble_sparse_system_z(int n, double fac_r, double fac_i, int *nnz,
                                  double *data_in, int *indices_in, int *indptr_in,
                                  doublecomplex *data_out, int *indices_out, int *indptr_out,
                                  int flag_mass, double *mass_diag):
    ## TODO: Add docstring
    cdef np.ndarray[double, mode="c", ndim=1] data_J_py = np.empty(nnz[0], dtype = np.double)
    cdef np.ndarray[int, mode="c", ndim=1] indices_J_py = np.empty(nnz[0], dtype = np.intc)
    cdef np.ndarray[int, mode="c", ndim=1] indptr_J_py = np.empty(n + 1, dtype = np.intc)

    c2py_d(data_J_py, data_in, nnz[0])
    c2py_i(indices_J_py, indices_in, nnz[0])
    c2py_i(indptr_J_py, indptr_in, n + 1)

    # Reconstruct Jacobian
    J = sps.csc_matrix((data_J_py, indices_J_py, indptr_J_py), shape = (n, n))

    cdef np.ndarray[double, mode="c", ndim=1] M_diag = np.empty(n, dtype = np.double)
    if flag_mass:
        c2py_d(M_diag, mass_diag, n)
        M = sps.diags(M_diag, offsets = 0, shape = J.shape, format = 'csc')
    else:
        M = sps.eye(*J.shape, k = 0, dtype = np.double, format = 'csc')

    A = (fac_r + 1j*fac_i) * M - J

    cdef np.ndarray[doublecomplex, mode="c", ndim=1] data_A_py = A.data
    cdef np.ndarray[int, mode="c", ndim=1] indices_A_py = A.indices.astype(np.intc)
    cdef np.ndarray[int, mode="c", ndim=1] indptr_A_py = A.indptr.astype(np.intc)

    py2c_z(data_out, data_A_py, len(A.data))
    py2c_i(indices_out, indices_A_py, len(A.indices))
    py2c_i(indptr_out, indptr_A_py, len(A.indptr))
    nnz[0] = A.nnz

    return 0

cpdef radau5(fcn_PY, doublereal x, np.ndarray y,
             doublereal xend, doublereal h__, np.ndarray rtol, np.ndarray atol,
             integer itol, jac_PY, integer ijac, integer mljac, integer mujac,
             mas_PY, integer imas, integer mlmas, integer mumas, solout_PY,
             integer iout, np.ndarray work, np.ndarray iwork):
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
            itol
                        - Switch for rtol and atol:
                          itol == 0: Both are scalars
                          itol == 1: Both are vectors
            jac_PY
                        - Jacobian function jac(x, y), where 'x' is time
            ijac
                        - Switch for Jacobian computation:
                          ijac == 0: C based finite differences
                          ijac == 1: Calls supplied 'jac_PY' function 
            mljac
                        - Switch for banded structure of Jacobian (used when solving the arising linear systems)
                          mljac == len(y): Full matrix Gauss-elimination
                          0 <= mljac < len(y): Size of non-zero lower diagonal bandwidth
            mujac
                        - Compare 'mljac', size of non-zero upper diagonal bandwidth, ignored if mljac == len(y)
            mas_PY
                        - Mass matrix function mas_PY(am)
            imas
                        - Switch for mass matrix:
                          imas == 0: Mass matrix is identity, 'mas_PY' never called
                          imas == 1: 'mas_PY' is used to determine mass matrix
            mlmas
                        - Switch for banded structure of Mass matrix, supposed to fulfill mlmas <= mljac
                          mlmas == len(y): Full matrix Gauss-elimination
                          0 <= mlmas < len(y): Size of non-zero lower diagonal bandwidth
            mumas
                        - Compare 'mumax', size of non-zero upper diagonal bandwidth, ignored if mlmas == len(y), supposed to fulfill mumas <= mujac
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
                        - Advanced tuning parameters of Radau solver, see radau_decsol.c for details
            iwork
                        - Advanced tuning parameters of Radau solver, see radau_decsol.c for details
        Returns::
            
            x
                        - Final time for which a solution has been computed, x == xend, if succesful
            y
                        - Final solution
            h__
                        - Prediced size of last accepted step
            iwork
                        - Statistics about number of function calls etc, see radau_decsol.c for details
            idid
                        - Return flag, see radau_decsol.c for details (1 == Successful computation)
    """
    # array lengthes, required for C call
    cdef integer n = len(y)
    cdef integer lwork = len(work)
    cdef integer liwork = len(iwork)
    
    # UNUSED: optional parameters used for communication between fcn, jac, mas, solout
    cdef doublereal rpar = 0
    cdef integer ipar = 0

    cdef integer idid = 1 ## "Successful compution"
    
    iwork_in = np.array(iwork, dtype = np.int32)
    cdef np.ndarray[double, mode="c", ndim=1] y_vec = y
    cdef np.ndarray[double, mode="c", ndim=1] rtol_vec = rtol
    cdef np.ndarray[double, mode="c", ndim=1] atol_vec = atol
    cdef np.ndarray[double, mode="c", ndim=1] work_vec = work
    cdef np.ndarray[integer, mode="c", ndim=1] iwork_vec = iwork_in
    
    radau5_c_py.radau5_c(n, callback_fcn, <void*>fcn_PY, &x, &y_vec[0], &xend,
                         &h__, &rtol_vec[0], &atol_vec[0], &itol, callback_jac, callback_jac_sparse, <void*> jac_PY,
                         &ijac, &mljac, &mujac, callback_mas, <void*> mas_PY, &imas, &mlmas, &mumas,
                         callback_solout, <void*>solout_PY, &iout, &work_vec[0], &lwork, &iwork_vec[0], &liwork, &rpar,
                         &ipar, &idid, assemble_sparse_system_d, assemble_sparse_system_z)
    
    return x, y, h__, np.array(iwork_in, dtype = np.int32), idid

cpdef contr5(integer i__, doublereal x, np.ndarray cont):
    """
        Python interface for calling the C based interpolation function using dense output
        Returns 'i'-th component at time 'x'. IMPORTANT: This function uses index 1 based notation.

            Parameters::

                i
                        - Which component to compute the solution for. IMPORTANT: starting index is 1, not 0
                x
                        - time-point at which the solution is requested. Needs to be within the time-interval defined by the last successful step.
                cont
                        - 'cont' input parameter to 'solout_PY' callback in 'radau5' function
            
            Returns::
                        - See function description

    """
    cdef integer lrc = len(cont)
    cdef np.ndarray[double, mode="c", ndim=1] cont_vec = cont
    return radau5_c_py.contr5_c(&i__, &x, &cont_vec[0], &lrc)
