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

from numpy cimport PyArray_DATA

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c(double* dest, object source, int dim):
    """
    Copy 1D numpy array data to (double *) C vector
    """
    cdef double* data
    if not (isinstance(source, np.ndarray) and source.flags.contiguous and source.dtype == np.float):
        source = np.ascontiguousarray(source, dtype=np.float)
    assert source.size >= dim, "The dimension of the vector is {} and not equal to the problem dimension {}. Please verify the output vectors from the min/max/nominal/evalute methods in the Problem class.".format(source.size, dim)
    data = <double*>PyArray_DATA(source)
    memcpy(dest, data, dim*sizeof(double))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c_matrix_flat_F(double* dest, object source, int nrow, int ncol):
    """
    Copy (square) 2D numpy array (order = c) to (double *) C matrix (with Fortran-style column major ordering)
    """
    cdef np.ndarray[double, ndim=2] source_np = np.array(source, copy=False, dtype = np.float)
    for i in range(ncol):
        for j in range(nrow):
            dest[j + i*nrow] = source_np[j][i]
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c2py(np.ndarray[double, ndim=1, mode='c'] dest, double* source, int dim):
    """
    Copy (double *) C vector to 1D numpy array
    """
    memcpy(PyArray_DATA(dest), source, dim*sizeof(double))
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c2py_mat_F(np.ndarray[double, ndim=2, mode='fortran'] dest, double* source, int dim):
    """
    Copy (double *) C matrix (Fotran-style column major ordering) to 2D numpy array
    """
    memcpy(PyArray_DATA(dest), source, dim*sizeof(double))

cdef int callback_fcn(integer* n, doublereal* x, doublereal* y_in, doublereal* y_out,
                      doublereal* rpar, integer* ipar, void* fcn_PY):
    """
    Internal callback function to enable call to Python based rhs function from C
    """
    cdef np.ndarray[double, ndim=1, mode="c"]y_py_in = np.empty(n[0], dtype = np.double)
    c2py(y_py_in, y_in, n[0])
    res = (<object>fcn_PY)(x[0], y_py_in)
    py2c(y_out, res[0], len(res[0]))
    ipar[0] = res[1][0]
    return 0

cdef int callback_jac(integer* n, doublereal* x, doublereal* y, doublereal* fjac,
                      integer* ldjac, doublereal* rpar, integer* ipar, void* jac_PY):
    """
    Internal callback function to enable call to Python based Jacobian function from C
    """
    cdef np.ndarray[double, ndim=1, mode="c"]y_py_in = np.empty(n[0], dtype = np.double)
    c2py(y_py_in, y, n[0])
    res = (<object>jac_PY)(x[0], y_py_in)
    py2c_matrix_flat_F(fjac, res, res.shape[0], res.shape[1])
    return 0

cdef int callback_mas(integer* n, doublereal* am, integer* lmas, doublereal* rpar,
                      integer* ipar, void* mas_PY):
    """
    Internal callback function to enable call to Python based mass matrix function from C
    """
    cdef np.ndarray[double, mode="fortran", ndim=2]am_py = np.empty((lmas[0], n[0]), order = 'F', dtype = np.double)
    c2py_mat_F(am_py, am, n[0]*lmas[0])
    res = (<object>mas_PY)(am_py)
    py2c_matrix_flat_F(am, res, res.shape[0], res.shape[1])
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
    c2py(y_py, y, nsolu[0])
    c2py(cont_py, cont, 4*nsolu[0])
    c2py(werr_py, werr, nsolu[0])

    irtrn[0] = (<object>solout_PY)(nrsol[0], xosol[0], xsol[0],
                                   y_py, cont_py, werr_py,
                                   lrc[0], irtrn[0])

    return irtrn[0]

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
    
    iwork_in = np.array(iwork, dtype = np.int64)
    cdef np.ndarray[double, mode="c", ndim=1] y_vec = y
    cdef np.ndarray[double, mode="c", ndim=1] rtol_vec = rtol
    cdef np.ndarray[double, mode="c", ndim=1] atol_vec = atol
    cdef np.ndarray[double, mode="c", ndim=1] work_vec = work
    cdef np.ndarray[integer, mode="c", ndim=1] iwork_vec = iwork_in
    
    radau5_c_py.radau5_c(&n, callback_fcn, <void*>fcn_PY, &x, &y_vec[0], &xend,
                         &h__, &rtol_vec[0], &atol_vec[0], &itol, callback_jac, <void*> jac_PY,
                         &ijac, &mljac, &mujac, callback_mas, <void*> mas_PY, &imas, &mlmas, &mumas,
                         callback_solout, <void*>solout_PY, &iout, &work_vec[0], &lwork, &iwork_vec[0], &liwork, &rpar,
                         &ipar, &idid)
    
    return x, y, h__, np.array(iwork_in, dtype = int), idid

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
