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
                      doublereal* rpar, integer* ipar, void* jac_PY):
    """
    Internal callback function to enable call to Python based Jacobian function from C
    """
    cdef np.ndarray[double, ndim=1, mode="c"]y_py_in = np.empty(n, dtype = np.double)
    c2py_d(y_py_in, y, n)
    res = (<object>jac_PY)(x[0], y_py_in)
    py2c_d_matrix_flat_F(fjac, res, res.shape[0], res.shape[1])
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

cdef class RadauSuperLUaux:
    """Auxiliary data structure to have the memory ownership of internal
    and auxiliary superLU data structures available on the highest possible level.
    """
    cdef SuperLU_aux_d* superLU_aux_struct_d # real
    cdef SuperLU_aux_z* superLU_aux_struct_z # complex
    ## data storage of sparse jacobian
    cdef double* jac_data
    cdef int* jac_indicies
    cdef int* jac_indptr
    
    cpdef int initialize(self, int nprocs, int n, int nnz):
        self.superLU_aux_struct_d = superlu_init_d(nprocs, n, nnz)
        self.superLU_aux_struct_z = superlu_init_z(nprocs, n, nnz)

        radau_sparse_aux_init(&self.jac_data, &self.jac_indicies, &self.jac_indptr, nnz, n)
        return 0

    cpdef int finalize(self):
        radau_sparse_aux_finalize(&self.jac_data, &self.jac_indicies, &self.jac_indptr)
        cdef int ret
        ret = superlu_finalize_d(self.superLU_aux_struct_d)
        if ret != 0:
            return ret
        ret = superlu_finalize_z(self.superLU_aux_struct_z)
        return ret

cdef int callback_jac_sparse(int n, double *x, double *y, int *nnz,
                             double * data, int *indices, int *indptr,
                             doublereal* rpar, integer* ipar,
                             void* jac_PY):
    """
    Internal callback function to enable call to Python based evaluation of sparse (csc) jacobians.

    Checks if jacobian includes the diagonal.
    If not, it adds corresponding zeroes in the CSC format.
    This is since the resulting system matrix will add values on the diagonal.
    """
    cdef np.ndarray[double, ndim=1, mode="c"]y_py = np.empty(n, dtype = np.double)
    c2py_d(y_py, y, n)

    J = (<object>jac_PY)(x[0], y_py)
    if not isinstance(J, sps.csc.csc_matrix):
        return -1

    if J.nnz > nnz[0]:
        return J.nnz

    cdef np.ndarray[double, mode="c", ndim=1] jac_data_py = J.data.astype(np.double)
    cdef np.ndarray[int, mode="c", ndim=1] jac_indices_py = J.indices.astype(np.intc)
    cdef np.ndarray[int, mode="c", ndim=1] jac_indptr_py = J.indptr.astype(np.intc)

    ## Check if any diagonal entries are missing in Jacobian

    missing = [] ## list of diagonal entries that are missing
    ## iterate over columns
    for i in range(n):
        indptr_base = jac_indptr_py[i]
        ## take slice of indices
        for j in range(jac_indptr_py[i+1] - jac_indptr_py[i]):
            if i == jac_indices_py[indptr_base + j]: ## diagonal index found in indices
                break
        else: ## loop done, diagonal index not found
            missing.append(i)

    if not missing: ## i.e., not empty
        ## copy data to output structures
        py2c_d(data, jac_data_py, len(J.data))
        py2c_i(indices, jac_indices_py, len(J.indices))
        py2c_i(indptr, jac_indptr_py, n + 1)
        nnz[0] = J.nnz
        return 0
    ## else: ## some diagonal entries are missing
    ## else skipped since cdefs not allowed in blocks
    ## CSC format of Jacobian needs to be updated

    cdef np.ndarray[double, mode="c", ndim=1] jac_data_py_new = np.empty(nnz[0] + len(missing), dtype = np.double)
    cdef np.ndarray[int, mode="c", ndim=1] jac_indices_py_new = np.empty(nnz[0] + len(missing), dtype = np.intc)
    cdef np.ndarray[int, mode="c", ndim=1] jac_indptr_py_new = np.empty(n + 1, dtype = np.intc)
    jac_indptr_py_new[0] = 0 ## always starts with 0

    cdef int idx_data_indices_next = 0
    cdef int indptr_incr = 0 ## increment of indptr from previous diagonal insertions
    cdef int intptr_base = 0

    ## iterate though all columns
    for i in range(n):
        intptr_base = jac_indptr_py[i] ## for current column
        if i in missing:
            indptr_incr += 1 ## increment indptr increment
            missing_added = False
            ## iterate through old column
            j = 0 ## required, since range(0) in next line would cause issue later on
            for j in range(jac_indptr_py[i+1] - jac_indptr_py[i]):
                ## insert diagonal entry if larger than current index
                if not missing_added:
                    if jac_indices_py[intptr_base + j] > i: ## correct position found
                        jac_indices_py_new[idx_data_indices_next] = i
                        jac_data_py_new[idx_data_indices_next] = 0
                        idx_data_indices_next += 1
                        missing_added = True
                ## add back old column index
                jac_indices_py_new[idx_data_indices_next] = jac_indices_py[intptr_base + j]
                jac_data_py_new[idx_data_indices_next] = jac_data_py[intptr_base + j]
                idx_data_indices_next += 1
            if not missing_added: ## column loop done, but diagonal entry not yet added
                jac_indices_py_new[idx_data_indices_next] = i
                jac_data_py_new[idx_data_indices_next] = 0
                idx_data_indices_next += 1
        else: ## column already contains diagonal entry
            ## simple loop to add old column
            for j in range(jac_indptr_py[i+1] - jac_indptr_py[i]):
                jac_indices_py_new[idx_data_indices_next] = jac_indices_py[intptr_base + j]
                jac_data_py_new[idx_data_indices_next] = jac_data_py[intptr_base + j]
                idx_data_indices_next += 1
        ## add old inptr entry, possibly increment by earlier added entries
        jac_indptr_py_new[i+1] = jac_indptr_py[i+1] + indptr_incr

    py2c_d(data, jac_data_py_new, len(J.data) + len(missing))
    py2c_i(indices, jac_indices_py_new, len(J.indices) + len(missing))
    py2c_i(indptr, jac_indptr_py_new, len(J.indptr))
    nnz[0] = J.nnz + len(missing)
    return 0

cpdef radau5(fcn_PY, doublereal x, np.ndarray y,
             doublereal xend, doublereal h__, np.ndarray rtol, np.ndarray atol,
             integer itol, jac_PY, integer ijac, integer mljac, integer mujac,
             mas_PY, integer imas, integer mlmas, integer mumas, solout_PY,
             integer iout, np.ndarray work, np.ndarray iwork,
             RadauSuperLUaux aux_class):
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
            aux_class
                        - instance of RadauSuperLUaux, needs to be initialized via aux_class.initialize for SPARSE solver
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

    if iwork[10]: ## sparse
        radau5_c_py.radau5_c(n, callback_fcn, <void*>fcn_PY, &x, &y_vec[0], &xend,
                            &h__, &rtol_vec[0], &atol_vec[0], &itol, callback_jac, callback_jac_sparse, <void*> jac_PY,
                            &ijac, &mljac, &mujac, &imas, &mlmas, &mumas,
                            callback_solout, <void*>solout_PY, &iout, &work_vec[0], &lwork, &iwork_vec[0], &liwork, &rpar,
                            &ipar, &idid,
                            aux_class.jac_data, aux_class.jac_indicies, aux_class.jac_indptr,
                            aux_class.superLU_aux_struct_d, aux_class.superLU_aux_struct_z)
    else: ## Dense
        radau5_c_py.radau5_c(n, callback_fcn, <void*>fcn_PY, &x, &y_vec[0], &xend,
                            &h__, &rtol_vec[0], &atol_vec[0], &itol, callback_jac, callback_jac_sparse, <void*> jac_PY,
                            &ijac, &mljac, &mujac, &imas, &mlmas, &mumas,
                            callback_solout, <void*>solout_PY, &iout, &work_vec[0], &lwork, &iwork_vec[0], &liwork, &rpar,
                            &ipar, &idid,
                            NULL, NULL, NULL, NULL, NULL)
    
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
