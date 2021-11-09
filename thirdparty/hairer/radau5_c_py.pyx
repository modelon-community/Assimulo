#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Copyright (C) 2018-2021 Modelon AB, all rights reserved.
""" 

cimport radau5_c_py
cimport numpy as np
cimport cython
import numpy as np
from cython.view cimport array as cvarray

from numpy cimport PyArray_DATA

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void py2c(double* dest, object source, int dim):
    cdef double* data
    if not (isinstance(source, np.ndarray) and source.flags.contiguous and source.dtype == np.float):
        source = np.ascontiguousarray(source, dtype=np.float)
    assert source.size >= dim, "The dimension of the vector is {} and not equal to the problem dimension {}. Please verify the output vectors from the min/max/nominal/evalute methods in the Problem class.".format(source.size, dim)
    data = <double*>PyArray_DATA(source)
    memcpy(dest, data, dim*sizeof(double))
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c2py(np.ndarray[double, ndim=1,mode='c'] dest, double* source, int dim):
    memcpy(dest.data, source, dim*sizeof(double))
    
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void c2py_mat(np.ndarray[double, ndim=2,mode='c'] dest, double* source, int dim):
    memcpy(dest.data, source, dim*sizeof(double))

cdef int callback_fcn(integer* n, doublereal* x, doublereal* y_in, doublereal* y_out,
                      doublereal* rpar, integer* ipar, void* fcn_PY):
    cdef np.ndarray[double,mode="c"]y_py_in = np.zeros(n[0])
    c2py(y_py_in, y_in, n[0])
    res = (<object>fcn_PY)(x[0], y_py_in)
    py2c(y_out, res[0], res[0].size)
    ipar[0] = res[1][0]
    return 0

cdef int callback_jac(integer* n, doublereal* x, doublereal* y, doublereal* fjac,
                      integer* ldjac, doublereal* rpar, integer* ipar, void* jac_PY):
    cdef np.ndarray[double,mode="c"]y_py = np.zeros(n[0])
    c2py(y_py, y, n[0])
    res = (<object>jac_PY)(x[0], y_py)
    res = res.flatten()
    py2c(fjac, res, res.size)
    return 0

cdef int callback_mas(integer* n, doublereal* am, integer* lmas, doublereal* rpar,
                      integer* ipar, void* mas_PY):
    cdef np.ndarray[double,mode="c",ndim=2]am_py = np.zeros((lmas[0], n[0]))
    c2py_mat(am_py, am, n[0]*lmas[0])
    res = (<object>mas_PY)(am_py)
    res = res.flatten()
    py2c(am, res, res.size)
    return 0

cdef int callback_solout(integer* nrsol, doublereal* xosol, doublereal* xsol, doublereal* y,
                         doublereal* cont, doublereal* werr, integer* lrc, integer* nsolu,
                         doublereal* rpar, integer* ipar, integer* irtrn, void* solout_PY):
    cdef double[:] y_py = cvarray(shape=(nsolu[0],), itemsize=sizeof(double), format="d")
    cdef double[:] cont_py = cvarray(shape=(4*nsolu[0],), itemsize=sizeof(double), format="d")
    cdef double[:] werr_py = cvarray(shape=(nsolu[0],), itemsize=sizeof(double), format="d")
    c2py(np.asarray(y_py), y, nsolu[0])
    c2py(np.asarray(cont_py), cont, 4*nsolu[0])
    c2py(np.asarray(werr_py), cont, nsolu[0])

    irtrn[0] = (<object>solout_PY)(nrsol[0], xosol[0], xsol[0],
                                   np.asarray(y_py), np.asarray(cont_py), np.asarray(werr_py),
                                   lrc[0], irtrn[0])
    return irtrn[0]

cpdef radau5(fcn_PY, doublereal x, np.ndarray y,
             doublereal xend, doublereal h__, np.ndarray rtol, np.ndarray atol,
             integer itol, jac_PY, integer ijac, integer mljac, integer mujac,
             mas_PY, integer imas, integer mlmas, integer mumas, solout_PY,
             integer iout, np.ndarray work, np.ndarray iwork):
    # array lengthes, required for C call
    cdef integer n = len(y)
    cdef integer lwork = len(work)
    cdef integer liwork = len(iwork)
    
    # UNUSED: optional parameters used for communication between fcn, jac, mas, solout
    cdef doublereal rpar = 0
    cdef integer ipar = 0

    cdef integer idid = 1 ## "Successful compution"
    
    cdef np.ndarray[double,mode="c"] y_vec = y
    cdef np.ndarray[double,mode="c"] rtol_vec = rtol
    cdef np.ndarray[double,mode="c"] atol_vec = atol
    cdef np.ndarray[double,mode="c"] work_vec = work
    cdef np.ndarray[int,mode="c"] iwork_vec = iwork
    
    radau5_c_py.radau5_c(&n, callback_fcn, <void*>fcn_PY, &x, &y_vec[0], &xend,
                         &h__, &rtol_vec[0], &rtol_vec[0], &itol, callback_jac, <void*> jac_PY,
                         &ijac, &mljac, &mujac, callback_mas, <void*> mas_PY, &imas, &mlmas, &mumas,
                         callback_solout, <void*>solout_PY, &iout, &work_vec[0], &lwork, &iwork_vec[0], &liwork, &rpar,
                         &ipar, &idid)
    return x, y, h__, iwork, idid

cpdef contr5(integer i__, doublereal x, np.ndarray cont):
    cdef np.ndarray[double,mode="c"] cont_vec = cont
    cdef integer lrc = len(cont)
    return radau5_c_py.contr5_c(&i__, &x, &cont_vec[0], &lrc)
