#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Copyright (C) 2018-2021 Modelon AB, all rights reserved.
""" 

a=1 ## seemingly important TODO: Find better fix
    
## TODO: Possibly rename files a little
cimport radau5_c_py ## .pxd
cimport numpy as np
import numpy as np


"""
TODO:
    1. Make sure callback functions are correctly passed
    2. Make sure output parameters are handled correctly
    3. Make sure the change of IWORK to dtype = np.intc doesn't mess anything up in Fortran
"""

## TODO: y and f are arrays
cdef int fcn(integer n, doublereal x, doublereal y, doublereal f, doublereal rpar, integer ipar):
    f=y
    return 0

## TODO: y is array, dfy is matrix
cdef int jac(integer n, doublereal x, doublereal y, doublereal dfy, integer ldfy, doublereal rpar, integer ipar):
    dfy=y
    return 0

## TODO: am is matrix
cdef int mas(integer n, doublereal am, integer lmas, doublereal rpar, integer ipar):
    am=1
    return 0

## TODO: y, cont are arrays
cdef int solout(integer nr, doublereal xold, doublereal x, doublereal y,
                doublereal cont, integer lrc, integer n, doublereal rpar,
                integer ipar, integer irtrn):
    cont=0
    return 0

cpdef radau5(f_rhs, doublereal x, np.ndarray y,
             doublereal xend, doublereal h__, np.ndarray rtol, np.ndarray atol,
             integer itol, f_jac, integer ijac, integer mljac, integer mujac,
             f_mas, integer imas, integer mlmas, integer mumas, f_solout,
             integer iout, np.ndarray work, np.ndarray iwork):
    ## TODO: define these outside the function call?
    cdef integer n = len(y)
    cdef integer lwork = len(work)
    cdef integer liwork = len(iwork)
    
    ## TODO: in fortran these are referenced as dimension(1), should they be actual arrays?
    cdef doublereal rpar = 0 ## TODO: which value to choose?
    cdef integer ipar = 0 ## TODO: which value to choose?
    
    cdef integer idid = 0 ## TODO: Formally output value
    
    ## Array inputs which require appropriate conversion
    cdef np.ndarray[double,mode="c"] y_vec = y
    cdef np.ndarray[double,mode="c"] rtol_vec = rtol
    cdef np.ndarray[double,mode="c"] atol_vec = atol
    cdef np.ndarray[double,mode="c"] work_vec = work
    cdef np.ndarray[int,mode="c"] iwork_vec = iwork
    
    return radau5_c_py.radau5_c(&n, &fcn, &x, &y_vec[0], &xend, &h__, &rtol_vec[0], &rtol_vec[0], &itol, &jac,
                                &ijac, &mljac, &mujac, &mas, &imas, &mlmas, &mumas,
                                &solout, &iout, &work_vec[0], &lwork, &iwork_vec[0], &liwork, &rpar,
                                &ipar, &idid)

cpdef contr5(integer i__, doublereal x, np.ndarray cont):
    cdef np.ndarray[double,mode="c"] cont_vec = cont
    cdef integer lrc = len(cont)
    return radau5_c_py.contr5_c(&i__, &x, &cont_vec[0], &lrc)
