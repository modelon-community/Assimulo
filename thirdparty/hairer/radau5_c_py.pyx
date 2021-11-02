#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Copyright (C) 2018-2021 Modelon AB, all rights reserved.
""" 

a=1 ## seemingly important TODO: Find better fix
    
## TODO: Possibly rename files a little
cimport radau5_c_py ## .pxd


"""
TODO:
    1. Figure out if libf2c.a was properly included
    2. See if I can find the actual forwarding code from fortran, check makefile log?
    3. Use actual sensible functions
    4. array stuff?
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


## TODO: inputs for actual functions currently only dummies
cpdef radau5(integer n, f_rhs, doublereal x, doublereal y,
             doublereal xend, doublereal h__, doublereal rtol, doublereal atol,
             integer itol, f_jac, integer ijac, integer mljac, integer mujac,
             f_mas, integer imas, integer mlmas, integer mumas, f_solout,
             integer iout, doublereal work, integer lwork, integer iwork,
             integer liwork, doublereal rpar, integer ipar, integer idid):
    return radau5_c_py.radau5_c(&n, &fcn, &x, &y, &xend, &h__, &rtol, &atol, &itol, &jac,
                                &ijac, &mljac, &mujac, &mas, &imas, &mlmas, &mumas,
                                &solout, &iout, &work, &lwork, &iwork, &liwork, &rpar,
                                &ipar, &idid)

cpdef contr5(integer i__, doublereal x, doublereal cont, integer lrc):
    return radau5_c_py.contr5_c(&i__, &x, &cont, &lrc)