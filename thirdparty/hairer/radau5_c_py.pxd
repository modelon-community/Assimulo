#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Copyright (C) 2018-2021 Modelon AB, all rights reserved.
"""

## .pxd acts as header file for .pyx file

cdef extern from "f2c.h":
    ctypedef int integer
    ctypedef double doublereal
    
## FP = Function Pointer
ctypedef int (*FP_f)(integer, doublereal, doublereal, doublereal, doublereal, integer)
ctypedef int (*FP_jac)(integer, doublereal, doublereal, doublereal, integer, doublereal, integer)
ctypedef int (*FP_mas)(integer, doublereal, integer, doublereal, integer)
ctypedef int (*FP_solout)(integer, doublereal, doublereal, doublereal, doublereal,
                          integer, integer, doublereal, integer, integer)
    
cdef extern from "radau_decsol_c.h":
    ## TODO: remove various input parameters here and instead infer them from the others
    ## e.g. n = len(y)
    ## See .pyf file for reference, try to get a signature identical to the fotran version
    int radau5_c(integer *n, FP_f fcn, doublereal *x, doublereal *y,
                 doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *atol,
                 integer *itol, FP_jac jac, integer *ijac, integer *mljac, integer *mujac,
                 FP_mas mas, integer *imas, integer *mlmas, integer *mumas, FP_solout solout,
                 integer *iout, doublereal *work, integer *lwork, integer *iwork,
                 integer *liwork, doublereal *rpar, integer *ipar, integer *idid)

    doublereal contr5_c(integer *i__, doublereal *x, doublereal *cont, integer *lrc)