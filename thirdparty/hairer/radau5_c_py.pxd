#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Copyright (C) 2018-2021 Modelon AB, all rights reserved.
"""

cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)

cdef extern from "f2c.h":
    ctypedef int integer
    ctypedef double doublereal
    
## FunctionPointer_CallBack
ctypedef int (*FP_CB_f)(integer*, doublereal*, doublereal*, doublereal*,
                        doublereal*, integer*, void*)
ctypedef int (*FP_CB_jac)(integer*, doublereal*, doublereal*, doublereal*,
                          integer*, doublereal*, integer*, void*)
ctypedef int (*FP_CB_mas)(integer*, doublereal*, integer*, doublereal*,
                          integer*, void*)
ctypedef int (*FP_CB_solout)(integer*, doublereal*, doublereal*, doublereal*,
                             doublereal*, doublereal*, integer*, integer*,
                             doublereal*, integer*, integer*, void*)

cdef extern from "radau_decsol_c.h":
    int radau5_c(integer*, FP_CB_f, void*, doublereal*, doublereal*,
                 doublereal*, doublereal*, doublereal*, doublereal*,
                 integer*, FP_CB_jac, void*, integer*, integer*, integer*,
                 FP_CB_mas, void*, integer*, integer*, integer*, FP_CB_solout,
                 void*, integer*, doublereal*, integer*, integer*, integer*,
                 doublereal*, integer*, integer*)

    doublereal contr5_c(integer*, doublereal*, doublereal*, integer*)