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

cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)

cdef extern from "radau5_superlu_double.h":
    ctypedef struct SuperLU_aux_d
    SuperLU_aux_d* superlu_init_d(int, int, int)
    int superlu_finalize_d(SuperLU_aux_d*)

cdef extern from "radau5_superlu_complex.h":
    ctypedef struct SuperLU_aux_z
    SuperLU_aux_z* superlu_init_z(int, int, int)
    int superlu_finalize_z(SuperLU_aux_z*)

cdef extern from "radau_decsol_c.h":
    ctypedef struct Radau_SuperLU_aux

    int RADAU_CALLBACK_OK
    int RADAU_CALLBACK_ERROR_RECOVERABLE
    int RADAU_CALLBACK_ERROR_NONRECOVERABLE
    int RADAU_CALLBACK_ERROR_INVALID_JAC_FORMAT
    int RADAU_CALLBACK_ERROR_INVALID_NNZ

    ## FunctionPointer_CallBack
    ctypedef int (*FP_CB_f)(int, double*, double*, double*, void*)
    ctypedef int (*FP_CB_jac)(int, double*, double*, double*, void*)
    ctypedef int (*FP_CB_solout)(int*, double*, double*, double*,
                                double*, double*, int*, int*,
                                int*, void*)
    ctypedef int (*FP_CB_jac_sparse)(int, double*, double*, int*, double*, int*, int*, void*)

    int radau5_c(int, FP_CB_f, void*,
			 double*, double*, double*, double*,
			 double*, double*, int*,
			 FP_CB_jac, FP_CB_jac_sparse, void*, int*, int,
			 FP_CB_solout, void*, int*,
			 double*, int*, int*, int*, int*,
			 Radau_SuperLU_aux*)

    double contr5_c(int*, double*, double*, int*)

    Radau_SuperLU_aux* radau_superlu_aux_setup(int, int, int, int*);
    int radau_superlu_aux_finalize(Radau_SuperLU_aux*);
