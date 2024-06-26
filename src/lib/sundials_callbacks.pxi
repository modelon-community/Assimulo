#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
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

from numpy cimport PyArray_DATA

#=================
# Module functions
#=================

cdef N_Vector N_VNewEmpty_Euclidean(long int n) noexcept:
    IF SUNDIALS_VERSION >= (6,0,0):
        cdef SUNDIALS.SUNContext ctx = NULL
        IF SUNDIALS_VERSION >= (7,0,0):
            cdef SUNDIALS.SUNComm comm = 0
        ELSE:
            cdef void* comm = NULL
        SUNDIALS.SUNContext_Create(comm, &ctx)
        cdef N_Vector v = N_VNew_Serial(n, ctx)
    ELSE:
        cdef N_Vector v = N_VNew_Serial(n)
    v.ops.nvwrmsnorm = v.ops.nvwl2norm #Overwrite the WRMS norm to the 2-Norm
    return v

cdef inline N_Vector arr2nv(x) noexcept:
    x=np.array(x)
    cdef long int n = len(x)
    cdef np.ndarray[realtype, ndim=1,mode='c'] ndx=x
    cdef void* data_ptr=PyArray_DATA(ndx)
    IF SUNDIALS_VERSION >= (6,0,0):
        cdef SUNDIALS.SUNContext ctx = NULL
        IF SUNDIALS_VERSION >= (7,0,0):
            cdef SUNDIALS.SUNComm comm = 0
        ELSE:
            cdef void* comm = NULL
        SUNDIALS.SUNContext_Create(comm, &ctx)
        cdef N_Vector v = N_VNew_Serial(n, ctx)
    ELSE:
        cdef N_Vector v = N_VNew_Serial(n)
    memcpy((<N_VectorContent_Serial>v.content).data, data_ptr, n*sizeof(realtype))
    return v

cdef inline N_Vector arr2nv_euclidean(x) noexcept:
    x=np.array(x)
    cdef long int n = len(x)
    cdef np.ndarray[realtype, ndim=1,mode='c'] ndx=x
    cdef void* data_ptr=PyArray_DATA(ndx)
    cdef N_Vector v=N_VNewEmpty_Euclidean(n)
    memcpy((<N_VectorContent_Serial>v.content).data, data_ptr, n*sizeof(realtype))
    return v
    
cdef inline void arr2nv_inplace(x, N_Vector out) noexcept:
    x=np.array(x)
    cdef long int n = len(x)
    cdef np.ndarray[realtype, ndim=1,mode='c'] ndx=x
    cdef void* data_ptr=PyArray_DATA(ndx)
    memcpy((<N_VectorContent_Serial>out.content).data, data_ptr, n*sizeof(realtype))
    
cdef inline np.ndarray nv2arr(N_Vector v):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    cdef np.ndarray[realtype, ndim=1, mode='c'] x=np.empty(n)
    memcpy(PyArray_DATA(x), v_data, n*sizeof(realtype))
    return x
    
cdef inline void nv2arr_inplace(N_Vector v, np.ndarray o) noexcept:
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    memcpy(PyArray_DATA(o), v_data, n*sizeof(realtype))
    
cdef inline void nv2mat_inplace(int Ns, N_Vector *v, np.ndarray o) noexcept:
    cdef long int i,j, Nf
    for i in range(Ns):
        Nf = (<N_VectorContent_Serial>v[i].content).length
        for j in range(Nf):
            o[j,i] = (<N_VectorContent_Serial>v[i].content).data[j]

cdef inline realtype2arr(realtype *data, int n):
    """Create new numpy array from realtype*"""
    cdef np.ndarray[realtype, ndim=1, mode='c'] x=np.empty(n)
    memcpy(PyArray_DATA(x), data, n*sizeof(realtype))
    return x
