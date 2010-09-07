#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


from __future__ import division


include "sundials_kinsol_core.pxd" # Includes fcts from other header files
include "sundials_kinsol_core.pxi" # Includes the constants (textual include)


cdef int kin_res(N_Vector xv, N_Vector fval, void *problem_data):
    """
    Residual fct called by KINSOL
    """
    cdef:
        int i
        ProblemData pData = <ProblemData>problem_data
        ndarray[realtype, ndim = 1, mode = 'c'] rhs # Container holding return from user fct
        
        ndarray x = nv2arr(xv)
        realtype* resptr = (<N_VectorContent_Serial>fval.content).data

    try:
        rhs = (<object>pData.RHS)(x)

        for i in range(pData.dim):
            resptr[i] = rhs[i]

        return KIN_SUCCESS
    except:
        return KIN_REC_ERR

class KINError(Exception):
    """
    Kinsol exception
    """
    
    msg = {KIN_MEM_FAIL : 'Memory allocation failed.',
                KIN_MEM_NULL : 'KINSOL not properly created.',
                KIN_ILL_INPUT: 'There was an input error.',
                KIN_NO_MALLOC: 'Memory not allocated, call KINSOL_init(...)',
                KIN_LINESEARCH_NONCONV: 'Linesearch could not find a step big enough.',
                KIN_MAXITER_REACHED: 'Max number of iterations reached',
                KIN_MXNEWT_5X_EXCEEDED: 'Max size of newton step exceeded 5 times.',
                KIN_LINESEARCH_BCFAIL: 'Linesearch beta-test failed, probably because of too poor progress.',
                KIN_LINSOLV_NO_RECOVERY: 'psolve encountered an error, but preconditioner is current',
                KIN_LINIT_FAIL: 'Init of linear solver failed.',
                KIN_LSETUP_FAIL: 'pset (preconditioner setup fct) encountered an error',
                KIN_LSOLVE_FAIL: 'Error in either psolve or in linear solver routine.',
                KIN_SYSFUNC_FAIL: 'Call to RHS failed.',
                KIN_FIRST_SYSFUNC_ERR: 'Call to RHS failed on first call',
                KIN_REPTD_SYSFUNC_ERR: 'Call to RHS failed multiple times.'}
    def __init__(self, value):
        self.value = value

    def __str__(self):
        try:
            return repr(self.msg[self.value])
        except KeyError:
            return repr('KINSOL failed with flag %s.'%(self.value))

cdef class KINSOL_wrap:

    cdef:
        void* solver 
        ProblemData pData # A struct containing problem data
        #ProblemData *ppData # Pointer to above mentioned struct
        N_Vector x_cur, x_scale, f_scale

    def __cinit__(self):
        self.pData = ProblemData() #Create a new problem struct
        #self.ppData = &self.pData

    cpdef KINSOL_set_problem_info(self, RHS, dim):
        """
        Sets the problem info to the desired struct
        """

        # Sets residual function
        self.pData.RHS = <void*>RHS
        self.pData.dim = dim

    def KINSOL_init(self,RHS,x0,dim):
        """
        Initializes solver
        """        
        cdef int flag
        # Set problem info
        self.KINSOL_set_problem_info(RHS,dim)

        # Create initial guess from the supplied numpy array
        self.x_cur = arr2nv(x0)

        if self.solver == NULL: # solver runs for the first time
            
            # Create solver
            self.solver = KINCreate()
            if self.solver == NULL:
                raise KINError(KIN_MEM_FAIL)
            
            # Allocate internal memory
            flag = KINInit(self.solver,kin_res,self.x_cur)
            if flag < 0:
                raise KINError(flag)

            # Link to linear solver, for the moment not specified by user
            # but this will eventually be implemented
            flag = KINDense(self.solver,self.pData.dim)
            if flag < 0:
                raise KINError(flag)
        else:
            pass

        # Link the solver to the supplied problem data
        flag = KINSetUserData(self.solver,<void*>self.pData)
        if flag < 0:
            raise KINError(flag)

        # Create scaling vectors filled with ones
        # since the problem is assumed to be scaled
        self.x_scale = arr2nv(np.ones(self.pData.dim))
        self.f_scale = arr2nv(np.ones(self.pData.dim))
    
    def KINSOL_solve(self):
        """
        Function that should be called after a call to KINSOL_init
        solves the function supplied as RHS
        """
        
        cdef int flag
        
        flag = KINSol(<void*>self.solver,self.x_cur,KIN_LINESEARCH,self.x_scale,self.f_scale)
        if flag < 0:
            print "Error in solve, flag: ", flag
            raise KINError(flag)

        return nv2arr(self.x_cur)
        
