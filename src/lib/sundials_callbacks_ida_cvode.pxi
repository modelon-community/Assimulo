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

import cython


cdef int cv_rhs(realtype t, N_Vector yv, N_Vector yvdot, void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.f to the Sundials
    right-hand-side function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray y = pData.work_y
    cdef realtype* resptr=(<N_VectorContent_Serial>yvdot.content).data
    cdef int i
    
    nv2arr_inplace(yv, y)
    
    if pData.dimSens>0: #Sensitivity activated
        p = realtype2arr(pData.p,pData.dimSens)
        try:
            if pData.sw != NULL:
                rhs = (<object>pData.RHS)(t,y,sw=<list>pData.sw, p=p)
            else:
                rhs = (<object>pData.RHS)(t,y,p)
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)
        
    else: #No sensitivity
        try:
            if pData.sw != NULL:
                rhs = (<object>pData.RHS)(t,y,<list>pData.sw)
            else:
                rhs = (<object>pData.RHS)(t,y)
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)
    
    for i in range(pData.dim):
        resptr[i] = rhs[i]
    
    return CV_SUCCESS
            
cdef int cv_sens_rhs_all(int Ns, realtype t, N_Vector yv, N_Vector yvdot,
                         N_Vector *yvS, N_Vector *yvSdot, void *problem_data, 
                         N_Vector tmp1, N_Vector tmp2):
    
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray y = pData.work_y
    cdef N.ndarray s = pData.work_ys
    cdef realtype* resptr
    cdef int i, j
    
    nv2arr_inplace(yv, y)
    nv2mat_inplace(Ns, yvS, s)
    p = realtype2arr(pData.p,pData.dimSens)
    
    try:
        if pData.sw != NULL:
            sens_rhs = (<object>pData.RHS_SENS_ALL)(t,y,s,p,<list>pData.sw)
        else:
            sens_rhs = (<object>pData.RHS_SENS_ALL)(t,y,s,p)
        
        for i in range(Ns):
            resptr=(<N_VectorContent_Serial>yvSdot[i].content).data
            for j in range(pData.dim):
                resptr[j] = sens_rhs[j,i]
        
        return CV_SUCCESS
    except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
        return CV_REC_ERR
    except:
        traceback.print_exc()
        return CV_UNREC_RHSFUNC_ERR 


IF SUNDIALS_VERSION >= (3,0,0):
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int cv_jac_sparse(realtype t, N_Vector yv, N_Vector fy, SUNMatrix Jac,
                    void *problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
        """
        This method is used to connect the Assimulo.Problem.jac to the Sundials
        Sparse Jacobian function.
        """
        cdef ProblemData pData = <ProblemData>problem_data
        cdef SUNMatrixContent_Sparse Jacobian = <SUNMatrixContent_Sparse>Jac.content
        cdef N.ndarray y = pData.work_y
        cdef int i
        cdef sunindextype nnz = Jacobian.NNZ
        cdef int ret_nnz
        cdef sunindextype dim = Jacobian.N
        cdef realtype* data = Jacobian.data
        cdef sunindextype* rowvals = Jacobian.rowvals[0]
        cdef sunindextype* colptrs = Jacobian.colptrs[0]
        
        nv2arr_inplace(yv, y)

        try:
            if pData.dimSens > 0: #Sensitivity activated
                p = realtype2arr(pData.p,pData.dimSens)
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(t,y,p=p,sw=<list>pData.sw)
                else:
                    jac=(<object>pData.JAC)(t,y,p=p)
            else:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(t,y,sw=<list>pData.sw)
                else:
                    jac=(<object>pData.JAC)(t,y)
                
            if not isinstance(jac, sparse.csc.csc_matrix):
                jac = sparse.csc.csc_matrix(jac)
                raise AssimuloException("The Jacobian must be stored on Scipy's CSC format.")
            ret_nnz = jac.nnz
            if ret_nnz > nnz:
                raise AssimuloException("The Jacobian has more entries than supplied to the problem class via 'jac_nnz'")    

            for i in range(min(ret_nnz,nnz)):
                data[i]    = jac.data[i]
                rowvals[i] = jac.indices[i]
            for i in range(dim+1):
                colptrs[i] = jac.indptr[i]
            
            return CVDLS_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
        except:
            traceback.print_exc()
            return CVDLS_JACFUNC_UNRECVR
ELSE:
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef int cv_jac_sparse(realtype t, N_Vector yv, N_Vector fy, SlsMat Jacobian,
                    void *problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
        """
        This method is used to connect the Assimulo.Problem.jac to the Sundials
        Sparse Jacobian function.
        """
        cdef ProblemData pData = <ProblemData>problem_data
        cdef N.ndarray y = pData.work_y
        cdef int i
        cdef int nnz = Jacobian.NNZ
        cdef int ret_nnz
        cdef int dim = Jacobian.N
        cdef realtype* data = Jacobian.data
        
        IF SUNDIALS_VERSION >= (2,6,3):
            cdef int* rowvals = Jacobian.rowvals[0]
            cdef int* colptrs = Jacobian.colptrs[0]
        ELSE:
            cdef int* rowvals = Jacobian.rowvals
            cdef int* colptrs = Jacobian.colptrs
        
        nv2arr_inplace(yv, y)
        """
            realtype *data;
            int *rowvals;
            int *colptrs;
        """
        try:
            if pData.dimSens > 0: #Sensitivity activated
                p = realtype2arr(pData.p,pData.dimSens)
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(t,y,p=p,sw=<list>pData.sw)
                else:
                    jac=(<object>pData.JAC)(t,y,p=p)
            else:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(t,y,sw=<list>pData.sw)
                else:
                    jac=(<object>pData.JAC)(t,y)
                
            if not isinstance(jac, sparse.csc.csc_matrix):
                jac = sparse.csc.csc_matrix(jac)
                raise AssimuloException("The Jacobian must be stored on Scipy's CSC format.")
            ret_nnz = jac.nnz
            if ret_nnz > nnz:
                raise AssimuloException("The Jacobian has more entries than supplied to the problem class via 'jac_nnz'")    

            for i in range(min(ret_nnz,nnz)):
                data[i]    = jac.data[i]
                rowvals[i] = jac.indices[i]
            for i in range(dim+1):
                colptrs[i] = jac.indptr[i]
            
            return CVDLS_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
        except:
            traceback.print_exc()
            return CVDLS_JACFUNC_UNRECVR


IF SUNDIALS_VERSION >= (3,0,0):
    cdef int cv_jac(realtype t, N_Vector yv, N_Vector fy, SUNMatrix Jac, 
                void *problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
        """
        This method is used to connect the Assimulo.Problem.jac to the Sundials
        Jacobian function.
        """
        cdef SUNMatrixContent_Dense Jacobian = <SUNMatrixContent_Dense>Jac.content
        cdef ProblemData pData = <ProblemData>problem_data
        cdef realtype* col_i=Jacobian.cols[0]
        cdef N.ndarray y = pData.work_y
        cdef int i,j, Neq = pData.dim
        
        nv2arr_inplace(yv, y)

        if pData.dimSens>0: #Sensitivity activated
            p = realtype2arr(pData.p,pData.dimSens)
            try:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(t,y,sw=<list>pData.sw,p=p)
                else:
                    jac=(<object>pData.JAC)(t,y,p)
            except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
                return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
            except:
                traceback.print_exc()
                return CVDLS_JACFUNC_UNRECVR
        else:
            try:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(t,y,sw=<list>pData.sw)
                else:
                    jac=(<object>pData.JAC)(t,y)
            except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
                return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
            except:
                traceback.print_exc()
                return CVDLS_JACFUNC_UNRECVR
        
        if isinstance(jac, sparse.csc.csc_matrix):
            for j in range(Neq):
                col_i = Jacobian.cols[j]
                for i in range(jac.indptr[j], jac.indptr[j+1]):
                    col_i[jac.indices[i]] = jac.data[i]
        else:
            for i in range(Neq):
                col_i = Jacobian.cols[i]
                for j in range(Neq):
                    col_i[j] = jac[j,i]
        
        return CVDLS_SUCCESS
ELSE:
    cdef int cv_jac(long int Neq, realtype t, N_Vector yv, N_Vector fy, DlsMat Jacobian, 
                    void *problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
        """
        This method is used to connect the Assimulo.Problem.jac to the Sundials
        Jacobian function.
        """
        cdef ProblemData pData = <ProblemData>problem_data
        cdef realtype* col_i=DENSE_COL(Jacobian,0)
        cdef N.ndarray y = pData.work_y
        cdef int i,j
        
        nv2arr_inplace(yv, y)

        if pData.dimSens>0: #Sensitivity activated
            p = realtype2arr(pData.p,pData.dimSens)
            try:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(t,y,sw=<list>pData.sw,p=p)
                else:
                    jac=(<object>pData.JAC)(t,y,p)
            except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
                return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
            except:
                traceback.print_exc()
                return CVDLS_JACFUNC_UNRECVR
        else:
            try:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(t,y,sw=<list>pData.sw)
                else:
                    jac=(<object>pData.JAC)(t,y)
            except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
                return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
            except:
                traceback.print_exc()
                return CVDLS_JACFUNC_UNRECVR
                
        if isinstance(jac, sparse.csc.csc_matrix):
            for j in range(Neq):
                col_i = DENSE_COL(Jacobian, j)
                for i in range(jac.indptr[j], jac.indptr[j+1]):
                    col_i[jac.indices[i]] = jac.data[i]
        else:
            for i in range(Neq):
                col_i = DENSE_COL(Jacobian, i)
                for j in range(Neq):
                    col_i[j] = jac[j,i]
        
        return CVDLS_SUCCESS
        
        
cdef int cv_jacv(N_Vector vv, N_Vector Jv, realtype t, N_Vector yv, N_Vector fyv,
				    void *problem_data, N_Vector tmp):
    """
    This method is used to connect the Assimulo.Problem.jacv to the Sundials
    Jacobian times vector function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray y  = nv2arr(yv)
    cdef N.ndarray v  = nv2arr(vv)
    cdef N.ndarray fy = nv2arr(fyv)
    cdef int i
    
    cdef realtype* jacvptr=(<N_VectorContent_Serial>Jv.content).data
    
    if pData.dimSens>0: #Sensitivity activated
        p = realtype2arr(pData.p,pData.dimSens)
        try:
            if pData.sw != NULL:
                jacv = (<object>pData.JACV)(t,y,fy,v,sw=<list>pData.sw,p=p)
            else:
                jacv = (<object>pData.JACV)(t,y,fy,v,p=p)
            
            for i in range(pData.dim):
                jacvptr[i] = jacv[i]
            
            return SPGMR_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return SPGMR_ATIMES_FAIL_REC
        except:
            traceback.print_exc()
            return SPGMR_PSOLVE_FAIL_UNREC 
    else:
        try:
            if pData.sw != NULL:
                jacv = (<object>pData.JACV)(t,y,fy,v,sw=<list>pData.sw)
            else:
                jacv = (<object>pData.JACV)(t,y,fy,v)
            
            for i in range(pData.dim):
                jacvptr[i] = jacv[i]
            
            return SPGMR_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return SPGMR_ATIMES_FAIL_REC
        except:
            traceback.print_exc()
            return SPGMR_PSOLVE_FAIL_UNREC

IF SUNDIALS_VERSION >= (3,0,0):
    cdef int cv_prec_setup(realtype t, N_Vector yy, N_Vector fyy,
                      bint jok, bint *jcurPtr,
                      realtype gamma, void *problem_data):
        """
        For information see CVODES documentation 4.6.9
        """
        cdef ProblemData pData = <ProblemData>problem_data
        cdef N.ndarray y   = nv2arr(yy)
        cdef N.ndarray fy  = nv2arr(fyy)
        cdef object ret
        
        try:
            ret = (<object>pData.PREC_SETUP)(t,y,fy,jok,gamma,pData.PREC_DATA)
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)
        
        jcurPtr[0] = 1 if ret[0] else 0
        pData.PREC_DATA = ret[1]
        
        return CVSPILS_SUCCESS

    cdef int cv_prec_solve(realtype t, N_Vector yy, N_Vector fyy,
                      N_Vector rr, N_Vector z,
                      realtype gamma, realtype delta,
                      int lr, void *problem_data):
        """
        For information see CVODES documentation 4.6.8
        """
        cdef ProblemData pData = <ProblemData>problem_data
        cdef N.ndarray y   = nv2arr(yy)
        cdef N.ndarray r   = nv2arr(rr)
        cdef N.ndarray fy  = nv2arr(fyy)
        cdef realtype* zptr=(<N_VectorContent_Serial>z.content).data
        cdef int i

        try:
            zres = (<object>pData.PREC_SOLVE)(t,y,fy,r,gamma,delta,pData.PREC_DATA)
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)
                    
        for i in range(pData.dim):
            zptr[i] = zres[i]
        
        return CVSPILS_SUCCESS
ELSE:
    cdef int cv_prec_setup(realtype t, N_Vector yy, N_Vector fyy,
                      bint jok, bint *jcurPtr,
                      realtype gamma, void *problem_data,
                      N_Vector tmp1, N_Vector tmp2,
                      N_Vector tmp3):
        """
        For information see CVODES documentation 4.6.9
        """
        cdef ProblemData pData = <ProblemData>problem_data
        cdef N.ndarray y   = nv2arr(yy)
        cdef N.ndarray fy  = nv2arr(fyy)
        cdef object ret
        
        try:
            ret = (<object>pData.PREC_SETUP)(t,y,fy,jok,gamma,pData.PREC_DATA)
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)
        
        jcurPtr[0] = 1 if ret[0] else 0
        pData.PREC_DATA = ret[1]
        
        return CVSPILS_SUCCESS

    cdef int cv_prec_solve(realtype t, N_Vector yy, N_Vector fyy,
                      N_Vector rr, N_Vector z,
                      realtype gamma, realtype delta,
                      int lr, void *problem_data, N_Vector tmp):
        """
        For information see CVODES documentation 4.6.8
        """
        cdef ProblemData pData = <ProblemData>problem_data
        cdef N.ndarray y   = nv2arr(yy)
        cdef N.ndarray r   = nv2arr(rr)
        cdef N.ndarray fy  = nv2arr(fyy)
        cdef realtype* zptr=(<N_VectorContent_Serial>z.content).data
        cdef int i

        try:
            zres = (<object>pData.PREC_SOLVE)(t,y,fy,r,gamma,delta,pData.PREC_DATA)
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)
                    
        for i in range(pData.dim):
            zptr[i] = zres[i]
        
        return CVSPILS_SUCCESS

"""
cdef int cv_prec(realtype t, N Vector yv, N Vector fyv, 
         N Vector rv, N Vector z, realtype gamma, realtype delta, int lr, 
         void *problem_data, N Vector tmp):
    
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray y  = nv2arr(yv)
    cdef N.ndarray fy = nv2arr(fyv)
    cdef N.ndarray r  = nv2arr(rv)
    cdef int i
    
    cdef realtype* zptr=(<N_VectorContent_Serial>z.content).data
    
    try:
    
        zres = (<object>pData.PREC)(t,y,fy,r,...)
    
        for i in range(pData.dim):
            zptr[i] = zres[i]
        
        return SPGMR_SUCCESS
    except:
        return SPGMR_PSOLVE_FAIL_UNREC
"""

cdef int cv_root(realtype t, N_Vector yv, realtype *gout,  void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.state_events to the Sundials
    Root-finding function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray y = pData.work_y
    cdef int i
    
    nv2arr_inplace(yv, y)
    
    try:
        if pData.sw != NULL:
            root=(<object>pData.ROOT)(t,y,<list>pData.sw) #Call to the Python root function 
        else:
            root=(<object>pData.ROOT)(t,y,None) #Call to the Python root function
            
        #memcpy(gout,<realtype*>root.data,pData.memSizeRoot) #Copy data from the return to the output
        for i in range(pData.dimRoot):
            gout[i]=root[i]
    
        return CV_SUCCESS
    except:
        return CV_RTFUNC_FAIL  # Unrecoverable Error

cdef int ida_res(realtype t, N_Vector yv, N_Vector yvdot, N_Vector residual, void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.f to the Sundials
    residual function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray[realtype, ndim=1, mode='c'] res #Used for return from the user function
    cdef N.ndarray y = pData.work_y
    cdef N.ndarray yd = pData.work_yd
    cdef realtype* resptr=(<N_VectorContent_Serial>residual.content).data
    cdef int i
    
    nv2arr_inplace(yv, y)
    nv2arr_inplace(yvdot, yd)
    
    if pData.dimSens!=0: #SENSITIVITY 
        p = realtype2arr(pData.p,pData.dimSens)
        try:
            if pData.sw != NULL:
                res=(<object>pData.RHS)(t,y,yd,sw=<list>pData.sw,p=p)  # call to the python residual function
            else:
                res=(<object>pData.RHS)(t,y,yd,p)
            
            #memcpy((<N_VectorContent_Serial>residual.content).data,<realtype*>res.data,pData.memSize)
            for i in range(pData.dim):
                resptr[i] = res[i]

            return IDA_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return IDA_REC_ERR # recoverable error (see Sundials description)
        except:
            traceback.print_exc()
            return IDA_RES_FAIL
    else: #NO SENSITIVITY
        try:
            if pData.sw != NULL:
                res=(<object>pData.RHS)(t,y,yd,<list>pData.sw)  #Call to the Python residual function
            else:
                res=(<object>pData.RHS)(t,y,yd)
                #res = (<object>pData.RHS)(t,y,yd)
            
            #memcpy((<N_VectorContent_Serial>residual.content).data,<realtype*>res.data,pData.memSize)
            for i in range(pData.dim):
                resptr[i] = res[i]
            
            return IDA_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return IDA_REC_ERR # recoverable error (see Sundials description)
        except:
            traceback.print_exc()
            return IDA_RES_FAIL

IF SUNDIALS_VERSION >= (3,0,0):
    cdef int ida_jac(realtype t, realtype c, N_Vector yv, N_Vector yvdot, N_Vector residual, SUNMatrix Jac,
                 void *problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
        """
        This method is used to connect the Assimulo.Problem.jac to the Sundials
        Jacobian function.
        """
        cdef SUNMatrixContent_Dense Jacobian = <SUNMatrixContent_Dense>Jac.content
        cdef ProblemData pData = <ProblemData>problem_data
        cdef N.ndarray[realtype, ndim=2, mode='c'] jac #Used for return from the user function
        cdef realtype* col_i=Jacobian.cols[0]
        cdef N.ndarray y = pData.work_y
        cdef N.ndarray yd = pData.work_yd
        cdef int i,j, Neq = pData.dim
        
        nv2arr_inplace(yv, y)
        nv2arr_inplace(yvdot, yd)
        
        if pData.dimSens!=0: #SENSITIVITY 
            p = realtype2arr(pData.p,pData.dimSens)
            try:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(c,t,y,yd,sw=<list>pData.sw,p=p)  # call to the python residual function
                else:
                    jac=(<object>pData.JAC)(c,t,y,yd,p=p)
                
                for i in range(Neq):
                    col_i = Jacobian.cols[i]
                    for j in range(Neq):
                        col_i[j] = jac[j,i]
                return IDADLS_SUCCESS
            except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
                return IDADLS_JACFUNC_RECVR #Recoverable Error
            except:
                traceback.print_exc()
                return IDADLS_JACFUNC_UNRECVR
        else:
            try:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(c,t,y,yd,<list>pData.sw)  # call to the python residual function
                else:
                    jac=(<object>pData.JAC)(c,t,y,yd)
                
                for i in range(Neq):
                    col_i = Jacobian.cols[i]
                    for j in range(Neq):
                        col_i[j] = jac[j,i]
                return IDADLS_SUCCESS
            except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
                return IDADLS_JACFUNC_RECVR #Recoverable Error
            except:
                traceback.print_exc()
                return IDADLS_JACFUNC_UNRECVR
ELSE:
    cdef int ida_jac(long int Neq, realtype t, realtype c, N_Vector yv, N_Vector yvdot, N_Vector residual, DlsMat Jacobian,
                 void* problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
        """
        This method is used to connect the Assimulo.Problem.jac to the Sundials
        Jacobian function.
        """
        cdef ProblemData pData = <ProblemData>problem_data
        cdef N.ndarray[realtype, ndim=2, mode='c'] jac #Used for return from the user function
        cdef realtype* col_i=DENSE_COL(Jacobian,0)
        cdef N.ndarray y = pData.work_y
        cdef N.ndarray yd = pData.work_yd
        cdef int i,j
        
        nv2arr_inplace(yv, y)
        nv2arr_inplace(yvdot, yd)
        
        if pData.dimSens!=0: #SENSITIVITY 
            p = realtype2arr(pData.p,pData.dimSens)
            try:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(c,t,y,yd,sw=<list>pData.sw,p=p)  # call to the python residual function
                else:
                    jac=(<object>pData.JAC)(c,t,y,yd,p=p)
                
                for i in range(Neq):
                    col_i = DENSE_COL(Jacobian, i)
                    for j in range(Neq):
                        col_i[j] = jac[j,i]
                return IDADLS_SUCCESS
            except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
                return IDADLS_JACFUNC_RECVR #Recoverable Error
            except:
                traceback.print_exc()
                return IDADLS_JACFUNC_UNRECVR
        else:
            try:
                if pData.sw != NULL:
                    jac=(<object>pData.JAC)(c,t,y,yd,<list>pData.sw)  # call to the python residual function
                else:
                    jac=(<object>pData.JAC)(c,t,y,yd)
                
                for i in range(Neq):
                    col_i = DENSE_COL(Jacobian, i)
                    for j in range(Neq):
                        col_i[j] = jac[j,i]
                return IDADLS_SUCCESS
            except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
                return IDADLS_JACFUNC_RECVR #Recoverable Error
            except:
                traceback.print_exc()
                return IDADLS_JACFUNC_UNRECVR
            

cdef int ida_root(realtype t, N_Vector yv, N_Vector yvdot, realtype *gout,  void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.state_events to the Sundials
    root function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray[realtype, ndim=1, mode='c'] root #Used for return from the user function
    cdef N.ndarray y = pData.work_y
    cdef N.ndarray yd = pData.work_yd
    cdef int i
    
    nv2arr_inplace(yv, y)
    nv2arr_inplace(yvdot, yd)
    
    try:
        if pData.sw != NULL:
            root=(<object>pData.ROOT)(t,y,yd,<list>pData.sw)  #Call to the Python root function
        else:
            root=(<object>pData.ROOT)(t,y,yd,None)  #Call to the Python root function
    
        #memcpy(gout,<realtype*>root.data,pData.memSizeRoot) #Copy data from the return to the output
        for i in range(pData.dimRoot):
            gout[i]=root[i]
        
        return IDA_SUCCESS
    except:
        return IDA_RTFUNC_FAIL  # Unrecoverable Error

cdef int ida_jacv(realtype t, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector vv, N_Vector Jv, realtype cj,
				    void *problem_data, N_Vector tmp1, N_Vector tmp2):
    """
    This method is used to connect the Assimulo.Problem.jacv to the Sundials
    Jacobian times vector function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray y  = nv2arr(yy)
    cdef N.ndarray yd = nv2arr(yp)
    cdef N.ndarray v  = nv2arr(vv)
    cdef N.ndarray res = nv2arr(rr)
    cdef int i
    
    cdef realtype* jacvptr=(<N_VectorContent_Serial>Jv.content).data
    
    if pData.dimSens>0: #Sensitivity activated
        p = realtype2arr(pData.p,pData.dimSens)
        try:
            if pData.sw != NULL:
                jacv = (<object>pData.JACV)(t,y,yd,res,v,cj,sw=<list>pData.sw,p=p)
            else:
                jacv = (<object>pData.JACV)(t,y,yd,res,v,cj,p=p)
        
            for i in range(pData.dim):
                jacvptr[i] = jacv[i]
            
            return SPGMR_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return SPGMR_ATIMES_FAIL_REC
        except:
            traceback.print_exc()
            return SPGMR_PSOLVE_FAIL_UNREC 
    else:
        try:
            if pData.sw != NULL:
                jacv = (<object>pData.JACV)(t,y,yd,res,v,cj,sw=<list>pData.sw)
            else:
                jacv = (<object>pData.JACV)(t,y,yd,res,v,cj)
            
            for i in range(pData.dim):
                jacvptr[i] = jacv[i]
            
            return SPGMR_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return SPGMR_ATIMES_FAIL_REC
        except:
            traceback.print_exc()
            return SPGMR_PSOLVE_FAIL_UNREC

# Error handling callback functions
# =================================

cdef void cv_err(int error_code, const char *module, const char *function, char *msg, void *problem_data):
    """
    This method overrides the default handling of error messages.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    
    if error_code > 0 and pData.verbose > 0: #Warning
        print '[CVode Warning]', msg
    
    if pData.verbose > 2: #Verbosity is greater than NORMAL, print warnings and errors
        if error_code < 0: #Error
            print '[CVode Error]', msg
            
cdef void ida_err(int error_code, const char *module, const char *function, char *msg, void *problem_data):
    """
    This method overrides the default handling of error messages.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    
    if error_code > 0 and pData.verbose > 0: #Warning
        print '[IDA Warning]', msg
    
    if pData.verbose > 2: #Verbosity is greater than NORMAL, print warnings and errors
        if error_code < 0: #Error
            print '[IDA Error]', msg


cdef class ProblemData:
    cdef:
        void *RHS          #Should store the residual or the right-hand-side
        void *RHS_SENS_ALL #Should store the sensitivty equation of all parameters
        void *ROOT         #Should store the root function
        void *JAC          #Should store the jacobian
        void *JACV         #Should store the jacobian times a vector
        void *SENS         #Should store the sensitivity function
        void *PREC_SOLVE   #Should store the preconditioner solve function
        void *PREC_SETUP   #Should store the preconditioner setup function
        void *y            #Temporary storage for the states
        void *yd           #Temporary storage for the derivatives
        void *sw           #Storage for the switches
        realtype *p            #Storage for the parameters
        realtype *pbar
        int dim            #Dimension of the problem
        int dimRoot        #Dimension of the roots
        int dimSens        #Dimension of the parameters (For sensitivity)
        int memSize        #dim*sizeof(realtype) used when copying memory
        int memSizeRoot    #dimRoot*sizeof(realtype) used when copying memory
        int memSizeJac     #dim*dim*sizeof(realtype) used when copying memory
        int verbose        #Defines the verbosity
        object PREC_DATA   #Arbitrary data from the preconditioner
        N.ndarray work_y
        N.ndarray work_yd
        N.ndarray work_ys
        
    cdef create_work_arrays(self):
        self.work_y = N.empty(self.dim)
        self.work_yd = N.empty(self.dim)
        self.work_ys = N.empty((self.dim, self.dimSens))
        
