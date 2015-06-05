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
    #cdef ndarray[realtype, ndim=1, mode='c'] rhs #Used for return from the user function
    #(<ndarray>pData.y).data =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
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
                
            #memcpy((<N_VectorContent_Serial>yvdot.content).data,<realtype*>rhs.data,pData.memSize)
            for i in range(pData.dim):
                resptr[i] = rhs[i]
            
            return CV_SUCCESS
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)
        
    else: #No sensitivity
        try:
            if pData.sw != NULL:
                rhs = (<object>pData.RHS)(t,y,<list>pData.sw)
            else:
                rhs = (<object>pData.RHS)(t,y)
            
            #memcpy((<N_VectorContent_Serial>yvdot.content).data,<realtype*>rhs.data,pData.memSize)
            for i in range(pData.dim):
                resptr[i] = rhs[i]
            
            return CV_SUCCESS
        except:
            return CV_REC_ERR #Recoverable Error (See Sundials description)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int cv_jac_sparse(realtype t, N_Vector yv, N_Vector fy, SlsMat Jacobian,
                void *problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    This method is used to connect the Assimulo.Problem.jac to the Sundials
    Sparse Jacobian function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray y = nv2arr(yv)
    cdef int i
    cdef int nnz = Jacobian.NNZ
    cdef int ret_nnz
    cdef int dim = Jacobian.N
    cdef realtype* data = Jacobian.data
    cdef int* rowvals = Jacobian.rowvals
    cdef int* colptrs = Jacobian.colptrs

    """
        realtype *data;
        int *rowvals;
        int *colptrs;
    """
    if pData.dimSens>0: #Sensitivity activated
        raise Exception("Not Suppported!")
    else:
        try:
            if pData.sw != NULL:
                jac=(<object>pData.JAC)(t,y,sw=<list>pData.sw)
            else:
                jac=(<object>pData.JAC)(t,y)
                
            if not isinstance(jac, sparse.csc.csc_matrix):
                jac = sparse.csc.csc_matrix(jac)
                raise AssimuloException("The Jacobian must be stored on Scipy's CSC format.")
            ret_nnz = jac.nnz
            if ret_nnz> nnz:
                raise AssimuloException("The Jacobian has more entries than supplied to the problem class via 'jac_nnz'")    
                
            for i in range(min(ret_nnz,nnz)):
                data[i]    = jac.data[i]
                rowvals[i] = jac.indices[i]
            for i in range(dim+1):
                colptrs[i] = jac.indptr[i]
            
            return CVDLS_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError):
            return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
        except:
            traceback.print_exc()
            return CVDLS_JACFUNC_UNRECVR

cdef int cv_jac(int Neq, realtype t, N_Vector yv, N_Vector fy, DlsMat Jacobian, 
                void *problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    This method is used to connect the Assimulo.Problem.jac to the Sundials
    Jacobian function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    #cdef ndarray[realtype, ndim=2, mode='c'] jac #Used for return from the user function
    cdef realtype* col_i=DENSE_COL(Jacobian,0)
    #(<ndarray>pData.y).data =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef int i,j

    if pData.dimSens>0: #Sensitivity activated
        p = realtype2arr(pData.p,pData.dimSens)
        try:
            if pData.sw != NULL:
                jac=(<object>pData.JAC)(t,y,sw=<list>pData.sw,p=p)
            else:
                jac=(<object>pData.JAC)(t,y,p)
                
            for i in range(Neq):
                col_i = DENSE_COL(Jacobian, i)
                for j in range(Neq):
                    col_i[j] = jac[j,i]

            return CVDLS_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError):
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
    
            for i in range(Neq):
                col_i = DENSE_COL(Jacobian, i)
                for j in range(Neq):
                    col_i[j] = jac[j,i]

            return CVDLS_SUCCESS
        except(N.linalg.LinAlgError,ZeroDivisionError):
            return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
        except:
            traceback.print_exc()
            return CVDLS_JACFUNC_UNRECVR
        
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
        except(N.linalg.LinAlgError,ZeroDivisionError):
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
        except(N.linalg.LinAlgError,ZeroDivisionError):
            return SPGMR_ATIMES_FAIL_REC
        except:
            traceback.print_exc()
            return SPGMR_PSOLVE_FAIL_UNREC
            
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
    #cdef ndarray[realtype, ndim=1, mode='c'] root #Used for return from the user function
    #(<ndarray>pData.y).data =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #cdef N.ndarray y = nv2arr(yv)
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
    #(<ndarray>pData.y).data  =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #(<ndarray>pData.yd).data =  <realtype*>((<N_VectorContent_Serial>yvdot.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef N.ndarray yd = nv2arr(yvdot)
    cdef realtype* resptr=(<N_VectorContent_Serial>residual.content).data
    cdef int i
    
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
        except(N.linalg.LinAlgError,ZeroDivisionError):
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
        except(N.linalg.LinAlgError,ZeroDivisionError):
            return IDA_REC_ERR # recoverable error (see Sundials description)
        except:
            traceback.print_exc()
            return IDA_RES_FAIL
            
cdef int ida_jac(int Neq, realtype t, realtype c, N_Vector yv, N_Vector yvdot, N_Vector residual, DlsMat Jacobian,
                 void* problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    This method is used to connect the Assimulo.Problem.jac to the Sundials
    Jacobian function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef N.ndarray[realtype, ndim=2, mode='c'] jac #Used for return from the user function
    cdef realtype* col_i=DENSE_COL(Jacobian,0)
    #(<ndarray>pData.y).data  =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #(<ndarray>pData.yd).data =  <realtype*>((<N_VectorContent_Serial>yvdot.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef N.ndarray yd = nv2arr(yvdot)
    cdef int i,j
    
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
        except(N.linalg.LinAlgError,ZeroDivisionError):
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
        except(N.linalg.LinAlgError,ZeroDivisionError):
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
    #(<ndarray>pData.y).data  =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #(<ndarray>pData.yd).data =  <realtype*>((<N_VectorContent_Serial>yvdot.content).data)
    cdef N.ndarray y = nv2arr(yv)
    cdef N.ndarray yd = nv2arr(yvdot)
    cdef int i
    
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
        except(N.linalg.LinAlgError,ZeroDivisionError):
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
        except(N.linalg.LinAlgError,ZeroDivisionError):
            return SPGMR_ATIMES_FAIL_REC
        except:
            traceback.print_exc()
            return SPGMR_PSOLVE_FAIL_UNREC
    

cdef int kin_jac(int Neq, N_Vector xv, N_Vector fval, DlsMat Jacobian, 
                void *problem_data, N_Vector tmp1, N_Vector tmp2):
    """
    This method is used to connect the assimulo.Problem.jac to the Sundials
    Jacobian function.
    """
    cdef ProblemDataEquationSolver pData = <ProblemDataEquationSolver>problem_data
    cdef realtype* col_i=DENSE_COL(Jacobian,0)
    cdef N.ndarray x = nv2arr(xv)
    cdef int i,j
    
    try:
        jac=(<object>pData.JAC)(x)

        for i in range(Neq):
            col_i = DENSE_COL(Jacobian, i)
            for j in range(Neq):
                col_i[j] = jac[j,i]

        return KINDLS_SUCCESS
    except:
        return KINDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
        
cdef int kin_jacv(N_Vector vv, N_Vector Jv, N_Vector vx, bint new_u,
            void *problem_data):
    cdef ProblemDataEquationSolver pData = <ProblemDataEquationSolver>problem_data
    cdef N.ndarray x  = nv2arr(vx)
    cdef N.ndarray v  = nv2arr(vv)
    cdef int i
    
    cdef realtype* jacvptr=(<N_VectorContent_Serial>Jv.content).data

    try:
        jacv = (<object>pData.JACV)(x,v)
        
        for i in range(pData.dim):
            jacvptr[i] = jacv[i]
        
        return SPGMR_SUCCESS
    except(N.linalg.LinAlgError,ZeroDivisionError, AssimuloRecoverableError):
        return SPGMR_ATIMES_FAIL_REC
    except:
        traceback.print_exc()
        return SPGMR_PSOLVE_FAIL_UNREC 
    
cdef int kin_res(N_Vector xv, N_Vector fval, void *problem_data):
    """
    Residual fct called by KINSOL
    """
    cdef ProblemDataEquationSolver pData = <ProblemDataEquationSolver>problem_data
    cdef N.ndarray x = nv2arr(xv)
    cdef realtype* resptr = (<N_VectorContent_Serial>fval.content).data
    cdef int i

    try:
        res = (<object>pData.RES)(x)

        for i in range(pData.dim):
            resptr[i] = res[i]

        return KIN_SUCCESS
    except(N.linalg.LinAlgError,ZeroDivisionError, AssimuloRecoverableError):
        return KIN_REC_ERR
    except:
        traceback.print_exc()
        return KIN_SYSFUNC_FAIL

cdef int kin_prec_solve(N_Vector u, N_Vector uscaleN, N_Vector fval, 
         N_Vector fscaleN, N_Vector v, void *problem_data, N_Vector tmp):
    """
    Preconditioning solve function
    
        Pz = r
        
        v on input == r
        v on output == z
    """
    cdef ProblemDataEquationSolver pData = <ProblemDataEquationSolver>problem_data
    
    cdef N.ndarray fscale  = nv2arr(fscaleN)
    cdef N.ndarray uscale  = nv2arr(uscaleN)
    cdef N.ndarray r       = nv2arr(v)
    cdef realtype* zptr=(<N_VectorContent_Serial>v.content).data
    
    try:
        zres = (<object>pData.PREC_SOLVE)(r)
    except(N.linalg.LinAlgError,ZeroDivisionError, AssimuloRecoverableError):
        return KIN_REC_ERR
    except:
        traceback.print_exc()
        return KIN_SYSFUNC_FAIL
                
    for i in range(pData.dim):
        zptr[i] = zres[i]
    
    return KIN_SUCCESS
    
cdef int kin_prec_setup(N_Vector uN, N_Vector uscaleN, N_Vector fvalN, 
         N_Vector fscaleN, void *problem_data, N_Vector tmp1, N_Vector tmp2):
    """
    Preconditioning setup function
    """
    cdef ProblemDataEquationSolver pData = <ProblemDataEquationSolver>problem_data
    
    cdef N.ndarray fscale  = nv2arr(fscaleN)
    cdef N.ndarray uscale  = nv2arr(uscaleN)
    cdef N.ndarray u       = nv2arr(uN)
    cdef N.ndarray fval    = nv2arr(fvalN)
    
    try:
        (<object>pData.PREC_SETUP)(u, fval, uscale, fscale)
    except(N.linalg.LinAlgError,ZeroDivisionError, AssimuloRecoverableError):
        return KIN_REC_ERR
    except:
        traceback.print_exc()
        return KIN_SYSFUNC_FAIL
    
    return KIN_SUCCESS
    

cdef void kin_err(int err_code, char *module, char *function, char *msg, void *eh_data):
    cdef ProblemDataEquationSolver pData = <ProblemDataEquationSolver>eh_data
    
    if err_code > 0: #Warning
        category = 1
    elif err_code < 0: #Error
        category = -1
    else:
        category = 0
    
    print "Error occured in <function: %s>."%function
    print "<message: %s>"%msg
    #print "<functionNorm: %g, scaledStepLength: %g, tolerance: %g>"%(fnorm, snorm, pData.TOL)


cdef void kin_info(char *module, char *function, char *msg, void *eh_data):
    cdef ProblemDataEquationSolver pData = <ProblemDataEquationSolver>eh_data
    cdef int flag
    cdef realtype fnorm
    
    if str(function) == "KINSol" and "fnorm" in str(msg):
        #fnorm = float(msg.split("fnorm = ")[-1].strip())
        flag = SUNDIALS.KINGetFuncNorm(pData.KIN_MEM, &fnorm)
        pData.nl_fnorm.append(fnorm)
        
    pData.log.append([module, function, msg])
    
    #print "KinsolInfo <calling_function:%s>"%function
    #print "<message: %s>"%msg
    """
    # Get the number of iterations
    KINGetNumNonlinSolvIters(kin_mem, &nniters)
    
    
    /* Only output an iteration under certain conditions:
     *  1. nle_solver_log > 2
     *  2. The calling function is either KINSolInit or KINSol
     *  3. The message string starts with "nni"
     *
     *  This approach gives one printout per iteration
    
    
    if ("KINSolInit" in function or "KINSol" in function) and "nni" in msg:
        print "<iteration_index:%d>"%nniters
        print "ivs", N_VGetArrayPointer(kin_mem->kin_uu), block->n);
        print "<scaled_residual_norm:%E>", kin_mem->kin_fnorm);
        print "residuals", 
            realtype* f = N_VGetArrayPointer(kin_mem->kin_fval);
            f[i]*residual_scaling_factors[i]
    """

# Error handling callback functions
# =================================

cdef int cv_err(int error_code, char *module, char *function, char *msg, void *problem_data):
    """
    This method overrides the default handling of error messages.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    
    if error_code > 0 and pData.verbose > 0: #Warning
        print '[CVode Warning]', msg
    
    if pData.verbose > 2: #Verbosity is greater than NORMAL, print warnings and errors
        if error_code < 0: #Error
            print '[CVode Error]', msg
            
cdef int ida_err(int error_code, char *module, char *function, char *msg, void *problem_data):
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
        
    cdef create_work_arrays(self):
        self.work_y = N.empty(self.dim)
        self.work_yd = N.empty(self.dim)
        

cdef class ProblemDataEquationSolver:
    cdef:
        void *RES          # Residual
        void *JAC          # Jacobian
        void *JACV
        void *PREC_SOLVE
        void *PREC_SETUP
        int dim            # Dimension of the problem
        void *KIN_MEM      # Kinsol memory
        list nl_fnorm      # The norm of the residual at each nonlinear iteration (if the verbosity is set high enough)
        list l_fnorm
        list log

#=================
# Module functions
#=================

cdef inline N_Vector arr2nv(x):
    x=N.array(x)
    cdef long int n = len(x)
    cdef N.ndarray[realtype, ndim=1,mode='c'] ndx=x
    cdef void* data_ptr=PyArray_DATA(ndx)
    cdef N_Vector v=N_VNew_Serial(n)
    memcpy((<N_VectorContent_Serial>v.content).data, data_ptr, n*sizeof(realtype))
    return v
    
cdef inline void arr2nv_inplace(x, N_Vector out):
    x=N.array(x)
    cdef long int n = len(x)
    cdef N.ndarray[realtype, ndim=1,mode='c'] ndx=x
    cdef void* data_ptr=PyArray_DATA(ndx)
    memcpy((<N_VectorContent_Serial>out.content).data, data_ptr, n*sizeof(realtype))
    
cdef inline N.ndarray nv2arr(N_Vector v):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    cdef N.ndarray[realtype, ndim=1, mode='c'] x=N.empty(n)
    memcpy(x.data, v_data, n*sizeof(realtype))
    return x
    
cdef inline void nv2arr_inplace(N_Vector v, N.ndarray o):
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    memcpy(o.data, v_data, n*sizeof(realtype))

cdef inline realtype2arr(realtype *data, int n):
    """Create new numpy array from realtype*"""
    cdef N.ndarray[realtype, ndim=1, mode='c'] x=N.empty(n)
    memcpy(x.data, data, n*sizeof(realtype))
    return x
