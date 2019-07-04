#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2017 Modelon AB
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

IF SUNDIALS_VERSION >= (3,0,0):
    cdef int kin_jac(N_Vector xv, N_Vector fval, SUNMatrix Jac, 
                    void *problem_data, N_Vector tmp1, N_Vector tmp2):
        """
        This method is used to connect the assimulo.Problem.jac to the Sundials
        Jacobian function.
        """
        cdef SUNMatrixContent_Dense Jacobian = <SUNMatrixContent_Dense>Jac.content
        cdef ProblemDataEquationSolver pData = <ProblemDataEquationSolver>problem_data
        cdef realtype* col_i=Jacobian.cols[0]
        cdef N.ndarray x = nv2arr(xv)
        cdef int i,j, Neq = pData.dim
        
        try:
            jac=(<object>pData.JAC)(x)

            for i in range(Neq):
                col_i = Jacobian.cols[i]
                for j in range(Neq):
                    col_i[j] = jac[j,i]

            return KINDLS_SUCCESS
        except:
            return KINDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
ELSE:
    cdef int kin_jac(long int Neq, N_Vector xv, N_Vector fval, DlsMat Jacobian, 
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
            
cdef int kin_jacv(N_Vector vv, N_Vector Jv, N_Vector vx, int* new_u,
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
    except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
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
    except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
        return KIN_REC_ERR
    except:
        traceback.print_exc()
        return KIN_SYSFUNC_FAIL


IF SUNDIALS_VERSION >= (3,0,0):
    cdef int kin_prec_solve(N_Vector u, N_Vector uscaleN, N_Vector fval, 
             N_Vector fscaleN, N_Vector v, void *problem_data):
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
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return KIN_REC_ERR
        except:
            traceback.print_exc()
            return KIN_SYSFUNC_FAIL
                    
        for i in range(pData.dim):
            zptr[i] = zres[i]
        
        return KIN_SUCCESS
        
    cdef int kin_prec_setup(N_Vector uN, N_Vector uscaleN, N_Vector fvalN, 
             N_Vector fscaleN, void *problem_data):
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
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return KIN_REC_ERR
        except:
            traceback.print_exc()
            return KIN_SYSFUNC_FAIL
        
        return KIN_SUCCESS
        

ELSE:
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
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
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
        except(N.linalg.LinAlgError,ZeroDivisionError,AssimuloRecoverableError):
            return KIN_REC_ERR
        except:
            traceback.print_exc()
            return KIN_SYSFUNC_FAIL
        
        return KIN_SUCCESS
        

cdef void kin_err(int err_code, const char *module, const char *function, char *msg, void *eh_data):
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


cdef void kin_info(const char *module, const char *function, char *msg, void *eh_data):
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
