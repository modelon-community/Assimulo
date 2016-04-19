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

import numpy as N 
cimport numpy as N
from numpy cimport PyArray_DATA

N.import_array()

import numpy.linalg
import traceback 
 
from assimulo.exception import * 
from assimulo.algebraic cimport Algebraic

cimport sundials_includes as SUNDIALS

#Various C includes transfered to namespace
from sundials_includes cimport N_Vector, realtype, N_VectorContent_Serial, DENSE_COL
from sundials_includes cimport memcpy, N_VNew_Serial, DlsMat, SlsMat
from sundials_includes cimport malloc, free, realtype, N_VCloneVectorArray_Serial
from sundials_includes cimport N_VConst_Serial, N_VDestroy_Serial

include "constants.pxi" #Includes the constants (textual include)
include "../lib/sundials_constants.pxi" #Sundials related constants
include "../lib/sundials_callbacks.pxi"

cdef class KINSOL(Algebraic):
    """
    This class provides a connection to the Sundials 
    (https://computation.llnl.gov/casc/sundials/main.html) solver 
    KINSOL.
    
    .. math::
    
        0 = F(y)
    
    """
    cdef void* kinsol_mem
    cdef ProblemDataEquationSolver pData #A struct containing information about the problem
    cdef N_Vector y_temp, y_scale, f_scale
    cdef public double _eps
    
    cdef object pt_fcn, pt_jac, pt_jacv, pt_prec_setup, pt_prec_solve
    cdef object _added_linear_solver
    
    def __init__(self, problem):
        Algebraic.__init__(self, problem) #Calls the base class

        self.pData = ProblemDataEquationSolver()
        
        self._eps = N.finfo('double').eps
        self._added_linear_solver = False
        
        #Populate the ProblemData
        self.set_problem_data()
        
        #Solver options
        self.options["ftol"] = self._eps**(1.0/3.0)
        self.options["stol"] = self._eps**(2.0/3.0)
        self.options["strategy"] = KIN_LINESEARCH
        self.options["y_scale"] = N.array([1.0]*self.problem_info["dim"])
        self.options["f_scale"] = N.array([1.0]*self.problem_info["dim"]) 
        #self.options["y_nominal"] = N.array([1.0]*self.problem_info["dim"])
        #self.options["y_min"] = N.array([MIN_VALUE]*self.problem_info["dim"])
        #self.options["y_max"] = N.array([MAX_VALUE]*self.problem_info["dim"])
        self.options["linear_solver"] = "DENSE"
        self.options["max_iter"] = 200 #Maximum number of nonlinear iterations
        self.options["no_initial_setup"] = False #Specifies wheter or not a call to the setup function should be made
        self.options["max_solves_between_setup_calls"] = 10 #Specifies the maximum number of allowed calls to the solve function between setup calls
        self.options["max_newton_step"] = 0.0 #Specifies the max allowed newton step
        self.options["no_min_epsilon"] = False #Specifies wheter the scaled linear residual is bounded from below
        self.options["max_beta_fails"] = 10
        self.options["max_krylov"] = 0
        
        #Statistics
        self.statistics["nfevals"]    = 0 #Function evaluations
        self.statistics["nniters"]    = 0 #Nonlinear iterations
        self.statistics["njevals"]    = 0 #Jacobian evaluations
        self.statistics["nfevalsLS"]  = 0 #Function evaluations due to Jac
        self.statistics["nbacktr"]    = 0 #The function KINGetNumBacktrackOps returns the number of backtrack operations (step length adjustments) performed by the line search algorithm.
        self.statistics["nbcfails"]   = 0 #The function KINGetNumBetaCondFails returns the number of β-condition failures.
        self.statistics["nliters"]    = 0
        self.statistics["nlcfails"]   = 0
        
        #Initialize Kinsol
        self.initialize_kinsol()
                
    def __dealloc__(self):
        
        if self.y_temp != NULL:
            #Deallocate N_Vector
            N_VDestroy_Serial(self.y_temp)
        
        if self.kinsol_mem != NULL:
            #Free Memory
            SUNDIALS.KINFree(&self.kinsol_mem)
        
    def update_variable_scaling(self, value="Automatic"):
        """
        Updates the variable scaling with either 
        """
        if isinstance(value, str) and value.upper() == "AUTOMATIC":
            for i in range(self.problem_info["dim"]):
                if self.options["y_nominal"]:
                    self.options["y_scale"][i] = self.options["y_nominal"][i]
                elif self.y[i] != 0.0:
                    self.options["y_scale"][i] = self.y[i]
                elif self.options["y_max"] and self.options["y_min"]:
                    self.options["y_scale"][i] = max(abs(self.options["y_max"][i]+self.options["y_min"][i])/2.0, 0.01*self.options["y_max"][i])
                elif self.options["y_max"]:
                    self.options["y_scale"][i] = max(1.0, abs(self.options["y_max"][i]))
                elif self.options["y_min"]:
                    self.options["y_scale"][i] = max(1.0, abs(self.options["y_min"][i]))
        else:
            self.options["y_scale"] = N.array([value]) if isinstance(value, float) or isinstance(value, int) else N.array(value)
            
        arr2nv_inplace(self.options["y_scale"], self.y_scale)
    
    def update_residual_scaling(self, value="Automatic"):
        """
        Updates the residual scaling.
        """
        if isinstance(value, str) and value.upper() == "AUTOMATIC":
            pass
        else:
            self.options["f_scale"] = N.array([value]) if isinstance(value, float) or isinstance(value, int) else N.array(value)
            
        arr2nv_inplace(self.options["f_scale"], self.f_scale)
            
    cdef set_problem_data(self):
        
        #Sets the f function
        self.pt_fcn = self.problem.res
        self.pData.RES = <void*>self.pt_fcn
        self.pData.dim = self.problem_info["dim"]
        self.pData.nl_fnorm = []
        self.pData.l_fnorm = []
        self.pData.log = []
        
        if self.problem_info["jac_fcn"] is True: #Sets the jacobian
            self.pt_jac = self.problem.jac
            self.pData.JAC = <void*>self.pt_jac
        
        if self.problem_info["jacv_fcn"] is True: #Sets the jacobian*vector function
            self.pt_jacv = self.problem.jacv
            self.pData.JACV = <void*>self.pt_jacv
            
        if self.problem_info["prec_solve"] is True: #Sets the preconditioner solve function
            self.pt_prec_solve = self.problem.prec_solve
            self.pData.PREC_SOLVE = <void*>self.pt_prec_solve
            
        if self.problem_info["prec_setup"] is True: #Sets the preconditioner setup function
            self.pt_prec_setup = self.problem.prec_setup
            self.pData.PREC_SETUP = <void*>self.pt_prec_setup
            
    cdef initialize_kinsol(self):
        cdef int flag #Used for return
        
        self.y_temp  = arr2nv(self.y)
        self.y_scale = arr2nv([1.0]*self.problem_info["dim"])
        self.f_scale = arr2nv([1.0]*self.problem_info["dim"])
   
        if self.kinsol_mem == NULL: #The solver is not initialized
            
            #Create the solver
            self.kinsol_mem = SUNDIALS.KINCreate()
            if self.kinsol_mem == NULL:
                raise KINSOLError(KIN_MEM_NULL)
            self.pData.KIN_MEM = self.kinsol_mem
            
            #Specify the residual and the initial conditions to the solver
            flag = SUNDIALS.KINInit(self.kinsol_mem, kin_res, self.y_temp)
            if flag < 0:
                raise KINSOLError(flag)
            
            #Specify the error handling
            flag = SUNDIALS.KINSetErrHandlerFn(self.kinsol_mem, kin_err, <void*>self.pData)
            if flag < 0:
                raise KINSOLError(flag) 
                
            #Specify the handling of info messages
            flag = SUNDIALS.KINSetInfoHandlerFn(self.kinsol_mem, kin_info, <void*>self.pData);
            if flag < 0:
                raise KINSOLError(flag)
            
        else: #The solver needs not to be reinitialized
            pass
            
        #Set the user data
        flag = SUNDIALS.KINSetUserData(self.kinsol_mem, <void*>self.pData)
        if flag < 0:
            raise KINSOLError(flag)
            
    cpdef add_linear_solver(self):
        if self.options["linear_solver"] == "DENSE":
            flag = SUNDIALS.KINDense(self.kinsol_mem, self.problem_info["dim"])
            if flag < 0:
                raise KINSOLError(flag)
            
            if self.problem_info["jac_fcn"]:
                flag = SUNDIALS.KINDlsSetDenseJacFn(self.kinsol_mem, kin_jac);
                if flag < 0:
                    raise KINSOLError(flag)
        elif self.options["linear_solver"] == "SPGMR":
            flag = SUNDIALS.KINSpgmr(self.kinsol_mem, self.options["max_krylov"])
            if flag < 0:
                raise KINSOLError(flag)
            
            if self.problem_info["jacv_fcn"]:    
                flag = SUNDIALS.KINSpilsSetJacTimesVecFn(self.kinsol_mem, kin_jacv)
                if flag < 0:
                    raise KINSOLError(flag)
            
            if self.problem_info["prec_setup"] or self.problem_info["prec_solve"]:
                if not self.problem_info["prec_setup"]:
                    flag = SUNDIALS.KINSpilsSetPreconditioner(self.kinsol_mem, NULL,kin_prec_solve)
                    if flag < 0:
                        raise KINSOLError(flag)
                elif not self.problem_info["prec_solve"]:
                    flag = SUNDIALS.KINSpilsSetPreconditioner(self.kinsol_mem, kin_prec_setup, NULL)
                    if flag < 0:
                        raise KINSOLError(flag)
                else:  
                    flag = SUNDIALS.KINSpilsSetPreconditioner(self.kinsol_mem, kin_prec_setup, kin_prec_solve)
                    if flag < 0:
                        raise KINSOLError(flag)
                
                    
        else:
            raise KINSOLError(-100)
        
        self._added_linear_solver = True
            
    cpdef update_options(self):
        """
        Updates the simulation options.
        """
        cdef int flag
        
        #Update scaling
        arr2nv_inplace(self.options["y_scale"], self.y_scale)
        arr2nv_inplace(self.options["f_scale"], self.f_scale)
        
        flag = SUNDIALS.KINSetFuncNormTol(self.kinsol_mem, self.options["ftol"]);
        if flag < 0:
            raise KINSOLError(flag)
        
        flag = SUNDIALS.KINSetScaledStepTol(self.kinsol_mem, self.options["stol"]);
        if flag < 0:
            raise KINSOLError(flag)
            
        flag = SUNDIALS.KINSetMaxBetaFails(self.kinsol_mem, self.options["max_beta_fails"]); #The function KINSetMaxBetaFails specifies the maximum number of β-condition failures in the linesearch algorithm.
        if flag < 0:
            raise KINSOLError(flag)
            
        flag = SUNDIALS.KINSetMaxNewtonStep(self.kinsol_mem, self.options["max_newton_step"]); #The function KINSetMaxNewtonStep specifies the maximum allowable scaled length of the Newton step.
        if flag < 0:
            raise KINSOLError(flag)
            
        flag = SUNDIALS.KINSetNoMinEps(self.kinsol_mem, self.options["no_min_epsilon"]); #The function KINSetNoMinEps specifies a flag that controls whether or not the value of ǫ, the scaled linear residual tolerance, is bounded from below.
        if flag < 0:
            raise KINSOLError(flag)
            
        #flag = SUNDIALS.KINSetEtaParams(self.kinsol_mem, egamma, ealpha); #The function KINSetEtaParams specifies the parameters γ and α in the formula for η, in the case etachoice = KIN ETACHOICE2.
        #if flag < 0:
        #    raise KINSOLError(flag)
            
        flag = SUNDIALS.KINSetMaxSetupCalls(self.kinsol_mem, self.options["max_solves_between_setup_calls"]); #The function KINSetMaxSetupCalls specifies the maximum number of nonlinear iterations that can be performed between calls to the preconditioner setup function.
        if flag < 0:
            raise KINSOLError(flag)
            
        flag = SUNDIALS.KINSetNoInitSetup(self.kinsol_mem, self.options["no_initial_setup"]); #The function KINSetNoInitSetup specifies whether an initial call to the preconditioner setup function should be made or not.
        if flag < 0:
            raise KINSOLError(flag)
            
        flag = SUNDIALS.KINSetNumMaxIters(self.kinsol_mem, self.options["max_iter"]); #The function KINSetNumMaxIters specifies the maximum number of nonlinear iterations allowed.
        if flag < 0:
            raise KINSOLError(flag)
        
        flag = SUNDIALS.KINSetPrintLevel(self.kinsol_mem, self._get_print_level()); #The function KINSetPrintLevel specifies the level of verbosity of the output.
        if flag < 0:
            raise KINSOLError(flag)
    
    def _get_print_level(self):
        """
        Converts Assimulos verbosity level to Kinsol print level. Kinsol
        print level::
        
            0 no information displayed.
            
            1 for each nonlinear iteration display the following information: the scaled
            Euclidean ℓ2 norm of the system function evaluated at the current iterate,
            the scaled norm of the Newton step (only if using KIN NONE), and the
            number of function evaluations performed so far.
            
            2 display level 1 output and the following values for each iteration:
            kF(u)kDF (only for KIN NONE).
            kF(u)kDF ,∞ (for KIN NONE and KIN LINESEARCH).
            
            3 display level 2 output plus additional values used by the global strategy
            (only if using KIN LINESEARCH), and statistical information for the linear
            solver.
        """
        if self.verbosity >= NORMAL:
            return 0
        elif self.verbosity >= LOUD:
            return 2
        else:
            return 3
    
    cpdef _solve(self, y0=None):
        """
        Solves the system.
        """
        if y0 is not None:
            arr2nv_inplace(y0, self.y_temp)
        else:
            arr2nv_inplace(self.y, self.y_temp)
            
        if not self._added_linear_solver:
            self.add_linear_solver()
            
        #Update the solver options
        self.update_options()
        
        flag = SUNDIALS.KINSol(self.kinsol_mem, self.y_temp, self.options["strategy"], self.y_scale, self.f_scale)
        self.y = nv2arr(self.y_temp)
        
        if flag < 0:
            raise KINSOLError(flag)
        if flag == KIN_STEP_LT_STPTOL:
            print 'Scaled step length too small. Either an approximate solution or a local minimum is reached. Check value of residual.'
        
        self.store_statistics()
        
        return nv2arr(self.y_temp)
    
    def get_last_flag(self):
        """
        Returns the last flag reported by Kinsol.
        """
        cdef int flag = 0, lsflag = 0
        
        flag = SUNDIALS.KINDlsGetLastFlag(self.kinsol_mem, &lsflag)
        if flag < 0:
            raise KINSOLError(flag)
            
        return lsflag
    
    def store_statistics(self):
        cdef int flag
        cdef long int nfevalsLS, njevals, nbacktr, nbcfails, nniters
        cdef long int nliters, nlcfails, npevals, npsolves
        cdef long int nfevals
        
        flag = SUNDIALS.KINGetNumFuncEvals(self.kinsol_mem, &nfevals)
        if flag < 0:
            raise KINSOLError(flag)
        self.statistics["nfevals"] = nfevals
            
        flag = SUNDIALS.KINGetNumNonlinSolvIters(self.kinsol_mem, &nniters)
        if flag < 0:
            raise KINSOLError(flag)
        self.statistics["nniters"] = nniters
            
        flag = SUNDIALS.KINGetNumBacktrackOps(self.kinsol_mem, &nbacktr) #The function KINGetNumBacktrackOps returns the number of backtrack operations (step length adjustments) performed by the line search algorithm.
        if flag < 0:
            raise KINSOLError(flag)
        self.statistics["nbacktr"] = nbacktr
            
        flag = SUNDIALS.KINGetNumBetaCondFails(self.kinsol_mem, &nbcfails) #The function KINGetNumBetaCondFails returns the number of β-condition failures.
        if flag < 0:
            raise KINSOLError(flag)
        self.statistics["nbcfails"] = nbcfails
        
        if self.options["linear_solver"] == "SPGMR":
            
            flag = SUNDIALS.KINSpilsGetNumLinIters(self.kinsol_mem, &nliters)
            if flag < 0:
                raise KINSOLError(flag)
            self.statistics["nliters"] = nliters
            
            flag = SUNDIALS.KINSpilsGetNumConvFails(self.kinsol_mem, &nlcfails)
            if flag < 0:
                raise KINSOLError(flag)
            self.statistics["nlcfails"] = nlcfails
            
            flag = SUNDIALS.KINSpilsGetNumPrecEvals(self.kinsol_mem, &npevals)
            if flag < 0:
                raise KINSOLError(flag)
            self.statistics["npevals"] = npevals
            
            flag = SUNDIALS.KINSpilsGetNumPrecSolves(self.kinsol_mem, &npsolves)
            if flag < 0:
                raise KINSOLError(flag)
            self.statistics["npsolves"] = npsolves
                
            flag = SUNDIALS.KINSpilsGetNumJtimesEvals(self.kinsol_mem, &njevals)
            if flag < 0:
                raise KINSOLError(flag)
            self.statistics["njevals"] = njevals
            
            flag = SUNDIALS.KINSpilsGetNumFuncEvals(self.kinsol_mem, &nfevalsLS)
            if flag < 0:
                raise KINSOLError(flag)
            self.statistics["nfevalsLS"] = nfevalsLS
            
        elif self.options["linear_solver"] == "DENSE":
        
            flag = SUNDIALS.KINDlsGetNumJacEvals(self.kinsol_mem, &njevals) #The function KINDlsGetNumJacEvals returns the number of calls to the dense Jacobian approximation function.
            if flag < 0:
                raise KINSOLError(flag)
            self.statistics["njevals"] = njevals
                
            flag = SUNDIALS.KINDlsGetNumFuncEvals(self.kinsol_mem, &nfevalsLS) #The function KINDlsGetNumFuncEvals returns the number of calls to the user system function used to compute the difference quotient approximation to the dense or banded Jacobian.
            if flag < 0:
                raise KINSOLError(flag)
            self.statistics["nfevalsLS"] = nfevalsLS
        
            
    def print_statistics(self, verbose=NORMAL):
        """
        Should print the statistics.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
                     
        self.log_message(' Number of function evaluations              : '+ str(self.statistics["nfevals"]),   verbose)
        self.log_message(' Number of Nonlinear Iterations              : '+ str(self.statistics["nniters"]),   verbose)
        self.log_message(' Number of Backtrack Operations (Linesearch) : '+ str(self.statistics["nbacktr"]),   verbose) #The function KINGetNumBacktrackOps returns the number of backtrack operations (step length adjustments) performed by the line search algorithm.
        self.log_message(' Number of Beta-condition Failures           : '+ str(self.statistics["nbcfails"]),  verbose) #The function KINGetNumBetaCondFails returns the number of β-condition failures.
        
        if self.options["linear_solver"] == "SPGMR":
            self.log_message(' Number of Jacobian*Vector Evaluations       : '+ str(self.statistics["njevals"]),   verbose)
            self.log_message(' Number of F-Eval During Jac*Vec-Eval        : '+ str(self.statistics["nfevalsLS"]), verbose)
            self.log_message(' Number of Linear Iterations                 : '+ str(self.statistics["nliters"]), verbose)
            self.log_message(' Number of Linear Convergence Failures       : '+ str(self.statistics["nlcfails"]), verbose)
        elif self.options["linear_solver"] == "DENSE":
            self.log_message(' Number of Jacobian evaluations              : '+ str(self.statistics["njevals"]),   verbose)
            self.log_message(' Number of F-eval during Jac-eval            : '+ str(self.statistics["nfevalsLS"]), verbose)
        
    
        self.log_message('\nSolver options:\n',                                     verbose)
        self.log_message(' Solver                  : Kinsol',                       verbose)
        self.log_message(' Linear Solver           : ' + str(self.options["linear_solver"]),                       verbose)
        self.log_message(' Globalization Strategy  : ' + ("LINESEARCH" if self.options["strategy"]==1 else "NONE"),                       verbose)
        self.log_message(' Function Tolerances     : ' + str(self.options["ftol"]),  verbose)
        self.log_message(' Step Tolerances         : ' + str(self.options["stol"]),  verbose)
        self.log_message(' Variable Scaling        : ' + str(self.options["y_scale"]),  verbose)
        self.log_message(' Function Scaling        : ' + str(self.options["f_scale"]),  verbose)
        self.log_message('',                                                        verbose)
        
    def _set_ftol_method(self,ftol):
        self.options["ftol"] = ftol
        
    def _get_ftol_method(self):
        """
        Specifies the scalar used as a stopping tolerance on the scaled
        maximum norm of the residual.
        """
        return self.options["ftol"]
        
    ftol = property(_get_ftol_method,_set_ftol_method)
    
    def _set_stol_method(self,stol):
        self.options["stol"] = stol
        
    def _get_stol_method(self):
        """
        Specifies the scalar used as a stopping tolerance on the
        minimum scaled step length.
        """
        return self.options["stol"]
        
    stol = property(_get_stol_method,_set_stol_method)
    
    def _set_max_iter_method(self, max_iter):
        self.options["max_iter"] = max_iter
        
    def _get_max_iter_method(self):
        """
        Specifies the maximum number of nonlinear iterations.
        """
        return self.options["max_iter"]
        
    max_iter = property(_get_max_iter_method,_set_max_iter_method)
    
    def _set_no_initial_setup_method(self, no_initial_setup):
        self.options["no_initial_setup"] = no_initial_setup
        
    def _get_no_initial_setup_method(self):
        """
        Specifies whether or not an initial call to the preconditioner
        setup call should be made.
        """
        return self.options["no_initial_setup"]
        
    no_initial_setup = property(_get_no_initial_setup_method,_set_no_initial_setup_method)
    
    def _set_max_solves_between_setup_calls_method(self, max_solves_between_setup_calls):
        self.options["max_solves_between_setup_calls"] = max_solves_between_setup_calls
        
    def _get_max_solves_between_setup_calls_method(self):
        """
        Specifies the maximum number of allowed solve calls inbetween
        setup calls to the linear solver.
        """
        return self.options["max_solves_between_setup_calls"]
        
    max_solves_between_setup_calls = property(_get_max_solves_between_setup_calls_method,_set_max_solves_between_setup_calls_method)
    
    def _set_max_newton_step_method(self, max_newton_step):
        self.options["max_newton_step"] = max_newton_step
        
    def _get_max_newton_step_method(self):
        """
        Specifies the maximum allowed scaled length of the Newton step.
        """
        return self.options["max_newton_step"]
        
    max_newton_step = property(_get_max_newton_step_method,_set_max_newton_step_method)
    
    def _set_no_min_epsilon_method(self, no_min_epsilon):
        self.options["no_min_epsilon"] = no_min_epsilon
        
    def _get_no_min_epsilon_method(self):
        """
        Specifies if the scaled linear residual tolerance is bounded
        from below
        """
        return self.options["no_min_epsilon"]
        
    no_min_epsilon = property(_get_no_min_epsilon_method,_set_no_min_epsilon_method)
    
    def _set_max_beta_fails_method(self, max_beta_fails):
        self.options["max_beta_fails"] = max_beta_fails
        
    def _get_max_beta_fails_method(self):
        """
        Specifies the maximum number of beta condition fails in the 
        linesearch algorithm.
        """
        return self.options["max_beta_fails"]
        
    max_beta_fails = property(_get_max_beta_fails_method,_set_max_beta_fails_method)
    
    def _set_linear_solver(self, lsolver):
        if lsolver.upper() == "DENSE" or lsolver.upper() == "SPGMR":
            self.options["linear_solver"] = lsolver.upper()
        else:
            raise Exception('The linear solver must be either "DENSE" or "SPGMR".')
        
    def _get_linear_solver(self):
        """
        Specifies the linear solver to be used.
        
            Parameters::
            
                linearsolver
                        - Default 'DENSE'. Can also be 'SPGMR'.
        """
        return self.options["linear_solver"]
    
    linear_solver = property(_get_linear_solver, _set_linear_solver)
    
    def _set_globalization_strategy(self, lsolver):
        if lsolver.upper() == "LINESEARCH":
            self.options["strategy"] = KIN_LINSEARCH
        elif lsolver.upper() == "NONE":
            self.options["strategy"] = KIN_NONE
        else:
            raise Exception('The linear solver must be either "LINESEARCH" or "NONE".')
        
    def _get_globalization_strategy(self):
        """
        Specifies the globalization strategy to be used.
        
            Parameters::
            
                linearsolver
                        - Default 'LINSEARCH'. Can also be 'NONE'.
        """
        return self.options["strategy"]
    
    globalization_strategy = property(_get_globalization_strategy, _set_globalization_strategy)
    
    def _set_max_krylov(self, max_krylov):
        try:
            self.options["max_krylov"] = int(max_krylov)
        except:
            raise Exception("Maximum number of krylov dimension should be an integer.")
        if self.options["max_krylov"] < 0:
            raise Exception("Maximum number of krylov dimension should be an positive integer.")
            
    def _get_max_krylov(self):
        """
        Specifies the maximum number of dimensions for the krylov subspace to be used.
        
            Parameters::
            
                    maxkrylov
                            - A positive integer.
                            - Default 0
            
            Returns::
            
                The current value of maxkrylov.
                
        See SUNDIALS documentation 'CVSpgmr'
        """
        return self.options["max_krylov"]
    
    max_dim_krylov_subspace = property(_get_max_krylov, _set_max_krylov)
    
    def get_residual_norm_nonlinear_iterations(self): 
        return self.pData.nl_fnorm
        
    def get_residual_norm_linear_iterations(self):
        return self.pData.l_fnorm
        
    def get_log(self):
        return self.pData.log
    
class KINSOLError(Exception):
    """  
    Defines the KINSOLError and provides the textual error message.
    """ 
    msg = { KIN_MEM_FAIL : 'Memory allocation failed.',
            KIN_MEM_NULL : 'KINSOL was not properly created.',
            KIN_ILL_INPUT: 'There was an input error.',
            KIN_NO_MALLOC: 'Memory not allocated, call KINSOL_init(...)',
            KIN_LINESEARCH_NONCONV: 'Linesearch could not find a step big enough.',
            KIN_MAXITER_REACHED: 'Max number of iterations reached',
            KIN_MXNEWT_5X_EXCEEDED: 'Max size of Newton step exceeded 5 times.',
            KIN_LINESEARCH_BCFAIL: 'Linesearch beta-test failed, probably because of too poor progress.',
            KIN_LINSOLV_NO_RECOVERY: 'psolve encountered an error, but preconditioner is current',
            KIN_LINIT_FAIL: 'Init of linear solver failed.',
            KIN_LSETUP_FAIL: 'pset (preconditioner setup fct) encountered an error',
            KIN_LSOLVE_FAIL: 'Error in either psolve or in linear solver routine.',
            KIN_SYSFUNC_FAIL: 'Call to RES failed.',
            KIN_FIRST_SYSFUNC_ERR: 'Call to RHS failed on first call',
            KIN_REPTD_SYSFUNC_ERR: 'Call to RHS failed multiple times.',
            KIN_STEP_LT_STPTOL: 'Scaled step length too small. Either an approximate solution or a local minimum is reached. Check value of residual.'}
    
    def __init__(self, value):
        self.value = value
        
    def __str__(self): 
        try:
            return repr(self.msg[self.value])    
        except KeyError:
            return repr('Sundials failed with flag %s.'%(self.value))

