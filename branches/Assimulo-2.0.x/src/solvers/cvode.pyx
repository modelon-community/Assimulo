#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2011 Modelon AB
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

cimport numpy as N
import numpy as N
import traceback 

from assimulo.exception import *

from assimulo.explicit_ode cimport Explicit_ODE

cimport sundials_includes as Sun

#Various C includes transfered to namespace
from sundials_includes cimport N_Vector, realtype, N_VectorContent_Serial, DENSE_COL
from sundials_includes cimport memcpy, N_VNew_Serial, DlsMat, PyArray_DATA, import_array
from sundials_includes cimport malloc, free, realtype

include "constants.pxi" #Includes the constants (textual include)
include "sundials_constants.pxi" #Sundials related constants
include "sundials_callbacks.pxi"

cdef class CVode(Explicit_ODE):
    """
    CVode.
    """
    cdef void* cvode_mem
    cdef ProblemData pData      #A struct containing information about the problem
    cdef N_Vector yTemp, ydTemp
    cdef object f
    cdef dict statistics

    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.options["atol"] = 1.0e-6        #The absolute tolerance
        self.options["rtol"] = 1.0e-6        #The relative tolerance
        self.options["maxh"] = 0.0             #Maximum step-size
        self.options["inith"] = 0.0
        self.options["suppress_alg"] = False
        self.options["suppress_sens"] = False 
        self.options["maxord"] = 5        #Maximum order allowed
        self.options["usejac"] = True
        self.options["usesens"] = True
        self.options["maxsteps"] = 10000
        self.options["maxkrylov"] = 5
        self.options["precond"] = PREC_NONE
        self.options["linear_solver"] = "DENSE"
        self.options["iter"] = "Newton"
        self.options["discr"] = "BDF"
        
        #Statistics
        self.statistics = {}
        self.statistics["nfevals"]    = 0 #Function evaluations
        self.statistics["nsteps"]     = 0 #Number of steps
        self.statistics["netfails"]   = 0 #Number of error test failures
        self.statistics["nlinsetups"] = 0
        self.statistics["nncfails"]   = 0 #Nonlinear fails
        self.statistics["nniters"]    = 0 #Nonlinear iterations
        self.statistics["ngevals"]    = 0 #Root evaluations
        self.statistics["njevals"]    = 0 #Jacobian evaluations
        self.statistics["nfevalsLS"]  = 0 #Function evaluations due to Jac
        self.statistics["njvevals"]   = 0 #Number of Jacobian*Vector eval
        
        #Internal temporary result vector
        #self.yd1 = N.array([0.0]*len(problem.y0))
        
        #RHS-Function
        self.f = problem.rhs_internal
        
        #Solver support
        self.supports["one_step_mode"] = True
        self.supports["interpolated_output"] = True
        
        self.pData = ProblemData()
        
        #Populated the ProblemData
        self.set_problem_data()
    
    cdef set_problem_data(self):
        
        #Sets the residual or rhs
        self.pData.RHS = <void*>self.problem.f
        self.pData.dim = self.problem_info["dim"]
        self.pData.memSize = self.pData.dim*sizeof(realtype)
        
        #Set the ndarray to the problem struct
        #self.yTemp   = N.zeros(self.pData.dim, dtype=N.float, order='c')
        #self.ydTemp  = N.zeros(self.pData.dim, dtype=N.float, order='c')
        #self.pData.y  = <void*>self.yTemp
        #self.pData.yd = <void*>self.ydTemp
        
        if self.problem_info["state_events"] is True: #Sets the root function
            self.pData.ROOT = <void*>self.problem.state_events
            self.pData.dimRoot = self.problem_info["dimRoot"]
            self.pData.memSizeRoot = self.pData.dimRoot*sizeof(realtype)
    
        if self.problem_info["jac_fcn"] is True: #Sets the jacobian
            self.pData.JAC = <void*>self.problem.jac
            self.pData.memSizeJac = self.pData.dim*self.pData.dim*sizeof(realtype)
        
        if self.problem_info["jacv_fcn"] is True: #Sets the jacobian times vector
            self.pData.JACV = <void*>self.problem.jacv
            
        if self.problem_info["sens_fcn"] is True: #Sets the sensitivity function
            self.pData.SENS = <void*>self.problem.sens
        
        if self.problem_info["dimSens"] > 0: #Sensitivity parameters (does not need the sensitivity function)
            self.pData.dimSens = self.problem_info["dimSens"]   
            self.pData.p = <realtype*> malloc(self.problem_info["dimSens"]*sizeof(realtype))
        else:
            self.pData.dimSens = 0
    
    cdef initialize_cvode(self):
        cdef int flag #Used for return

        self.yTemp = arr2nv(self.y_cur)
        
        #Updates the switches
        if self.problem_info["switches"]:
            #self.switches = sw0
            #self.pData.sw = <void*>self.switches
            pass
            
        if self.cvode_mem == NULL: #The solver is not initialized
            
            #Create the solver
            self.cvode_mem= Sun.CVodeCreate(CV_BDF if self.options["discr"] == "BDF" else CV_ADAMS, CV_NEWTON if self.options["iter"] == "Newton" else CV_FUNCTIONAL)
            if self.cvode_mem == NULL:
                raise CVodeError(CV_MEM_FAIL)
            
            #Specify the residual and the initial conditions to the solver
            flag = Sun.CVodeInit(self.cvode_mem, cv_rhs, self.t_cur, self.yTemp)
            if flag < 0:
                raise CVodeError(flag, self.t_cur)
                
            #Specify the root function to the solver
            if self.problem_info["state_events"]:
                flag = Sun.CVodeRootInit(self.cvode_mem, self.pData.dimRoot, cv_root)
                if flag < 0:
                    raise CVodeError(flag, self.t_cur)
                    
            #Specify the error handling
            flag = Sun.CVodeSetErrHandlerFn(self.cvode_mem, cv_err, <void*>self.pData)
            if flag < 0:
                raise CVodeError(flag, self.t_cur)
            
        else: #The solver needs to be reinitialized
            #Reinitialize
            flag = Sun.CVodeReInit(self.cvode_mem, self.t_cur, self.yTemp)
            if flag < 0:
                raise CVodeError(flag, self.t_cur)
        
        #Set the user data
        flag = Sun.CVodeSetUserData(self.cvode_mem, <void*>self.pData)
        if flag < 0:
            raise CVodeError(flag, self.t_cur)
            
    
    cpdef N.ndarray interpolate(self,double t,int k = 0):
        """
        Calls the internal CVodeGetDky for the interpolated values at time t.
        t must be within the last internal step. k is the derivative of y which
        can be from zero to the current order.
        """
        cdef flag
        cdef N_Vector dky=N_VNew_Serial(self.pData.dim)
        
        flag = Sun.CVodeGetDky(self.cvode_mem, t, k, dky)
        
        if flag < 0:
            raise CVodeError(flag, t)
        
        return nv2arr(dky)
        
    cpdef N.ndarray interpolate_sensitivity(self, realtype t, int k, int i=-1):
        """
        This method calls the internal method CVodeGetSensDky which computes the k-th derivatives
        of the interpolating polynomials for the sensitivity variables at time t.
        
            Parameters::
                    
                    t
                        - Specifies the time at which sensitivity information is requested. The time
                          t must fall within the interval defined by the last successful step taken
                          by CVodeS.
                    
                    k   
                        - The order of derivatives.
                        
                    i
                        - Specifies the sensitivity derivative vector to be returned (0<=i<=Ns)
                        
            Return::
            
                    A matrix containing the Ns vectors or a vector if i is specified.
        """
        cdef N_Vector dkyS=N_VNew_Serial(self.pData.dimSens)
        cdef int flag
        
        if i==-1:
            
            matrix = []
            
            for x in range(self.pData.dimSens):
                flag = Sun.CVodeGetSensDky1(self.cvode_mem, t, k, x, dkyS)
                if flag<0:
                    raise CVodeError(flag, t)
                
                matrix += [nv2arr(dkyS)]
            
            return N.array(matrix)
        else:
            flag = Sun.CVodeGetSensDky1(self.cvode_mem, t, k, i, dkyS)
            if flag <0:
                raise CVodeError(flag, t)
            
            return nv2arr(dkyS)
        
    cpdef step(self,double t,N.ndarray y,double tf,dict opts):
        cdef int flag
        cdef N_Vector yout
        cdef double tret = 0.0
        cdef double tr
        cdef N.ndarray yr
        
        yout = arr2nv(y)
        
        #Get options
        initialize  = opts["initialize"]
        output_list = opts["output_list"]        
        
        #Initialize?
        if initialize:
            self.initialize_cvode()
            self.initialize_options()
        
        #Set stop time
        flag = Sun.CVodeSetStopTime(self.cvode_mem, tf)
        if flag < 0:
            raise CVodeError(flag, t)
        
        #Integration loop
        flag = Sun.CVode(self.cvode_mem,tf,yout,&tret,CV_ONE_STEP)
        if flag < 0:
            raise CVodeError(flag, tret)
            
        #Store results
        tr = tret
        yr = nv2arr(yout)
        
        if flag == CV_ROOT_RETURN: #Found a root
            flag = ID_EVENT #Convert to Assimulo flags
            
        if flag == CV_TSTOP_RETURN: #Reached tf
            flag = ID_COMPLETE
            self.store_statistics()
            
                
        return flag, tr, yr
    
    cpdef integrate(self,double t,N.ndarray[ndim=1, dtype=realtype] y,double tf,dict opts):
        cdef int flag
        cdef N_Vector yout
        cdef double tret = 0.0 
        cdef list tr = [], yr = []
        
        yout = arr2nv(y)
        
        #Get options
        initialize  = opts["initialize"]
        output_list = opts["output_list"]        
        
        #Initialize? 
        if initialize:
            self.initialize_cvode() 
            self.initialize_options()
        
        #Set stop time
        flag = Sun.CVodeSetStopTime(self.cvode_mem, tf)
        if flag < 0:
            raise CVodeError(flag, t)
         
        #Integration loop
        while True:
            #print self.y_cur, yout
            flag = Sun.CVode(self.cvode_mem,tf,yout,&tret,CV_ONE_STEP)
            if flag < 0:
                raise CVodeError(flag, tret)
            
            #Store results
            tr.append(tret)
            yr.append(nv2arr(yout))
            
            if flag == CV_ROOT_RETURN: #Found a root
                flag = ID_EVENT #Convert to Assimulo flags
                break
            if flag == CV_TSTOP_RETURN: #Reached tf
                flag = ID_COMPLETE
                self.store_statistics()
                break
                
        return flag, tr, yr
    
    cdef state_event_info(self):
        """
        Returns the event info.
        """
        cdef int* c_info
        cdef flag
        
        # Allocate memory for the event_info_ vector and initialize to zeros
        c_info = <int*> malloc(self.pData.dimRoot*sizeof(int))
        for k in range(self.pData.dimRoot):
            c_info[k] = 0
        
        # Fetch data on which root functions that became zero and store in class
        flag = Sun.CVodeGetRootInfo(self.cvode_mem, c_info)
        if flag < 0:
            raise CVodeError(flag)
        
        #event_info = PyArray_SimpleNew(1,&self.pData.dimRoot,NPY_INT)
        event_info = [0]*self.pData.dimRoot
        
        for k in range(self.pData.dimRoot):
            event_info[k] = c_info[k]
        
        # Remember to deallocate
        free(c_info)
        
        return event_info
        
    cpdef initialize_options(self):
        """
        Updates the simulation options.
        """
        cdef flag
        
        #Choose a linear solver if and only if NEWTON is choosen
        if self.options["linear_solver"] == 'DENSE' and self.options["iter"] == "Newton":
            #Specify the use of the internal dense linear algebra functions.
            flag = Sun.CVDense(self.cvode_mem, self.pData.dim)
            if flag < 0:
                raise CVodeError(flag)
                
            #Specify the jacobian to the solver
            if self.pData.JAC != NULL and self.options["usejac"]:
                flag = Sun.CVDlsSetDenseJacFn(self.cvode_mem, cv_jac)
                if flag < 0:
                    raise CVodeError(flag)
            else:
                flag = Sun.CVDlsSetDenseJacFn(self.cvode_mem, NULL)
                if flag < 0:
                    raise CVodeError(flag)
                    
        elif self.options["linear_solver"] == 'SPGMR' and self.options["iter"] == "Newton":
            #Specify the use of CVSPGMR linear solver.
            flag = Sun.CVSpgmr(self.cvode_mem, self.options["precond"], self.options["max_krylov"])
            if flag < 0:
                raise CVodeError(flag)
                
            #Specify the jacobian times vector function
            if self.pData.JACV != NULL and self.options["usejac"]:
                flag = Sun.CVSpilsSetJacTimesVecFn(self.cvode_mem, cv_jacv)
                if flag < 0:
                    raise CVodeError(flag)
            else:
                flag = Sun.CVSpilsSetJacTimesVecFn(self.cvode_mem, NULL)
                if flag < 0:
                    raise CVodeError(flag)
        else: #Functional Iteration choosen.
            pass #raise CVodeError(100,t0) #Unknown error message

        #Maximum order
        flag = Sun.CVodeSetMaxOrd(self.cvode_mem, int(self.options["maxord"]))
        if flag < 0:
            raise CVodeError(flag)
            
        #Initial step
        flag = Sun.CVodeSetInitStep(self.cvode_mem, self.options["inith"])
        if flag < 0:
            raise CVodeError(flag)
        
        #Maximum step
        flag = Sun.CVodeSetMaxStep(self.cvode_mem, self.options["maxh"])
        if flag < 0:
            raise CVodeError(flag)
            
        #Maximum Number of steps
        flag = Sun.CVodeSetMaxNumSteps(self.cvode_mem, self.options["maxsteps"])
        if flag < 0:
            raise CVodeError(flag)
        
        #Tolerances
        flag = Sun.CVodeSVtolerances(self.cvode_mem, self.options["rtol"], arr2nv(self.options["atol"]))
        if flag < 0:
            raise CVodeError(flag)
    
    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
    def _get_atol(self):
        """
        Defines the absolute tolerance(s) that is to be used by the solver.
        Can be set differently for each variable.
        
            Parameters::
            
                atol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float or a numpy vector
                          of floats.
                        
                            Example:
                                atol = [1.0e-4, 1.0e-6]
        
        See SUNDIALS IDA documentation 4.5.2 for more details.
        """
        return self.options["atol"]
    
    atol=property(_get_atol,_set_atol)
    
    cdef void store_statistics(self):
        """
        Retrieves and stores the statistics.
        """
        cdef long int nsteps = 0, njevals = 0, ngevals = 0, netfails = 0, nniters = 0, nncfails = 0
        cdef long int nSniters = 0, nSncfails = 0, nfevalsLS = 0, njvevals = 0, nfevals = 0
        cdef long int nfSevals = 0,nfevalsS = 0,nSetfails = 0,nlinsetupsS = 0, nlinsetups = 0
        cdef int qlast = 0, qcur = 0
        cdef realtype hinused = 0.0, hlast = 0.0, hcur = 0.0, tcur = 0.0

        if self.options["linear_solver"] == "SPGMR":
            flag = Sun.CVSpilsGetNumJtimesEvals(self.cvode_mem, &njvevals) #Number of jac*vector
            flag = Sun.CVSpilsGetNumRhsEvals(self.cvode_mem, &nfevalsLS) #Number of rhs due to jac*vector
        else:
            flag = Sun.CVDlsGetNumJacEvals(self.cvode_mem, &njevals) #Number of jac evals
            flag = Sun.CVDlsGetNumRhsEvals(self.cvode_mem, &nfevalsLS) #Number of res evals due to jac evals
            
        flag = Sun.CVodeGetNumGEvals(self.cvode_mem, &ngevals) #Number of root evals
        
        #Get all integrator statistics
        flag = Sun.CVodeGetIntegratorStats(self.cvode_mem, &nsteps, &nfevals, &nlinsetups, &netfails, &qlast,
                                       &qcur, &hinused, &hlast, &hcur, &tcur)
        
        flag = Sun.CVodeGetNonlinSolvStats(self.cvode_mem, &nniters, &nncfails) #Number of nonlinear iteration
                                                                            #Number of nonlinear conv failures
        
        self.statistics["nsteps"] = nsteps
        self.statistics["nfevals"] = nfevals
        self.statistics["netfails"] = netfails
        self.statistics["nniters"]  = nniters
        self.statistics["nncfails"] = nncfails
        self.statistics["nfevalsLS"] = nfevalsLS
        self.statistics["njvevals"] = njvevals
        self.statistics["njevals"] = njevals
        self.statistics["ngevals"] = ngevals
        
        #If sensitivity    
        if self.pData.dimSens > 0:
            flag = Sun.CVodeGetSensStats(self.cvode_mem, &nfSevals, &nfevalsS, &nSetfails, &nlinsetupsS)
            flag = Sun.CVodeGetSensNonlinSolvStats(self.cvode_mem, &nSniters, &nSncfails)
                
    def print_statistics(self, verbose=NORMAL):
        """
        Should print the statistics.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of Steps                          : '+str(self.statistics["nsteps"]),          verbose)               
        self.log_message(' Number of Function Evaluations           : '+str(self.statistics["nfevals"]),         verbose)
        if self.options["linear_solver"] == "SPGMR":
            self.log_message(' Number of Jacobian*Vector Evaluations    : ' + str(self.statistics["njvevals"]),  verbose)
            self.log_message(' Number of F-Evals During Jac*Vec-Evals   : ' + str(self.statistics["nfevalsLS"]), verbose)
        else:     
            self.log_message(' Number of Jacobian Evaluations           : '+ str(self.statistics["njevals"]),    verbose)
            self.log_message(' Number of F-Eval During Jac-Eval         : '+ str(self.statistics["nfevalsLS"]),  verbose)
        self.log_message(' Number of Root Evaluations               : '+ str(self.statistics["ngevals"]),        verbose)
        self.log_message(' Number of Error Test Failures            : '+ str(self.statistics["netfails"]),       verbose)
        self.log_message(' Number of Nonlinear Iterations           : '+ str(self.statistics["nniters"]),        verbose)
        self.log_message(' Number of Nonlinear Convergence Failures : '+ str(self.statistics["nncfails"]),       verbose)
        
        """
            if self.problem_info['dimSens'] > 0: #Senstivity calculations is on
                sens_stats = self.get_sensitivity_self.statistics[""]()
                
                self.log_message('\nSensitivity self.statistics[""]:\n'
                self.log_message(' Number of Sensitivity Calculations             :', sens_stats[0]
                self.log_message(' Number of F-Evals Due to Finite Approximation  :', sens_stats[1]
                self.log_message(' Number of Local Error Test Failures            :', sens_stats[2]
                self.log_message(' Number of Linear Setups                        :', sens_stats[3]
                self.log_message(' Number of Nonlinear iterations                 :', sens_stats[4]
                self.log_message(' Number of Nonlinear Convergance Failures       :', sens_stats[5]
                
                self.log_message('\nSensitivity options:\n'
                self.log_message(' Method                   : ' ,self.sensmethod
                self.log_message(' Difference quotient type : ' ,self.dqtype
                self.log_message(' Suppress Sens            : ' ,self.suppress_sens
            """
        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : CVode',                         verbose)
        self.log_message(' Linear Multistep Method : ' +self.options["discr"],       verbose)
        self.log_message(' Nonlinear Solver        : ' + self.options["iter"],       verbose)
        self.log_message(' Maxord                  : ' + str(self.options["maxord"]),verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self.options["atol"]),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)


class CVodeError(Exception):
    """
    Defines the CVodeError and provides the textual error message.
    """
    msg = { CV_TOO_MUCH_WORK     : 'The solver took max internal steps but could not reach tout.',
            CV_TOO_MUCH_ACC      : 'The solver could not satisfy the accuracy demanded by the user for some internal step.',
            CV_ERR_FAIL          : 'Error test failures occurred too many times during one internal time step or minimum step size was reached.',
            CV_CONV_FAIL         : 'Convergence test failures occurred too many times during one internal time step or minimum step size was reached.',
            CV_LINIT_FAIL        : 'The linear solvers initialization function failed.',
            CV_LSETUP_FAIL       : 'The linear solvers setup function failed in an unrecoverable manner.',
            CV_LSOLVE_FAIL       : 'The linear solvers solve function failed in an unrecoverable manner.',
            CV_RHSFUNC_FAIL      : 'The user-provided rhs function failed in an unrecoverable manner.',
            CV_FIRST_RHSFUNC_ERR : 'The right-hand side function failed at the first call.',
            CV_REPTD_RHSFUNC_ERR : 'The right-hand side function had repetead recoverable errors.',
            CV_UNREC_RHSFUNC_ERR : 'The right-hand side function had a recoverable error, but no recovery is possible.',
            CV_RTFUNC_FAIL       : 'The rootfinding function failed in an unrecoverable manner.',
            CV_MEM_FAIL          : 'A memory allocation failed.',
            CV_MEM_NULL          : 'The cvode_mem argument was NULL.',
            CV_ILL_INPUT         : 'One of the function inputs is illegal.',
            CV_NO_MALLOC         : 'The CVode memory block was not allocated by a call to CVodeMalloc.',
            CV_BAD_K             : 'The derivative order k is larger than the order used.',
            CV_BAD_T             : 'The time t is outside the last step taken.',
            CV_BAD_DKY           : 'The output derivative vector is NULL.',
            CV_TOO_CLOSE         : 'The output and initial times are too close to each other.'}
    
    def __init__(self, value, t = 0.0):
        self.value = value
        self.t = t
        
    def __str__(self): 
        try:
            return repr(self.msg[self.value]+' At time %f.'%self.t)    
        except KeyError:
            return repr('Sundials failed with flag %s. At time %f.'%(self.value, self.t))   



class IDAError(Exception):
    """
    Defines the IDAError and provides the textual error message.
    """
    msg = { IDA_TOO_MUCH_WORK    : 'The solver took max internal steps but could not reach tout.',
            IDA_TOO_MUCH_ACC     : 'The solver could not satisfy the accuracy demanded by the user for some internal step.',
            IDA_ERR_FAIL         : 'Error test failures occurred too many times during one internal time step or minimum step size was reached.',
            IDA_CONV_FAIL        : 'Convergence test failures occurred too many times during one internal time step or minimum step size was reached.',
            IDA_LINIT_FAIL       : 'The linear solvers initialization function failed.',
            IDA_LSETUP_FAIL      : 'The linear solvers setup function failed in an unrecoverable manner.',
            IDA_LSOLVE_FAIL      : 'The linear solvers solve function failed in an unrecoverable manner.',
            IDA_RES_FAIL         : 'The user-provided residual function failed in an unrecoverable manner.',
            IDA_REP_RES_FAIL     : 'The user-provided residual function repeatedly returned a recoverable error flag, but the solver was unable to recover.',
            IDA_RTFUNC_FAIL      : 'The rootfinding function failed in an unrecoverable manner.',
            IDA_CONSTR_FAIL      : 'The inequality constraints were violated and the solver was unable to recover.',
            IDA_FIRST_RES_FAIL   : 'The user-provided residual function failed recoverable on the first call.',
            IDA_LINESEARCH_FAIL  : 'The line search failed.',
            IDA_NO_RECOVERY      : 'The residual function, linear solver setup function or linear solver solve function had a recoverable failure. But IDACalcIC could not recover.',
            IDA_MEM_NULL         : 'The ida_mem argument was NULL.',
            IDA_MEM_FAIL         : 'A memory allocation failed.',
            IDA_ILL_INPUT        : 'One of the function inputs is illegal.',
            IDA_NO_MALLOC        : 'The IDA memory was not allocated by a call to IDAInit.',
            IDA_BAD_EWT          : 'Zero value of some error weight component.',
            IDA_BAD_K            : 'The k-th derivative is not available.',
            IDA_BAD_T            : 'The time t is outside the last step taken.',
            IDA_BAD_DKY          : 'The vector argument where derivative should be stored is NULL.',
            IDA_SRES_FAIL        : 'The user-provided sensitivity residual function failed in an unrecoverable manner.',
            IDA_REP_SRES_ERR     : 'The user-provided sensitivity residual function repeatedly returned a recoverable error flag, but the solver was unable to recover.',
            IDA_BAD_IS           : 'The sensitivity identifier is not valid.'}
    
    def __init__(self, value, t = 0.0):
        self.value = value
        self.t = t
        
    def __str__(self): 
        try:
            return repr(self.msg[self.value]+' At time %f.'%self.t)    
        except KeyError:
            return repr('Sundials failed with flag %s. At time %f.'%(self.value, self.t))
