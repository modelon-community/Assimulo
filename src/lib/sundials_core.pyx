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

"""
Cython Wrapper for interfacing Python with CVode and IDA (Sundials Version 2.4.0)
Claus Fuhrer,        Lund University        
Christian Andersson, Lund University        

see also Jon Olav Vik: 
http://codespeak.net/pipermail/cython-dev/2009-June/005947.html

"""
from __future__ import division
#from sundials_core cimport * sundials_core
include "sundials_core.pxd" #Should really be "from sundials_core.pxd cimport *"
                            #But the setup.py script does not seem to work
                            #if this is the case.
include "sundials_core.pxi" #Includes the constants (textual include)


#=====================================
#This section contains the callback functions used for connecting Sundials
#to Assimulo.Problem.
#=====================================

cdef int cv_rhs(realtype t, N_Vector yv, N_Vector yvdot, void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.f to the Sundials
    right-hand-side function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    #cdef ndarray[realtype, ndim=1, mode='c'] rhs #Used for return from the user function
    #(<ndarray>pData.y).data =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    cdef ndarray y = nv2arr(yv)
    cdef realtype* resptr=(<N_VectorContent_Serial>yvdot.content).data
    
    if pData.dimSens>0: #Sensitivity activated
        p = realtype2arr(pData.p,pData.dimSens)
        try:
            if pData.sw != NULL:
                rhs = (<object>pData.RHS)(t,y,<list>pData.sw, p=p)
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
    cdef ndarray y = nv2arr(yv)
    
    try:
        if pData.sw != NULL:
            #jac=(<object>pData.JAC)(t,(<ndarray>pData.y),<list>pData.sw)
            jac=(<object>pData.JAC)(t,y,<list>pData.sw)
        else:
            jac=(<object>pData.JAC)(t,y)
        
        #This needs further investigations:
        #memcpy(Jacobian.data,<realtype*>jac.data, pData.memSizeJac)
        
        for i in range(Neq):
            col_i = DENSE_COL(Jacobian, i)
            for j in range(Neq):
                col_i[j] = jac[j,i]

        return CVDLS_SUCCESS
    except:
        return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)
        
cdef int cv_jacv(N_Vector vv, N_Vector Jv, realtype t, N_Vector yv, N_Vector fyv,
				    void *problem_data, N_Vector tmp):
    """
    This method is used to connect the Assimulo.Problem.jacv to the Sundials
    Jacobian times vector function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef ndarray y  = nv2arr(yv)
    cdef ndarray v  = nv2arr(vv)
    cdef ndarray fy = nv2arr(fyv)
    
    cdef realtype* jacvptr=(<N_VectorContent_Serial>Jv.content).data
    
    try:
        jacv = (<object>pData.JACV)(t,y,fy,v)
        
        for i in range(pData.dim):
                jacvptr[i] = jacv[i]
        
        return SPGMR_SUCCESS
    except:
        return SPGMR_ATIMES_FAIL_REC
    

cdef int cv_root(realtype t, N_Vector yv, realtype *gout,  void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.state_events to the Sundials
    Root-finding function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    #cdef ndarray[realtype, ndim=1, mode='c'] root #Used for return from the user function
    #(<ndarray>pData.y).data =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    cdef ndarray y = nv2arr(yv)
    
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
    cdef ndarray[realtype, ndim=1, mode='c'] res #Used for return from the user function
    #(<ndarray>pData.y).data  =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #(<ndarray>pData.yd).data =  <realtype*>((<N_VectorContent_Serial>yvdot.content).data)
    cdef ndarray y = nv2arr(yv)
    cdef ndarray yd = nv2arr(yvdot)
    cdef realtype* resptr=(<N_VectorContent_Serial>residual.content).data
     
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
        except(LinAlgError,ZeroDivisionError):
            return IDA_REC_ERR # recoverable error (see Sundials description)
        except:
            print "Unexpected error, probably due to a programing error in rhs/res function:\n"
            traceback.print_exc()
            return -1
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
        except(LinAlgError,ZeroDivisionError):
            return IDA_REC_ERR # recoverable error (see Sundials description)
        except:
            print "Unexpected error, probably due to a programing error in rhs/res function:\n"
            traceback.print_exc()
            return -1
            
cdef int ida_jac(int Neq, realtype t, realtype c, N_Vector yv, N_Vector yvdot, N_Vector residual, DlsMat Jacobian,
                 void* problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    This method is used to connect the Assimulo.Problem.jac to the Sundials
    Jacobian function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef ndarray[realtype, ndim=2, mode='c'] jac #Used for return from the user function
    cdef realtype* col_i=DENSE_COL(Jacobian,0)
    #(<ndarray>pData.y).data  =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #(<ndarray>pData.yd).data =  <realtype*>((<N_VectorContent_Serial>yvdot.content).data)
    cdef ndarray y = nv2arr(yv)
    cdef ndarray yd = nv2arr(yvdot)
    
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
    except: 
        return IDADLS_JACFUNC_RECVR #Recoverable Error

cdef int ida_root(realtype t, N_Vector yv, N_Vector yvdot, realtype *gout,  void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.state_events to the Sundials
    root function.
    """
    cdef ProblemData pData = <ProblemData>problem_data
    cdef ndarray[realtype, ndim=1, mode='c'] root #Used for return from the user function
    #(<ndarray>pData.y).data  =  <realtype*>((<N_VectorContent_Serial>yv.content).data)
    #(<ndarray>pData.yd).data =  <realtype*>((<N_VectorContent_Serial>yvdot.content).data)
    cdef ndarray y = nv2arr(yv)
    cdef ndarray yd = nv2arr(yvdot)
     
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

#=====================
# CLASS IMPLEMENTATION
#=====================

# Error handling
# ============== 
 
class SundialsError(Exception):
    """
    Sundials exception
    """
    def __init__(self, value, t = 0.0):
        self.value = value
        self.t = t
        
    def __str__(self): 
        try:
            return repr(self.msg[self.value]+' At time %f.'%self.t)    
        except KeyError:
            return repr('Sundials failed with flag %s. At time %f.'%(self.value, self.t))

class IDAError(SundialsError):
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
    pass  
    
class CVodeError(SundialsError):
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
    pass    

        
# Solver classes
#===============

cdef class Sundials:
    """
    Base class for wrapping the Sundials solvers to Python.
    """
    cdef void* solver           #Contains the solver (IDA or CVode)
    cdef ProblemData pData      #A struct containing information about the problem
    #cdef ProblemData *ppData    #A pointer to the problem data
    cdef ndarray y_nd           #Storage for the states
    cdef ndarray yd_nd          #Storage for the derivatives
    cdef ndarray p_nd           #Storage for the parameters
    cdef realtype* pp_nd        #Pointer to pData.p.data
    cdef list sol               #List for storing the solution
    cdef N_Vector y_cur, yd_cur #Store the states and state derivatives
    cdef public switches           #Store the switches
    cdef public booleantype usejac #Option for turning on and off the jacobian
    cdef public int maxord         #Defines the maximum order allowed.
    cdef public realtype inith     #The initial step that is to be used.
    cdef public ndarray atol       #The absolute tolerance
    cdef public realtype rtol      #The relative tolerance
    cdef public realtype maxh      #The maximum step-size
    cdef public solver_stats       #Stores the solver statistics
    cdef public solver_sens_stats  #Stores the solver sensitivity statistics
    cdef public booleantype suppress_alg  #Suppress algebraic components or not
    cdef public booleantype sensToggleOff #Should the sensitivity calculations be turned off?
    cdef public booleantype errconS     #Error control strategy
    cdef public int maxcorS             #Maximum number of nonlinear solver iterations (sens) per step
    cdef public int DQtype              #Difference quotient type
    cdef public realtype DQrhomax       #Selection parameter (sens)
    cdef booleantype _flag_active_sens  #Flag used for determine if reinit sens is required
    cdef public int ism                 #Sensitivity solution method
    cdef public object p                #The parameter information
    cdef public object pbar 
    cdef N_Vector *ySO                  #The sensitivity start matrix
    cdef public yS0                     #The sensitivity Provided start matrix
    cdef public npy_intp nbrRoot        #The number of root functions
    cdef public list p_result           #The sensitivity result matrix
    cdef public booleantype save_detailed_info #Save detailed information about the solver
    cdef public linear_solver #The linear solver used
    cdef public int pretype #Specifies the preconditioner type
    cdef public int max_krylov #Maximum number of krylov dimensions
    cdef public statistics
    
    def __cinit__(self):
        self.pData = ProblemData() #Create a new problem struct
        #self.ppData = &self.pData
        
        #Default values
        self.usejac = True
        self.maxord = 12
        self.inith = 0.0
        self.maxh = 0.0
        self.rtol = 1.e-6
        self.solver_stats = [0,0,0,0,0,0,0,0]
        self.solver_sens_stats = [0,0,0,0,0,0]
        self.suppress_alg = False
        self.yS0 = None
        self.save_detailed_info = False
        self.linear_solver = 'DENSE'
        self.pretype = PREC_NONE
        self.max_krylov = 5
        self.statistics = {"JVEC":0, "RHSJVEC":0}
        
        #Default values (Sensitivity)
        self.sensToggleOff = False #Toggle the sensitivity calculations off
        self.errconS  = False #Specifies whether sensitivity variables are included in the error control mechanism
        self.maxcorS  = 3 #Maximum number of nonlinear solver iterations for sensitivity variables per step.
        self.DQtype   = CV_CENTERED #(IDA_CENTERED) #Specifies the difference quotient type
        self.DQrhomax = 0.0 #Positive value of the selction parameter used in deciding switching
        self.ism      = CV_STAGGERED #(IDA_STAGGERED) #The sensitivity solution method
        self._flag_active_sens = False #The sensitivity is not activated
        
    cpdef set_problem_info(self, RHS, dim, ROOT = None, dimRoot = None, JAC = None, SENS = None, dimSens = None, JACV = None):
        """
        Sets the problem information to the problem struct.
        """
        #Sets the residual or rhs
        self.pData.RHS = <void*>RHS
        self.pData.dim = dim
        self.pData.memSize = dim*sizeof(realtype)
        
        #Default verbosity
        self.pData.verbose = 0
        
        #Set the ndarray to the problem struct
        self.y_nd   = np.zeros(dim, dtype=float, order='c')
        self.yd_nd  = np.zeros(dim, dtype=float, order='c')
        self.pData.y  = <void*>self.y_nd
        self.pData.yd = <void*>self.yd_nd
        
        if ROOT != None: #Sets the root function
            self.pData.ROOT = <void*>ROOT
            self.pData.dimRoot = dimRoot
            self.nbrRoot = dimRoot
            self.pData.memSizeRoot = dimRoot*sizeof(realtype)
        else:
            self.nbrRoot = 0
    
        if JAC != None: #Sets the jacobian
            self.pData.JAC = <void*>JAC
            self.pData.memSizeJac = dim*dim*sizeof(realtype)
        
        if JACV != None: #Sets the jacobian times vector
            self.pData.JACV = <void*>JACV
            
        if SENS != None: #Sets the sensitivity function
            self.pData.SENS = <void*>SENS
        
        if dimSens != None: #Sensitivity parameters (does not need the sensitivity function)
            self.pData.dimSens = dimSens
            self.pData.p = <realtype*> malloc(dimSens*sizeof(realtype))
            #self.p_nd = np.empty(dimSens, dtype=float, order='c')
            #self.pData.p = <void*>self.p_nd
            #self.pp_nd = <realtype*>(<ndarray>self.pData.p).data
        else:
            self.pData.dimSens = 0
    
    cdef inline void add_sens_point(self, realtype t):
        """
        Store a sensitivity solution point in the solution list.
        """
        if self.pData.dimSens > 0:
            self.p_result += [self.interpolate_sensitivity(t, 0)]
    
    def __dealloc__(self):
        """Free allocated data."""
        self.y_nd = None
        self.pData.y = NULL
        self.pData.yd = NULL
        if self.pData.dimSens > 0:
            free(self.pData.p) #Free the allocated data


cdef class CVode_wrap(Sundials):
    """Class to wrap CVode"""
    cdef:
        public int discr, iter, dim, _ordersum,_count_output
        public long int max_steps
        public realtype event_time
        public realtype t0
        public ndarray event_info
        public dict stats
        public dict detailed_info
        public booleantype jacobian, store_cont, comp_step, store_state
        public booleantype sim_complete
    method=['Adams','BDF']
    iteration=['Fixed Point','Newton']
    
    def __init__(self,dim):
        self.comp_step = False
        self.dim=dim
        self.discr=1
        self.iter=1
        self.store_cont = False
        self.store_state = False
        self.sim_complete = False
        
        #Default values
        self.atol = np.ones(dim)*1.e-6 #Default tolerance (absolute)
        
    def cvinit(self,t0,u, maxsteps, verbosity, switches = None):
        cdef flag
        if self.sim_complete:
            CVodeFree(&self.solver)
        self.store_state = True
        self.sim_complete = False
        self._ordersum=self._count_output=0 # initialize ordersum and output count for avarage order
        
        self.update_solver(t0, u, switches) #Update the solver
        self.update_options(verbosity, maxsteps)
        
        #Are there sensitivities to be calculated?
        if self.pData.dimSens > 0:
            self.update_sens_options(t0,u)
    
    cdef update_options(self, verbosity, maxsteps):
        """
        Updates the simulation options.
        """
        cdef flag
        
        #Verbosity
        self.pData.verbose = verbosity

        #Choose a linear solver if and only if NEWTON is choosen
        if self.linear_solver == 'DENSE' and self.iter == CV_NEWTON:
            #Specify the use of the internal dense linear algebra functions.
            flag = CVDense(self.solver, self.pData.dim)
            if flag < 0:
                raise CVodeError(flag)
                
            #Specify the jacobian to the solver
            if self.pData.JAC != NULL and self.usejac:
                flag = CVDlsSetDenseJacFn(self.solver, cv_jac)
                if flag < 0:
                    raise CVodeError(flag)
            else:
                flag = CVDlsSetDenseJacFn(self.solver, NULL)
                if flag < 0:
                    raise CVodeError(flag)
        elif self.linear_solver == 'SPGMR' and self.iter == CV_NEWTON:
            #Specify the use of CVSPGMR linear solver.
            flag = CVSpgmr(self.solver, self.pretype, self.max_krylov)
            if flag < 0:
                raise CVodeError(flag)
                
            #Specify the jacobian times vector function
            if self.pData.JACV != NULL and self.usejac:
                flag = CVSpilsSetJacTimesVecFn(self.solver, cv_jacv)
                if flag < 0:
                    raise CVodeError(flag)
            else:
                flag = CVSpilsSetJacTimesVecFn(self.solver, NULL)
                if flag < 0:
                    raise CVodeError(flag)
        else: #Functional Iteration choosen.
            pass #raise CVodeError(100,t0) #Unknown error message

        #Maximum order
        flag = CVodeSetMaxOrd(self.solver, self.maxord)
        if flag < 0:
            raise CVodeError(flag)
            
        #Initial step
        flag = CVodeSetInitStep(self.solver, self.inith)
        if flag < 0:
            raise CVodeError(flag)
            
        #Maximum step
        flag = CVodeSetMaxStep(self.solver, self.maxh)
        if flag < 0:
            raise CVodeError(flag)
            
        #Maximum Number of steps
        flag = CVodeSetMaxNumSteps(self.solver, maxsteps)
        if flag < 0:
            raise CVodeError(flag)
            
        flag = CVodeSVtolerances(self.solver, self.rtol, arr2nv(self.atol))
        if flag < 0:
            raise CVodeError(flag)
    
    cdef update_sens_options(self,t0,y):
        """
        Updates the sensitivity options.
        """
        cdef int flag
        cdef realtype *pbar
        cdef realtype ZERO = 0.0
        cdef realtype ONE = 1.0
        self.ySO  = N_VCloneVectorArray_Serial(self.pData.dimSens, arr2nv(y))
        
        #Filling the start vectors
        for i in range(self.pData.dimSens):
             N_VConst_Serial(ZERO,  self.ySO[i]);
             if self.yS0 != None:
                for j in range(self.pData.dim):
                    (<N_VectorContent_Serial>self.ySO[i].content).data[j] = self.yS0[i,j]
        
        if self._flag_active_sens:
            flag = CVodeSensReInit(self.solver, self.ism, self.ySO)
        else:
            flag = CVodeSensInit(self.solver, self.pData.dimSens, self.ism, NULL, self.ySO)
            self._flag_active_sens = True
        if flag < 0:
            raise CVodeError(flag, t0)
            
        #Sets the parameters to the userdata object.
        for i in range(self.pData.dimSens):
            self.pData.p[i] = self.p[i]
        
        #Sets the pbar.
        pbar = <realtype*> malloc(self.pData.dimSens*sizeof(realtype))
        if self.pbar != None:
            for i in range(self.pData.dimSens):
                pbar[i] = self.pbar[i]
        else:
            for i in range(self.pData.dimSens):
                pbar[i] = 1.0
        
        #Problem parameter information
        flag = CVodeSetSensParams(self.solver, self.pData.p, pbar, NULL)
        if flag < 0:
            raise CVodeError(flag, t0)
        
        free(pbar) #Free the allocated memory
        
        #Difference quotient strategy
        flag = CVodeSetSensDQMethod(self.solver, self.DQtype, self.DQrhomax)
        if flag < 0:
            raise CVodeError(flag, t0)
        
        #Maximum number of nonlinear iterations
        flag = CVodeSetSensMaxNonlinIters(self.solver, self.maxcorS)
        if flag < 0:
            raise CVodeError(flag, t0)
        
        #Specify the error control strategy
        flag = CVodeSetSensErrCon(self.solver, self.errconS)
        if flag < 0:
            raise CVodeError(flag, t0)
        
        #Estimate the sensitivity
        flag = CVodeSensEEtolerances(self.solver)
        if flag<0:
            raise CVodeError(flag, t0)
        
        #Should the sensitivities be calculated this time around?
        if self.sensToggleOff:
            flag = CVodeSensToggleOff(self.solver)
            
            if flag<0:
                raise CVodeError(flag, t0)
    
    cdef update_solver(self, t0, y0, sw0 = None):
        """
        Create or reinitiate the solver.
        """
        cdef int flag #Used for return

        self.y_cur = arr2nv(y0)
        
        #Updates the switches
        if sw0 != None:
            self.switches = sw0
            self.pData.sw = <void*>self.switches
        
        if self.solver == NULL: #The solver is not initialized
            
            self.solver=CVodeCreate(self.discr, self.iter) #Create solver
            if self.solver == NULL:
                raise CVodeError(CV_MEM_FAIL)
            
            #Specify the residual and the initial conditions to the solver
            flag = CVodeInit(self.solver, cv_rhs, t0, self.y_cur)
            if flag < 0:
                raise CVodeError(flag, t0)
                
            #Specify the root function to the solver
            if self.pData.ROOT != NULL:
                flag = CVodeRootInit(self.solver, self.pData.dimRoot, cv_root)
                if flag < 0:
                    raise CVodeError(flag,t0)
                    
            #Specify the error handling
            flag = CVodeSetErrHandlerFn(self.solver, cv_err, <void*>self.pData)
            if flag < 0:
                raise CVodeError(flag, t0)
            
        else: #The solver needs to be reinitialized
            #Reinitialize
            flag = CVodeReInit(self.solver, t0, self.y_cur)
            if flag < 0:
                raise CVodeError(flag, t0)
        
        #Set the user data
        flag = CVodeSetUserData(self.solver, <void*>self.pData)
        if flag < 0:
            raise CVodeError(flag, t0)
    
    cpdef interpolate(self, t, k):
        """
        Calls the internal CVodeGetDky for the interpolated values at time t.
        t must be within the last internal step. k is the derivative of y which
        can be from zero to the current order.
        """
        cdef flag
        cdef N_Vector dky=N_VNew_Serial(self.dim)
        
        flag = CVodeGetDky(self.solver, t, k, dky)
        
        if flag < 0:
            raise CVodeError(flag, t)
        
        return nv2arr(dky)
    
    cpdef interpolate_sensitivity(self,realtype t, int k, int i=-1):
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
                flag = CVodeGetSensDky1(self.solver, t, k, x, dkyS)
                if flag<0:
                    raise CVodeError(flag, t)
                
                matrix += [nv2arr(dkyS)]
            
            return np.array(matrix)
        else:
            flag = CVodeGetSensDky1(self.solver, t, k, i, dkyS)
            if flag <0:
                raise CVodeError(flag, t)
            
            return nv2arr(dkyS)
    
    cpdef store_statistics(self):
        """
        Retrieves and stores the statistics.
        """
        cdef long int nsteps, nrevals, njevals, nrevalsLS, ngevals, netfails, nniters, nncfails
        cdef long int nSniters, nSncfails, nfevalsLS, njvevals
        cdef long int nfSevals,nfevalsS,nSetfails,nlinsetupsS
        
        if self.store_state:
            flag = CVodeGetNumSteps(self.solver, &nsteps) #Number of steps
            flag = CVodeGetNumRhsEvals(self.solver, &nrevals) #Number of function evals
            if self.iter == 1:
                njevals = 0
                nrevalsLS = 0
            elif self.linear_solver == "SPGMR":
                flag = CVSpilsGetNumJtimesEvals(self.solver, &njvevals) #Number of jac*vector
                flag = CVSpilsGetNumRhsEvals(self.solver, &nfevalsLS) #Number of rhs due to jac*vector
                self.statistics["JVEC"] += njvevals
                self.statistics["RHSJVEC"] += nfevalsLS
            else:
                flag = CVDlsGetNumJacEvals(self.solver, &njevals) #Number of jac evals
                flag = CVDlsGetNumRhsEvals(self.solver, &nrevalsLS) #Number of res evals due to jac evals
            flag = CVodeGetNumGEvals(self.solver, &ngevals) #Number of root evals
            flag = CVodeGetNumErrTestFails(self.solver, &netfails) #Number of local error test failures
            flag = CVodeGetNumNonlinSolvIters(self.solver, &nniters) #Number of nonlinear iteration
            flag = CVodeGetNumNonlinSolvConvFails(self.solver, &nncfails) #Number of nonlinear conv failures
            
            stats_values = [nsteps, nrevals, njevals, nrevalsLS, ngevals, netfails, nniters, nncfails]
            
            for i in range(len(stats_values)):
                self.solver_stats[i] += stats_values[i]
            
            if self.pData.dimSens > 0:
                
                flag = CVodeGetSensStats(self.solver, &nfSevals, &nfevalsS, &nSetfails, &nlinsetupsS)
                flag = CVodeGetSensNonlinSolvStats(self.solver, &nSniters, &nSncfails)
                
                stats_values = [nfSevals, nfevalsS, nSetfails, nlinsetupsS,nSniters, nSncfails]
                
                for x in range(len(stats_values)):
                    self.solver_sens_stats[x] += stats_values[x]
                
        self.store_state = False
    
    cdef save_event_info(self,tret):
        """
        Store the discontinuity information in event_info and event_time
        """
        cdef int* event_info_
        cdef flag
        # Allocate memory for the event_info_ vector and initialize to zeros
        event_info_ = <int*> malloc(self.nbrRoot*sizeof(int))
        for k in range(self.nbrRoot):
            event_info_[k] = 0
        
        # Fetch data on which root functions that became zero and store in class
        flag = CVodeGetRootInfo(self.solver, event_info_)
        if flag < 0:
            raise CVodeError(flag,tret)
        
        self.event_info =  PyArray_SimpleNew(1,&self.nbrRoot ,NPY_INT)
        for k in range(self.nbrRoot):
            self.event_info[k] = event_info_[k]
        
        # Remember to deallocate
        free(event_info_)
        self.event_time=tret

    cdef inline void add_sol_point(self,realtype t, N_Vector y):
        """
        Store a solution point in the solution list.
        """
        self.sol.append((np.array(t), nv2arr(y)))
        
    cpdef run(self,realtype t0,realtype tf,float dt):
        """
        Runs the simulation.
        """
        cdef realtype tret           #Return time (not neceeserily tout)
        cdef realtype tout           #Communication time
        cdef int i,itask,nt
        cdef realtype hinused,hlast,hcur,tcur
        cdef int qlast, qcurrent
        cdef flag, solveFlag
        cdef N_Vector errW
        
        flag = CVodeSetStopTime(self.solver, tf) #Set the stop time
        if flag < 0:
            raise CVodeError(flag, t0)
        
        self.sol=[] #Reset the solution list
        self.p_result = [] #Reset the sensitivity solution list
        tret=t0
        errW = arr2nv(np.empty(self.pData.dim))
        
        
        if dt > 0.0:
            nt = int(math.ceil((tf-t0)/dt))
            for i in xrange(1, nt+1):
                tout=t0+i*dt
                
                solveFlag=CVode(self.solver,tout,self.y_cur,&tret,CV_NORMAL)
                if solveFlag < 0:
                    raise CVodeError(solveFlag, tret)
                
                self.add_sol_point(tret, self.y_cur)
                self.add_sens_point(tret)
                
                if solveFlag == CV_ROOT_RETURN: #Found a root
                    self.save_event_info(tret)
                    break

                flag = CVodeGetLastOrder(self.solver, &qlast)
                self._count_output+=1
                self._ordersum+=qlast
                avar=float(self._ordersum)/self._count_output
                if solveFlag == 1:
                    break
                if i == nt:
                    solveFlag=1
                if self.store_cont:
                    break
                    
        else: # one step mode
            if self.detailed_info == None:
                self.detailed_info = {}
                self.detailed_info['qlast'] = []
                self.detailed_info['qcurrent'] = []
                self.detailed_info['hlast'] = []
                self.detailed_info['errorWeights'] = []
            while tret < tf:
                solveFlag=CVode(self.solver,tf,self.y_cur,&tret,CV_ONE_STEP)
                if solveFlag < 0:
                    raise CVodeError(solveFlag, tret)
                
                self.add_sol_point(tret, self.y_cur)
                self.add_sens_point(tret)
                
                if solveFlag == CV_ROOT_RETURN: #Found a root
                    self.save_event_info(tret)
                    break
                
                if self.save_detailed_info:
                    flag = CVodeGetLastOrder(self.solver, &qlast)
                    flag = CVodeGetCurrentOrder(self.solver, &qcurrent)
                    flag = CVodeGetLastStep(self.solver, &hlast)
                    flag = CVodeGetEstLocalErrors(self.solver, errW)
                    self.detailed_info['qlast'].append(qlast)
                    self.detailed_info['qcurrent'].append(qcurrent)
                    self.detailed_info['hlast'].append(hlast)
                    self.detailed_info['errorWeights'].append(nv2arr(errW))
                self._count_output+=1
                self._ordersum+=qlast
                avar=float(self._ordersum)/self._count_output
                if self.comp_step:
                    break
                if self.store_cont:
                    break
            else:
                solveFlag = 1
                
        if solveFlag >= 1:
            self.store_statistics()
        if solveFlag >= 1 and solveFlag != CV_ROOT_RETURN:
            self.sim_complete = True

        return self.sol
    
    def __dealloc__(self):
        pass
        #if self.solver != NULL:
        #    CVodeFree(&self.solver)


cdef class IDA_wrap(Sundials):
    """Class to wrap Sundials IDA"""
    cdef:
        #void* comp_step_method
        public int dim, _ordersum,_count_output
        public realtype event_time
        public realtype t0
        public ndarray algvar,event_info
        public dict detailed_info
        public booleantype jacobian, store_cont,store_state,comp_step
        public booleantype sens_activated
        public int icopt,Ns
        public booleantype sim_complete
        N_Vector temp_nvector
        N_Vector *ydSO
    def __init__(self,dim):
        
        self.dim=dim
        self.store_cont = False
        self.store_state = False
        self.comp_step = False
        self.sim_complete = False
        
        #Default values (No sensitivities)
        self.algvar = np.ones(dim)
        self.atol = np.ones(dim)*1.e-6 #Default tolerance (absolute)
        
    def idinit(self,t0,u,ud, maxsteps, verbosity, switches = None):
        cdef flag
        self.t0 = t0
        self.store_state = True
        self.sim_complete = False
        self._ordersum=self._count_output=0 # initialize ordersum and output count for avarage order
        
        self.update_solver(t0, u, ud, switches)
        self.update_options(verbosity, maxsteps)
        
        #Are there sensitivities to be calculated?
        if self.pData.dimSens > 0:
            self.update_sens_options(t0,u,ud)
    
    cdef update_options(self, verbosity, maxsteps):
        """
        Updates the simulation options.
        """
        cdef flag
        
        #Verbosity
        self.pData.verbose = verbosity
        
        #Maximum order
        flag =IDASetMaxOrd(self.solver, self.maxord)
        if flag < 0:
            raise IDAError(flag)
            
        #Initial step
        flag = IDASetInitStep(self.solver, self.inith)
        if flag < 0:
            raise IDAError(flag)
            
        #Maximum step
        flag = IDASetMaxStep(self.solver, self.maxh)
        if flag < 0:
            raise IDAError(flag)
            
        #Maximum Number of steps
        flag = IDASetMaxNumSteps(self.solver, maxsteps)
        if flag < 0:
            raise IDAError(flag)
            
        #Maximum Number of error test failures
        flag = IDASetMaxErrTestFails(self.solver, 20)
        if flag < 0:
            raise IDAError(flag) 
        
        #Set the algebraic components and the differential
        flag = IDASetId(self.solver, arr2nv(self.algvar))
        if flag < 0:
            raise IDAError(flag)
        
        #Suppress algebraic components on the error test
        flag = IDASetSuppressAlg(self.solver, self.suppress_alg)
        if flag < 0:
            raise IDAError(flag)
            
        #Set the sensitivity
        flag = IDASVtolerances(self.solver, self.rtol, arr2nv(self.atol))
        if flag < 0:
            raise IDAError(flag)
    
    cdef update_solver(self, t0, y0, yd0, sw0 = None):
        """
        Create or reinitiate the solver.
        """
        cdef int flag #Used for return

        self.y_cur  = arr2nv(y0)
        self.yd_cur = arr2nv(yd0)
        
        #Updates the switches
        if sw0 != None:
            self.switches = sw0
            self.pData.sw = <void*>self.switches
        
        if self.solver == NULL: #The solver is not initialized
        
            self.solver=IDACreate() #Create solver
            if self.solver == NULL:
                raise IDAError(IDA_MEM_FAIL)
            
            #Specify the residual and the initial conditions to the solver
            flag = IDAInit(self.solver, ida_res, t0, self.y_cur, self.yd_cur)
            if flag < 0:
                raise IDAError(flag, t0)
                
            #Specify the use of the internal dense linear algebra functions.
            flag = IDADense(self.solver, self.pData.dim)
            if flag < 0:
                raise IDAError(flag, t0)
            
            #Specify the root function to the solver
            if self.pData.ROOT != NULL:
                flag = IDARootInit(self.solver, self.pData.dimRoot, ida_root)
                
                if flag < 0:
                    raise IDAError(flag,t0)
            
            #Specify the error handling
            flag = IDASetErrHandlerFn(self.solver, ida_err, <void*>self.pData)
            if flag < 0:
                raise IDAError(flag, t0)
            
        else: #The solver needs to be reinitialized
            
            #Reinitialize
            flag = IDAReInit(self.solver, t0, self.y_cur, self.yd_cur)
            if flag < 0:
                raise IDAError(flag, t0)
        
        #Specify the jacobian to the solver
        if self.pData.JAC != NULL and self.usejac:
            flag = IDADlsSetDenseJacFn(self.solver, ida_jac)
            if flag < 0:
                raise IDAError(flag,t0)
        else:
            flag = IDADlsSetDenseJacFn(self.solver, NULL)
            if flag < 0:
                raise IDAError(flag,t0)
        
        #Set the user data
        flag = IDASetUserData(self.solver, <void*>self.pData)
        if flag < 0:
            raise IDAError(flag, t0)
    
    def interpolate(self, t, k):
        """
        Calls the internal IDAGetDky for the interpolated values at time t.
        t must be within the last internal step. k is the derivative of y which
        can be from zero to the current order.
        """
        cdef flag
        cdef N_Vector dky=N_VNew_Serial(self.dim)
        
        flag = IDAGetDky(self.solver, t, k, dky)
        
        if flag < 0:
            raise IDAError(flag, t)
        
        return nv2arr(dky)
    
    cdef update_sens_options(self,t0,y,yd):
        """
        Sets the sensitivity information.
        """
        #Create the initial matrices
        self.ySO  = N_VCloneVectorArray_Serial(self.pData.dimSens, self.y_cur)
        self.ydSO = N_VCloneVectorArray_Serial(self.pData.dimSens, self.yd_cur)
        cdef realtype ZERO = 0.0
        cdef realtype *pbar
        
        #Filling the start vectors
        for i in range(self.pData.dimSens):
             N_VConst_Serial(ZERO,  self.ySO[i]);
             N_VConst_Serial(ZERO, self.ydSO[i]);
             if self.yS0 != None:
                for j in range(self.pData.dim):
                    (<N_VectorContent_Serial>self.ySO[i].content).data[j] = self.yS0[i,j]

        if self._flag_active_sens:
            flag = IDASensReInit(self.solver, self.ism, self.ySO, self.ydSO)
        else:
            flag = IDASensInit(self.solver, self.pData.dimSens, self.ism, NULL, self.ySO, self.ydSO)
            self._flag_active_sens = True
        if flag<0:
            raise IDAError(flag, t0)    
            
         
        #Sets the parameters to the userdata object.
        for i in range(self.pData.dimSens):
            self.pData.p[i] = self.p[i]
        
        #Sets the pbar.
        pbar = <realtype*> malloc(self.pData.dimSens*sizeof(realtype))
        if self.pbar != None:
            for i in range(self.pData.dimSens):
                pbar[i] = self.pbar[i]
        else:
            for i in range(self.pData.dimSens):
                pbar[i] = 1.0
        
        #Specify problem parameter information for sensitivity calculations
        flag = IDASetSensParams(self.solver, self.pData.p, pbar, NULL)
        
        if self.pbar != None:
            free(pbar) #Free the allocated space.
            
        if flag<0:
            raise IDAError(flag, t0)
        
        #Specify the difference quotient strategy
        flag = IDASetSensDQMethod(self.solver, self.DQtype, self.DQrhomax)
        if flag<0:
            raise IDAError(flag, t0)
        
        #Specify the error control strategy
        flag = IDASetSensErrCon(self.solver, self.errconS)
        if flag<0:
            raise IDAError(flag, t0)
        
        #Specify the maximum number of nonlinear solver iterations
        flag = IDASetSensMaxNonlinIters(self.solver, self.maxcorS)
        if flag<0:
            raise IDAError(flag, t0)
        
        #Estimate the sensitivity  ----SHOULD BE IMPROVED with IDASensSVTolerances ...
        flag = IDASensEEtolerances(self.solver)
        if flag<0:
            raise IDAError(flag, t0)
        
        #Should the sensitivities be calculated this time around?
        if self.sensToggleOff:
            flag = IDASensToggleOff(self.solver)
            if flag<0:
                raise IDAError(flag, t0)
    
    def interpolate_sensitivity(self,realtype t, int k, int i=-1):
        """
        This method calls the internal method IDAGetSensDky which computes the k-th derivatives
        of the interpolating polynomials for the sensitivity variables at time t.
        
            Parameters::
                    
                    t
                        - Specifies the time at which sensitivity information is requested. The time
                          t must fall within the interval defined by the last successful step taken
                          by IDAS.
                    
                    k   
                        - The order of derivatives.
                        
                    i
                        - Specifies the sensitivity derivative vector to be returned (0<=i<=Ns)
                        
            Return::
            
                    A matrix containing the Ns vectors or a vector if i is specified.
        """
        cdef N_Vector dkyS=N_VNew_Serial(self.dim)
        cdef flag
        
        if i==-1:
            
            matrix = []
            
            for x in range(self.pData.dimSens):
                flag = IDAGetSensDky1(self.solver, t, k, x, dkyS)
                
                if flag<0:
                    raise IDAError(flag, t)
                
                matrix += [nv2arr(dkyS)]
            
            return np.array(matrix)
        else:
            flag = IDAGetSensDky1(self.solver, t, k, i, dkyS)
            
            if flag <0:
                raise IDAError(flag, t)
            
            return nv2arr(dkyS)

        
    def store_statistics(self):
        """
        Retrieves and stores the statistics.
        """
        cdef long int nsteps, nrevals,njevals,nrevalsLS,ngevals,netfails,nniters,nncfails
        cdef long int nSniters, nSncfails
        cdef long int nfSevals,nfevalsS,nSetfails,nlinsetupsS
        
        if self.store_state:
            flag = IDAGetNumSteps(self.solver, &nsteps) #Number of steps
            flag = IDAGetNumResEvals(self.solver, &nrevals) #Number of res evals
            flag = IDADlsGetNumJacEvals(self.solver, &njevals) #Number of jac evals
            flag = IDADlsGetNumResEvals(self.solver, &nrevalsLS) #Number of res evals due to jac evals
            flag = IDAGetNumGEvals(self.solver, &ngevals) #Number of root evals
            flag = IDAGetNumErrTestFails(self.solver, &netfails) #Number of local error test failures
            flag = IDAGetNumNonlinSolvIters(self.solver, &nniters) #Number of nonlinear iteration
            flag = IDAGetNumNonlinSolvConvFails(self.solver, &nncfails) #Number of nonlinear conv failures
            
            stats_values = [nsteps, nrevals, njevals, nrevalsLS, ngevals, netfails, nniters, nncfails]
            
            for i in range(len(stats_values)):
                self.solver_stats[i] += stats_values[i]
            
            if self.pData.dimSens > 0:
                
                flag = IDAGetSensStats(self.solver, &nfSevals, &nfevalsS, &nSetfails, &nlinsetupsS)
                flag = IDAGetSensNonlinSolvStats(self.solver, &nSniters, &nSncfails)
                
                stats_values = [nfSevals, nfevalsS, nSetfails, nlinsetupsS,nSniters, nSncfails]
                
                for x in range(len(stats_values)):
                    self.solver_sens_stats[x] += stats_values[x]
            
            
        self.store_state = False
        
    def calc_IC(self,method, direction, lsoff):
        """
        This calculates the initial conditions with the built in SUNDIALS
        solver. See IDA Documentation section 4.5.4.
        """
        cdef realtype tout1
        if method == 'IDA_Y_INIT':
            self.icopt = IDA_Y_INIT
        if method == 'IDA_YA_YDP_INIT':
            self.icopt = IDA_YA_YDP_INIT 
        
        tout1 = self.t0+direction #tout1 is needed for the solver to determine the direction of the integration
        if self.solver == NULL:
            raise Exception, "IDA must be initialized"
        
        flag = IDASetLineSearchOffIC(self.solver, lsoff)
        
        flag = IDACalcIC(self.solver, self.icopt, tout1)
        
        if flag == IDA_SUCCESS: #Gets the calculated values
            flag = IDAGetConsistentIC(self.solver, self.y_cur, self.yd_cur)
        
        return [flag, nv2arr(self.y_cur), nv2arr(self.yd_cur)]
        
    
    cdef save_event_info(self,tret):
        """
        Store the discontinuity information in event_info and event_time
        """
        cdef int* event_info_
        cdef flag
        
        # Allocate memory for the event_info_ vector and initialize to zeros
        event_info_ = <int*> malloc(self.nbrRoot*sizeof(int))
        for k in range(self.nbrRoot):
            event_info_[k] = 0
            # Fetch data on which root functions that became zero and store in class
        flag = IDAGetRootInfo(self.solver, event_info_)
        if flag < 0:
            raise IDAError(flag, tret)
            
        self.event_info =  PyArray_SimpleNew(1,&self.nbrRoot ,NPY_INT)
        for k in range(self.nbrRoot):
            self.event_info[k] = event_info_[k]
        # Remember to deallocat
        free(event_info_)
        self.event_time=tret
    
    cdef inline void add_sol_point(self,realtype t, N_Vector y, N_Vector yd):
        """
        Store a solution point in the solution list.
        """
        self.sol.append((np.array(t), nv2arr(y), nv2arr(yd)))
    
    cpdef run(self,realtype t0,realtype tf,float dt):
        #cdef realtype dt             # time increment
        cdef realtype tret           # return time (not neceeserily tout)
        cdef realtype tout           # communication time
        cdef int i,itask, nt
        cdef flag, solveFlag
        #cdef realtype hinused,hlast,hcur,tcur
        #cdef long int nsteps, fevals, nlinsetups, netfails
        cdef int  qlast, qcurrent
        cdef realtype hlast
        cdef N_Vector errW
        flag = IDASetStopTime(self.solver, tf)
        self.sol=[]
        self.p_result = []
        tret=t0
        errW = arr2nv(np.empty(self.pData.dim))
        
        if dt > 0.0:
            nt = int(math.ceil((tf-t0)/dt))
            for i in xrange(1, nt+1):
                tout=t0+i*dt
                flag=0
                flags=0
                solveFlag=IDASolve(self.solver,tout,&tret, self.y_cur, self.yd_cur,IDA_NORMAL)
                if solveFlag < 0:
                    raise IDAError(solveFlag, tret)
                
                self.add_sol_point(tret, self.y_cur, self.yd_cur)
                self.add_sens_point(tret)
                
                if solveFlag == IDA_ROOT_RETURN: #Found a root
                    self.save_event_info(tret)
                    break

                flag = IDAGetLastOrder(self.solver, &qlast)

                if solveFlag == 1:
                    break
                if i == nt:
                    solveFlag =1
                if self.store_cont:
                    break
        else: # one step mode
            if self.detailed_info == None:
                self.detailed_info = {}
                self.detailed_info['qlast'] = []
                self.detailed_info['qcurrent'] = []
                self.detailed_info['hlast'] = []
                self.detailed_info['errorWeights'] = []
            while tret < tf:
                
                solveFlag=IDASolve(self.solver,tf,&tret, self.y_cur,self.yd_cur,IDA_ONE_STEP)
                if solveFlag < 0:
                    raise IDAError(solveFlag, tret)

                self.add_sol_point(tret, self.y_cur, self.yd_cur)
                self.add_sens_point(tret)
                
                if solveFlag == IDA_ROOT_RETURN: #Found a root
                    self.save_event_info(tret)
                    break
                
                if self.save_detailed_info:
                    flag = IDAGetLastOrder(self.solver, &qlast)
                    flag = IDAGetCurrentOrder(self.solver, &qcurrent)
                    flag = IDAGetLastStep(self.solver, &hlast)
                    flag = IDAGetEstLocalErrors(self.solver, errW)
                    self.detailed_info['qlast'].append(qlast)
                    self.detailed_info['qcurrent'].append(qcurrent)
                    self.detailed_info['hlast'].append(hlast)
                    self.detailed_info['errorWeights'].append(nv2arr(errW))

                if self.store_cont:
                    break
                if self.comp_step:
                    break
            else:
                solveFlag=1
        if solveFlag >= 1:
            self.store_statistics()
        if solveFlag >= 1 and solveFlag != IDA_ROOT_RETURN:
            self.sim_complete = True
        # Free memory
        #IDAFree(&self.solver)
        #N_VDestroy_Serial(self.curr_state)
        #N_VDestroy_Serial(self.curr_deriv)
        
        return self.sol                
            


        
    
            
        
