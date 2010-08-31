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
#from sundials_core cimport *
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
    Wraps  Python rhs-callback function to obtain CVode required interface
    see also ctypedef statement above
    """
    cdef ProblemData *pData = <ProblemData*>problem_data
    y=nv2arr(yv)
    cdef realtype* ydotptr=(<N_VectorContent_Serial>yvdot.content).data
    cdef long int n=(<N_VectorContent_Serial>yv.content).length

    try:
        if pData.sw != NULL:
            ydot=(<object>pData.RHS)(t,y,<list>pData.sw) #Call the Python rhs function
        else:
            ydot=(<object>pData.RHS)(t,y) #Call the Python rhs function
        for i in range(n):
            ydotptr[i]=ydot[i]
        return CV_SUCCESS
    except:
        return CV_REC_ERR #Recoverable Error (See Sundials description)

cdef int cv_jac(int Neq, realtype t, N_Vector yv, N_Vector fy, DlsMat Jac, 
                void *problem_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    Wraps Python jacobian-callback function to obtain CV required interface.
    """
    cdef ProblemData *pData = <ProblemData*>problem_data
    y = nv2arr(yv)

    cdef realtype* col_i=DENSE_COL(Jac,0)
    try:
        if pData.sw != NULL:
            jacobian=(<object>pData.JAC)(t,y,<list>pData.sw)  # call to the python residual function
        else:
            jacobian=(<object>pData.JAC)(t,y)
        
        for i in range(Neq):
            col_i = DENSE_COL(Jac, i)
            for j in range(Neq):
                col_i[j] = jacobian[j,i]
        return CVDLS_SUCCESS
    except:
        return CVDLS_JACFUNC_RECVR #Recoverable Error (See Sundials description)

cdef int cv_root(realtype t, N_Vector yv, realtype *gout,  void* problem_data):
    """
    This method is used to connect the Assimulo.Problem.state_events to the Sundials
    Root-finding function.
    """
    cdef ProblemData *pData = <ProblemData*>problem_data
    cdef ndarray[realtype, ndim=1, mode='c'] root #Used for return from the user function
    (<ndarray>pData.y).data =  <char*>((<N_VectorContent_Serial>yv.content).data)

    try:
        if pData.sw != NULL:
            root=(<object>pData.ROOT)(t,(<ndarray>pData.y),<list>pData.sw) #Call to the Python root function 
        else:
            root=(<object>pData.ROOT)(t,(<ndarray>pData.y),None) #Call to the Python root function
            
        memcpy(gout,<realtype*>root.data,pData.memSizeRoot) #Copy data from the return to the output
    
        return CV_SUCCESS
    except:
        return CV_RTFUNC_FAIL  # Unrecoverable Error
        
cdef int ida_res(realtype t, N_Vector yv, N_Vector yvdot, N_Vector residual, void* user_data):
    """
    Wraps  Python res-callback function to obtain IDA required interface
    see also ctypedef statement above
    """
    y=nv2arr(yv)
    yd=nv2arr(yvdot)
    try:
        switch=(<object> (<ProblemData*> user_data).data)[IDA_ROOT_IND][IDA_SW_IND]
    except:
        switch=False
    cdef realtype* resptr=(<N_VectorContent_Serial>residual.content).data
    cdef long int n=(<N_VectorContent_Serial>yv.content).length
    
    nparam = (<object> (<ProblemData*> user_data).data)[IDA_RES_IND][IDA_SENS_IND]
    if nparam!=0: #SENSITIVITY
        p = realtype2arr((<ProblemData*> user_data).params,nparam)
        try:
            if switch:
                res=(<object> (<ProblemData*> user_data).data)[IDA_RES_IND][IDA_RESF_IND](t,y,yd,sw=switch,p=p)  # call to the python residual function
            else:
                res=(<object> (<ProblemData*> user_data).data)[IDA_RES_IND][IDA_RESF_IND](t,y,yd,p)
            for i in range(n):
                resptr[i]=res[i]
            return 0
        except(LinAlgError,ZeroDivisionError):
            return 1 # recoverable error (see Sundials description)
        except:
            print "Unexpected error, probably due to a programing error in rhs/res function:\n"
            traceback.print_exc()
            return -1
    else: #NO SENSITIVITY
        try:
            if switch:
                res=(<object> (<ProblemData*> user_data).data)[IDA_RES_IND][IDA_RESF_IND](t,y,yd,switch)  # call to the python residual function
            else:
                res=(<object> (<ProblemData*> user_data).data)[IDA_RES_IND][IDA_RESF_IND](t,y,yd)
            for i in range(n):
                resptr[i]=res[i]
            return 0
        except(LinAlgError,ZeroDivisionError):
            return 1 # recoverable error (see Sundials description)
        except:
            print "Unexpected error, probably due to a programing error in rhs/res function:\n"
            traceback.print_exc()
            return -1
            
cdef int ida_jac(int Neq, realtype t, realtype c, N_Vector yv, N_Vector yvdot, N_Vector residual, DlsMat Jac,
                 void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    Wraps Python jacobian-callback function to obtain IDA required interface.
    """
    y = nv2arr(yv)
    yd = nv2arr(yvdot)
    try:
        switch=(<object> (<ProblemData*> user_data).data)[IDA_ROOT_IND][IDA_SW_IND]
    except:
        switch=False
    cdef realtype* col_i=DENSE_COL(Jac,0)
    try:
        if switch:
            jacobian=(<object> (<ProblemData*> user_data).data)[IDA_RES_IND][IDA_JAC_IND](c,t,y,yd,switch)  # call to the python residual function
        else:
            jacobian=(<object> (<ProblemData*> user_data).data)[IDA_RES_IND][IDA_JAC_IND](c,t,y,yd)
        
        for i in range(Neq):
            col_i = DENSE_COL(Jac, i)
            for j in range(Neq):
                col_i[j] = jacobian[j,i]
        return 0
    except: #None recoverable
        return -1

cdef int ida_root(realtype t, N_Vector yv, N_Vector yvdot, realtype *gout,  void* user_data):
    """
    Wraps  Python root-callback function to obtain IDA required interface
    see also ctypedef statement above
    """
    y=nv2arr(yv)
    yd=nv2arr(yvdot)
    try:
        switch=(<object> (<ProblemData*> user_data).data)[IDA_ROOT_IND][IDA_SW_IND]
    except:
        switch=False
    try:
        rootf=(<object> (<ProblemData*> user_data).data)[IDA_ROOT_IND][IDA_ROOTF_IND](t,y,yd,switch)  # call to the python root function 
        rootf = np.asarray(rootf).reshape(-1) # Make sure we get a vector
        for i in range(rootf.shape[0]):
            gout[i]=rootf[i]
        return 0
    except:
        return 1 # generates an error of type IDA_RTFUNC_FAIL
    

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
    cdef void* solver
    cdef ProblemData pData #A struct containing information about the problem
    cdef ProblemData *ppData #A pointer to the problem data
    cdef ndarray y_nd
    
    def __cinit__(self):
        self.pData = ProblemData() 
        self.ppData = &self.pData
    
    cpdef set_problem_info(self, RHS, dim, ROOT = None, dimRoot = None, JAC = None, SENS = None):
        """
        Sets the problem information to the problem struct.
        """
        #Sets the residual or rhs
        self.pData.RHS = <void*>RHS
        self.pData.dim = dim
        self.pData.memSize = dim*sizeof(realtype)
        
        if ROOT != None: #Sets the root function
            self.pData.ROOT = <void*>ROOT
            self.pData.dimRoot = dimRoot
            self.pData.memSizeRoot = dimRoot*sizeof(realtype)
    
        if JAC != None: #Sets the jacobian
            self.pData.JAC = <void*>JAC
            self.pData.memSizeJac = dim*dim*sizeof(realtype)
            
        if SENS != None: #Sets the sensitivity function
            self.pData.SENS = <void*>SENS


cdef class CVode_wrap(Sundials):
    """Class to wrap CVode"""
    cdef:
        public int discr, iter, dim, maxord, _ordersum,_count_output, max_h
        public long int max_steps
        public realtype abstol,reltol,event_time
        public realtype t0
        public ndarray abstol_ar,event_info
        public dict stats
        public dict detailed_info
        public booleantype jacobian, store_cont, comp_step, store_state
        public npy_intp num_state_events
        public booleantype sim_complete
        public switches
        #void* comp_step_method
        N_Vector curr_state, y_cur
        N_Vector temp_nvector
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
    def cvinit(self,t0,user_data,u,maxord, max_steps, init_step, switches = None):
        cdef flag
        self.store_state = True
        self.sim_complete = False
        self.max_steps = max_steps
        self._ordersum=self._count_output=0 # initialize ordersum and output count for avarage order
        
        self.update_solver(t0, u, switches) #Update the solver
        
        if self.abstol_ar[0] > 0:
            flag = CVodeSVtolerances(self.solver, self.reltol, arr2nv(self.abstol_ar))
        else:
            flag = CVodeSStolerances(self.solver, self.reltol, self.abstol)
        if flag!=CV_SUCCESS:
                raise Exception,"CVode Tolerance Initialization  Error"
        if maxord:
            flag=CVodeSetMaxOrd(self.solver, maxord)
        flag = CVodeSetMaxNumSteps(self.solver, self.max_steps)
        flag = CVodeSetInitStep(self.solver, init_step)
        
        try:
            flag= CVodeSetMaxStep(self.solver, self.max_h)
        except AttributeError:
            pass
    
    cdef update_solver(self, t0, y0, sw0 = None):
        """
        Create or reinitiate the solver.
        """
        cdef int flag #Used for return
        cdef ndarray[realtype, ndim=1, mode='c'] yy = y0
        
        self.y_nd = yy
        self.y_cur = arr2nv(y0)
        
        #Set the ndarray to the problem struct
        self.pData.y = <void*>self.y_nd
        
        #Updates the switches
        if sw0 != None:
            self.switches = sw0
            self.pData.sw = <void*>self.switches
        
        if self.solver == NULL: #The solver is not initialized
        
            self.solver=CVodeCreate(self.discr, self.iter) #Create solver
            if self.solver == NULL:
                CVodeError(CV_MEM_FAIL)
            
            #Specify the residual and the initial conditions to the solver
            flag = CVodeInit(self.solver, cv_rhs, t0, self.y_cur)
            if flag < 0:
                CVodeError(flag, t0)
                
            #Specify the use of the internal dense linear algebra functions.
            flag = CVDense(self.solver, self.pData.dim)
            if flag < 0:
                CVodeError(flag, t0)
            
            #Specify the root function to the solver
            if self.pData.ROOT != NULL:
                flag = CVodeRootInit(self.solver, self.pData.dimRoot, cv_root)
                
                if flag < 0:
                    CVodeError(flag,t0)
            
            #Specify the jacobian to the solver
            if self.pData.JAC != NULL:
                flag = CVDlsSetDenseJacFn(self.solver, cv_jac)
                if flag < 0:
                    CVodeError(flag,t0)
            
        else: #The solver needs to be reinitialized
            
            #Reinitialize
            flag = CVodeReInit(self.solver, t0, self.y_cur)
            if flag < 0:
                CVodeError(flag, t0)
        
        #Set the user data
        flag = CVodeSetUserData(self.solver, <void*>self.ppData)
        if flag < 0:
            CVodeError(flag, t0)
    
    def interpolate(self, t, k):
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
    
    def store_statistics(self):
        """
        Retrieves and stores the statistics.
        """
        cdef long int nsteps, nrevals, njevals, nrevalsLS, ngevals, netfails, nniters, nncfails
        
        if self.store_state:
            flag = CVodeGetNumSteps(self.solver, &nsteps) #Number of steps
            flag = CVodeGetNumRhsEvals(self.solver, &nrevals) #Number of function evals
            if self.iter == 1:
                njevals = 0
                nrevalsLS = 0
            else:
                flag = CVDlsGetNumJacEvals(self.solver, &njevals) #Number of jac evals
                flag = CVDlsGetNumRhsEvals(self.solver, &nrevalsLS) #Number of res evals due to jac evals            
            flag = CVodeGetNumGEvals(self.solver, &ngevals) #Number of root evals
            flag = CVodeGetNumErrTestFails(self.solver, &netfails) #Number of local error test failures
            flag = CVodeGetNumNonlinSolvIters(self.solver, &nniters) #Number of nonlinear iteration
            flag = CVodeGetNumNonlinSolvConvFails(self.solver, &nncfails) #Number of nonlinear conv failures
            
            stats_values = [nsteps, nrevals, njevals, nrevalsLS, ngevals, netfails, nniters, nncfails]
            stats_text = ['Number of Steps                          ',
                          'Number of Function Evaluations           ',
                          'Number of Jacobian Evaluations           ',
                          'Number of F-Eval During Jac-Eval         ',
                          'Number of Root Evaluations               ',
                          'Number of Error Test Failures            ',
                          'Number of Nonlinear Iterations           ',
                          'Number of Nonlinear Convergence Failures ']
            
            if self.stats != None:
                for x in range(len(stats_text)):
                    self.stats[stats_text[x]] += stats_values[x]
            else:
                self.stats = {}
                for x in range(len(stats_text)):
                    self.stats[stats_text[x]] = stats_values[x]
                    
        self.store_state = False
    
    def treat_disc(self,flag,tret):
        cdef int* event_info_
        """
        Treats solver returns with root_found flag
        """
        if flag == CV_ROOT_RETURN:
            # Allocate memory for the event_info_ vector and initialize to zeros
            event_info_ = <int*> malloc(self.num_state_events*sizeof(int))
            for k in range(self.num_state_events):
                event_info_[k] = 0
                # Fetch data on which root functions that became zero and store in class
            flag = CVodeGetRootInfo(self.solver, event_info_)
            self.event_info =  PyArray_SimpleNew(1,&self.num_state_events ,NPY_INT)
            for k in range(self.num_state_events):
                self.event_info[k] = event_info_[k]
            # Remember to deallocate
            free(event_info_)
            self.event_time=tret
            return True
        else:
            return False

        
    def run(self,realtype t0,realtype tf,float dt):
        #cdef realtype dt             # time increment
        cdef realtype tret           # return time (not neceeserily tout)
        cdef realtype tout           # communication time
        cdef int i,itask,nt
        cdef realtype hinused,hlast,hcur,tcur
        #cdef long int nsteps, fevals, nlinsetups, netfails
        cdef int  qlast, qcurrent
        flag = CVodeSetStopTime(self.solver, tf)
        sol=[]
        tret=t0
        if dt > 0.0:
            nt = int(math.ceil((tf-t0)/dt))
            for i in xrange(1, nt+1):
                tout=t0+i*dt
                flag=0
                flags=CVode(self.solver,tout,self.y_cur,&tret,CV_NORMAL)
                if flags<0 and flags!=CV_TSTOP_RETURN:
                    raise CVodeError(flags, tret)
                sol.append((np.array(tret),nv2arr(self.y_cur)))
                flag = CVodeGetLastOrder(self.solver, &qlast)
                self._count_output+=1
                self._ordersum+=qlast
                avar=float(self._ordersum)/self._count_output
                if flags == 1:
                    break
                if self.treat_disc(flags,tret):
                    break
                if i == nt:
                    flags=1
                if self.store_cont:
                    break
        else: # one step mode
            if self.detailed_info == None:
                self.detailed_info = {}
                self.detailed_info['qlast'] = []
                self.detailed_info['qcurrent'] = []
            while tret < tf:
                flag=0
                flags=CVode(self.solver,tf,self.y_cur,&tret,CV_ONE_STEP)
                if flags<0 and flags!=CV_TSTOP_RETURN:
                    raise CVodeError(flags, tret)
                sol.append((np.array(tret),nv2arr(self.y_cur)))
                flag = CVodeGetLastOrder(self.solver, &qlast)
                flag = CVodeGetCurrentOrder(self.solver, &qcurrent)
                self.detailed_info['qlast'].append(qlast)
                self.detailed_info['qcurrent'].append(qcurrent)
                self._count_output+=1
                self._ordersum+=qlast
                avar=float(self._ordersum)/self._count_output
                if self.treat_disc(flags,tret):
                    break
                if self.comp_step:
                    break
                if self.store_cont:
                    break
            else:
                flags = 1
                
        if flags >= 1:
            self.sim_complete = True
            self.store_statistics()

        return sol
    
    def __dealloc__(self):
        if self.solver != NULL:
            CVodeFree(&self.solver)


cdef class IDA_wrap(Sundials):
    """Class to wrap Sundials IDA"""
    cdef:
        void* mem
        #void* comp_step_method
        public int dim, maxord, _ordersum,_count_output, max_h
        public long int max_steps
        public realtype abstol,reltol,event_time
        public realtype t0, DQrhomax
        public ndarray abstol_ar,algvar,event_info
        public dict stats
        public list sens_stats
        public dict detailed_info
        public booleantype suppress_alg,jacobian, store_cont,store_state,comp_step
        public booleantype sens_activated
        public int icopt, nbr_params, DQtype, maxcorS, ism, Ns
        public npy_intp num_state_events
        public booleantype sim_complete, errconS, sensToggleOff
        N_Vector curr_state
        N_Vector curr_deriv
        N_Vector temp_nvector
        N_Vector *ySO, *ydSO
        ProblemData *uDataT
        cdef ProblemData tempStruct
        public object p, pbar
    def __init__(self,dim):
        
        self.dim=dim
        self.store_cont = False
        self.store_state = False
        self.comp_step = False
        self.sim_complete = False
        self.nbr_params = 0
        
        #Default values
        self.DQtype   = IDA_CENTERED #Specifies the difference quotient type
        self.DQrhomax = 0.0 #Positive value of the selction parameter used in deciding switching
        self.errconS  = False #Specifies whether sensitivity variables are included in the error control mechanism
        self.maxcorS  = 3 #Maximum number of nonlinear solver iterations for sensitivity variables per step.
        self.sens_activated = False #The sensitivities are not allocated
        self.ism = IDA_STAGGERED #The corrector step for the sensitivity variables takes place at the same time for all sensitivity equations
        self.sensToggleOff = False #Toggle the sensitivity calculations off
        self.uDataT = &self.tempStruct
        
    def __del__(self):
        free(self.uDataT.params)
        
    def idinit(self,t0,user_data,u,ud,maxord, max_steps, init_step, max_h):
        cdef flag
        self.uDataT.data = <void*>user_data
        self.t0 = t0
        self.store_state = True
        self.sim_complete = False
        self.curr_state=arr2nv(u)
        self.curr_deriv=arr2nv(ud)
        self.max_steps = max_steps
        self._ordersum=self._count_output=0 # initialize ordersum and output count for avarage order
        if self.mem == NULL:
            # Newinitialization
            self.mem=IDACreate()
            if self.mem == NULL:
                raise Exception, 'IDA: Memory allocation failed'
            flag=IDAInit(self.mem, ida_res, t0, self.curr_state, self.curr_deriv)
            if flag!=IDA_SUCCESS:
                    raise Exception,"IDA Initialization Error"
            if self.num_state_events>0: 
                flag = IDARootInit(self.mem, self.num_state_events, ida_root)
                if flag!=IDA_SUCCESS:
                    raise Exception,"IDA root-finding initialization error"
        else:
            flag = IDAReInit(self.mem, t0, self.curr_state, self.curr_deriv)
        if self.abstol_ar[0] > 0:
            flag = IDASVtolerances(self.mem, self.reltol, arr2nv(self.abstol_ar))
        else:
            flag = IDASStolerances(self.mem, self.reltol, self.abstol)
        if flag!=IDA_SUCCESS:
                raise Exception,"IDA Tolerance Initialization  Error"
        if maxord:
            flag=IDASetMaxOrd(self.mem, maxord)
        flag = IDASetMaxNumSteps(self.mem, self.max_steps)
        flag = IDASetInitStep(self.mem, init_step)
        flag = IDASetMaxStep(self.mem, max_h)
        flag = IDADense(self.mem, self.dim)
        if self.jacobian:
            flag = IDADlsSetDenseJacFn(self.mem, ida_jac)
        #flag = IDASetUserData(self.mem, <void*> user_data)
        flag = IDASetUserData(self.mem, self.uDataT)
        flag = IDASetId(self.mem, arr2nv(self.algvar))
        flag = IDASetSuppressAlg(self.mem, self.suppress_alg)
        
        #Are there sensitivities to be calculated?
        if self.nbr_params > 0:
            self.set_sensitivity_options(t0,user_data,u,ud)
        
    
    def interpolate(self, t, k):
        """
        Calls the internal IDAGetDky for the interpolated values at time t.
        t must be within the last internal step. k is the derivative of y which
        can be from zero to the current order.
        """
        cdef flag
        cdef N_Vector temp=N_VNew_Serial(self.dim)
        
        flag = IDAGetDky(self.mem, t, k, temp)
        
        if flag < 0:
            raise IDAError(flag, t)
        
        return nv2arr(temp)
    
    def set_sensitivity_options(self,t0,user_data,y,yd):
        """
        Sets the sensitivity information.
        """
        #Create the initial matrices
        self.ySO  = N_VCloneVectorArray_Serial(self.nbr_params, arr2nv(y))
        self.ydSO = N_VCloneVectorArray_Serial(self.nbr_params, arr2nv(yd))
        cdef IDASensResFn empty_p
        cdef realtype ZERO = 0.0
        cdef realtype *pbar
        cdef int *plist
        
        #Filling the start vectors
        for i in range(self.nbr_params):
             N_VConst_Serial(ZERO,  self.ySO[i]);
             N_VConst_Serial(ZERO, self.ydSO[i]);

        if self.sens_activated:
            flag = IDASensReInit(self.mem, self.ism, self.ySO, self.ydSO)
            
            if flag<0:
                raise IDAError(flag, t0)
        else:
            flag = IDASensInit(self.mem, self.nbr_params, self.ism, empty_p, self.ySO, self.ydSO)
            
            if flag<0:
                raise IDAError(flag, t0)
            
            self.sens_activated = True
        
        #Sets the parameters to the userdata object.
        self.uDataT.params = <realtype*> malloc(self.nbr_params*sizeof(realtype))
        for i in range(self.nbr_params):
            self.uDataT.params[i] = self.p[i]
        
        #Sets the pbar.
        pbar = <realtype*> malloc(self.nbr_params*sizeof(realtype))
        if self.pbar != None:
            for i in range(self.nbr_params):
                pbar[i] = self.pbar[i]
        else:
            for i in range(self.nbr_params):
                pbar[i] = 1.0
        
        #Specify problem parameter information for sensitivity calculations
        flag = IDASetSensParams(self.mem, self.uDataT.params, pbar, plist)
        
        if self.pbar != None:
            free(pbar) #Free the allocated space.
        
        if flag<0:
            raise IDAError(flag, t0)
        
        #Specify the difference quotient strategy
        flag = IDASetSensDQMethod(self.mem, self.DQtype, self.DQrhomax)
        
        if flag<0:
            raise IDAError(flag, t0)
        
        #Specify the error control strategy
        flag = IDASetSensErrCon(self.mem, self.errconS)
        
        if flag<0:
            raise IDAError(flag, t0)
        
        #Specify the maximum number of nonlinear solver iterations
        flag = IDASetSensMaxNonlinIters(self.mem, self.maxcorS)
        
        if flag<0:
            raise IDAError(flag, t0)
        
        #Estimate the sensitivity  ----SHOULD BE IMPROVED with IDASensSVTolerances ...
        flag = IDASensEEtolerances(self.mem)
        
        if flag<0:
            raise IDAError(flag, t0)
        
        #Should the sensitivities be calculated this time around?
        if self.sensToggleOff:
            flag = IDASensToggleOff(self.mem)
            
            if flag<0:
                raise IDAError(flag, t0)
    
    def interpolate_sensitivity(self,realtype t, int k, int i=-1):
        """
        This method class the internal method IDAGetSensDky which computes the k-th derivatives
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
            
            for x in range(self.nbr_params):
                flag = IDAGetSensDky1(self.mem, t, k, x, dkyS)
                
                if flag<0:
                    raise IDAError(flag, t)
                
                matrix += [nv2arr(dkyS)]
            
            return np.array(matrix)
        else:
            flag = IDAGetSensDky1(self.mem, t, k, i, dkyS)
            
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
            flag = IDAGetNumSteps(self.mem, &nsteps) #Number of steps
            flag = IDAGetNumResEvals(self.mem, &nrevals) #Number of res evals
            flag = IDADlsGetNumJacEvals(self.mem, &njevals) #Number of jac evals
            flag = IDADlsGetNumResEvals(self.mem, &nrevalsLS) #Number of res evals due to jac evals
            flag = IDAGetNumGEvals(self.mem, &ngevals) #Number of root evals
            flag = IDAGetNumErrTestFails(self.mem, &netfails) #Number of local error test failures
            flag = IDAGetNumNonlinSolvIters(self.mem, &nniters) #Number of nonlinear iteration
            flag = IDAGetNumNonlinSolvConvFails(self.mem, &nncfails) #Number of nonlinear conv failures
            
            stats_values = [nsteps, nrevals, njevals, nrevalsLS, ngevals, netfails, nniters, nncfails]
            stats_text = ['Number of Steps                          ',
                          'Number of Function Evaluations           ',
                          'Number of Jacobian Evaluations           ',
                          'Number of F-Eval During Jac-Eval         ',
                          'Number of Root Evaluations               ',
                          'Number of Error Test Failures            ',
                          'Number of Nonlinear Iterations           ',
                          'Number of Nonlinear Convergence Failures ']
            if self.stats != None:
                for x in range(len(stats_text)):
                    self.stats[stats_text[x]] += stats_values[x]
            else:
                self.stats = {}
                for x in range(len(stats_text)):
                    self.stats[stats_text[x]] = stats_values[x]
            
            if self.nbr_params > 0:
                
                flag = IDAGetSensStats(self.mem, &nfSevals, &nfevalsS, &nSetfails, &nlinsetupsS)
                flag = IDAGetSensNonlinSolvStats(self.mem, &nSniters, &nSncfails)
                
                stats_values = [nfSevals, nfevalsS, nSetfails, nlinsetupsS,nSniters, nSncfails]
                
                if self.sens_stats != None:
                    for x in range(len(stats_values)):
                        self.sens_stats[x] += stats_values[x]
                else:
                    self.sens_stats = []
                    for x in range(len(stats_values)):
                        self.sens_stats += [stats_values[x]]
            
            
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
        if self.mem == NULL:
            raise Exception, "IDA must be initialized"
        
        flag = IDASetLineSearchOffIC(self.mem, lsoff)
        
        flag = IDACalcIC(self.mem, self.icopt, tout1)
        
        if flag == IDA_SUCCESS: #Gets the calculated values
            flag = IDAGetConsistentIC(self.mem, self.curr_state, self.curr_deriv)
        
        return [flag, nv2arr(self.curr_state), nv2arr(self.curr_deriv)]
        
    
    def treat_disc(self,flag,tret):
        cdef int* event_info_
        """
        Treats solver returns with root_found flag
        """
        if flag == IDA_ROOT_RETURN:
            # Allocate memory for the event_info_ vector and initialize to zeros
            event_info_ = <int*> malloc(self.num_state_events*sizeof(int))
            for k in range(self.num_state_events):
                event_info_[k] = 0
                # Fetch data on which root functions that became zero and store in class
            flag = IDAGetRootInfo(self.mem, event_info_)
            self.event_info =  PyArray_SimpleNew(1,&self.num_state_events ,NPY_INT)
            for k in range(self.num_state_events):
                self.event_info[k] = event_info_[k]

            # Remember to deallocate
            free(event_info_)
            self.event_time=tret
            return True
        else:
            return False
   
    def run(self,realtype t0,realtype tf,float dt):
        #cdef realtype dt             # time increment
        cdef realtype tret           # return time (not neceeserily tout)
        cdef realtype tout           # communication time
        cdef int i,itask, nt
        #cdef realtype hinused,hlast,hcur,tcur
        #cdef long int nsteps, fevals, nlinsetups, netfails
        cdef int  qlast, qcurrent
        flag = IDASetStopTime(self.mem, tf)
        sol=[]
        tret=t0
        if dt > 0.0:
            nt = int(math.ceil((tf-t0)/dt))
            for i in xrange(1, nt+1):
                tout=t0+i*dt
                flag=0
                flags=0
                flags=IDASolve(self.mem,tout,&tret, self.curr_state, self.curr_deriv,IDA_NORMAL)
                if flags<0 and flags!=IDA_TSTOP_RETURN:
                    raise IDAError(flags, tret)
                sol.append((np.array(tret),nv2arr(self.curr_state),nv2arr(self.curr_deriv)))
                flag = IDAGetLastOrder(self.mem, &qlast)
                self._count_output+=1
                self._ordersum+=qlast
                avar=float(self._ordersum)/self._count_output
                if flags == 1:
                    break
                if self.treat_disc(flags,tret):
                    break
                if i == nt:
                    flags =1
                if self.store_cont:
                    break
        else: # one step mode
            if self.detailed_info == None:
                self.detailed_info = {}
                self.detailed_info['qlast'] = []
                self.detailed_info['qcurrent'] = []
                
            while tret < tf:
                flag=0
                flags=0
                flags=IDASolve(self.mem,tf,&tret, self.curr_state,self.curr_deriv,IDA_ONE_STEP)
                if flags<0 and flags!=IDA_TSTOP_RETURN:
                    raise IDAError(flags, tret)
                sol.append((np.array(tret),nv2arr(self.curr_state),nv2arr(self.curr_deriv)))
                flag = IDAGetLastOrder(self.mem, &qlast)
                flag = IDAGetCurrentOrder(self.mem, &qcurrent)
                self.detailed_info['qlast'].append(qlast)
                self.detailed_info['qcurrent'].append(qcurrent)
                self._count_output+=1
                self._ordersum+=qlast
                avar=float(self._ordersum)/self._count_output
                if self.treat_disc(flags,tret):
                    break
                if self.store_cont:
                    break
                if self.comp_step:
                    break
            else:
                flags=1
        if flags >= 1:
            self.sim_complete = True
            self.store_statistics()
        # Free memory
        #IDAFree(&self.mem)
        #N_VDestroy_Serial(self.curr_state)
        #N_VDestroy_Serial(self.curr_deriv)
        
        return sol                
            


        
    
            
        
