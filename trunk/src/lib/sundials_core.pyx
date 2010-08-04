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
Claus Fuhrer,        Lund University        October 2009
Christian Andersson, Lund University        Februari 2010

see also Jon Olav Vik: 
http://codespeak.net/pipermail/cython-dev/2009-June/005947.html

"""
from __future__ import division
import numpy as np
import math
from numpy cimport ndarray, NPY_DOUBLE, npy_intp, NPY_INT

# ==============================================
# external definitions from numpy headers
# ==============================================
cdef extern from "numpy/arrayobject.h":
    cdef object PyArray_SimpleNew(int nd, npy_intp* dims, int typenum)
    cdef object PyArray_SimpleNewFromData(int nd, npy_intp *dims,
                                           int typenum, void *data)
    void import_array() 
    void* PyArray_GetPtr(ndarray aobj, npy_intp* ind)
    void *PyArray_DATA(ndarray aobj)
    
cdef extern from "Python.h":
    cdef void Py_INCREF( object )
# ==============================================
#  C headers
# ==============================================
cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)
cdef extern from "stdlib.h":
    void *malloc(int size)
    void free(void *ptr)
    
# ==============================================
#  external definitions from Sundial headers
# ==============================================

cdef extern from "sundials/sundials_types.h":
    ctypedef double realtype
    ctypedef bint booleantype # should be bool instead of bint, but there is a bug in Cython
    # This bug is fixed in http://trac.cython.org/cython_trac/ticket/227

cdef extern from "sundials/sundials_nvector.h":
    cdef struct _generic_N_Vector:
        void* content
    ctypedef _generic_N_Vector* N_Vector
    N_Vector N_VNew_Serial(long int vec_length)
    void N_VDestroy_Serial(N_Vector v)
    void N_VPrint_Serial(N_Vector v)
    
    
cdef extern from "nvector/nvector_serial.h":

    cdef struct _N_VectorContent_Serial:
        long int length
        realtype* data
    ctypedef _N_VectorContent_Serial* N_VectorContent_Serial
    cdef N_Vector N_VMake_Serial(long int vec_length, realtype *v_data)

cdef extern from "sundials/sundials_direct.h":
    cdef struct _DlsMat:
        int type
        int M
        int N
        int ldim
        int mu
        int ml
        int s_mu
        realtype *data
        int ldata
        realtype **cols
    ctypedef _DlsMat* DlsMat
    cdef realtype* DENSE_COL(DlsMat A, int j)

cdef extern from "cvode/cvode.h":
    void* CVodeCreate(int lmm, int iter)
    ctypedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *f_data)
    int CVodeInit(void* cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)
    int CVodeReInit(void* cvode_mem, realtype t0, N_Vector y0)
    int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol)
    int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol)
    int CVodeSetStopTime(void* cvode_mem, realtype tstop)
    int CVodeGetIntegratorStats(void* cvode_mem, long int *nsteps, long int *nfevals,
                                long int *nlinsetups, long int *netfails, int *qlast, int *qcur,
                                realtype *hinused, realtype *hlast, realtype *hcur, realtype *tcur)
    int CVodeSetMaxOrd(void * cvode_mem, int maxord)
    int CVodeSetMaxNumSteps(void * cvode_mem, long int mxsteps)
    int CVodeSetMaxStep(void* cvode_mem, realtype hmax)
    int CVodeSetInitStep(void * cvode_mem, realtype hin)
    void CVodeFree(void **cvode_mem)
    int CVode(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, 
        int itask)
    int CVodeSetUserData(void *cvode_mem,void *user_data)
    int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky)
    # functions for discontinuity handling
    ctypedef int (*CVRootFn)(realtype tt, N_Vector yy, realtype *gout, void *user_data)
    int CVodeRootDirection(void *cvode_mem, int *rootdir)
    int CVodeSetNoInactiveRootWarn(void *cvode_mem)
    int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g)
    int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
    #Functions for retrieving statistics
    int CVodeGetLastOrder(void * cvode_mem,int *qlast)
    int CVodeGetCurrentOrder(void * cvode_mem,int *qcurrent)
    int CVodeGetNumSteps(void *cvode_mem, long int *nsteps) #Number of steps
    int CVodeGetNumRhsEvals(void *cvode_mem, long int *nrevals) #Number of function evals
    int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals) #Number of jac evals
    int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nrevalsLS) #Number of res evals due to jac evals
    int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals) #Number of root evals
    int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails) #Number of local error test failures
    int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters) #Number of nonlinear iteration
    int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails) #Number of nonlinear conv failures
cdef extern from "cvode/cvode_dense.h":
    int CVDense(void *cvode_mem, long int N)
    ctypedef int (*CVDlsDenseJacFn)(int N, realtype t, N_Vector y, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn djac)
    
cdef extern from "idas/idas.h":
    ctypedef int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
    void* IDACreate()
    int IDAInit(void* ida_mem, IDAResFn res, realtype t0, N_Vector y0, N_Vector yp0)
    int IDAReInit(void* ida_mem, realtype t0, N_Vector y0, N_Vector yp0)
    int IDASetStopTime(void* ida_mem, realtype tstop)
    int IDASetMaxNumSteps(void * cvode_mem, long int mxsteps)
    int IDASetMaxOrd(void * cvode_mem, int maxord)
    int IDASetMaxStep(void* ida_mem, realtype hmax)
    void IDAFree(void **cvode_mem)
    int IDAGetIntegratorStats(void* ida_mem,long int  *nsteps, long int *nrevals, 
                            long int *nlinsetups, long int *netfails, int *klast, 
                            int *kcur, realtype *hinused, realtype *hlast, 
                            realtype *hcur, realtype *tcur)
    int IDASolve(void* ida_mem, realtype tout,realtype  *tret, N_Vector yret, 
                            N_Vector ypret, int itask)
    int IDASetUserData(void *ida_mem,void *user_data)
    int IDASetInitStep(void *ida_mem, realtype hin)
    int IDAGetDky(void *ida_mem, realtype t, int k, N_Vector dky)
    # functions to control the error test
    int IDASStolerances(void *ida_mem, realtype reltol, realtype abstol)
    int IDASVtolerances(void *ida_mem, realtype reltol, N_Vector abstol)
    int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg)
    int IDASetId(void *ida_mem, N_Vector id)
    # functions for discontinuity handling
    ctypedef int (*IDARootFn)(realtype tt, N_Vector yy, N_Vector yp, realtype *gout, void *user_data)
    int IDASetRootDirection(void *ida_mem, int *rootdir)
    int IDASetNoInactiveRootWarn(void *ida_mem)
    int IDARootInit(void *ida_mem, int nrtfn, IDARootFn g)
    int IDAGetRootInfo(void *ida_mem, int *rootsfound)
    int IDACalcIC(void *ida_men, int icopt, realtype tout1)
    int IDAGetConsistentIC(void *ida_mem, N_Vector y0, N_Vector yp0)
    int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff)
    #Functions for retrieving statistics
    int IDAGetLastOrder(void *ida_mem,int *qlast) #Last order used
    int IDAGetCurrentOrder(void *ida_mem,int *qcurrent) #Order that is about to be tried
    int IDAGetNumSteps(void *ida_mem, long int *nsteps) #Number of steps
    int IDAGetNumResEvals(void *ida_mem, long int *nrevals) #Number of res evals
    int IDADlsGetNumJacEvals(void *ida_mem, long int *njevals) #Number of jac evals
    int IDADlsGetNumResEvals(void *ida_mem, long int *nrevalsLS) #Number of res evals due to jac evals
    int IDAGetNumGEvals(void *ida_mem, long int *ngevals) #Number of root evals
    int IDAGetNumErrTestFails(void *ida_mem, long int *netfails) #Number of local error test failures
    int IDAGetNumNonlinSolvIters(void *ida_mem, long int *nniters) #Number of nonlinear iteration
    int IDAGetNumNonlinSolvConvFails(void *ida_mem, long int *nncfails) #Number of nonlinear conv failures
    
    #Start Sensitivities
    #---------------------------
    ctypedef int (*IDASensResFn)(int Ns, realtype t, N_Vector yy, N_Vector yp, N_Vector *yS, N_Vector *ypS, N_Vector *resvalS, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int IDASensInit(void *ida_mem, int Ns, int ism, IDASensResFn resS, N_Vector *ySO, N_Vector *ypSO)
    int IDASensReInit(void *ida_mem, int ism, N_Vector *ySO, N_Vector *ypSO)
    
    #Options
    int IDASensToggleOff(void *ida_mem)
    int IDASensSStolerances(void *ida_mem, realtype reltolS, realtype *abstolS)
    int IDASensSVtolerances(void *ida_mem, realtype reltolS, N_Vector *abstolS)
    int IDAEEtolerances(void *ida_mem)
    
    #Results
    int IDAGetSens(void *ida_mem, realtype tret, N_Vector *yS)
    int IDAGetSensDky(void *ida_mem, realtype t, int k, N_Vector *dkyS)
    int IDAGetSensDky1(void *ida_mem, realtype t, int k, int i, N_Vector dkyS)
    
    #Options (optional)
    int IDASetSensParams(void *ida_mem, realtype *p, realtype *pbar, int *plist)
    int IDASetSensDQMethod(void *ida_mem, int DQtype, realtype DQrhomax)
    int IDASetSensErrCon(void *ida_mem, booleantype errconS)
    int IDASetSensMaxNonlinIters(void *ida_mem, int maxcorS)
    
    #Statistics
    int IDAGetSensNumResEvals(void *ida_mem, long int nfSevals)
    int IDAGetNumResEvalsSens(void *ida_mem, long int nfevalsS)
    int IDAGetSensNumErrTestFails(void *ida_mem, long int nSetfails)
    int IDAGetSensNumLinSolvSetups(void *ida_mem, long int nlinsetupsS)
    int IDAGetSensStats(void *ida_mem, long int nfSevals, long int nfevalsS, long int nSetfails, long int nlinsetupsS)
    int IDAGetSensNumNonlinSolvIters(void *ida_mem, long int nSniters)
    int IDAGetSeonsNumNonlinSolvConvFails(void *ida_mem, long int nSncfails)
    int IDAGetSensNonlinSolvStats(void *ida_mem, long int nSniters, long int nSncfails)
    
    #End Sensitivities
    #--------------------------

cdef extern from "idas/idas_dense.h":
    int IDADense(void *ida_mem, long int N)
    ctypedef int (*IDADlsDenseJacFn)(int Neq, realtype tt, realtype cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int IDADlsSetDenseJacFn(void *ida_mem, IDADlsDenseJacFn djac)

#===============================================================
# Constants
#===============================================================
#a) CVODE in
DEF CV_RHS_IND        = 0   # Index to user data rhs handling
DEF CV_RHSF_IND       = 0   # Index to user data rhs
DEF CV_JAC_IND        = 1   # Index to user data jacobian
DEF CV_ROOT_IND       = 1   # Index to user data root handling
DEF CV_ROOTF_IND      = 0   # Index to user data root function
DEF CV_SW_IND         = 1   # Index to user data root switches
DEF CV_ADAMS = 1
DEF CV_BDF   = 2
DEF CV_FUNCTIONAL = 1
DEF CV_NEWTON     = 2
DEF CV_SS = 1
DEF CV_SV = 2
DEF CV_WF = 3
DEF CV_NORMAL         = 1
DEF CV_ONE_STEP       = 2
DEF CV_NORMAL_TSTOP   = 3
DEF CV_ONE_STEP_TSTOP = 4
#b) CVODE out
DEF CV_SUCCESS = 0
DEF CV_TSTOP_RETURN = 1
DEF CV_ROOT_RETURN = 2
DEF CV_ROOT_RETURN    = 2   # CVSolve succeeded and found one or more roots.
DEF CV_RTFUNC_FAIL    = -10 # The rootfinding function failed in an unrecoverable manner.

#c) IDA in
DEF IDA_NORMAL         = 1   # Solver returns at specified output time.
DEF IDA_ONE_STEP       = 2   # Solver returns after each successful step.
DEF IDA_RES_IND        = 0   # Index to user data residual handling
DEF IDA_RESF_IND       = 0   # Index to user data residual
DEF IDA_JAC_IND        = 1   # Index to user data jacobian
DEF IDA_ROOT_IND       = 1   # Index to user data root handling
DEF IDA_ROOTF_IND      = 0   # Index to user data root function
DEF IDA_SW_IND         = 1   # Index to user data root switches
DEF IDA_YA_YDP_INIT    = 1   # See IDA Documentation 4.5.4
DEF IDA_Y_INIT         = 2   # See IDA Documentation 4.5.4
DEF IDA_SIMULTANEOUS   = 1   # Simultaneous corrector forward sensitivity method.
DEF IDA_STAGGERED      = 2   # Staggered corrector forward sensitivity method.
DEF IDA_CENTERED       = 1   # Central difference quotient approximation (2nd order) of the sensitivity RHS.
DEF IDA_FORWARD        = 2   # Forward difference quotient approximation (1st order) of the sensitivity RHS.
#d) IDA out
DEF IDA_SUCCESS        = 0   # Successful function return.   
DEF IDA_TSTOP_RETURN   = 1   # IDASolve succeeded by reaching the specified stopping point.
DEF IDA_ROOT_RETURN    = 2   # IDASolve succeeded and found one or more roots.
DEF IDA_RTFUNC_FAIL    = -10 # The rootfinding function failed in an unrecoverable manner.

# ===============================================================
#  Module level functions
# ===============================================================
cdef N_Vector arr2nv(x):
    """Create new N_Vector from numpy array"""
    x=np.array(x)
    cdef long int n = len(x)

    cdef ndarray[double, ndim=1,mode='c'] ndx=x
    import_array()
    cdef void* data_ptr=PyArray_DATA(ndx)
    cdef N_Vector v=N_VNew_Serial(n)
    memcpy((<N_VectorContent_Serial>v.content).data, data_ptr, n*sizeof(double))
    return v
    
cdef nv2arr(N_Vector v):
    """Create new numpy array from N_Vector"""
    cdef long int n = (<N_VectorContent_Serial>v.content).length
    cdef realtype* v_data = (<N_VectorContent_Serial>v.content).data
    cdef long int i
    cdef npy_intp dims = <npy_intp>n
    import_array()
    cdef ndarray x=np.empty(n)
    #cdef ndarray x = PyArray_SimpleNewFromData(1, &dims, NPY_DOUBLE,v_data)
    memcpy(x.data, v_data, n*sizeof(double))
#    x = np.empty(n)
#    for i in range(n):
#        x[i] = v_data[i]
    return x

cdef int sundials_error(int flag, int solver, float time) except -1:
    """
    Raises an Sundials error and prints suitable error message.
    """
    if solver == 1: #IDA
        err_mess = {-1: 'The solver took max internal steps but could not reach tout.',
                    -2: 'The solver could not satisfy the accuracy demanded by the user for some internal step.',
                    -3: 'Error test failures occurred too many times during one internal time step or minimum step size was reached.',
                    -4: 'Convergence test failures occurred too many times during one internal time step or minimum step size was reached.',
                    -5: 'The linear solvers initialization function failed.',
                    -6: 'The linear solvers setup function failed in an unrecoverable manner.',
                    -7: 'The linear solvers solve function failed in an unrecoverable manner.',
                    -8: 'The user-provided residual function failed in an unrecoverable manner.',
                    -9: 'The user-provided residual function repeatedly returned a recoverable error flag, but the solver was unable to recover.',
                    -10: 'The rootfinding function failed in an unrecoverable manner.',
                    -11: 'The inequality constraints were violated and the solver was unable to recover.',
                    -12: 'The user-provided residual function failed recoverable on the first call.',
                    -13: 'The line search failed.',
                    -14: 'The residual function, linear solver setup function or linear solver solve function had a recoverable failure. But IDACalcIC could not recover.',
                    -20: 'The ida_mem argument was NULL.',
                    -21: 'A memory allocation failed.',
                    -22: 'One of the function inputs is illegal.',
                    -23: 'The IDA memory was not allocated by a call to IDAInit.',
                    -24: 'Zero value of some error weight component.',
                    -25: 'The k-th derivative is not available.',
                    -26: 'The time t is outside the last step taken.',
                    -27: 'The vector argument where derivative should be stored is NULL.',
                    -41: 'The user-provided sensitivity residual function failed in an unrecoverable manner.',
                    -42: 'The user-provided sensitivity residual function repeatedly returned a recoverable error flag, but the solver was unable to recover.',
                    -43: 'The sensitivity identifier is not valid.'}
    else:           #CVode
        err_mess = {-1: 'The solver took max internal steps but could not reach tout.',
                    -2: 'The solver could not satisfy the accuracy demanded by the user for some internal step.',
                    -3: 'Error test failures occurred too many times during one internal time step or minimum step size was reached.',
                    -4: 'Convergence test failures occurred too many times during one internal time step or minimum step size was reached.',
                    -5: 'The linear solvers initialization function failed.',
                    -6: 'The linear solvers setup function failed in an unrecoverable manner.',
                    -7: 'The linear solvers solve function failed in an unrecoverable manner.',
                    -8: 'The user-provided rhs function failed in an unrecoverable manner.',
                    -9: 'The right-hand side function failed at the first call.',
                    -10: 'The right-hand side function had repetead recoverable errors.',
                    -11: 'The right-hand side function had a recoverable error, but no recovery is possible.',
                    -12: 'The rootfinding function failed in an unrecoverable manner.',
                    -20: 'A memory allocation failed.',
                    -21: 'The cvode_mem argument was NULL.',
                    -22: 'One of the function inputs is illegal.',
                    -23: 'The CVode memory block was not allocated by a call to CVodeMalloc.',
                    -24: 'The derivative order k is larger than the order used.',
                    -25: 'The time t is outside the last step taken.',
                    -26: 'The output derivative vector is NULL.',
                    -27: 'The output and initial times are too close to each other.'}
    try:    
        raise SundialsError, err_mess[flag]+' At time %f'%time
    except KeyError:
        raise SundialsError, 'Sundials failed with flag %s. At time %f'%(flag,time)

cdef int cv_rhs(realtype t, N_Vector yv, N_Vector yvdot, void* user_data):
    """
    Wraps  Python rhs-callback function to obtain CVode required interface
    see also ctypedef statement above
    """
    y=nv2arr(yv)
    cdef realtype* ydotptr=(<N_VectorContent_Serial>yvdot.content).data
    cdef long int n=(<N_VectorContent_Serial>yv.content).length
    try:
        switch=(<object> user_data)[CV_ROOT_IND][CV_SW_IND]
    except:
        switch=False
        pass
    try:
        if switch:
            ydot=(<object> user_data)[CV_RHS_IND][CV_RHSF_IND](t,y,switch)  # call to the python rhs function
        else:
            ydot=(<object> user_data)[CV_RHS_IND][CV_RHSF_IND](t,y)
        for i in range(n):
            ydotptr[i]=ydot[i]
        return 0
    except:
        return 1 # recoverable error (see Sundials description)

cdef int cv_jac(int Neq, realtype t, N_Vector yv, N_Vector fy, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    Wraps Python jacobian-callback function to obtain CV required interface.
    """
    y = nv2arr(yv)
    try:
        switch=(<object> user_data)[CV_ROOT_IND][CV_SW_IND]
    except:
        switch=False
    cdef realtype* col_i=DENSE_COL(Jac,0)
    try:
        if switch:
            jacobian=(<object> user_data)[CV_RHS_IND][CV_JAC_IND](t,y,switch)  # call to the python residual function
        else:
            jacobian=(<object> user_data)[CV_RHS_IND][CV_JAC_IND](t,y)
        
        for i in range(Neq):
            col_i = DENSE_COL(Jac, i)
            for j in range(Neq):
                col_i[j] = jacobian[j,i]
        return 0
    except: #None recoverable
        return -1

cdef int cv_root(realtype t, N_Vector yv, realtype *gout,  void* user_data):
    """
    Wraps  Python root-callback function to obtain CV required interface
    see also ctypedef statement above
    """
    y=nv2arr(yv)
    try:
        switch=(<object> user_data)[CV_ROOT_IND][CV_SW_IND]
    except:
        switch=False
    try:
        rootf=(<object> user_data)[CV_ROOT_IND][CV_ROOTF_IND](t,y,switch)  # call to the python root function 
        rootf = np.asarray(rootf).reshape(-1) # Make sure we get a vector
        for i in range(rootf.shape[0]):
            gout[i]=rootf[i]
        return 0
    except:
        return 1 # generates an error of type IDA_RTFUNC_FAIL
cdef int ida_res(realtype t, N_Vector yv, N_Vector yvdot, N_Vector residual, void* user_data):
    """
    Wraps  Python res-callback function to obtain IDA required interface
    see also ctypedef statement above
    """
    y=nv2arr(yv)
    yd=nv2arr(yvdot)
    try:
        switch=(<object> user_data)[IDA_ROOT_IND][IDA_SW_IND]
    except:
        switch=False
    cdef realtype* resptr=(<N_VectorContent_Serial>residual.content).data
    cdef long int n=(<N_VectorContent_Serial>yv.content).length
    try:
        if switch:
            res=(<object> user_data)[IDA_RES_IND][IDA_RESF_IND](t,y,yd,switch)  # call to the python residual function
        else:
            res=(<object> user_data)[IDA_RES_IND][IDA_RESF_IND](t,y,yd)
        for i in range(n):
            resptr[i]=res[i]
        return 0
    except:
        return 1 # recoverable error (see Sundials description)
        
cdef int ida_jac(int Neq, realtype t, realtype c, N_Vector yv, N_Vector yvdot, N_Vector residual, DlsMat Jac,
                 void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3):
    """
    Wraps Python jacobian-callback function to obtain IDA required interface.
    """
    y = nv2arr(yv)
    yd = nv2arr(yvdot)
    try:
        switch=(<object> user_data)[IDA_ROOT_IND][IDA_SW_IND]
    except:
        switch=False
    cdef realtype* col_i=DENSE_COL(Jac,0)
    try:
        if switch:
            jacobian=(<object> user_data)[IDA_RES_IND][IDA_JAC_IND](c,t,y,yd,switch)  # call to the python residual function
        else:
            jacobian=(<object> user_data)[IDA_RES_IND][IDA_JAC_IND](c,t,y,yd)
        
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
        switch=(<object> user_data)[IDA_ROOT_IND][IDA_SW_IND]
    except:
        switch=False
    try:
        rootf=(<object> user_data)[IDA_ROOT_IND][IDA_ROOTF_IND](t,y,yd,switch)  # call to the python root function 
        rootf = np.asarray(rootf).reshape(-1) # Make sure we get a vector
        for i in range(rootf.shape[0]):
            gout[i]=rootf[i]
        return 0
    except:
        return 1 # generates an error of type IDA_RTFUNC_FAIL
        
#cdef int completed_step(void* user_data):
#    return (<object> user_data)[0]()

class SundialsError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def eval_rhs(t,y,func):
    """
    Evaluates specific CVode rhs function for test purposes only
    """
    cdef N_Vector yv=arr2nv(y)
    cdef N_Vector yvdot=arr2nv(np.zeros(len(y)))
    cv_rhs(t,yv,yvdot,<void *>func)
    return nv2arr(yvdot)

def eval_res(t,y,yd,func):
    """
    Evaluates specific IDA res function for test purposes only
    """
    cdef N_Vector yv=arr2nv(y)
    cdef N_Vector yvd=arr2nv(yd)
    cdef N_Vector resid=arr2nv(np.zeros(len(y)))
    flag=ida_res(t,yv,yvd,resid,<void *>func)
    return nv2arr(resid),flag
        
# =====================================================================
#  Wrapper Class definition
# =====================================================================
cdef class CVode_wrap:
    """Class to wrap CVode"""
    cdef:
        void* mem
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
        #void* comp_step_method
        N_Vector curr_state
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
    def cvinit(self,t0,user_data,u,maxord, max_steps, init_step):
        cdef flag
        self.curr_state=arr2nv(u)
        self.store_state = True
        self.max_steps = max_steps
        self._ordersum=self._count_output=0 # initialize ordersum and output count for avarage order
        if self.mem == NULL:
            # Newinitialization
            self.mem=CVodeCreate(self.discr, self.iter)
            if self.mem == NULL:
                raise Exception, 'CVodeCreate: Memory allocation failed'
            flag=CVodeInit(self.mem, cv_rhs, t0,self.curr_state)
            if flag!=CV_SUCCESS:
                raise Exception,"CVode Initialization Error"
            if self.num_state_events>0: 
                flag = CVodeRootInit(self.mem, self.num_state_events, cv_root)
                if flag!=CV_SUCCESS:
                    raise Exception,"CV root-finding initialization error"
        else:
            flag = CVodeReInit(self.mem, t0, self.curr_state)
        if self.abstol_ar[0] > 0:
            flag = CVodeSVtolerances(self.mem, self.reltol, arr2nv(self.abstol_ar))
        else:
            flag = CVodeSStolerances(self.mem, self.reltol, self.abstol)
        if flag!=CV_SUCCESS:
                raise Exception,"CVode Tolerance Initialization  Error"
        if maxord:
            flag=CVodeSetMaxOrd(self.mem, maxord)
        flag = CVodeSetMaxNumSteps(self.mem, self.max_steps)
        flag = CVodeSetInitStep(self.mem, init_step)
        flag=CVDense(self.mem, self.dim)
        if self.jacobian:
            flag = CVDlsSetDenseJacFn(self.mem, cv_jac)
        flag = CVodeSetUserData(self.mem, <void*>user_data)
        try:
            flag= CVodeSetMaxStep(self.mem, self.max_h)
        except AttributeError:
            pass
    
    #def set_completed_method(self,data):
    #    self.comp_step_method = <void*>data
    
    def interpolate(self, t, k):
        """
        Calls the internal CVodeGetDky for the interpolated values at time t.
        t must be within the last internal step. k is the derivative of y which
        can be from zero to the current order.
        """
        cdef N_Vector temp=N_VNew_Serial(self.dim)
        
        flag = CVodeGetDky(self.mem, t, k, temp)
        
        if flag < 0:
            sundials_error(flag,2,t)
        
        return nv2arr(temp)
    
    def store_statistics(self):
        """
        Retrieves and stores the statistics.
        """
        cdef long int nsteps, nrevals, njevals, nrevalsLS, ngevals, netfails, nniters, nncfails
        
        if self.store_state:
            flag = CVodeGetNumSteps(self.mem, &nsteps) #Number of steps
            flag = CVodeGetNumRhsEvals(self.mem, &nrevals) #Number of function evals
            if self.iter == 1:
                njevals = 0
                nrevalsLS = 0
            else:
                flag = CVDlsGetNumJacEvals(self.mem, &njevals) #Number of jac evals
                flag = CVDlsGetNumRhsEvals(self.mem, &nrevalsLS) #Number of res evals due to jac evals            
            flag = CVodeGetNumGEvals(self.mem, &ngevals) #Number of root evals
            flag = CVodeGetNumErrTestFails(self.mem, &netfails) #Number of local error test failures
            flag = CVodeGetNumNonlinSolvIters(self.mem, &nniters) #Number of nonlinear iteration
            flag = CVodeGetNumNonlinSolvConvFails(self.mem, &nncfails) #Number of nonlinear conv failures
            
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
            flag = CVodeGetRootInfo(self.mem, event_info_)
            self.event_info =  PyArray_SimpleNew(1,&self.num_state_events ,NPY_INT)
            for k in range(self.num_state_events):
                self.event_info[k] = event_info_[k]
            # Remember to deallocate
            free(event_info_)
            self.event_time=tret
            return True
        else:
            return False

        
    def run(self,float t0,float tf,float dt):
        #cdef realtype dt             # time increment
        cdef realtype tret           # return time (not neceeserily tout)
        cdef realtype tout           # communication time
        cdef int i,itask,nt
        cdef realtype hinused,hlast,hcur,tcur
        #cdef long int nsteps, fevals, nlinsetups, netfails
        cdef int  qlast, qcurrent
        flag = CVodeSetStopTime(self.mem, tf)
        sol=[]
        tret=t0
        if dt > 0.0:
            nt = int(math.ceil((tf-t0)/dt))
            for i in xrange(1, nt+1):
                tout=t0+i*dt
                flag=0
                flags=CVode(self.mem,tout,self.curr_state,&tret,CV_NORMAL)
                if flags<0 and flags!=CV_TSTOP_RETURN:
                    sundials_error(flags,2,tret)
                sol.append((np.array(tret),nv2arr(self.curr_state)))
                flag = CVodeGetLastOrder(self.mem, &qlast)
                self._count_output+=1
                self._ordersum+=qlast
                avar=float(self._ordersum)/self._count_output
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
                flags=CVode(self.mem,tf,self.curr_state,&tret,CV_ONE_STEP)
                if flags<0 and flags!=CV_TSTOP_RETURN:
                    sundials_error(flags,2,tret)
                sol.append((np.array(tret),nv2arr(self.curr_state)))
                flag = CVodeGetLastOrder(self.mem, &qlast)
                flag = CVodeGetCurrentOrder(self.mem, &qcurrent)
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
        # Free memory
        #CVodeFree(&self.mem)
        #N_VDestroy_Serial(self.curr_state)
        return sol
        
        
cdef class IDA_wrap:
    """Class to wrap Sundials IDA"""
    cdef:
        void* mem
        #void* comp_step_method
        public int dim, maxord, _ordersum,_count_output, max_h
        public long int max_steps
        public realtype abstol,reltol,event_time
        public realtype t0
        public ndarray abstol_ar,algvar,event_info
        public dict stats
        public dict detailed_info
        public booleantype suppress_alg,jacobian, store_cont,store_state,comp_step
        public int icopt, nbr_params
        public npy_intp num_state_events
        public booleantype sim_complete
        N_Vector curr_state
        N_Vector curr_deriv
        N_Vector temp_nvector
    def __init__(self,dim):
        self.dim=dim
        self.store_cont = False
        self.store_state = False
        self.comp_step = False
        self.sim_complete = False
        self.nbr_params = 0
    def idinit(self,t0,user_data,u,ud,maxord, max_steps, init_step, max_h):
        cdef flag
        self.t0 = t0
        self.store_state = True
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
        flag = IDASetUserData(self.mem, <void*> user_data)
        flag = IDASetId(self.mem, arr2nv(self.algvar))
        flag = IDASetSuppressAlg(self.mem, self.suppress_alg)
        
    #def set_completed_method(self,data):
    #    self.comp_step_method = <void*>data
    
    def interpolate(self, t, k):
        """
        Calls the internal IDAGetDky for the interpolated values at time t.
        t must be within the last internal step. k is the derivative of y which
        can be from zero to the current order.
        """
        cdef N_Vector temp=N_VNew_Serial(self.dim)
        
        flag = IDAGetDky(self.mem, t, k, temp)
        
        if flag < 0:
            sundials_error(flag,1,t)
        
        return nv2arr(temp)
        
    def get_sens_res(self,realtype t, int k, int i=-1):
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
        
        if i==-1:
            
            matrix = []
            
            for x in range(self.nbr_params):
                flag = IDAGetSensDky1(self.mem, t, k, x, dkyS)
                
                if flag<0:
                    sundials_error(flag,1,t)
                
                matrix += nv2arr(dkyS)
            
            return np.array(matrix)
        else:
            flag = IDAGetSensDky1(self.mem, t, k, i, dkyS)
            
            if flag <0:
                sundials_error(flag,1,t)
            
            return nv2arr(dkyS)

        
    def store_statistics(self):
        """
        Retrieves and stores the statistics.
        """
        cdef long int nsteps, nrevals,njevals,nrevalsLS,ngevals,netfails,nniters,nncfails

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
   
    def run(self,float t0,float tf,float dt):
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
                    sundials_error(flags,1,tret)
                sol.append((np.array(tret),nv2arr(self.curr_state),nv2arr(self.curr_deriv)))
                flag = IDAGetLastOrder(self.mem, &qlast)
                self._count_output+=1
                self._ordersum+=qlast
                avar=float(self._ordersum)/self._count_output
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
                    sundials_error(flags,1,tret)
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
            


        
    
            
        
