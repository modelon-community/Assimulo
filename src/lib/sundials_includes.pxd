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

"""
Cython Wrapper for interfacing Python with CVode and IDA (Sundials Version 2.4.0)
Claus Fuhrer,        Lund University        
Christian Andersson, Lund University        

see also Jon Olav Vik: 
http://codespeak.net/pipermail/cython-dev/2009-June/005947.html

"""
#import numpy as N
#cimport numpy as N

from numpy cimport NPY_DOUBLE, npy_intp, NPY_INT

#==============================================
#External definitions from Sundials headers
#==============================================

IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "sundials/sundials_context.h":
        ctypedef _SUNContext * SUNContext
        cdef struct _SUNContext:
            pass
        int SUNContext_Create(void* comm, SUNContext* ctx)

IF SUNDIALS_VERSION >= (7,0,0):
    cdef extern from "sundials/sundials_context.h":
        ctypedef int SUNErrCode
        ctypedef void (*SUNErrHandlerFn)(int line, const char* func, const char* file, const char* msg, SUNErrCode err_code, void* err_user_data, SUNContext sunctx)
        SUNErrCode SUNContext_PushErrHandler(SUNContext sunctx, SUNErrHandlerFn err_fn, void* err_user_data)

IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "sundials/sundials_types.h":
        ctypedef double sunrealtype
        ctypedef bint sunbooleantype
    ctypedef double realtype
    ctypedef bint booleantype
ELSE:
    cdef extern from "sundials/sundials_types.h":
        ctypedef double realtype
        ctypedef bint booleantype # should be bool instead of bint, but there is a bug in Cython

cdef extern from "sundials/sundials_nvector.h":
    ctypedef _generic_N_Vector* N_Vector
    
    cdef struct _generic_N_Vector_Ops:
        realtype    (*nvwrmsnorm)(N_Vector, N_Vector)
        realtype    (*nvwl2norm)(N_Vector, N_Vector)
    ctypedef _generic_N_Vector_Ops *N_Vector_Ops
    
    cdef struct _generic_N_Vector:
        void* content
        N_Vector_Ops ops

cdef extern from "nvector/nvector_serial.h":
    cdef struct _N_VectorContent_Serial:
        long int length
        booleantype own_data
        realtype* data
    ctypedef _N_VectorContent_Serial* N_VectorContent_Serial
    cdef N_Vector N_VMake_Serial(long int vec_length, realtype *v_data)
    void N_VSetArrayPointer_Serial(realtype *v_data, N_Vector v)
    void N_VConst_Serial(realtype c, N_Vector z)
    IF SUNDIALS_VERSION >= (6,0,0):
        N_Vector N_VNew_Serial(long int vec_length, SUNContext ctx)
        N_Vector *N_VCloneVectorArray(int count, N_Vector w)
        N_Vector *N_VCloneVectorArrayEmpty(int count, N_Vector w)
        void N_VDestroy(N_Vector v)
    ELSE:
        N_Vector N_VNew_Serial(long int vec_length)
        N_Vector *N_VCloneVectorArray_Serial(int count, N_Vector w)
        N_Vector *N_VCloneVectorArrayEmpty_Serial(int count, N_Vector w)
        void N_VDestroy_Serial(N_Vector v)
    void N_VPrint_Serial(N_Vector v)


IF SUNDIALS_VERSION >= (4,0,0):
    cdef extern from "sundials/sundials_nonlinearsolver.h":
        ctypedef _generic_SUNNonlinearSolver *SUNNonlinearSolver
        
        cdef struct _generic_SUNNonlinearSolver:
            pass
ELSE:
    #Dummy defines
    ctypedef void *SUNNonlinearSolver
        
IF SUNDIALS_VERSION >= (3,0,0):
    cdef extern from "sundials/sundials_types.h":
        IF SUNDIALS_VECTOR_SIZE == "64":
            ctypedef long int sunindextype
        ELSE:
            ctypedef int sunindextype
    cdef extern from "sundials/sundials_matrix.h":
        ctypedef _generic_SUNMatrix *SUNMatrix
        void SUNMatDestroy(SUNMatrix A)
        
        cdef struct _generic_SUNMatrix_Ops:
            SUNMatrix_ID (*getid)(SUNMatrix)
            SUNMatrix    (*clone)(SUNMatrix)
            void         (*destroy)(SUNMatrix)
            int          (*zero)(SUNMatrix)
            int          (*copy)(SUNMatrix, SUNMatrix)
            int          (*scaleadd)(realtype, SUNMatrix, SUNMatrix)
            int          (*scaleaddi)(realtype, SUNMatrix)
            int          (*matvec)(SUNMatrix, N_Vector, N_Vector)
            int          (*space)(SUNMatrix, long int*, long int*)

        cdef struct _generic_SUNMatrix:
            void *content
            _generic_SUNMatrix_Ops *ops
            
        cdef enum SUNMatrix_ID:
            SUNMATRIX_DENSE, 
            SUNMATRIX_BAND, 
            SUNMATRIX_SPARSE, 
            SUNMATRIX_CUSTOM
    
    cdef extern from "sundials/sundials_linearsolver.h":
        ctypedef _generic_SUNLinearSolver *SUNLinearSolver
        int SUNLinSolFree(SUNLinearSolver S)
        
        cdef struct _generic_SUNLinearSolver_Ops:
            SUNLinearSolver_Type (*gettype)(SUNLinearSolver)
            int                  (*setatimes)(SUNLinearSolver, void*, ATimesFn)
            int                  (*setpreconditioner)(SUNLinearSolver, void*, 
                                                    PSetupFn, PSolveFn)
            int                  (*setscalingvectors)(SUNLinearSolver,
                                                    N_Vector, N_Vector)
            int                  (*initialize)(SUNLinearSolver)
            int                  (*setup)(SUNLinearSolver, SUNMatrix)
            int                  (*solve)(SUNLinearSolver, SUNMatrix, N_Vector, 
                                        N_Vector, realtype)
            int                  (*numiters)(SUNLinearSolver)
            realtype             (*resnorm)(SUNLinearSolver)
            long int             (*lastflag)(SUNLinearSolver)
            int                  (*space)(SUNLinearSolver, long int*, long int*)
            N_Vector             (*resid)(SUNLinearSolver)
            int                  (*free)(SUNLinearSolver)
        
        cdef struct _generic_SUNLinearSolver:
            void *content
            _generic_SUNLinearSolver_Ops *ops
            
        cdef enum SUNLinearSolver_Type:
            SUNLINEARSOLVER_DIRECT,
            SUNLINEARSOLVER_ITERATIVE,
            SUNLINEARSOLVER_CUSTOM
    
    cdef extern from "sunmatrix/sunmatrix_dense.h":
        ctypedef _SUNMatrixContent_Dense *SUNMatrixContent_Dense
        cdef struct _SUNMatrixContent_Dense:
            sunindextype M
            sunindextype N
            realtype *data
            sunindextype ldata
            realtype **cols
        IF SUNDIALS_VERSION >= (6,0,0):
            SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N, SUNContext ctx)
        ELSE:
            SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N)
    cdef extern from "sunmatrix/sunmatrix_sparse.h":
        ctypedef _SUNMatrixContent_Sparse *SUNMatrixContent_Sparse
        cdef struct _SUNMatrixContent_Sparse:
            sunindextype M
            sunindextype N
            sunindextype NNZ
            sunindextype NP
            realtype *data
            int sparsetype
            sunindextype *indexvals
            sunindextype *indexptrs
            sunindextype **rowvals
            sunindextype **colptrs
            sunindextype **colvals
            sunindextype **rowptrs
        IF SUNDIALS_VERSION >= (6,0,0):
            SUNMatrix SUNSparseMatrix(sunindextype M, sunindextype N, sunindextype NNZ, int sparsetype, SUNContext ctx)
        ELSE:
            SUNMatrix SUNSparseMatrix(sunindextype M, sunindextype N, sunindextype NNZ, int sparsetype)
    cdef extern from "sunlinsol/sunlinsol_dense.h":
        IF SUNDIALS_VERSION >= (4,0,0):
            IF SUNDIALS_VERSION >= (6,0,0):
                SUNLinearSolver SUNLinSol_Dense(N_Vector y, SUNMatrix A, SUNContext ctx)
            ELSE:
                SUNLinearSolver SUNLinSol_Dense(N_Vector y, SUNMatrix A)
        ELSE:
            SUNLinearSolver SUNDenseLinearSolver(N_Vector y, SUNMatrix A)
    cdef extern from "sunlinsol/sunlinsol_spgmr.h":
        IF SUNDIALS_VERSION >= (4,0,0):
            IF SUNDIALS_VERSION >= (6,0,0):
                SUNLinearSolver SUNLinSol_SPGMR(N_Vector y, int pretype, int maxl, SUNContext ctx)
            ELSE:
                SUNLinearSolver SUNLinSol_SPGMR(N_Vector y, int pretype, int maxl)
        ELSE:
            SUNLinearSolver SUNSPGMR(N_Vector y, int pretype, int maxl)

ELSE: 
    #Dummy defines
    ctypedef void *SUNLinearSolver
    ctypedef void *SUNMatrix
    ctypedef void *SUNMatrixContent_Dense
    ctypedef void *SUNMatrixContent_Sparse
    ctypedef int sunindextype


#Struct for handling the Jacobian data
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

IF SUNDIALS_VERSION >= (5,0,0):
    pass
ELIF SUNDIALS_VERSION >= (2,6,3):
    cdef extern from "sundials/sundials_sparse.h":
        cdef struct _SlsMat:
            int M
            int N
            int NNZ
            int NP
            realtype *data
            int sparsetype
            int *indexvals
            int *indexptrs
            int **rowvals
            int **colptrs
            int **colvals
            int **rowptrs
        ctypedef _SlsMat* SlsMat
ELIF SUNDIALS_VERSION >= (2,6,0):
    cdef extern from "sundials/sundials_sparse.h":
        cdef struct _SlsMat:
            int M
            int N
            int NNZ
            realtype *data
            int *rowvals
            int *colptrs
        ctypedef _SlsMat* SlsMat
ELSE:
    cdef struct _SlsMat:
        int M
        int N
        int NNZ
        realtype *data
        int *rowvals
        int *colptrs
    ctypedef _SlsMat* SlsMat
    
#==============================================
# C headers
#==============================================
cdef extern from "string.h":
    void *memcpy(void *s1, void *s2, int n)
cdef extern from "stdlib.h":
    void *malloc(int size)
    void free(void *ptr)
    
#==============================================
#External definitions from Sundials headers
#==============================================

IF SUNDIALS_WITH_SUPERLU:
    cdef inline int with_superlu(): return 1
ELSE:
    cdef inline int with_superlu(): return 0

IF SUNDIALS_VERSION >= (4,0,0):
    cdef extern from "cvodes/cvodes.h":
        IF SUNDIALS_VERSION >= (6,0,0):
            void* CVodeCreate(int lmm, SUNContext ctx)
        ELSE:
            void* CVodeCreate(int lmm)

        int CVodeSetNonlinearSolver(void *cvode_mem, SUNNonlinearSolver NLS)
        int CVodeSetNonlinearSolverSensSim(void *cvode_mem, SUNNonlinearSolver NLS)
        int CVodeSetNonlinearSolverSensStg(void *cvode_mem, SUNNonlinearSolver NLS)
    
    cdef extern from "sunnonlinsol/sunnonlinsol_newton.h":
        IF SUNDIALS_VERSION >= (6,0,0):
            SUNNonlinearSolver SUNNonlinSol_Newton(N_Vector y, SUNContext ctx)
            SUNNonlinearSolver SUNNonlinSol_NewtonSens(int count, N_Vector y, SUNContext ctx)
        ELSE:
            SUNNonlinearSolver SUNNonlinSol_Newton(N_Vector y)
            SUNNonlinearSolver SUNNonlinSol_NewtonSens(int count, N_Vector y)
    cdef extern from "sunnonlinsol/sunnonlinsol_fixedpoint.h":
        IF SUNDIALS_VERSION >= (6,0,0):
            SUNNonlinearSolver SUNNonlinSol_FixedPoint(N_Vector y, int m, SUNContext ctx)
            SUNNonlinearSolver SUNNonlinSol_FixedPointSens(int count, N_Vector y, int m, SUNContext ctx)
        ELSE:
            SUNNonlinearSolver SUNNonlinSol_FixedPoint(N_Vector y, int m)
            SUNNonlinearSolver SUNNonlinSol_FixedPointSens(int count, N_Vector y, int m)
ELSE:
    cdef extern from "cvodes/cvodes.h":
        void* CVodeCreate(int lmm, int iter)
        
cdef extern from "cvodes/cvodes.h":
    ctypedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *f_data)
    int CVodeInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)
    int CVodeReInit(void *cvode_mem, realtype t0, N_Vector y0)
    void CVodeFree(void **cvode_mem)
    int CVode(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, int itask)
    
    #Functions for settings options
    int CVodeSetMaxOrd(void *cvode_mem, int maxord)
    int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps)
    int CVodeSetMaxStep(void   *cvode_mem, realtype hmax)
    int CVodeSetMinStep(void   *cvode_mem, realtype hmin)
    int CVodeSetInitStep(void  *cvode_mem, realtype hin)
    int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol)
    int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol)
    IF SUNDIALS_CVODE_RTOL_VEC:
        int CVodeVVtolerances(void *cvode_mem, N_Vector reltol, N_Vector abstol)
    int CVodeSetStopTime(void  *cvode_mem, realtype tstop)
    int CVodeSetUserData(void  *cvode_mem,void *user_data)
    int CVodeSetMaxConvFails(void *cvode_mem, int maxncf)
    int CVodeSetMaxErrTestFails(void *cvode_mem, int maxnef)
    int CVodeSetMaxNonlinIters(void *cvode_mem, int maxcor)
    
    #Functions for retrieving results
    int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky)

    #Functions for discontinuity handling
    ctypedef int (*CVRootFn)(realtype tt, N_Vector yy, realtype *gout, void *user_data)
    int CVodeRootDirection(void *cvode_mem, int *rootdir)
    int CVodeSetNoInactiveRootWarn(void *cvode_mem)
    int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g)
    int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
    
    #Functions for retrieving statistics
    int CVodeGetLastOrder(void * cvode_mem,int *qlast)
    int CVodeGetLastStep(void * cvode_mem, realtype *hlast)
    int CVodeGetCurrentOrder(void * cvode_mem,int *qcurrent)
    int CVodeGetActualInitStep(void * cvode_mem, realtype *hinused)
    int CVodeGetNumSteps(void *cvode_mem, long int *nsteps) #Number of steps
    int CVodeGetNumRhsEvals(void *cvode_mem, long int *nrevals) #Number of function evals
    IF SUNDIALS_VERSION >= (4,0,0):
        int CVodeGetNumJacEvals(void *cvode_mem, long int *njevals) #Number of jac evals
        int CVodeGetNumLinRhsEvals(void *cvode_mem, long int *nrevalsLS) #Number of res evals due to jac evals
    ELSE:
        int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals) #Number of jac evals
        int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nrevalsLS) #Number of res evals due to jac evals
    int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals) #Number of root evals
    int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails) #Number of local error test failures
    int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters) #Number of nonlinear iteration
    int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails) #Number of nonlinear conv failures
    int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters, long int *nncfails)
    int CVodeGetIntegratorStats(void* cvode_mem, long int *nsteps, long int *nfevals,
                                long int *nlinsetups, long int *netfails, int *qlast, int *qcur,
                                realtype *hinused, realtype *hlast, realtype *hcur, realtype *tcur)
    int CVodeGetNumStabLimOrderReds(void *cvode_mem, long int *nslred)
    
    #Sensitivity methods
    ctypedef int (*CVSensRhsFn)(int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS,
                                N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2)
    ctypedef int (*CVSensRhs1Fn)(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector *yS,
                                N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2)
    int CVodeSensInit(void *cvode_mem, int Ns, int ism, CVSensRhsFn fS, N_Vector *ySO)
    int CVodeSensInit1(void *cvode_mem, int Ns, int ism, CVSensRhs1Fn fS1, N_Vector *ySO)
    int CVodeSensReInit(void *cvode_mem, int ism, N_Vector *ySO)
    int CVodeSensFree(void *cvode_mem)
    int CVodeSensToggleOff(void *cvode_mem)
    int CVodeSensEEtolerances(void *cvode_mem)
    int CVodeGetSens(void *cvode_mem, realtype *tret, N_Vector *yS)
    int CVodeGetSensDky(void *cvode_mem, realtype t, int k, N_Vector *dkyS)
    int CVodeGetSens1(void *cvode_mem, realtype *tret, int iss, N_Vector yS)
    int CVodeGetSensDky1(void *cvode_mem, realtype t, int k, int iss, N_Vector dkyS)
    int CVodeSetSensParams(void *cvode_mem, realtype *p, realtype *pbar, int *plist)
    int CVodeSetSensDQMethod(void *cvode_mem, int DQtype, realtype DQrhomax)
    int CVodeSetSensErrCon(void *cvode_mem, booleantype errconS)
    int CVodeSetSensMaxNonlinIters(void *cvode_mem, int maxcorS)
    int CVodeSetStabLimDet(void *cvode_mem, booleantype stldet)
    
    
    
    #Statistics
    int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele)               #Estimated local errors
    int CVodeGetErrWeights(void *cvode_mem, N_Vector eweight)               #Estimated local errors
    int CVodeGetSensNumRhsEvals(void *cvode_mem, long int *nfSevals)
    int CVodeGetNumRhsEvalsSens(void *cvode_mem, long int *nfevalsS)
    int CVodeGetSensNumErrTestFails(void *cvode_mem, long int *nSetfails)
    int CVodeGetSensNumLinSolvSetups(void *cvode_mem, long int *nlinsetupsS)
    int CVodeGetSensStats(void *cvode_mem, long int *nfSevals, long int *nfevalsS,
                         long int *nSetfails, long int *nlinsetupsS)
    int CVodeGetSensErrWeights(void *cvode_mem, N_Vector *eSweight)
    int CVodeGetSensNumNonlinSolvIters(void *cvode_mem, long int *nSniters)
    int CVodeGetSensNumNonlinSolvConvFails(void *cvode_mem, long int *nSncfails)
    int CVodeGetSensNonlinSolvStats(void *cvode_mem, long int *nSniters, long int *nSncfails)
    int CVodeGetStgrSensNumNonlinSolvIters(void *cvode_mem, long int *nSTGR1niters)
    int CVodeGetStgrSensNumNonlinSolvConvFails(void *cvode_mem, long int *nSTGR1ncfails)

IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "cvodes/cvodes_ls.h":
        ctypedef int (*CVLsJacFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
        ctypedef int (*CVLsPrecSolveFn)(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, sunrealtype gamma, sunrealtype delta, int lr, void* user_data);
        ctypedef int (*CVLsJacTimesSetupFn)(realtype t, N_Vector y, N_Vector fy, void *user_data)
        ctypedef int (*CVLsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp)
ELSE:
    cdef extern from "cvodes/cvodes_spils.h":
        ctypedef int (*CVSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp)


IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "cvode/cvode_ls.h":
       int CVodeSetJacTimes(void *cvode_mem, CVLsJacTimesSetupFn jtsetup, CVLsJacTimesVecFn jtimes)
       int CVodeSetLinearSolver(void* cvode_mem, SUNLinearSolver LS, SUNMatrix A)
       int CVodeSetJacFn(void* cvode_mem, CVLsJacFn jac)
    IF SUNDIALS_WITH_SUPERLU:
        cdef extern from "sunlinsol/sunlinsol_superlumt.h":
            SUNLinearSolver SUNLinSol_SuperLUMT(N_Vector y, SUNMatrix A, int num_threads, SUNContext ctx)
    ELSE:
        cdef inline SUNLinearSolver SUNLinSol_SuperLUMT(N_Vector y, SUNMatrix A, int num_threads, SUNContext ctx): return NULL

    cdef inline int cv_spils_jtsetup_dummy(realtype t, N_Vector y, N_Vector fy, void *user_data): return 0
    cdef inline tuple version(): return (6,0,0)
ELIF SUNDIALS_VERSION >= (3,0,0):
    cdef extern from "cvodes/cvodes_direct.h":
        ctypedef int (*CVDlsDenseJacFn)(realtype t, N_Vector y, N_Vector fy, 
                       SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
        IF SUNDIALS_VERSION >= (4,0,0):
            int CVodeSetLinearSolver(void *cvode_mem, SUNLinearSolver LS, SUNMatrix A)
            int CVodeSetJacFn(void *cvode_mem, CVDlsDenseJacFn djac)
        ELSE:
            int CVDlsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS, SUNMatrix A)
            int CVDlsSetJacFn(void *cvode_mem, CVDlsDenseJacFn djac)
    cdef extern from "cvodes/cvodes_spils.h":
        ctypedef int (*CVSpilsJacTimesSetupFn)(realtype t, N_Vector y, N_Vector fy, void *user_data)
        IF SUNDIALS_VERSION >= (4,0,0):
            int CVodeSetJacTimes(void *cvode_mem, CVSpilsJacTimesSetupFn jtsetup, CVSpilsJacTimesVecFn jtimes)
        ELSE:
            int CVSpilsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS)
            int CVSpilsSetJacTimes(void *cvode_mem, CVSpilsJacTimesSetupFn jtsetup, CVSpilsJacTimesVecFn jtimes)
        ctypedef int (*CVSpilsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy,
				  booleantype jok, booleantype *jcurPtr, realtype gamma, void *user_data)
        ctypedef int (*CVSpilsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy,
				  N_Vector r, N_Vector z,
				  realtype gamma, realtype delta, int lr, void *user_data)

    IF SUNDIALS_WITH_SUPERLU:
        cdef extern from "sunlinsol/sunlinsol_superlumt.h":
            SUNLinearSolver SUNSuperLUMT(N_Vector y, SUNMatrix A, int num_threads)
    ELSE:
         cdef inline SUNLinearSolver SUNSuperLUMT(N_Vector y, SUNMatrix A, int num_threads): return NULL
    
    cdef inline int cv_spils_jtsetup_dummy(realtype t, N_Vector y, N_Vector fy, void *user_data): return 0
    cdef inline tuple version(): return (3,0,0)
ELSE:
    cdef extern from "cvodes/cvodes_dense.h":
        int CVDense(void *cvode_mem, long int n)
        ctypedef int (*CVDlsDenseJacFn)(long int n, realtype t, N_Vector y, N_Vector fy, 
                       DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
        int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn djac)

    cdef extern from "cvodes/cvodes_spgmr.h":
        int CVSpgmr(void *cvode_mem, int pretype, int max1)
    
    cdef extern from "cvodes/cvodes_spils.h":
        int CVSpilsSetJacTimesVecFn(void *cvode_mem,  CVSpilsJacTimesVecFn jtv)
        ctypedef int (*CVSpilsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy,
				  booleantype jok, booleantype *jcurPtr,
				  realtype gamma, void *user_data,
				  N_Vector tmp1, N_Vector tmp2,
				  N_Vector tmp3)
        ctypedef int (*CVSpilsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy,
				  N_Vector r, N_Vector z,
				  realtype gamma, realtype delta,
				  int lr, void *user_data, N_Vector tmp)
    
    IF SUNDIALS_VERSION >= (2,6,0):
        cdef extern from "cvodes/cvodes_sparse.h":
            ctypedef int (*CVSlsSparseJacFn)(realtype t, N_Vector y, N_Vector fy,
                                      SlsMat Jac, void *user_data, N_Vector tmp1,
                                        N_Vector tmp2, N_Vector tmp3)
            int CVSlsSetSparseJacFn(void *cvode_mem, CVSlsSparseJacFn jac)
            int CVSlsGetNumJacEvals(void *cvode_mem, long int *njevals)
        cdef inline tuple version(): return (2,6,0)
        IF SUNDIALS_WITH_SUPERLU:
            cdef extern from "cvodes/cvodes_superlumt.h":
                int CVSuperLUMT(void *cvode_mem, int numthreads, int n, int nnz)
        ELSE:
            cdef inline int CVSuperLUMT(void *cvode_mem, int numthreads, int n, int nnz): return -1
    ELSE:
        cdef inline int CVSuperLUMT(void *cvode_mem, int numthreads, int n, int nnz): return -1
        ctypedef int (*CVSlsSparseJacFn)(realtype t, N_Vector y, N_Vector fy,
                                  SlsMat Jac, void *user_data, N_Vector tmp1,
                                    N_Vector tmp2, N_Vector tmp3)
        cdef inline int CVSlsSetSparseJacFn(void *cvode_mem, CVSlsSparseJacFn jac): return -1
        cdef inline int CVSlsGetNumJacEvals(void *cvode_mem, long int *njevals): return -1
        cdef inline tuple version(): return (2,5,0)

IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "cvode/cvode_ls.h":
        ctypedef int (*CVLsPrecSetupFn)(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok, sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)
        int CVodeSetPreconditioner(void* cvode_mem, CVLsPrecSetupFn pset, CVLsPrecSolveFn psolve)
        int CVodeGetNumJtimesEvals(void *cvode_mem, long int *njvevals) #Number of jac*vector evals
        int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS) #Number of res evals due to jacÄvector evals
        int CVodeGetNumPrecEvals(void *cvode_mem, long int *npevals)
        int CVodeGetNumPrecSolves(void *cvode_mem, long int *npsolves)
ELSE:
    cdef extern from "cvodes/cvodes_spils.h":
        IF SUNDIALS_VERSION >= (4,0,0):
            int CVodeSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn psetup, CVSpilsPrecSolveFn psolve)
            int CVodeGetNumJtimesEvals(void *cvode_mem, long int *njvevals) #Number of jac*vector evals
            int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS) #Number of res evals due to jacÄvector evals
            int CVodeGetNumPrecEvals(void *cvode_mem, long int *npevals)
            int CVodeGetNumPrecSolves(void *cvode_mem, long int *npsolves)
        ELSE:
            int CVSpilsSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn psetup, CVSpilsPrecSolveFn psolve)
            int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals) #Number of jac*vector evals
            int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS) #Number of res evals due to jacÄvector evals
            int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals)
            int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves)

cdef extern from "idas/idas.h":
    ctypedef int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
    IF SUNDIALS_VERSION >= (6,0,0):
        void* IDACreate(SUNContext ctx)
    ELSE:
        void* IDACreate()
    int IDAInit(void* ida_mem, IDAResFn res, realtype t0, N_Vector y0, N_Vector yp0)
    int IDAReInit(void* ida_mem, realtype t0, N_Vector y0, N_Vector yp0)
    void IDAFree(void **ida_mem)
    int IDASolve(void* ida_mem, realtype tout,realtype  *tret, N_Vector yret, 
                            N_Vector ypret, int itask)
    
    #Functions for settings options
    int IDASStolerances(void *ida_mem, realtype reltol, realtype abstol)
    int IDASVtolerances(void *ida_mem, realtype reltol, N_Vector abstol)
    int IDASetSuppressAlg(void *ida_mem, booleantype suppressalg)
    int IDASetId(void *ida_mem, N_Vector id)
    int IDASetUserData(void *ida_mem,void *user_data)
    int IDASetInitStep(void *ida_mem, realtype hin)
    int IDASetStopTime(void *ida_mem, realtype tstop)
    int IDASetMaxErrTestFails(void *ida_mem, int maxnef)
    int IDASetMaxNumSteps(void *ida_mem, long int mxsteps)
    int IDASetMaxOrd(void *ida_mem, int maxord)
    int IDASetMaxStep(void* ida_mem, realtype hmax)
    
    #Functions for retrieving results
    int IDAGetDky(void *ida_mem, realtype t, int k, N_Vector dky)

    #Functions for discontinuity handling
    ctypedef int (*IDARootFn)(realtype tt, N_Vector yy, N_Vector yp, realtype *gout, void *user_data)
    int IDASetRootDirection(void *ida_mem, int *rootdir)
    int IDASetNoInactiveRootWarn(void *ida_mem)
    int IDARootInit(void *ida_mem, int nrtfn, IDARootFn g)
    int IDAGetRootInfo(void *ida_mem, int *rootsfound)
    int IDACalcIC(void *ida_men, int icopt, realtype tout1)
    int IDAGetConsistentIC(void *ida_mem, N_Vector y0, N_Vector yp0)
    int IDASetLineSearchOffIC(void *ida_mem, booleantype lsoff)
    
    #Functions for retrieving statistics
    int IDAGetEstLocalErrors(void *ida_mem, N_Vector ele)               #Estimated local errors
    int IDAGetErrWeights(void *ida_mem, N_Vector eweight)
    int IDAGetLastStep(void *ida_mem, realtype *hlast)
    int IDAGetLastOrder(void *ida_mem,int *qlast)                       #Last order used
    int IDAGetCurrentOrder(void *ida_mem,int *qcurrent)                 #Order that is about to be tried
    int IDAGetNumSteps(void *ida_mem, long int *nsteps)                 #Number of steps
    int IDAGetNumResEvals(void *ida_mem, long int *nrevals)             #Number of res evals
    IF SUNDIALS_VERSION >= (4,0,0):
        int IDAGetNumJacEvals(void *ida_mem, long int *njevals)          #Number of jac evals
        int IDAGetNumResEvals(void *ida_mem, long int *nrevalsLS)        #Number of res evals due to jac evals
    ELSE:
        int IDADlsGetNumJacEvals(void *ida_mem, long int *njevals)          #Number of jac evals
        int IDADlsGetNumResEvals(void *ida_mem, long int *nrevalsLS)        #Number of res evals due to jac evals
    int IDAGetNumGEvals(void *ida_mem, long int *ngevals)               #Number of root evals
    int IDAGetNumErrTestFails(void *ida_mem, long int *netfails)        #Number of local error test failures
    int IDAGetNumNonlinSolvIters(void *ida_mem, long int *nniters)      #Number of nonlinear iteration
    int IDAGetNumNonlinSolvConvFails(void *ida_mem, long int *nncfails) #Number of nonlinear conv failures
    int IDAGetNonlinSolvStats(void *ida_mem, long int *nniters, long int *nncfails)
    int IDAGetIntegratorStats(void* ida_mem,long int  *nsteps, long int *nrevals, 
                            long int *nlinsetups, long int *netfails, int *klast, 
                            int *kcur, realtype *hinused, realtype *hlast, 
                            realtype *hcur, realtype *tcur)
    
    #Start Sensitivities
    #===================
    ctypedef int (*IDASensResFn)(int Ns, realtype t, N_Vector yy, N_Vector yp, 
                                N_Vector *yS, N_Vector *ypS, N_Vector *resvalS, 
                                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    int IDASensInit(void *ida_mem, int Ns, int ism, IDASensResFn resS, N_Vector *ySO, N_Vector *ypSO)
    int IDASensReInit(void *ida_mem, int ism, N_Vector *ySO, N_Vector *ypSO)
    
    #Options
    int IDASensToggleOff(void *ida_mem)
    int IDASensSStolerances(void *ida_mem, realtype reltolS, realtype *abstolS)
    int IDASensSVtolerances(void *ida_mem, realtype reltolS, N_Vector *abstolS)
    int IDASensEEtolerances(void *ida_mem)
    
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
    int IDAGetSensStats(void *ida_mem, long int *nfSevals, long int *nfevalsS, 
                        long int *nSetfails, long int *nlinsetupsS)
    int IDAGetSensNumNonlinSolvIters(void *ida_mem, long int nSniters)
    int IDAGetSeonsNumNonlinSolvConvFails(void *ida_mem, long int nSncfails)
    int IDAGetSensNonlinSolvStats(void *ida_mem, long int *nSniters, long int *nSncfails)
    
    #End Sensitivities
    #=================

IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "ida/ida_ls.h":
        ctypedef int (*IDALsJacFn)(sunrealtype t, sunrealtype c_j, N_Vector y, N_Vector yp, N_Vector r, SUNMatrix Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
        ctypedef int (*IDALsJacTimesSetupFn)(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, sunrealtype c_j, void* user_data)
        ctypedef int (*IDALsJacTimesVecFn)(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector v, N_Vector Jv, sunrealtype c_j, void* user_data, N_Vector tmp1, N_Vector tmp2)
ELSE:
    cdef extern from "idas/idas_spils.h":
        ctypedef int (*IDASpilsJacTimesVecFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, 
                N_Vector v, N_Vector Jv, realtype cj, void *user_data,N_Vector tmp1, N_Vector tmp2)

IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "ida/ida_ls.h":
        int IDASetJacFn(void* ida_mem, IDALsJacFn jac)
        int IDASetLinearSolver(void* ida_mem, SUNLinearSolver LS, SUNMatrix A)
        int IDASetJacTimes(void* ida_mem, IDALsJacTimesSetupFn jtsetup, IDALsJacTimesVecFn jtimes)

    cdef inline int ida_spils_jtsetup_dummy(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, realtype c_j, void *user_data): return 0
ELIF SUNDIALS_VERSION >= (3,0,0):
    cdef extern from "idas/idas_direct.h":
        ctypedef int (*IDADlsDenseJacFn)(realtype tt, realtype cj, N_Vector yy, 
                       N_Vector yp, N_Vector rr, SUNMatrix Jac, void *user_data, 
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
        IF SUNDIALS_VERSION >= (4,0,0):
            int IDASetJacFn(void *ida_mem, IDADlsDenseJacFn djac)
            int IDASetLinearSolver(void *ida_mem, SUNLinearSolver LS, SUNMatrix A)
        ELSE:
            int IDADlsSetJacFn(void *ida_mem, IDADlsDenseJacFn djac)
            int IDADlsSetLinearSolver(void *ida_mem, SUNLinearSolver LS, SUNMatrix A)
    
    cdef extern from "idas/idas_spils.h":
        ctypedef int (*IDASpilsJacTimesSetupFn)(realtype tt, N_Vector yy,
                      N_Vector yp, N_Vector rr, realtype c_j, void *user_data)
        IF SUNDIALS_VERSION >= (4,0,0):
            int IDASetJacTimes(void *ida_mem,
                IDASpilsJacTimesSetupFn jtsetup, IDASpilsJacTimesVecFn jtimes)
        ELSE:
            int IDASpilsSetLinearSolver(void *ida_mem, SUNLinearSolver LS)
            int IDASpilsSetJacTimes(void *ida_mem,
                IDASpilsJacTimesSetupFn jtsetup, IDASpilsJacTimesVecFn jtimes)


                
    cdef inline int ida_spils_jtsetup_dummy(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, realtype c_j, void *user_data): return 0
ELSE:
    cdef extern from "idas/idas_dense.h":
        int IDADense(void *ida_mem, long int n)
        ctypedef int (*IDADlsDenseJacFn)(long int Neq, realtype tt, realtype cj, N_Vector yy, 
                       N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, 
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
        int IDADlsSetDenseJacFn(void *ida_mem, IDADlsDenseJacFn djac)
    
    cdef extern from "idas/idas_spgmr.h":
        int IDASpgmr(void *ida_mem, int max1)
        
    cdef extern from "idas/idas_spils.h":
        int IDASpilsSetJacTimesVecFn(void *ida_mem, IDASpilsJacTimesVecFn ida_jacv)

IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "ida/ida_ls.h": 
        int IDAGetNumJtimesEvals(void *ida_mem, long int *njvevals) #Number of jac*vector
        int IDAGetNumResEvals(void *ida_mem, long int *nfevalsLS) #Number of rhs due to jac*vector
ELSE:
    cdef extern from "idas/idas_spils.h":
        IF SUNDIALS_VERSION >= (4,0,0):
            int IDAGetNumJtimesEvals(void *ida_mem, long int *njvevals) #Number of jac*vector
            int IDAGetNumResEvals(void *ida_mem, long int *nfevalsLS) #Number of rhs due to jac*vector
        ELSE:
            int IDASpilsGetNumJtimesEvals(void *ida_mem, long int *njvevals) #Number of jac*vector
            int IDASpilsGetNumResEvals(void *ida_mem, long int *nfevalsLS) #Number of rhs due to jac*vector


####################
# KINSOL
####################


# KINSOL functions and routines
cdef extern from "kinsol/kinsol.h":
    # user defined functions
    ctypedef int (*KINSysFn)(N_Vector uu, N_Vector fval, void *user_data )
    ctypedef void (*KINInfoHandlerFn)(const char *module, const char *function, char *msg, void *user_data)
    # initialization routines
    IF SUNDIALS_VERSION >= (6,0,0):
        void *KINCreate(SUNContext ctx)
    ELSE:
        void *KINCreate()
    int KINInit(void *kinmem, KINSysFn func, N_Vector tmpl)

    # optional input spec. functions,
    # for specificationsdocumentation cf. kinsol.h line 218-449
    IF SUNDIALS_VERSION < (7,0,0):
        int KINSetInfoHandlerFn(void *kinmem, KINInfoHandlerFn ihfun, void *ih_data)
        int KINSetPrintLevel(void *kinmemm, int printfl)
    int KINSetUserData(void *kinmem, void *user_data)
    int KINSetNumMaxIters(void *kinmem, long int mxiter)
    int KINSetNoInitSetup(void *kinmem, booleantype noInitSetup)
    int KINSetNoResMon(void *kinmem, booleantype noNNIResMon)
    int KINSetMaxSetupCalls(void *kinmem, long int msbset)
    int KINSetMaxSubSetupCalls(void *kinmem, long int msbsetsub)
    int KINSetEtaForm(void *kinmem, int etachoice)
    int KINSetEtaConstValue(void *kinmem, realtype eta)
    int KINSetEtaParams(void *kinmem, realtype egamma, realtype ealpha)
    int KINSetResMonParams(void *kinmem, realtype omegamin, realtype omegamax)
    int KINSetResMonConstValue(void *kinmem, realtype omegaconst)
    int KINSetNoMinEps(void *kinmem, booleantype noMinEps)
    int KINSetMaxNewtonStep(void *kinmem, realtype mxnewtstep)
    int KINSetMaxBetaFails(void *kinmem, long int mxnbcf)
    int KINSetRelErrFunc(void *kinmem, realtype relfunc)
    int KINSetFuncNormTol(void *kinmem, realtype fnormtol)
    int KINSetScaledStepTol(void *kinmem, realtype scsteptol)
    int KINSetConstraints(void *kinmem, N_Vector constraints)
    int KINSetSysFunc(void *kinmem, KINSysFn func)

    # solver routine
    int KINSol(void *kinmem, N_Vector uu, int strategy, N_Vector u_scale, N_Vector f_scale)

    # optional output routines.
    # Documentation see kinsol.h line 670-735
    int KINGetWorkSpace(void *kinmem, long int *lenrw, long int *leniw)
    int KINGetNumNonlinSolvIters(void *kinmem, long int *nniters)
    int KINGetNumFuncEvals(void *kinmem, long int *nfevals)
    int KINGetNumBetaCondFails(void *kinmem, long int *nbcfails) 
    int KINGetNumBacktrackOps(void *kinmem, long int *nbacktr)
    int KINGetFuncNorm(void *kinmem, realtype *fnorm)
    int KINGetStepLength(void *kinmem, realtype *steplength)
    char *KINGetReturnFlagName(int flag)

    # fuction used to deallocate memory used by KINSOL
    void KINFree(void **kinmem)

# Functions for error handling
IF SUNDIALS_VERSION < (7,0,0):
    cdef extern from "kinsol/kinsol.h":
        ctypedef void (*KINErrHandlerFn)(int error_code, char *module, char *function, char *msg, void *user_data)
        int KINSetErrHandlerFn(void *kinmem, KINErrHandlerFn ehfun, void *eh_data)
    cdef extern from "cvodes/cvodes.h":
        ctypedef void (*CVErrHandlerFn)(int error_code, const char *module, const char *function, char *msg, void *eh_data)
        int CVodeSetErrHandlerFn(void *cvode_mem, CVErrHandlerFn ehfun, void* eh_data)
    cdef extern from "idas/idas.h":
        ctypedef void (*IDAErrHandlerFn)(int error_code, const char *module, const char *function, char *msg, void *eh_data) 
        int IDASetErrHandlerFn(void *ida_mem,IDAErrHandlerFn ehfun, void* eh_data)

IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "kinsol/kinsol_ls.h":
        ctypedef int (*KINLsJacFn)(N_Vector u, N_Vector fu, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2)
        int KINSetLinearSolver(void* kinmem, SUNLinearSolver LS, SUNMatrix A)
        int KINSetJacFn(void* kinmem, KINLsJacFn jac)
ELIF SUNDIALS_VERSION >= (3,0,0):
    cdef extern from "kinsol/kinsol_direct.h":
        ctypedef int (*KINDlsDenseJacFn)(N_Vector u, N_Vector fu, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2)
        IF SUNDIALS_VERSION < (4,0,0):
            int KINDlsSetLinearSolver(void *kinmem, SUNLinearSolver LS, SUNMatrix A)
            int KINDlsSetJacFn(void *kinmem, KINDlsDenseJacFn djac)
        ELSE:
            int KINSetLinearSolver(void *kinmem, SUNLinearSolver LS, SUNMatrix A)
            int KINSetJacFn(void *kinmem, KINDlsDenseJacFn djac)
    cdef extern from "kinsol/kinsol_spils.h":
        int KINSpilsSetLinearSolver(void *kinsol_mem, SUNLinearSolver LS)
        
        ctypedef int (*KINSpilsPrecSolveFn)(N_Vector u, N_Vector uscale,
                    N_Vector fval, N_Vector fscale, N_Vector v, void *problem_data)
        ctypedef int (*KINSpilsPrecSetupFn)(N_Vector u, N_Vector uscale,
                    N_Vector fval, N_Vector fscale, void *problem_data)
ELSE:
    # functions used for supplying jacobian, and receiving info from linear solver
    cdef extern from "kinsol/kinsol_direct.h":
        # user functions
        ctypedef int (*KINDlsDenseJacFn)(long int dim, N_Vector u, N_Vector fu, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2)
        
        # function used to link user functions to KINSOL
        int KINDlsSetDenseJacFn(void *kinmem, KINDlsDenseJacFn jac)
    
    cdef extern from "kinsol/kinsol_dense.h":
        int KINDense(void *kinmem, int dim)
    
    cdef extern from "kinsol/kinsol_spgmr.h":
        int KINSpgmr(void *kinmem, int maxl)
        
    cdef extern from "kinsol/kinsol_spils.h":
        ctypedef int (*KINSpilsPrecSolveFn)(N_Vector u, N_Vector uscale,
                    N_Vector fval, N_Vector fscale, N_Vector v, void *problem_data, N_Vector tmp)
        ctypedef int (*KINSpilsPrecSetupFn)(N_Vector u, N_Vector uscale,
                    N_Vector fval, N_Vector fscale, void *problem_data, N_Vector tmp1, N_Vector tmp2)

IF SUNDIALS_VERSION >= (6,0,0):
    cdef extern from "kinsol/kinsol_ls.h":
        ctypedef int (*KINLsPrecSetupFn)(N_Vector uu, N_Vector uscale, N_Vector fval, N_Vector fscale, void* user_data);
        ctypedef int (*KINLsPrecSolveFn)(N_Vector uu, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector vv, void* user_data)
        ctypedef int (*KINLsJacTimesVecFn)(N_Vector v, N_Vector Jv, N_Vector uu, sunbooleantype* new_uu, void* J_data)
        int KINGetLastLinFlag(void *kinmem, long int *flag)
        int KINGetNumJacEvals(void *kinmem, long int *njevalsB)
        int KINGetNumFuncEvals(void *kinmem, long int *nfevalsB)
        int KINSetJacTimesVecFn(void* kinmem, KINLsJacTimesVecFn jtv)
        int KINSetPreconditioner(void* kinmem, KINLsPrecSetupFn psetup, KINLsPrecSolveFn psolve)
        int KINGetNumLinIters(void *kinmem, long int *nliters)
        int KINGetNumLinConvFails(void *kinmem, long int *nlcfails)
        int KINGetNumPrecEvals(void *kinmem, long int *npevals)
        int KINGetNumPrecSolves(void *kinmem, long int *npsolves)
        int KINGetNumJtimesEvals(void *kinmem, long int *njevals)
        int KINGetNumFuncEvals(void *kinmem, long int *nfevalsLS)
ELSE:
    cdef extern from "kinsol/kinsol_direct.h":
        # optional output fcts for linear direct solver
        int KINDlsGetWorkSpace(void *kinmem, long int *lenrwB, long int *leniwB)
        IF SUNDIALS_VERSION >= (4,0,0):
            int KINGetLastLinFlag(void *kinmem, long int *flag)
            int KINGetNumJacEvals(void *kinmem, long int *njevalsB)
            int KINGetNumFuncEvals(void *kinmem, long int *nfevalsB)
        ELSE:
            int KINDlsGetLastFlag(void *kinmem, long int *flag)
            int KINDlsGetNumJacEvals(void *kinmem, long int *njevalsB)
            int KINDlsGetNumFuncEvals(void *kinmem, long int *nfevalsB)
        char *KINDlsGetReturnFlagName(int flag)

    cdef extern from "kinsol/kinsol_spils.h":
        ctypedef int (*KINSpilsJacTimesVecFn)(N_Vector vv, N_Vector Jv, N_Vector vx, int* new_u,
                    void *problem_data)
        IF SUNDIALS_VERSION >= (4,0,0):
            int KINSetJacTimesVecFn(void *kinmem, KINSpilsJacTimesVecFn jacv)
            int KINSetPreconditioner(void *kinmem, KINSpilsPrecSetupFn psetup, KINSpilsPrecSolveFn psolve)
            int KINGetNumLinIters(void *kinmem, long int *nliters)
            int KINGetNumLinConvFails(void *kinmem, long int *nlcfails)
            int KINGetNumPrecEvals(void *kinmem, long int *npevals)
            int KINGetNumPrecSolves(void *kinmem, long int *npsolves)
            int KINGetNumJtimesEvals(void *kinmem, long int *njevals)
            int KINGetNumFuncEvals(void *kinmem, long int *nfevalsLS)
        ELSE:
            int KINSpilsSetJacTimesVecFn(void *kinmem, KINSpilsJacTimesVecFn jacv)
            int KINSpilsSetPreconditioner(void *kinmem, KINSpilsPrecSetupFn psetup, KINSpilsPrecSolveFn psolve)
            int KINSpilsGetNumLinIters(void *kinmem, long int *nliters)
            int KINSpilsGetNumConvFails(void *kinmem, long int *nlcfails)
            int KINSpilsGetNumPrecEvals(void *kinmem, long int *npevals)
            int KINSpilsGetNumPrecSolves(void *kinmem, long int *npsolves)
            int KINSpilsGetNumJtimesEvals(void *kinmem, long int *njevals)
            int KINSpilsGetNumFuncEvals(void *kinmem, long int *nfevalsLS)

#=========================
# END SUNDIALS DEFINITIONS
#=========================
