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

from numpy cimport NPY_DOUBLE, npy_intp, NPY_INT

#==============================================
#External definitions from Sundials headers
#==============================================

IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "sundials/sundials_context.h":
        ctypedef _SUNContext * SUNContext
        cdef struct _SUNContext:
            pass
        IF SUNDIALS_VERSION_NR >= 700000:
            ctypedef int SUNComm
            int SUNContext_Create(SUNComm comm, SUNContext* ctx) noexcept
        ELSE:
            int SUNContext_Create(void* comm, SUNContext* ctx) noexcept

IF SUNDIALS_VERSION_NR >= 700000:
    cdef extern from "sundials/sundials_context.h":
        ctypedef int SUNErrCode
        ctypedef void (*SUNErrHandlerFn)(int line, const char* func, const char* file, const char* msg, SUNErrCode err_code, void* err_user_data, SUNContext sunctx) noexcept
        SUNErrCode SUNContext_PushErrHandler(SUNContext sunctx, SUNErrHandlerFn err_fn, void* err_user_data) noexcept

IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "sundials/sundials_types.h":
        ctypedef double sunrealtype
        ctypedef bint sunbooleantype
        IF SUNDIALS_VERSION_NR >= 700000:
            cdef int SUN_COMM_NULL
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
    cdef N_Vector N_VMake_Serial(long int vec_length, realtype *v_data) noexcept
    void N_VSetArrayPointer_Serial(realtype *v_data, N_Vector v) noexcept
    void N_VConst_Serial(realtype c, N_Vector z) noexcept
    IF SUNDIALS_VERSION_NR >= 600000:
        N_Vector N_VNew_Serial(long int vec_length, SUNContext ctx) noexcept
        N_Vector *N_VCloneVectorArray(int count, N_Vector w) noexcept
        N_Vector *N_VCloneVectorArrayEmpty(int count, N_Vector w) noexcept
        void N_VDestroy(N_Vector v) noexcept
    ELSE:
        N_Vector N_VNew_Serial(long int vec_length) noexcept
        N_Vector *N_VCloneVectorArray_Serial(int count, N_Vector w) noexcept
        N_Vector *N_VCloneVectorArrayEmpty_Serial(int count, N_Vector w) noexcept
        void N_VDestroy_Serial(N_Vector v) noexcept
    void N_VPrint_Serial(N_Vector v) noexcept


IF SUNDIALS_VERSION_NR >= 400000:
    cdef extern from "sundials/sundials_nonlinearsolver.h":
        ctypedef _generic_SUNNonlinearSolver *SUNNonlinearSolver
        
        cdef struct _generic_SUNNonlinearSolver:
            pass
ELSE:
    #Dummy defines
    ctypedef void *SUNNonlinearSolver
        
IF SUNDIALS_VERSION_NR >= 300000:
    cdef extern from "sundials/sundials_types.h":
        IF SUNDIALS_VECTOR_SIZE == "64":
            ctypedef long int sunindextype
        ELSE:
            ctypedef int sunindextype
    cdef extern from "sundials/sundials_matrix.h":
        ctypedef _generic_SUNMatrix *SUNMatrix
        void SUNMatDestroy(SUNMatrix A) noexcept
        
        cdef struct _generic_SUNMatrix_Ops:
            SUNMatrix_ID (*getid)(SUNMatrix) noexcept
            SUNMatrix    (*clone)(SUNMatrix) noexcept
            void         (*destroy)(SUNMatrix) noexcept
            int          (*zero)(SUNMatrix) noexcept
            int          (*copy)(SUNMatrix, SUNMatrix) noexcept
            int          (*scaleadd)(realtype, SUNMatrix, SUNMatrix) noexcept
            int          (*scaleaddi)(realtype, SUNMatrix) noexcept
            int          (*matvec)(SUNMatrix, N_Vector, N_Vector) noexcept
            int          (*space)(SUNMatrix, long int*, long int*) noexcept

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
        int SUNLinSolFree(SUNLinearSolver S) noexcept
        
        cdef struct _generic_SUNLinearSolver_Ops:
            SUNLinearSolver_Type (*gettype)(SUNLinearSolver) noexcept
            int                  (*setatimes)(SUNLinearSolver, void*, ATimesFn) noexcept
            int                  (*setpreconditioner)(SUNLinearSolver, void*, PSetupFn, PSolveFn) noexcept
            int                  (*setscalingvectors)(SUNLinearSolver, N_Vector, N_Vector) noexcept
            int                  (*initialize)(SUNLinearSolver) noexcept
            int                  (*setup)(SUNLinearSolver, SUNMatrix) noexcept
            int                  (*solve)(SUNLinearSolver, SUNMatrix, N_Vector, N_Vector, realtype) noexcept
            int                  (*numiters)(SUNLinearSolver) noexcept
            realtype             (*resnorm)(SUNLinearSolver) noexcept
            long int             (*lastflag)(SUNLinearSolver) noexcept
            int                  (*space)(SUNLinearSolver, long int*, long int*) noexcept
            N_Vector             (*resid)(SUNLinearSolver) noexcept
            int                  (*free)(SUNLinearSolver) noexcept
        
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
        IF SUNDIALS_VERSION_NR >= 600000:
            SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N, SUNContext ctx) noexcept
        ELSE:
            SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N) noexcept
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
        IF SUNDIALS_VERSION_NR >= 600000:
            SUNMatrix SUNSparseMatrix(sunindextype M, sunindextype N, sunindextype NNZ, int sparsetype, SUNContext ctx) noexcept
        ELSE:
            SUNMatrix SUNSparseMatrix(sunindextype M, sunindextype N, sunindextype NNZ, int sparsetype) noexcept
    cdef extern from "sunlinsol/sunlinsol_dense.h":
        IF SUNDIALS_VERSION_NR >= 400000:
            IF SUNDIALS_VERSION_NR >= 600000:
                SUNLinearSolver SUNLinSol_Dense(N_Vector y, SUNMatrix A, SUNContext ctx) noexcept
            ELSE:
                SUNLinearSolver SUNLinSol_Dense(N_Vector y, SUNMatrix A) noexcept
        ELSE:
            SUNLinearSolver SUNDenseLinearSolver(N_Vector y, SUNMatrix A) noexcept
    cdef extern from "sunlinsol/sunlinsol_spgmr.h":
        IF SUNDIALS_VERSION_NR >= 400000:
            IF SUNDIALS_VERSION_NR >= 600000:
                SUNLinearSolver SUNLinSol_SPGMR(N_Vector y, int pretype, int maxl, SUNContext ctx) noexcept
            ELSE:
                SUNLinearSolver SUNLinSol_SPGMR(N_Vector y, int pretype, int maxl) noexcept
        ELSE:
            SUNLinearSolver SUNSPGMR(N_Vector y, int pretype, int maxl) noexcept

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
    cdef realtype* DENSE_COL(DlsMat A, int j) noexcept

IF SUNDIALS_VERSION_NR >= 500000:
    pass
ELIF SUNDIALS_VERSION_NR >= 200603:
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
ELIF SUNDIALS_VERSION_NR >= 200600:
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
    cdef inline int with_superlu() noexcept: return 1
ELSE:
    cdef inline int with_superlu() noexcept: return 0

IF SUNDIALS_VERSION_NR >= 400000:
    cdef extern from "cvodes/cvodes.h":
        IF SUNDIALS_VERSION_NR >= 600000:
            void* CVodeCreate(int lmm, SUNContext ctx) noexcept
        ELSE:
            void* CVodeCreate(int lmm) noexcept

        int CVodeSetNonlinearSolver(void *cvode_mem, SUNNonlinearSolver NLS) noexcept
        int CVodeSetNonlinearSolverSensSim(void *cvode_mem, SUNNonlinearSolver NLS) noexcept
        int CVodeSetNonlinearSolverSensStg(void *cvode_mem, SUNNonlinearSolver NLS) noexcept
    
    cdef extern from "sunnonlinsol/sunnonlinsol_newton.h":
        IF SUNDIALS_VERSION_NR >= 600000:
            SUNNonlinearSolver SUNNonlinSol_Newton(N_Vector y, SUNContext ctx) noexcept
            SUNNonlinearSolver SUNNonlinSol_NewtonSens(int count, N_Vector y, SUNContext ctx) noexcept
        ELSE:
            SUNNonlinearSolver SUNNonlinSol_Newton(N_Vector y) noexcept
            SUNNonlinearSolver SUNNonlinSol_NewtonSens(int count, N_Vector y) noexcept
    cdef extern from "sunnonlinsol/sunnonlinsol_fixedpoint.h":
        IF SUNDIALS_VERSION_NR >= 600000:
            SUNNonlinearSolver SUNNonlinSol_FixedPoint(N_Vector y, int m, SUNContext ctx) noexcept
            SUNNonlinearSolver SUNNonlinSol_FixedPointSens(int count, N_Vector y, int m, SUNContext ctx) noexcept
        ELSE:
            SUNNonlinearSolver SUNNonlinSol_FixedPoint(N_Vector y, int m) noexcept
            SUNNonlinearSolver SUNNonlinSol_FixedPointSens(int count, N_Vector y, int m) noexcept
ELSE:
    cdef extern from "cvodes/cvodes.h":
        void* CVodeCreate(int lmm, int iter) noexcept
        
cdef extern from "cvodes/cvodes.h":
    ctypedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *f_data) noexcept
    int CVodeInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0) noexcept
    int CVodeReInit(void *cvode_mem, realtype t0, N_Vector y0) noexcept
    void CVodeFree(void **cvode_mem) noexcept
    int CVode(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, int itask) noexcept
    
    #Functions for settings options
    int CVodeSetMaxOrd(void *cvode_mem, int maxord) noexcept
    int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps) noexcept
    int CVodeSetMaxStep(void   *cvode_mem, realtype hmax) noexcept
    int CVodeSetMinStep(void   *cvode_mem, realtype hmin) noexcept
    int CVodeSetInitStep(void  *cvode_mem, realtype hin) noexcept
    int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol) noexcept
    int CVodeSVtolerances(void *cvode_mem, realtype reltol, N_Vector abstol) noexcept
    IF SUNDIALS_CVODE_RTOL_VEC:
        int CVodeVVtolerances(void *cvode_mem, N_Vector reltol, N_Vector abstol) noexcept
    int CVodeSetStopTime(void  *cvode_mem, realtype tstop) noexcept
    int CVodeSetUserData(void  *cvode_mem,void *user_data) noexcept
    int CVodeSetMaxConvFails(void *cvode_mem, int maxncf) noexcept
    int CVodeSetMaxErrTestFails(void *cvode_mem, int maxnef) noexcept
    int CVodeSetMaxNonlinIters(void *cvode_mem, int maxcor) noexcept
    
    #Functions for retrieving results
    int CVodeGetDky(void *cvode_mem, realtype t, int k, N_Vector dky) noexcept
    #Functions for discontinuity handling
    ctypedef int (*CVRootFn)(realtype tt, N_Vector yy, realtype *gout, void *user_data) noexcept
    int CVodeRootDirection(void *cvode_mem, int *rootdir) noexcept
    int CVodeSetNoInactiveRootWarn(void *cvode_mem) noexcept
    int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g) noexcept
    int CVodeGetRootInfo(void *cvode_mem, int *rootsfound) noexcept
    
    #Functions for retrieving statistics
    int CVodeGetLastOrder(void * cvode_mem,int *qlast) noexcept
    int CVodeGetLastStep(void * cvode_mem, realtype *hlast) noexcept
    int CVodeGetCurrentOrder(void * cvode_mem,int *qcurrent) noexcept
    int CVodeGetActualInitStep(void * cvode_mem, realtype *hinused) noexcept
    int CVodeGetNumSteps(void *cvode_mem, long int *nsteps) noexcept #Number of steps
    int CVodeGetNumRhsEvals(void *cvode_mem, long int *nrevals) noexcept #Number of function evals
    IF SUNDIALS_VERSION_NR >= 400000:
        int CVodeGetNumJacEvals(void *cvode_mem, long int *njevals) noexcept #Number of jac evals
        int CVodeGetNumLinRhsEvals(void *cvode_mem, long int *nrevalsLS) noexcept #Number of res evals due to jac evals
    ELSE:
        int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals) noexcept #Number of jac evals
        int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nrevalsLS) noexcept #Number of res evals due to jac evals
    int CVodeGetNumGEvals(void *cvode_mem, long int *ngevals) noexcept #Number of root evals
    int CVodeGetNumErrTestFails(void *cvode_mem, long int *netfails) noexcept#Number of local error test failures
    int CVodeGetNumNonlinSolvIters(void *cvode_mem, long int *nniters) noexcept #Number of nonlinear iteration
    int CVodeGetNumNonlinSolvConvFails(void *cvode_mem, long int *nncfails) noexcept #Number of nonlinear conv failures
    int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters, long int *nncfails) noexcept
    int CVodeGetIntegratorStats(void* cvode_mem, long int *nsteps, long int *nfevals,
                                long int *nlinsetups, long int *netfails, int *qlast, int *qcur,
                                realtype *hinused, realtype *hlast, realtype *hcur, realtype *tcur) noexcept
    int CVodeGetNumStabLimOrderReds(void *cvode_mem, long int *nslred) noexcept
    
    #Sensitivity methods
    ctypedef int (*CVSensRhsFn)(int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS,
                                N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) noexcept
    ctypedef int (*CVSensRhs1Fn)(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector *yS,
                                N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2) noexcept
    int CVodeSensInit(void *cvode_mem, int Ns, int ism, CVSensRhsFn fS, N_Vector *ySO) noexcept
    int CVodeSensInit1(void *cvode_mem, int Ns, int ism, CVSensRhs1Fn fS1, N_Vector *ySO) noexcept
    int CVodeSensReInit(void *cvode_mem, int ism, N_Vector *ySO) noexcept
    int CVodeSensFree(void *cvode_mem) noexcept
    int CVodeSensToggleOff(void *cvode_mem) noexcept
    int CVodeSensEEtolerances(void *cvode_mem) noexcept
    int CVodeGetSens(void *cvode_mem, realtype *tret, N_Vector *yS) noexcept
    int CVodeGetSensDky(void *cvode_mem, realtype t, int k, N_Vector *dkyS) noexcept
    int CVodeGetSens1(void *cvode_mem, realtype *tret, int iss, N_Vector yS) noexcept
    int CVodeGetSensDky1(void *cvode_mem, realtype t, int k, int iss, N_Vector dkyS) noexcept
    int CVodeSetSensParams(void *cvode_mem, realtype *p, realtype *pbar, int *plist) noexcept
    int CVodeSetSensDQMethod(void *cvode_mem, int DQtype, realtype DQrhomax) noexcept
    int CVodeSetSensErrCon(void *cvode_mem, booleantype errconS) noexcept
    int CVodeSetSensMaxNonlinIters(void *cvode_mem, int maxcorS) noexcept
    int CVodeSetStabLimDet(void *cvode_mem, booleantype stldet) noexcept
    
    
    
    #Statistics
    int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele) noexcept #Estimated local errors
    int CVodeGetErrWeights(void *cvode_mem, N_Vector eweight) noexcept #Estimated local errors
    int CVodeGetSensNumRhsEvals(void *cvode_mem, long int *nfSevals) noexcept
    int CVodeGetNumRhsEvalsSens(void *cvode_mem, long int *nfevalsS) noexcept
    int CVodeGetSensNumErrTestFails(void *cvode_mem, long int *nSetfails) noexcept
    int CVodeGetSensNumLinSolvSetups(void *cvode_mem, long int *nlinsetupsS) noexcept
    int CVodeGetSensStats(void *cvode_mem, long int *nfSevals, long int *nfevalsS,
                         long int *nSetfails, long int *nlinsetupsS) noexcept
    int CVodeGetSensErrWeights(void *cvode_mem, N_Vector *eSweight) noexcept
    int CVodeGetSensNumNonlinSolvIters(void *cvode_mem, long int *nSniters) noexcept
    int CVodeGetSensNumNonlinSolvConvFails(void *cvode_mem, long int *nSncfails) noexcept
    int CVodeGetSensNonlinSolvStats(void *cvode_mem, long int *nSniters, long int *nSncfails) noexcept
    int CVodeGetStgrSensNumNonlinSolvIters(void *cvode_mem, long int *nSTGR1niters) noexcept
    int CVodeGetStgrSensNumNonlinSolvConvFails(void *cvode_mem, long int *nSTGR1ncfails) noexcept

IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "cvodes/cvodes_ls.h":
        ctypedef int (*CVLsJacFn)(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) noexcept
        ctypedef int (*CVLsPrecSolveFn)(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, N_Vector z, sunrealtype gamma, sunrealtype delta, int lr, void* user_data) noexcept
        ctypedef int (*CVLsJacTimesSetupFn)(realtype t, N_Vector y, N_Vector fy, void *user_data) noexcept
        ctypedef int (*CVLsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp) noexcept
ELSE:
    cdef extern from "cvodes/cvodes_spils.h":
        ctypedef int (*CVSpilsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy, void *user_data, N_Vector tmp) noexcept


IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "cvode/cvode_ls.h":
       int CVodeSetJacTimes(void *cvode_mem, CVLsJacTimesSetupFn jtsetup, CVLsJacTimesVecFn jtimes) noexcept
       int CVodeSetLinearSolver(void* cvode_mem, SUNLinearSolver LS, SUNMatrix A) noexcept
       int CVodeSetJacFn(void* cvode_mem, CVLsJacFn jac) noexcept
    IF SUNDIALS_WITH_SUPERLU:
        cdef extern from "sunlinsol/sunlinsol_superlumt.h":
            SUNLinearSolver SUNLinSol_SuperLUMT(N_Vector y, SUNMatrix A, int num_threads, SUNContext ctx) noexcept
    ELSE:
        cdef inline SUNLinearSolver SUNLinSol_SuperLUMT(N_Vector y, SUNMatrix A, int num_threads, SUNContext ctx) noexcept: return NULL

    cdef inline int cv_spils_jtsetup_dummy(realtype t, N_Vector y, N_Vector fy, void *user_data) noexcept: return 0
    cdef inline tuple version(): return (6,0,0)
ELIF SUNDIALS_VERSION_NR >= 300000:
    cdef extern from "cvodes/cvodes_direct.h":
        ctypedef int (*CVDlsDenseJacFn)(realtype t, N_Vector y, N_Vector fy, 
                       SUNMatrix Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) noexcept
        IF SUNDIALS_VERSION_NR >= 400000:
            int CVodeSetLinearSolver(void *cvode_mem, SUNLinearSolver LS, SUNMatrix A) noexcept
            int CVodeSetJacFn(void *cvode_mem, CVDlsDenseJacFn djac) noexcept
        ELSE:
            int CVDlsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS, SUNMatrix A) noexcept
            int CVDlsSetJacFn(void *cvode_mem, CVDlsDenseJacFn djac) noexcept
    cdef extern from "cvodes/cvodes_spils.h":
        ctypedef int (*CVSpilsJacTimesSetupFn)(realtype t, N_Vector y, N_Vector fy, void *user_data) noexcept
        IF SUNDIALS_VERSION_NR >= 400000:
            int CVodeSetJacTimes(void *cvode_mem, CVSpilsJacTimesSetupFn jtsetup, CVSpilsJacTimesVecFn jtimes) noexcept
        ELSE:
            int CVSpilsSetLinearSolver(void *cvode_mem, SUNLinearSolver LS) noexcept
            int CVSpilsSetJacTimes(void *cvode_mem, CVSpilsJacTimesSetupFn jtsetup, CVSpilsJacTimesVecFn jtimes) noexcept
        ctypedef int (*CVSpilsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy,
				  booleantype jok, booleantype *jcurPtr, realtype gamma, void *user_data) noexcept
        ctypedef int (*CVSpilsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy,
				  N_Vector r, N_Vector z,
				  realtype gamma, realtype delta, int lr, void *user_data) noexcept

    IF SUNDIALS_WITH_SUPERLU:
        cdef extern from "sunlinsol/sunlinsol_superlumt.h":
            SUNLinearSolver SUNSuperLUMT(N_Vector y, SUNMatrix A, int num_threads) noexcept
    ELSE:
         cdef inline SUNLinearSolver SUNSuperLUMT(N_Vector y, SUNMatrix A, int num_threads) noexcept: return NULL
    
    cdef inline int cv_spils_jtsetup_dummy(realtype t, N_Vector y, N_Vector fy, void *user_data) noexcept: return 0
    cdef inline tuple version() noexcept: return (3,0,0)
ELSE:
    cdef extern from "cvodes/cvodes_dense.h":
        int CVDense(void *cvode_mem, long int n) noexcept
        ctypedef int (*CVDlsDenseJacFn)(long int n, realtype t, N_Vector y, N_Vector fy, 
                       DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) noexcept
        int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn djac) noexcept

    cdef extern from "cvodes/cvodes_spgmr.h":
        int CVSpgmr(void *cvode_mem, int pretype, int max1) noexcept
    
    cdef extern from "cvodes/cvodes_spils.h":
        int CVSpilsSetJacTimesVecFn(void *cvode_mem,  CVSpilsJacTimesVecFn jtv) noexcept
        ctypedef int (*CVSpilsPrecSetupFn)(realtype t, N_Vector y, N_Vector fy,
				  booleantype jok, booleantype *jcurPtr,
				  realtype gamma, void *user_data,
				  N_Vector tmp1, N_Vector tmp2,
				  N_Vector tmp3) noexcept
        ctypedef int (*CVSpilsPrecSolveFn)(realtype t, N_Vector y, N_Vector fy,
				  N_Vector r, N_Vector z,
				  realtype gamma, realtype delta,
				  int lr, void *user_data, N_Vector tmp) noexcept
    
    IF SUNDIALS_VERSION_NR >= 200600:
        cdef extern from "cvodes/cvodes_sparse.h":
            ctypedef int (*CVSlsSparseJacFn)(realtype t, N_Vector y, N_Vector fy,
                                      SlsMat Jac, void *user_data, N_Vector tmp1,
                                        N_Vector tmp2, N_Vector tmp3) noexcept
            int CVSlsSetSparseJacFn(void *cvode_mem, CVSlsSparseJacFn jac) noexcept
            int CVSlsGetNumJacEvals(void *cvode_mem, long int *njevals) noexcept
        cdef inline tuple version(): return (2,6,0)
        IF SUNDIALS_WITH_SUPERLU:
            cdef extern from "cvodes/cvodes_superlumt.h":
                int CVSuperLUMT(void *cvode_mem, int numthreads, int n, int nnz) noexcept
        ELSE:
            cdef inline int CVSuperLUMT(void *cvode_mem, int numthreads, int n, int nnz) noexcept: return -1
    ELSE:
        cdef inline int CVSuperLUMT(void *cvode_mem, int numthreads, int n, int nnz) noexcept: return -1
        ctypedef int (*CVSlsSparseJacFn)(realtype t, N_Vector y, N_Vector fy,
                                  SlsMat Jac, void *user_data, N_Vector tmp1,
                                    N_Vector tmp2, N_Vector tmp3) noexcept
        cdef inline int CVSlsSetSparseJacFn(void *cvode_mem, CVSlsSparseJacFn jac) noexcept: return -1
        cdef inline int CVSlsGetNumJacEvals(void *cvode_mem, long int *njevals) noexcept: return -1
        cdef inline tuple version() noexcept: return (2,5,0)

IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "cvode/cvode_ls.h":
        ctypedef int (*CVLsPrecSetupFn)(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok, sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data) noexcept
        int CVodeSetPreconditioner(void* cvode_mem, CVLsPrecSetupFn pset, CVLsPrecSolveFn psolve) noexcept
        int CVodeGetNumJtimesEvals(void *cvode_mem, long int *njvevals) noexcept #Number of jac*vector evals
        int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS) noexcept #Number of res evals due to jac*vector evals
        int CVodeGetNumPrecEvals(void *cvode_mem, long int *npevals) noexcept
        int CVodeGetNumPrecSolves(void *cvode_mem, long int *npsolves) noexcept
ELSE:
    cdef extern from "cvodes/cvodes_spils.h":
        IF SUNDIALS_VERSION_NR >= 400000:
            int CVodeSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn psetup, CVSpilsPrecSolveFn psolve) noexcept
            int CVodeGetNumJtimesEvals(void *cvode_mem, long int *njvevals) noexcept #Number of jac*vector evals
            int CVodeGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS) noexcept #Number of res evals due to jac*vector evals
            int CVodeGetNumPrecEvals(void *cvode_mem, long int *npevals) noexcept
            int CVodeGetNumPrecSolves(void *cvode_mem, long int *npsolves) noexcept
        ELSE:
            int CVSpilsSetPreconditioner(void *cvode_mem, CVSpilsPrecSetupFn psetup, CVSpilsPrecSolveFn psolve) noexcept
            int CVSpilsGetNumJtimesEvals(void *cvode_mem, long int *njvevals) noexcept #Number of jac*vector evals
            int CVSpilsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS) noexcept #Number of res evals due to jac*vector evals
            int CVSpilsGetNumPrecEvals(void *cvode_mem, long int *npevals) noexcept
            int CVSpilsGetNumPrecSolves(void *cvode_mem, long int *npsolves) noexcept

cdef extern from "idas/idas.h":
    ctypedef int (*IDAResFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data) noexcept
    IF SUNDIALS_VERSION_NR >= 600000:
        void* IDACreate(SUNContext ctx) noexcept
    ELSE:
        void* IDACreate() noexcept
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
    ctypedef int (*IDARootFn)(realtype tt, N_Vector yy, N_Vector yp, realtype *gout, void *user_data) noexcept
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
    IF SUNDIALS_VERSION_NR >= 400000:
        int IDAGetNumJacEvals(void *ida_mem, long int *njevals)          #Number of jac evals
        int IDAGetNumLinResEvals(void *ida_mem, long int *nrevalsLS)        #Number of res evals due to jac evals
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
                                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) noexcept
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

IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "ida/ida_ls.h":
        ctypedef int (*IDALsJacFn)(sunrealtype t, sunrealtype c_j, N_Vector y, N_Vector yp, N_Vector r, SUNMatrix Jac, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) noexcept
        ctypedef int (*IDALsJacTimesSetupFn)(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, sunrealtype c_j, void* user_data) noexcept
        ctypedef int (*IDALsJacTimesVecFn)(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector v, N_Vector Jv, sunrealtype c_j, void* user_data, N_Vector tmp1, N_Vector tmp2) noexcept
ELSE:
    cdef extern from "idas/idas_spils.h":
        ctypedef int (*IDASpilsJacTimesVecFn)(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, 
                N_Vector v, N_Vector Jv, realtype cj, void *user_data,N_Vector tmp1, N_Vector tmp2) noexcept

IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "ida/ida_ls.h":
        int IDASetJacFn(void* ida_mem, IDALsJacFn jac) noexcept
        int IDASetLinearSolver(void* ida_mem, SUNLinearSolver LS, SUNMatrix A) noexcept
        int IDASetJacTimes(void* ida_mem, IDALsJacTimesSetupFn jtsetup, IDALsJacTimesVecFn jtimes) noexcept

    cdef inline int ida_spils_jtsetup_dummy(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, realtype c_j, void *user_data) noexcept: return 0
ELIF SUNDIALS_VERSION_NR >= 300000:
    cdef extern from "idas/idas_direct.h":
        ctypedef int (*IDADlsDenseJacFn)(realtype tt, realtype cj, N_Vector yy, 
                       N_Vector yp, N_Vector rr, SUNMatrix Jac, void *user_data, 
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) noexcept
        IF SUNDIALS_VERSION_NR >= 400000:
            int IDASetJacFn(void *ida_mem, IDADlsDenseJacFn djac)
            int IDASetLinearSolver(void *ida_mem, SUNLinearSolver LS, SUNMatrix A)
        ELSE:
            int IDADlsSetJacFn(void *ida_mem, IDADlsDenseJacFn djac)
            int IDADlsSetLinearSolver(void *ida_mem, SUNLinearSolver LS, SUNMatrix A)
    
    cdef extern from "idas/idas_spils.h":
        ctypedef int (*IDASpilsJacTimesSetupFn)(realtype tt, N_Vector yy,
                      N_Vector yp, N_Vector rr, realtype c_j, void *user_data) noexcept
        IF SUNDIALS_VERSION_NR >= 400000:
            int IDASetJacTimes(void *ida_mem,
                IDASpilsJacTimesSetupFn jtsetup, IDASpilsJacTimesVecFn jtimes)
        ELSE:
            int IDASpilsSetLinearSolver(void *ida_mem, SUNLinearSolver LS)
            int IDASpilsSetJacTimes(void *ida_mem,
                IDASpilsJacTimesSetupFn jtsetup, IDASpilsJacTimesVecFn jtimes)


                
    cdef inline int ida_spils_jtsetup_dummy(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr, realtype c_j, void *user_data) noexcept: return 0
ELSE:
    cdef extern from "idas/idas_dense.h":
        int IDADense(void *ida_mem, long int n)
        ctypedef int (*IDADlsDenseJacFn)(long int Neq, realtype tt, realtype cj, N_Vector yy, 
                       N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, 
                       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) noexcept
        int IDADlsSetDenseJacFn(void *ida_mem, IDADlsDenseJacFn djac)
    
    cdef extern from "idas/idas_spgmr.h":
        int IDASpgmr(void *ida_mem, int max1)
        
    cdef extern from "idas/idas_spils.h":
        int IDASpilsSetJacTimesVecFn(void *ida_mem, IDASpilsJacTimesVecFn ida_jacv)

IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "ida/ida_ls.h": 
        int IDAGetNumJtimesEvals(void *ida_mem, long int *njvevals) #Number of jac*vector
        int IDAGetNumResEvals(void *ida_mem, long int *nfevalsLS) #Number of rhs due to jac*vector
ELSE:
    cdef extern from "idas/idas_spils.h":
        IF SUNDIALS_VERSION_NR >= 400000:
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
    ctypedef int (*KINSysFn)(N_Vector uu, N_Vector fval, void *user_data) noexcept
    ctypedef void (*KINInfoHandlerFn)(const char *module, const char *function, char *msg, void *user_data) noexcept
    # initialization routines
    IF SUNDIALS_VERSION_NR >= 600000:
        void *KINCreate(SUNContext ctx) noexcept
    ELSE:
        void *KINCreate() noexcept
    int KINInit(void *kinmem, KINSysFn func, N_Vector tmpl) noexcept

    # optional input spec. functions,
    # for specificationsdocumentation cf. kinsol.h line 218-449
    IF SUNDIALS_VERSION_NR < 700000:
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
IF SUNDIALS_VERSION_NR < 700000:
    cdef extern from "kinsol/kinsol.h":
        ctypedef void (*KINErrHandlerFn)(int error_code, char *module, char *function, char *msg, void *user_data) noexcept
        int KINSetErrHandlerFn(void *kinmem, KINErrHandlerFn ehfun, void *eh_data) noexcept
    cdef extern from "cvodes/cvodes.h":
        ctypedef void (*CVErrHandlerFn)(int error_code, const char *module, const char *function, char *msg, void *eh_data) noexcept
        int CVodeSetErrHandlerFn(void *cvode_mem, CVErrHandlerFn ehfun, void* eh_data) noexcept
    cdef extern from "idas/idas.h":
        ctypedef void (*IDAErrHandlerFn)(int error_code, const char *module, const char *function, char *msg, void *eh_data) noexcept
        int IDASetErrHandlerFn(void *ida_mem,IDAErrHandlerFn ehfun, void* eh_data) noexcept

IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "kinsol/kinsol_ls.h":
        ctypedef int (*KINLsJacFn)(N_Vector u, N_Vector fu, SUNMatrix J, void* user_data, N_Vector tmp1, N_Vector tmp2) noexcept
        int KINSetLinearSolver(void* kinmem, SUNLinearSolver LS, SUNMatrix A) noexcept
        int KINSetJacFn(void* kinmem, KINLsJacFn jac) noexcept
ELIF SUNDIALS_VERSION_NR >= 300000:
    cdef extern from "kinsol/kinsol_direct.h":
        ctypedef int (*KINDlsDenseJacFn)(N_Vector u, N_Vector fu, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2) noexcept
        IF SUNDIALS_VERSION_NR < 400000:
            int KINDlsSetLinearSolver(void *kinmem, SUNLinearSolver LS, SUNMatrix A)
            int KINDlsSetJacFn(void *kinmem, KINDlsDenseJacFn djac)
        ELSE:
            int KINSetLinearSolver(void *kinmem, SUNLinearSolver LS, SUNMatrix A)
            int KINSetJacFn(void *kinmem, KINDlsDenseJacFn djac)
    cdef extern from "kinsol/kinsol_spils.h":
        int KINSpilsSetLinearSolver(void *kinsol_mem, SUNLinearSolver LS)
        
        ctypedef int (*KINSpilsPrecSolveFn)(N_Vector u, N_Vector uscale,
                    N_Vector fval, N_Vector fscale, N_Vector v, void *problem_data) noexcept
        ctypedef int (*KINSpilsPrecSetupFn)(N_Vector u, N_Vector uscale,
                    N_Vector fval, N_Vector fscale, void *problem_data) noexcept
ELSE:
    # functions used for supplying jacobian, and receiving info from linear solver
    cdef extern from "kinsol/kinsol_direct.h":
        # user functions
        ctypedef int (*KINDlsDenseJacFn)(long int dim, N_Vector u, N_Vector fu, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2) noexcept
        
        # function used to link user functions to KINSOL
        int KINDlsSetDenseJacFn(void *kinmem, KINDlsDenseJacFn jac)
    
    cdef extern from "kinsol/kinsol_dense.h":
        int KINDense(void *kinmem, int dim)
    
    cdef extern from "kinsol/kinsol_spgmr.h":
        int KINSpgmr(void *kinmem, int maxl)
        
    cdef extern from "kinsol/kinsol_spils.h":
        ctypedef int (*KINSpilsPrecSolveFn)(N_Vector u, N_Vector uscale,
                    N_Vector fval, N_Vector fscale, N_Vector v, void *problem_data, N_Vector tmp) noexcept
        ctypedef int (*KINSpilsPrecSetupFn)(N_Vector u, N_Vector uscale,
                    N_Vector fval, N_Vector fscale, void *problem_data, N_Vector tmp1, N_Vector tmp2) noexcept

IF SUNDIALS_VERSION_NR >= 600000:
    cdef extern from "kinsol/kinsol_ls.h":
        ctypedef int (*KINLsPrecSetupFn)(N_Vector uu, N_Vector uscale, N_Vector fval, N_Vector fscale, void* user_data) noexcept
        ctypedef int (*KINLsPrecSolveFn)(N_Vector uu, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector vv, void* user_data) noexcept
        ctypedef int (*KINLsJacTimesVecFn)(N_Vector v, N_Vector Jv, N_Vector uu, sunbooleantype* new_uu, void* J_data) noexcept
        int KINGetLastLinFlag(void *kinmem, long int *flag) noexcept
        int KINGetNumJacEvals(void *kinmem, long int *njevalsB) noexcept
        int KINGetNumFuncEvals(void *kinmem, long int *nfevalsB) noexcept
        int KINSetJacTimesVecFn(void* kinmem, KINLsJacTimesVecFn jtv) noexcept
        int KINSetPreconditioner(void* kinmem, KINLsPrecSetupFn psetup, KINLsPrecSolveFn psolve) noexcept
        int KINGetNumLinIters(void *kinmem, long int *nliters) noexcept
        int KINGetNumLinConvFails(void *kinmem, long int *nlcfails) noexcept
        int KINGetNumPrecEvals(void *kinmem, long int *npevals) noexcept
        int KINGetNumPrecSolves(void *kinmem, long int *npsolves) noexcept
        int KINGetNumJtimesEvals(void *kinmem, long int *njevals) noexcept
        int KINGetNumFuncEvals(void *kinmem, long int *nfevalsLS) noexcept
ELSE:
    cdef extern from "kinsol/kinsol_direct.h":
        # optional output fcts for linear direct solver
        int KINDlsGetWorkSpace(void *kinmem, long int *lenrwB, long int *leniwB)
        IF SUNDIALS_VERSION_NR >= 400000:
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
                    void *problem_data) noexcept
        IF SUNDIALS_VERSION_NR >= 400000:
            int KINSetJacTimesVecFn(void *kinmem, KINSpilsJacTimesVecFn jacv) noexcept
            int KINSetPreconditioner(void *kinmem, KINSpilsPrecSetupFn psetup, KINSpilsPrecSolveFn psolve) noexcept
            int KINGetNumLinIters(void *kinmem, long int *nliters) noexcept
            int KINGetNumLinConvFails(void *kinmem, long int *nlcfails) noexcept
            int KINGetNumPrecEvals(void *kinmem, long int *npevals) noexcept
            int KINGetNumPrecSolves(void *kinmem, long int *npsolves) noexcept
            int KINGetNumJtimesEvals(void *kinmem, long int *njevals) noexcept
            int KINGetNumFuncEvals(void *kinmem, long int *nfevalsLS) noexcept
        ELSE:
            int KINSpilsSetJacTimesVecFn(void *kinmem, KINSpilsJacTimesVecFn jacv) noexcept
            int KINSpilsSetPreconditioner(void *kinmem, KINSpilsPrecSetupFn psetup, KINSpilsPrecSolveFn psolve) noexcept
            int KINSpilsGetNumLinIters(void *kinmem, long int *nliters) noexcept
            int KINSpilsGetNumConvFails(void *kinmem, long int *nlcfails) noexcept
            int KINSpilsGetNumPrecEvals(void *kinmem, long int *npevals) noexcept
            int KINSpilsGetNumPrecSolves(void *kinmem, long int *npsolves) noexcept
            int KINSpilsGetNumJtimesEvals(void *kinmem, long int *njevals) noexcept
            int KINSpilsGetNumFuncEvals(void *kinmem, long int *nfevalsLS) noexcept

#=========================
# END SUNDIALS DEFINITIONS
#=========================
