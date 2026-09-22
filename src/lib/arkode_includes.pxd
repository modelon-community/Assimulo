#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2026 Modelon AB
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
Declarations of the SUNDIALS ARKODE API used by assimulo.solvers.arkode.

Only the stepper-independent ARKode* functions of SUNDIALS >= 7.1 are declared
(the ARKStep*/ERKStep* calls are being folded into them upstream); the
stepper-specific creation and Butcher-table selection calls are the exception.
The vector, matrix, linear-solver and context types come from sundials_includes.
"""

from sundials_includes cimport realtype, N_Vector, SUNMatrix, SUNLinearSolver, SUNContext, sunindextype

cdef extern from "sundials/sundials_context.h":
    int SUNContext_Free(SUNContext* ctx) noexcept

# a custom SUNLinearSolver (the block-diagonal direct solver of ARKODE.linear_solver = "BLOCK")
cdef extern from "sundials/sundials_linearsolver.h":
    SUNLinearSolver SUNLinSolNewEmpty(SUNContext ctx) noexcept
    void SUNLinSolFreeEmpty(SUNLinearSolver S) noexcept

# SUNDIALS' own dense LU (what SUNLinSol_Dense uses), on column-pointer matrices
cdef extern from "sundials/sundials_direct.h":
    realtype** SUNDlsMat_newDenseMat(sunindextype m, sunindextype n) noexcept
    void SUNDlsMat_destroyMat(realtype** a) noexcept
    sunindextype* SUNDlsMat_newIndexArray(sunindextype n) noexcept
    void SUNDlsMat_destroyArray(void* v) noexcept

cdef extern from "sundials/sundials_dense.h":
    sunindextype SUNDlsMat_denseGETRF(realtype** a, sunindextype m, sunindextype n, sunindextype* p) noexcept
    void SUNDlsMat_denseGETRS(realtype** a, sunindextype n, sunindextype* p, realtype* b) noexcept

cdef extern from "arkode/arkode.h":
    # itask
    enum: ARK_NORMAL
    enum: ARK_ONE_STEP
    # return flags
    enum: ARK_SUCCESS
    enum: ARK_TSTOP_RETURN
    enum: ARK_ROOT_RETURN
    enum: ARK_WARNING
    enum: ARK_TOO_MUCH_WORK
    enum: ARK_TOO_MUCH_ACC
    enum: ARK_ERR_FAILURE
    enum: ARK_CONV_FAILURE
    enum: ARK_LINIT_FAIL
    enum: ARK_LSETUP_FAIL
    enum: ARK_LSOLVE_FAIL
    enum: ARK_RHSFUNC_FAIL
    enum: ARK_FIRST_RHSFUNC_ERR
    enum: ARK_REPTD_RHSFUNC_ERR
    enum: ARK_UNREC_RHSFUNC_ERR
    enum: ARK_RTFUNC_FAIL
    enum: ARK_MEM_FAIL
    enum: ARK_MEM_NULL
    enum: ARK_ILL_INPUT
    enum: ARK_TOO_CLOSE
    enum: ARK_INTERP_FAIL
    enum: ARK_STEPPER_UNSUPPORTED
    # interpolation modules
    enum: ARK_INTERP_HERMITE
    enum: ARK_INTERP_LAGRANGE

    ctypedef int (*ARKRhsFn)(realtype t, N_Vector y, N_Vector ydot, void* user_data) noexcept
    ctypedef int (*ARKRootFn)(realtype t, N_Vector y, realtype* gout, void* user_data) noexcept

    void ARKodeFree(void** arkode_mem) noexcept
    int ARKodeReset(void* arkode_mem, realtype tR, N_Vector yR) noexcept
    int ARKodeSStolerances(void* arkode_mem, realtype reltol, realtype abstol) noexcept
    int ARKodeSVtolerances(void* arkode_mem, realtype reltol, N_Vector abstol) noexcept
    int ARKodeRootInit(void* arkode_mem, int nrtfn, ARKRootFn g) noexcept
    int ARKodeSetRootDirection(void* arkode_mem, int* rootdir) noexcept
    int ARKodeSetNoInactiveRootWarn(void* arkode_mem) noexcept
    int ARKodeSetOrder(void* arkode_mem, int maxord) noexcept
    int ARKodeSetInterpolantType(void* arkode_mem, int itype) noexcept
    int ARKodeSetInterpolantDegree(void* arkode_mem, int degree) noexcept
    int ARKodeSetMaxNumSteps(void* arkode_mem, long int mxsteps) noexcept
    int ARKodeSetInterpolateStopTime(void* arkode_mem, int interp) noexcept
    int ARKodeSetStopTime(void* arkode_mem, realtype tstop) noexcept
    int ARKodeClearStopTime(void* arkode_mem) noexcept
    int ARKodeSetFixedStep(void* arkode_mem, realtype hfixed) noexcept
    int ARKodeSetUseCompensatedSums(void* arkode_mem, int onoff) noexcept
    int ARKodeSetUserData(void* arkode_mem, void* user_data) noexcept
    int ARKodeSetLinear(void* arkode_mem, int timedepend) noexcept
    int ARKodeSetNonlinear(void* arkode_mem) noexcept
    int ARKodeSetDeduceImplicitRhs(void* arkode_mem, int deduce) noexcept
    int ARKodeSetLSetupFrequency(void* arkode_mem, int msbp) noexcept
    int ARKodeSetDeltaGammaMax(void* arkode_mem, realtype dgmax) noexcept
    int ARKodeSetPredictorMethod(void* arkode_mem, int method) noexcept
    int ARKodeSetMaxNonlinIters(void* arkode_mem, int maxcor) noexcept
    int ARKodeSetMaxConvFails(void* arkode_mem, int maxncf) noexcept
    int ARKodeSetNonlinConvCoef(void* arkode_mem, realtype nlscoef) noexcept
    int ARKodeSetAdaptControllerByName(void* arkode_mem, const char* cname) noexcept
    int ARKodeSetSafetyFactor(void* arkode_mem, realtype safety) noexcept
    int ARKodeSetMaxGrowth(void* arkode_mem, realtype mx_growth) noexcept
    int ARKodeSetMaxFirstGrowth(void* arkode_mem, realtype etamx1) noexcept
    int ARKodeSetMaxEFailGrowth(void* arkode_mem, realtype etamxf) noexcept
    int ARKodeSetMaxCFailGrowth(void* arkode_mem, realtype etacf) noexcept
    int ARKodeSetMaxErrTestFails(void* arkode_mem, int maxnef) noexcept
    int ARKodeSetMaxHnilWarns(void* arkode_mem, int mxhnil) noexcept
    int ARKodeSetInitStep(void* arkode_mem, realtype hin) noexcept
    int ARKodeSetMinStep(void* arkode_mem, realtype hmin) noexcept
    int ARKodeSetMaxStep(void* arkode_mem, realtype hmax) noexcept
    int ARKodeEvolve(void* arkode_mem, realtype tout, N_Vector yout, realtype* tret, int itask) noexcept
    int ARKodeGetDky(void* arkode_mem, realtype t, int k, N_Vector dky) noexcept
    int ARKodeGetNumRhsEvals(void* arkode_mem, int partition_index, long int* num_rhs_evals) noexcept
    int ARKodeGetNumStepAttempts(void* arkode_mem, long int* step_attempts) noexcept
    int ARKodeGetNumSteps(void* arkode_mem, long int* nsteps) noexcept
    int ARKodeGetLastStep(void* arkode_mem, realtype* hlast) noexcept
    int ARKodeGetCurrentStep(void* arkode_mem, realtype* hcur) noexcept
    int ARKodeGetCurrentTime(void* arkode_mem, realtype* tcur) noexcept
    int ARKodeGetNumGEvals(void* arkode_mem, long int* ngevals) noexcept
    int ARKodeGetRootInfo(void* arkode_mem, int* rootsfound) noexcept
    int ARKodeGetNumErrTestFails(void* arkode_mem, long int* netfails) noexcept
    int ARKodeGetNumLinSolvSetups(void* arkode_mem, long int* nlinsetups) noexcept
    int ARKodeGetNumNonlinSolvIters(void* arkode_mem, long int* nniters) noexcept
    int ARKodeGetNumNonlinSolvConvFails(void* arkode_mem, long int* nnfails) noexcept
    int ARKodeGetNumStepSolveFails(void* arkode_mem, long int* nncfails) noexcept
    int ARKodeGetNumJacEvals(void* arkode_mem, long int* njevals) noexcept
    int ARKodeGetNumLinIters(void* arkode_mem, long int* nliters) noexcept
    int ARKodeGetNumJtimesEvals(void* arkode_mem, long int* njvevals) noexcept
    int ARKodeGetNumLinRhsEvals(void* arkode_mem, long int* nfevalsLS) noexcept
    int ARKodeGetEstLocalErrors(void* arkode_mem, N_Vector ele) noexcept
    int ARKodeGetErrWeights(void* arkode_mem, N_Vector eweight) noexcept

cdef extern from "arkode/arkode_ls.h":
    ctypedef int (*ARKLsJacFn)(realtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) noexcept
    ctypedef int (*ARKLsJacTimesSetupFn)(realtype t, N_Vector y, N_Vector fy, void* user_data) noexcept
    ctypedef int (*ARKLsJacTimesVecFn)(N_Vector v, N_Vector Jv, realtype t, N_Vector y, N_Vector fy,
                                       void* user_data, N_Vector tmp) noexcept
    int ARKodeSetLinearSolver(void* arkode_mem, SUNLinearSolver LS, SUNMatrix A) noexcept
    int ARKodeSetJacEvalFrequency(void* arkode_mem, long int msbj) noexcept
    int ARKodeSetJacFn(void* arkode_mem, ARKLsJacFn jac) noexcept
    int ARKodeSetJacTimes(void* arkode_mem, ARKLsJacTimesSetupFn jtsetup, ARKLsJacTimesVecFn jtimes) noexcept

cdef extern from "arkode/arkode_arkstep.h":
    void* ARKStepCreate(ARKRhsFn fe, ARKRhsFn fi, realtype t0, N_Vector y0, SUNContext sunctx) noexcept
    int ARKStepReInit(void* arkode_mem, ARKRhsFn fe, ARKRhsFn fi, realtype t0, N_Vector y0) noexcept
    int ARKStepSetTableName(void* arkode_mem, const char* itable, const char* etable) noexcept

cdef extern from "arkode/arkode_sprkstep.h":
    void* SPRKStepCreate(ARKRhsFn f1, ARKRhsFn f2, realtype t0, N_Vector y0, SUNContext sunctx) noexcept
    int SPRKStepSetMethodName(void* arkode_mem, const char* method) noexcept

cdef extern from "arkode/arkode_erkstep.h":
    void* ERKStepCreate(ARKRhsFn f, realtype t0, N_Vector y0, SUNContext sunctx) noexcept
    int ERKStepReInit(void* arkode_mem, ARKRhsFn f, realtype t0, N_Vector y0) noexcept
    int ERKStepSetTableName(void* arkode_mem, const char* etable) noexcept

cdef extern from "arkode/arkode_butcher.h":
    cdef struct ARKodeButcherTableMem:
        int q
        int p
        int stages
    ctypedef ARKodeButcherTableMem* ARKodeButcherTable
    void ARKodeButcherTable_Free(ARKodeButcherTable B) noexcept

cdef extern from "arkode/arkode_butcher_dirk.h":
    ARKodeButcherTable ARKodeButcherTable_LoadDIRKByName(const char* imethod) noexcept

cdef extern from "arkode/arkode_butcher_erk.h":
    ARKodeButcherTable ARKodeButcherTable_LoadERKByName(const char* emethod) noexcept
