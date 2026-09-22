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
ARKODE (SUNDIALS) as an Assimulo solver: adaptive one-step Runge-Kutta methods,
implicit (DIRK/ESDIRK), explicit, additive implicit-explicit (IMEX) or symplectic,
with ARKODE's rootfinding for state events.
Requires SUNDIALS >= 7.1 (the stepper-independent ARKode* API).
"""

import numpy as np
cimport numpy as np

from assimulo.exception import AssimuloException
from assimulo.explicit_ode cimport Explicit_ODE
from assimulo.support import set_type_shape_array

cimport sundials_includes as SUNDIALS
cimport arkode_includes as ARK
from arkode_includes cimport (SUNLinSolNewEmpty, SUNLinSolFreeEmpty, SUNDlsMat_newDenseMat, SUNDlsMat_destroyMat,
                              SUNDlsMat_newIndexArray, SUNDlsMat_destroyArray, SUNDlsMat_denseGETRF,
                              SUNDlsMat_denseGETRS)
from arkode_includes cimport (ARK_NORMAL, ARK_ONE_STEP, ARK_SUCCESS, ARK_TSTOP_RETURN, ARK_ROOT_RETURN,
                              ARK_TOO_MUCH_WORK, ARK_TOO_MUCH_ACC, ARK_ERR_FAILURE, ARK_CONV_FAILURE,
                              ARK_LINIT_FAIL, ARK_LSETUP_FAIL, ARK_LSOLVE_FAIL, ARK_RHSFUNC_FAIL,
                              ARK_FIRST_RHSFUNC_ERR, ARK_REPTD_RHSFUNC_ERR, ARK_UNREC_RHSFUNC_ERR,
                              ARK_RTFUNC_FAIL, ARK_MEM_FAIL, ARK_MEM_NULL, ARK_ILL_INPUT, ARK_TOO_CLOSE,
                              ARK_INTERP_FAIL, ARK_STEPPER_UNSUPPORTED)

from sundials_includes cimport N_Vector, realtype, N_VectorContent_Serial, DENSE_COL, sunindextype
from sundials_includes cimport memcpy, N_VNew_Serial, DlsMat, SUNMatrix, SUNMatrixContent_Dense, SUNMatrixContent_Sparse
from sundials_includes cimport malloc, free, N_VConst_Serial
from sundials_includes cimport N_VCloneVectorArray, N_VDestroy
from sundials_includes cimport SUNLinearSolver, SUNMatrix, N_VectorContent_Serial

include "constants.pxi"                          # Assimulo's ID_* flags
include "../lib/sundials_constants.pxi"
include "../lib/sundials_callbacks.pxi"          # arr2nv, nv2arr
include "../lib/sundials_callbacks_ida_cvode.pxi" # cv_rhs, cv_jac, cv_jacv, cv_root, cv_err, ProblemData


cdef inline int ark_rhs_part(object fn, realtype t, N_Vector yv, N_Vector yvdot, ProblemData pData) noexcept:
    """One part of an additive rhs (`fn` is problem.rhs_implicit or problem.rhs_explicit),
    with cv_rhs's calling convention and failure codes."""
    cdef np.ndarray y = pData.work_y
    cdef realtype* resptr = (<N_VectorContent_Serial>yvdot.content).data
    cdef int i
    nv2arr_inplace(yv, y)
    try:
        if pData.sw != NULL:
            rhs = fn(t, y, <list>pData.sw)
        else:
            rhs = fn(t, y)
    except Exception:
        pData.nrhsfails += 1
        return CV_REC_ERR
    except BaseException:
        return CV_UNREC_RHSFUNC_ERR
    for i in range(pData.dim):
        resptr[i] = rhs[i]
    return CV_SUCCESS

cdef int ark_rhs_implicit(realtype t, N_Vector yv, N_Vector yvdot, void* problem_data) noexcept:
    """ARKStep's fi: problem.rhs_implicit (imex)."""
    cdef ProblemData pData = <ProblemData>problem_data
    return ark_rhs_part(<object>pData.RHS_I, t, yv, yvdot, pData)

cdef int ark_rhs_explicit(realtype t, N_Vector yv, N_Vector yvdot, void* problem_data) noexcept:
    """ARKStep's fe: problem.rhs_explicit (imex)."""
    cdef ProblemData pData = <ProblemData>problem_data
    return ark_rhs_part(<object>pData.RHS_E, t, yv, yvdot, pData)


# ---------------------------------------------------------------------------------------------
# The block-diagonal direct linear solver (ARKODE.linear_solver = "BLOCK")
#
# ARKLS forms the Newton matrix A = I - gamma*J in a dense SUNMatrix and calls setup/solve on
# it. When the implicit part's states fall into blocks whose derivatives do not depend on each
# other's states (several fast subsystems that only talk to each other through the explicit
# part: a vehicle body with its suspensions), the rows of A are, after ordering,
#     identity rows for the states outside the blocks (their implicit rhs is zero),
#     one dense block A_BB per block, plus the columns A_BS coupling it to the outside states.
# The solve is then x_S = b_S and, per block, A_BB x_B = b_B - A_BS x_S: the dense LU's result
# to rounding, at sum(n_b^3) instead of n^3. The structure is the problem's claim
# (problem.implicit_blocks); every setup checks it against the matrix it is given.

cdef class BlockLU:
    """Content of the BLOCK SUNLinearSolver: the block index sets and their LU factors."""
    cdef int n, nblocks, nout, check
    cdef double check_tol
    cdef sunindextype* bstart      # block b holds bidx[bstart[b] : bstart[b+1]]
    cdef sunindextype* bidx        # the state indices of all blocks, block by block
    cdef sunindextype* sidx        # the state indices outside every block (the "outside" set S)
    cdef sunindextype* block_of    # state -> block number, -1 outside
    cdef realtype*** lu            # per block: SUNDlsMat_newDenseMat(n_b, n_b), LU-factored in setup
    cdef sunindextype** piv        # per block: the pivots
    cdef realtype* work            # max block size
    cdef public object blocks      # the index arrays as given (validated), for comparisons
    cdef public object message     # the structure check's complaint, None if none
    cdef public long int last_flag
    cdef public long int nsetups, nsolves

    def __cinit__(self):
        self.bstart = NULL; self.bidx = NULL; self.sidx = NULL; self.block_of = NULL
        self.lu = NULL; self.piv = NULL; self.work = NULL
        self.nblocks = 0

    def __init__(self, blocks, int n, int cover_all, int check, double check_tol):
        cdef int b, k, i, nb, maxnb = 0
        used = np.zeros(n, dtype=bool)
        arrays = []
        for b, blk in enumerate(blocks):
            a = np.asarray(blk)
            if a.dtype == bool:
                if len(a) != n:
                    raise AssimuloException("ARKODE BLOCK: a boolean block mask must have one entry per state (block %d)." % b)
                a = np.flatnonzero(a)
            a = np.unique(a.astype(np.int64))
            if len(a) == 0:
                raise AssimuloException("ARKODE BLOCK: block %d is empty." % b)
            if a[0] < 0 or a[-1] >= n:
                raise AssimuloException("ARKODE BLOCK: block %d has an index outside 0..%d." % (b, n - 1))
            if used[a].any():
                raise AssimuloException("ARKODE BLOCK: block %d overlaps an earlier block." % b)
            used[a] = True
            arrays.append(a)
        if not arrays:
            raise AssimuloException("ARKODE BLOCK: 'implicit_blocks' is empty.")
        if cover_all and not used.all():
            raise AssimuloException("ARKODE BLOCK: with method = 'implicit' the blocks must cover every state "
                                    "(%d of %d are outside); states outside the blocks are only allowed "
                                    "for method = 'imex', where their implicit rhs is zero." % ((~used).sum(), n))
        self.blocks = arrays
        self.n = n
        self.nblocks = len(arrays)
        self.check = check
        self.check_tol = check_tol
        self.message = None
        self.last_flag = 0
        self.nsetups = 0
        self.nsolves = 0
        self.bstart = SUNDlsMat_newIndexArray(self.nblocks + 1)
        self.bidx = SUNDlsMat_newIndexArray(max(1, int(used.sum())))
        outside = np.flatnonzero(~used)
        self.nout = len(outside)
        self.sidx = SUNDlsMat_newIndexArray(max(1, self.nout))
        self.block_of = SUNDlsMat_newIndexArray(n)
        self.lu = <realtype***>malloc(self.nblocks * sizeof(realtype**))
        self.piv = <sunindextype**>malloc(self.nblocks * sizeof(sunindextype*))
        for i in range(n):
            self.block_of[i] = -1
        for i in range(self.nout):
            self.sidx[i] = outside[i]
        k = 0
        for b in range(self.nblocks):
            self.bstart[b] = k
            nb = len(arrays[b])
            for i in range(nb):
                self.bidx[k] = arrays[b][i]
                self.block_of[arrays[b][i]] = b
                k += 1
            self.lu[b] = SUNDlsMat_newDenseMat(nb, nb)
            self.piv[b] = SUNDlsMat_newIndexArray(nb)
            if nb > maxnb: maxnb = nb
        self.bstart[self.nblocks] = k
        self.work = <realtype*>malloc(maxnb * sizeof(realtype))

    def __dealloc__(self):
        cdef int b
        if self.lu != NULL:
            for b in range(self.nblocks):
                if self.lu[b] != NULL: SUNDlsMat_destroyMat(self.lu[b])
            free(self.lu)
        if self.piv != NULL:
            for b in range(self.nblocks):
                if self.piv[b] != NULL: SUNDlsMat_destroyArray(self.piv[b])
            free(self.piv)
        if self.bstart != NULL: SUNDlsMat_destroyArray(self.bstart)
        if self.bidx != NULL: SUNDlsMat_destroyArray(self.bidx)
        if self.sidx != NULL: SUNDlsMat_destroyArray(self.sidx)
        if self.block_of != NULL: SUNDlsMat_destroyArray(self.block_of)
        if self.work != NULL: free(self.work)

    def sizes(self):
        return [len(a) for a in self.blocks]

    cdef int setup(self, SUNMatrix A) noexcept:
        """Copies each diagonal block out of the dense A = I - gamma*J and LU-factors it; with
        `check`, first verifies that A has the declared structure. Returns 0, +1 for a singular
        block (recoverable: ARKODE retries with a smaller step), -1 for a structure violation."""
        cdef SUNDIALS.SUNMatrixContent_Dense c = <SUNDIALS.SUNMatrixContent_Dense>A.content
        cdef realtype* d = c.data
        cdef sunindextype ld = c.M
        cdef int b, i, j, nb, k0, ii, jj
        cdef sunindextype ier
        cdef double tol, aii, dev, worst = 0.0
        cdef int wi = -1, wj = -1
        self.nsetups += 1
        if self.check:
            # rows outside the blocks are identity rows; a block's row is zero in the columns of
            # the other blocks (its columns in the outside set S may be anything: A_BS)
            for i in range(self.n):
                b = self.block_of[i]
                aii = d[i * ld + i]
                tol = self.check_tol * (aii if aii > 1.0 else (-aii if aii < -1.0 else 1.0))
                if b < 0:
                    dev = aii - 1.0 if aii > 1.0 else 1.0 - aii
                    if dev > tol and dev > worst:
                        worst = dev; wi = i; wj = i
                for j in range(self.n):
                    if j == i or not (b < 0 or (self.block_of[j] >= 0 and self.block_of[j] != b)):
                        continue
                    dev = d[j * ld + i]
                    if dev < 0: dev = -dev
                    if dev > tol and dev > worst:
                        worst = dev; wi = i; wj = j
            if wi >= 0:
                if self.block_of[wi] < 0:
                    self.message = ("the Newton matrix I - gamma*J has the entry (%d, %d) = %g in the row of "
                                    "state %d, which is outside every block and must therefore have a zero "
                                    "implicit rhs" % (wi, wj, d[wj * ld + wi], wi))
                else:
                    self.message = ("the Newton matrix I - gamma*J has the entry (%d, %d) = %g coupling block %d "
                                    "(state %d) to block %d (state %d), which 'implicit_blocks' declares "
                                    "independent" % (wi, wj, d[wj * ld + wi], self.block_of[wi], wi,
                                                     self.block_of[wj], wj))
                self.last_flag = -1
                return -1
        for b in range(self.nblocks):
            k0 = self.bstart[b]
            nb = self.bstart[b + 1] - k0
            for jj in range(nb):
                j = self.bidx[k0 + jj]
                for ii in range(nb):
                    self.lu[b][jj][ii] = d[j * ld + self.bidx[k0 + ii]]
            ier = SUNDlsMat_denseGETRF(self.lu[b], nb, nb, self.piv[b])
            if ier > 0:
                self.last_flag = ier
                return 1
        self.last_flag = 0
        return 0

    cdef int solve(self, SUNMatrix A, N_Vector x, N_Vector b) noexcept:
        cdef SUNDIALS.SUNMatrixContent_Dense c = <SUNDIALS.SUNMatrixContent_Dense>A.content
        cdef realtype* d = c.data
        cdef sunindextype ld = c.M
        cdef realtype* xd = (<N_VectorContent_Serial>x.content).data
        cdef realtype* bd = (<N_VectorContent_Serial>b.content).data
        cdef int blk, i, k, s, nb, k0
        cdef double acc
        self.nsolves += 1
        for k in range(self.nout):
            i = self.sidx[k]
            xd[i] = bd[i]
        for blk in range(self.nblocks):
            k0 = self.bstart[blk]
            nb = self.bstart[blk + 1] - k0
            for k in range(nb):
                i = self.bidx[k0 + k]
                acc = bd[i]
                for s in range(self.nout):
                    acc -= d[self.sidx[s] * ld + i] * xd[self.sidx[s]]
                self.work[k] = acc
            SUNDlsMat_denseGETRS(self.lu[blk], nb, self.piv[blk], self.work)
            for k in range(nb):
                xd[self.bidx[k0 + k]] = self.work[k]
        self.last_flag = 0
        return 0


cdef SUNDIALS.SUNLinearSolver_Type block_ls_gettype(SUNLinearSolver S) noexcept:
    return SUNDIALS.SUNLINEARSOLVER_DIRECT

cdef int block_ls_setup(SUNLinearSolver S, SUNMatrix A) noexcept:
    return (<BlockLU>S.content).setup(A)

cdef int block_ls_solve(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol) noexcept:
    return (<BlockLU>S.content).solve(A, x, b)

cdef sunindextype block_ls_lastflag(SUNLinearSolver S) noexcept:
    return (<BlockLU>S.content).last_flag

cdef int block_ls_free(SUNLinearSolver S) noexcept:
    # the content is a Python object owned by the ARKODE instance
    S.content = NULL
    SUNLinSolFreeEmpty(S)
    return 0

cdef SUNLinearSolver block_linear_solver(BlockLU content, SUNDIALS.SUNContext ctx) noexcept:
    cdef SUNLinearSolver S = SUNLinSolNewEmpty(ctx)
    if S == NULL:
        return NULL
    S.ops.gettype = block_ls_gettype
    S.ops.setup = block_ls_setup
    S.ops.solve = block_ls_solve
    S.ops.lastflag = block_ls_lastflag
    S.ops.free = block_ls_free
    S.content = <void*>content
    return S


cdef int sprk_f1(realtype t, N_Vector yv, N_Vector yvdot, void* problem_data) noexcept:
    """SPRKStep's f1: the momentum-state part of the rhs (the position part zeroed)."""
    cdef ProblemData pData = <ProblemData>problem_data
    cdef int flag = cv_rhs(t, yv, yvdot, problem_data)
    cdef np.ndarray[double, ndim=1, mode="c"] mask = pData.sprk_qmask
    cdef realtype* d = (<N_VectorContent_Serial>yvdot.content).data
    cdef int i
    if flag == 0:
        for i in range(pData.dim):
            if mask[i] != 0.0:
                d[i] = 0.0
    return flag

cdef int sprk_f2(realtype t, N_Vector yv, N_Vector yvdot, void* problem_data) noexcept:
    """SPRKStep's f2: the position-state part of the rhs (the momentum part zeroed)."""
    cdef ProblemData pData = <ProblemData>problem_data
    cdef int flag = cv_rhs(t, yv, yvdot, problem_data)
    cdef np.ndarray[double, ndim=1, mode="c"] mask = pData.sprk_qmask
    cdef realtype* d = (<N_VectorContent_Serial>yvdot.content).data
    cdef int i
    if flag == 0:
        for i in range(pData.dim):
            if mask[i] == 0.0:
                d[i] = 0.0
    return flag



cdef class ARKODE(Explicit_ODE):
    r"""
    SUNDIALS ARKODE: adaptive one-step Runge-Kutta methods for

    .. math::

        \dot{y} = f(t,y), \quad y(t_0) = y_0.

    ``method = "implicit"`` (default) runs ARKStep with a diagonally implicit RK
    table (orders 2-5, L-stable ones at every order), Newton iteration and a direct
    dense linear solver with the problem's Jacobian if ``usejac``;
    ``method = "explicit"`` runs an explicit RK table (orders 2-9) without any
    linear algebra; ``method = "imex"`` runs ARKStep's additive Runge-Kutta pairs
    (orders 2-5) on a split rhs ``f = rhs_implicit + rhs_explicit`` that the
    problem provides as two methods with the signature of ``rhs``, each returning
    a full-length vector (zeros outside its part): Newton iteration on the implicit
    part only, and ``problem.jac`` (if ``usejac``) is then the Jacobian of
    ``rhs_implicit``. When the implicit part falls into blocks of states whose
    derivatives do not depend on each other's states (fast subsystems that only
    talk to each other through the explicit part), ``problem.implicit_blocks`` (a
    list of index arrays) with ``linear_solver = "BLOCK"`` factors the Newton
    matrix block by block instead of as one dense matrix; every setup checks the
    matrix against the declared structure (``block_check``) and stops the run
    with a message naming the offending entry if it does not hold.
    ``method = "symplectic"`` runs SPRKStep, a symplectic
    partitioned RK method (orders 1-6, 8, 10) for separable Hamiltonian systems
    at a fixed step ``fixed_h``: the states listed in ``q_states`` are the
    positions (their derivatives are the momenta), the rest the momenta; the
    problem's rhs is evaluated twice per stage and masked. The table is chosen by ``order`` (ARKODE's default table of that
    order) or by name with ``table`` (e.g. ``"ARKODE_TRBDF2_3_3_2"``,
    ``"ARKODE_ESDIRK324L2SA_4_2_3"``, ``"ARKODE_DORMAND_PRINCE_7_4_5"``).

    Being one-step methods, they restart at full order after every event; state
    events are located by ARKODE's rootfinding (``external_event_detection = False``)
    or by Assimulo's locator on the dense output, time events are stop times.
    Only the stepper-independent ARKode* API of SUNDIALS >= 7.1 is used.

    On a rhs with a jump the solution rides on (a switch without an event indicator)
    Newton fails on most stages of a many-stage method whatever the Jacobian, while a
    2-stage method gets through: when more than ``fallback_conv_fail_rate`` of the last
    ``fallback_window`` step attempts failed by nonlinear convergence, the run continues
    with ``fallback_table`` (default TRBDF2) from that point, once per simulation, logged
    at normal verbosity and shown in the statistics; ``fallback_time`` records when.
    """
    cdef void* ark_mem
    cdef ProblemData pData
    cdef N_Vector yTemp, nv_atol
    cdef SUNDIALS.SUNContext ctx
    cdef SUNDIALS.SUNMatrix sun_matrix
    cdef SUNDIALS.SUNLinearSolver sun_linearsolver
    cdef object pt_root, pt_fcn, pt_jac, pt_jacv
    cdef public object event_func
    cdef public np.ndarray g_old
    cdef object _created_method     # the 'method' the memory block was created for
    cdef object pt_rhs_i, pt_rhs_e  # problem.rhs_implicit / rhs_explicit (imex), None if absent
    cdef BlockLU _block_lu           # the BLOCK linear solver's content while it is attached
    cdef dict _last_counters        # ARKODE's cumulative counters at the last store_statistics
    cdef int _fresh_memory          # the memory block was just created: the method-defining options are still to be set
    cdef object _active_table       # the table in use: options['table'] or, after a fallback, options['fallback_table']
    cdef long int _fb_attempts0, _fb_fails0, _fb_rhsfails0   # counters at the start of the current failure-rate window
    cdef int _consecutive_rhsfails   # restarts after a refused stage point since the last accepted step
    cdef long int _rhsfails_seen     # pData.nrhsfails at the last check: a retry needs a new failure
    cdef public object fallback_time  # time of the fallback to the 2-stage table, None if it did not happen
    cdef public object fallback_return_time   # time of the switch back to the user's table after a hard failure of the fallback

    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem)

        self.pData = ProblemData()
        self.ark_mem = NULL
        self.sun_matrix = NULL
        self.sun_linearsolver = NULL
        self.yTemp = NULL
        self.nv_atol = NULL
        self.ctx = NULL
        self._created_method = None
        self._last_counters = {}
        self._fresh_memory = 0
        self._active_table = None
        self._fb_attempts0 = 0
        self._fb_fails0 = 0
        self._fb_rhsfails0 = 0
        self.fallback_time = None
        self.fallback_return_time = None
        SUNDIALS.SUNContext_Create(SUNDIALS.SUN_COMM_NULL, &self.ctx)

        self.set_problem_data()

        # Solver options (names shared with CVode where the meaning is the same)
        self.options["atol"] = np.array([1.0e-6] * self.problem_info["dim"])
        self.options["rtol"] = 1.0e-6
        self.options["maxh"] = 0.0            # 0: no maximum step
        self.options["minh"] = 0.0
        self.options["inith"] = 0.0           # 0: ARKODE estimates the first step
        self.options["maxsteps"] = 10000
        self.options["usejac"] = True if (self.problem_info["jac_fcn"] or self.problem_info["jacv_fcn"]) else False
        self.options["linear_solver"] = "DENSE"   # or SPGMR (Jacobian-vector products), or BLOCK (problem.implicit_blocks)
        self.options["block_check"] = True        # BLOCK: verify the declared structure at every setup
        self.options["block_check_tol"] = 1.0e-10 # BLOCK: an entry above this (relative to max(1, |A_ii|)) is a violation
        self.options["maxkrylov"] = 5
        self.options["external_event_detection"] = False   # ARKODE rootfinding by default
        # ARKODE-specific
        self.options["method"] = "implicit"   # "implicit" (DIRK), "explicit" (ERK), "imex" (ARK pair) or "symplectic" (SPRK)
        self.options["q_states"] = None       # symplectic: indices (or a boolean mask) of the position states
        self.options["fixed_h"] = 0.0         # symplectic: the fixed step size (SPRKStep has no adaptivity)
        self.options["compensated_sums"] = False   # symplectic: Kahan-compensated stage sums
        self.options["order"] = 4             # method order when 'table' is None
        self.options["table"] = None          # Butcher table name, overrides 'order'
        self.options["predictor"] = 2         # ARKODE predictor 0-5; 0 (y_n) needs 2.4x the steps on Van der Pol, 2 (variable-order) halves the Newton failures of 1 on the data-center FMU
        self.options["nonlin_conv_coef"] = 0.1
        self.options["max_nonlin_iters"] = 4  # ARKODE default 3; 4 saves 60 % of the Jacobians on Van der Pol
        self.options["deduce_implicit_rhs"] = True   # stage derivatives from the stage equations
        self.options["interpolant_degree"] = 3       # -1: ARKODE's default for the method (its order); 3: the cubic Hermite interpolant -- the quintic of order 5 costs rhs evaluations at every output point and reports spurious events under rootfinding, orders <= 4 are unchanged
        self.options["maxncf"] = 10           # max convergence failures per step
        self.options["lsetup_frequency"] = 20     # steps between linear-solver setups (ARKODE default 20)
        self.options["jac_eval_frequency"] = 51   # setups between Jacobian evaluations (ARKODE default 51)
        self.options["delta_gamma_max"] = 0.05    # relative change of gamma that forces a new setup (ARKODE default 0.2; 0.05 saves 10-15 % rhs and most error-test failures on the FMUs)
        self.options["maxnef"] = 7            # max error test failures per step
        self.options["restart_h"] = "estimate"   # after an event: "estimate" a new first step or "keep" the last
        # step-size adaptivity bounds (ARKODE defaults: safety 0.96, max_growth 20, max_first_growth 1e4,
        # max_efail_growth 0.3, max_cfail_growth 0.25)
        self.options["safety"] = 0.96
        self.options["max_growth"] = 20.0
        self.options["max_first_growth"] = 10000.0
        self.options["max_efail_growth"] = 0.3
        self.options["max_cfail_growth"] = 0.25
        self.options["report_continuously"] = False
        # Fallback for a discontinuous rhs (a sliding mode, e.g. an all-or-nothing phase change
        # without an event indicator): Newton on the stages of a many-stage method fails on most
        # steps there whatever the Jacobian, while a 2-stage method rides it. When more than
        # 'fallback_conv_fail_rate' of the last 'fallback_window' step attempts failed by nonlinear
        # convergence (steps failed by a recoverable rhs failure not counted: those are the
        # model's, not the method's), the memory is recreated with 'fallback_table' at the current point and the
        # run continues with it (once per simulation, logged; None disables). If the fallback
        # table then fails hard (repeated error-test or convergence failures), the run switches
        # back to the user's table at that point, once.
        self.options["fallback_table"] = "ARKODE_TRBDF2_3_3_2"
        self.options["fallback_conv_fail_rate"] = 0.25
        self.options["fallback_window"] = 50
        # A recoverable rhs failure (the FMU refusing a trial point) at a directly evaluated
        # stage -- every explicit stage of an ERK or ARK pair, the first stage of an ESDIRK --
        # is ARK_UNREC_RHSFUNC_ERR to ARKODE: it retries only the failures inside Newton. As
        # Radau5/TRBDF2, the run is then restarted at the last accepted point with the step
        # reduced by 'fail_factor', up to 'fail_max' consecutive times.
        self.options["fail_factor"] = 0.25
        self.options["fail_max"] = 40

        self.statistics.add_key("nstepattempts", "Number of step attempts")
        self.statistics.add_key("nconvfails", "Number of steps failed by nonlinear convergence")
        self.statistics.add_key("nfcns_implicit", "Number of implicit rhs evaluations (imex: part of nfcns)")
        self.statistics.add_key("nrhsfails", "Number of recoverable rhs failures")

        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = True
        self.supports["state_events"] = True

    def __dealloc__(self):
        if self.yTemp != NULL:
            N_VDestroy(self.yTemp)
        if self.nv_atol != NULL:
            N_VDestroy(self.nv_atol)
        if self.ark_mem != NULL:
            ARK.ARKodeFree(&self.ark_mem)
        if self.sun_matrix != NULL:
            SUNDIALS.SUNMatDestroy(self.sun_matrix)
        if self.sun_linearsolver != NULL:
            SUNDIALS.SUNLinSolFree(self.sun_linearsolver)
        if self.ctx != NULL:
            ARK.SUNContext_Free(&self.ctx)

    # ------------------------------------------------------------------ set-up
    cdef set_problem_data(self):
        self.pt_fcn = self.problem.rhs
        self.pData.RHS = <void*>self.pt_fcn
        self.pt_rhs_i = getattr(self.problem, "rhs_implicit", None)
        self.pt_rhs_e = getattr(self.problem, "rhs_explicit", None)
        if self.pt_rhs_i is not None:
            self.pData.RHS_I = <void*>self.pt_rhs_i
        if self.pt_rhs_e is not None:
            self.pData.RHS_E = <void*>self.pt_rhs_e
        self.pData.dim = self.problem_info["dim"]
        self.pData.memSize = self.pData.dim * sizeof(realtype)
        if self.problem_info["state_events"] is True:
            self.pt_root = self.problem.state_events
            self.pData.ROOT = <void*>self.pt_root
            self.pData.dimRoot = self.problem_info["dimRoot"]
            self.pData.memSizeRoot = self.pData.dimRoot * sizeof(realtype)
        if self.problem_info["jac_fcn"] is True:
            self.pt_jac = self.problem.jac
            self.pData.JAC = <void*>self.pt_jac
            self.pData.memSizeJac = self.pData.dim * self.pData.dim * sizeof(realtype)
        if self.problem_info["jacv_fcn"] is True:
            self.pt_jacv = self.problem.jacv
            self.pData.JACV = <void*>self.pt_jacv
        self.pData.verbose = 2
        self.pData.create_work_arrays()

    cpdef initialize(self):
        self.statistics.reset()
        if self._active_table != self.options["table"] and self.ark_mem != NULL:
            self._free_memory()          # a previous simulation fell back: start again from the user's table
        self._active_table = self.options["table"]
        self.fallback_time = None
        self.fallback_return_time = None
        self._fb_attempts0 = 0
        self._fb_fails0 = 0
        self._fb_rhsfails0 = 0
        self.pData.nrhsfails = 0
        self._consecutive_rhsfails = 0
        self._rhsfails_seen = 0
        self.initialize_arkode()

    cdef _free_memory(self):
        if self.ark_mem != NULL:
            ARK.ARKodeFree(&self.ark_mem)
            self.ark_mem = NULL
        self._last_counters = {}
        self._free_linear_solver()

    cdef _free_linear_solver(self):
        if self.sun_linearsolver != NULL:
            SUNDIALS.SUNLinSolFree(self.sun_linearsolver)
            self.sun_linearsolver = NULL
        self._block_lu = None
        if self.sun_matrix != NULL:
            SUNDIALS.SUNMatDestroy(self.sun_matrix)
            self.sun_matrix = NULL

    cdef initialize_arkode(self):
        """Creates the ARKODE memory on the first call, resets it to (t, y) afterwards."""
        cdef int flag
        method = self.options["method"]
        if method not in ("implicit", "explicit", "imex", "symplectic"):
            raise ARKODEError(ARK_ILL_INPUT, self.t, "'method' must be 'implicit', 'explicit', 'imex' or 'symplectic'")
        if method == "symplectic":
            self.pData.sprk_qmask = self._q_mask()
        if method == "imex":
            # taken from the problem here, not only in __init__: PyFMI sets the partition on
            # the problem after creating the solver
            self.pt_rhs_i = getattr(self.problem, "rhs_implicit", None)
            self.pt_rhs_e = getattr(self.problem, "rhs_explicit", None)
            for name, fn in (("rhs_implicit", self.pt_rhs_i), ("rhs_explicit", self.pt_rhs_e)):
                if fn is None:
                    raise AssimuloException("ARKODE imex: the problem must define '%s(t, y)' (a full-length vector, "
                                            "zero outside its part; rhs_implicit + rhs_explicit = rhs)." % name)
            self.pData.RHS_I = <void*>self.pt_rhs_i
            self.pData.RHS_E = <void*>self.pt_rhs_e
        if method in ("implicit", "imex") and self.options["linear_solver"] == "BLOCK":
            blocks = getattr(self.problem, "implicit_blocks", None)
            if blocks is None:
                raise AssimuloException("ARKODE BLOCK: the problem must define 'implicit_blocks' (a list of index "
                                        "arrays of states whose implicit derivatives do not depend on each "
                                        "other's states).")
            if self._block_lu is not None and not _same_blocks(self._block_lu.blocks, blocks):
                self._free_memory()          # the structure changed: a new linear solver on a new memory block

        if self.yTemp != NULL:
            N_VDestroy(self.yTemp)
        self.yTemp = arr2nv(self.y, <void*>self.ctx)
        self.pData.verbose = 2 if self.verbosity <= NORMAL else 3

        if self.problem_info["switches"]:
            self.pData.sw = <void*>self.sw

        if self.ark_mem != NULL and self._created_method != method:
            self._free_memory()

        if self.ark_mem == NULL:
            if method == "implicit":
                self.ark_mem = ARK.ARKStepCreate(NULL, cv_rhs, self.t, self.yTemp, self.ctx)
            elif method == "explicit":
                self.ark_mem = ARK.ARKStepCreate(cv_rhs, NULL, self.t, self.yTemp, self.ctx)
            elif method == "imex":
                self.ark_mem = ARK.ARKStepCreate(ark_rhs_explicit, ark_rhs_implicit, self.t, self.yTemp, self.ctx)
            else:
                self.ark_mem = ARK.SPRKStepCreate(sprk_f1, sprk_f2, self.t, self.yTemp, self.ctx)
            if self.ark_mem == NULL:
                raise ARKODEError(ARK_MEM_FAIL, self.t)
            self._created_method = method
            self._fresh_memory = 1

            if self.problem_info["state_events"]:
                if self.options["external_event_detection"]:
                    flag = ARK.ARKodeRootInit(self.ark_mem, 0, cv_root)
                else:
                    flag = ARK.ARKodeRootInit(self.ark_mem, self.pData.dimRoot, cv_root)
                if flag < 0:
                    raise ARKODEError(flag, self.t)
            flag = SUNDIALS.SUNContext_PushErrHandler(self.ctx, cv_err, <void*>self.pData)
            if flag < 0:
                raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetUserData(self.ark_mem, <void*>self.pData)
            if flag < 0:
                raise ARKODEError(flag, self.t)
        else:
            # after an event: keep the method and its statistics, move to (t, y)
            flag = ARK.ARKodeReset(self.ark_mem, self.t, self.yTemp)
            if flag < 0:
                raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetUserData(self.ark_mem, <void*>self.pData)
            if flag < 0:
                raise ARKODEError(flag, self.t)

    def _block_message(self):
        """The BLOCK solver's structure complaint, if its last setup refused the matrix."""
        if self._block_lu is not None and self._block_lu.message is not None:
            return "BLOCK: " + self._block_lu.message + "."
        return None

    def _q_mask(self):
        """The position-state mask (1.0 / 0.0 per state) from the 'q_states' option."""
        q = self.options["q_states"]
        if q is None:
            raise AssimuloException("ARKODE symplectic: 'q_states' (the indices of the position states) must be set.")
        q = np.asarray(q)
        mask = np.zeros(self.pData.dim)
        if q.dtype == bool:
            if len(q) != self.pData.dim:
                raise AssimuloException("ARKODE symplectic: a boolean 'q_states' must have one entry per state.")
            mask[q] = 1.0
        else:
            mask[q.astype(int)] = 1.0
        if mask.sum() == 0 or mask.sum() == self.pData.dim:
            raise AssimuloException("ARKODE symplectic: 'q_states' must name some but not all states.")
        return np.ascontiguousarray(mask)

    cdef int _switch_table(self, double t, np.ndarray y, double tf, table) except -1:
        """Recreates the memory block with `table` at (t, y), keeping the last step size."""
        cdef double hlast = 0.0
        cdef int flag
        ARK.ARKodeGetLastStep(self.ark_mem, &hlast)
        self.store_statistics(ARK_TSTOP_RETURN)
        self.t = t
        self.y = y
        self._free_memory()
        self._active_table = table
        self.initialize_arkode()
        self.initialize_options()
        if hlast > 0.0:
            flag = ARK.ARKodeSetInitStep(self.ark_mem, hlast)   # the controller's last step, not a fresh estimate
            if flag < 0: raise ARKODEError(flag, t)
        flag = ARK.ARKodeSetStopTime(self.ark_mem, tf)
        if flag < 0: raise ARKODEError(flag, t)
        return 0

    cdef int _fallback_check(self, double t, np.ndarray y, double tf) except -1:
        """Failure-rate monitor (see the 'fallback_*' options). Returns 1 after switching the
        memory block to the fallback table at (t, y), else 0."""
        cdef long int nattempts = 0, nfails = 0
        table = self.options["fallback_table"]
        if table is None or self.fallback_time is not None or self.options["method"] != "implicit" \
                or table == self._active_table:
            return 0
        ARK.ARKodeGetNumStepAttempts(self.ark_mem, &nattempts)
        ARK.ARKodeGetNumStepSolveFails(self.ark_mem, &nfails)
        window = int(self.options["fallback_window"])
        if nattempts - self._fb_attempts0 < window:
            return 0
        # ARKODE counts a step failed by a recoverable rhs failure at a stage (an exception in
        # the problem's rhs: the FMU's own solver not converging at a trial point) together with
        # the steps Newton could not converge; only the latter say something about the method,
        # so the rhs failures of the window are taken out (at most one per failed step)
        nrhs = min(self.pData.nrhsfails - self._fb_rhsfails0, nfails - self._fb_fails0)
        rate = (nfails - self._fb_fails0 - nrhs) / float(nattempts - self._fb_attempts0)
        self._fb_attempts0 = nattempts
        self._fb_fails0 = nfails
        self._fb_rhsfails0 = self.pData.nrhsfails
        if rate <= float(self.options["fallback_conv_fail_rate"]):
            return 0
        self.log_message("ARKODE: %.0f %% of the last %d step attempts failed by nonlinear convergence (rhs failures "
                         "excluded) at t = %g; continuing with the table %s" % (100 * rate, window, t, table), NORMAL)
        self.fallback_time = t
        self._switch_table(t, y, tf, table)
        return 1

    cdef int _retry_after_failure(self, int flag, double t, double tf) except -1:
        """After a hard failure of the fallback table (repeated error-test or convergence
        failures), switch back to the user's table at the failure point, once. Returns 1 if
        the integration can go on, else 0."""
        cdef double tcur = t
        cdef N_Vector ycur
        cdef np.ndarray y
        if flag not in (ARK_ERR_FAILURE, ARK_CONV_FAILURE) or self.fallback_time is None \
                or self.fallback_return_time is not None:
            return 0
        ARK.ARKodeGetCurrentTime(self.ark_mem, &tcur)
        ycur = N_VNew_Serial(self.pData.dim, self.ctx)
        ARK.ARKodeGetDky(self.ark_mem, tcur, 0, ycur)
        y = nv2arr(ycur)
        N_VDestroy(ycur)
        self.log_message("ARKODE: the table %s failed at t = %g (flag %d); continuing with the table %s"
                         % (self._active_table, tcur, flag, self.options["table"] if self.options["table"]
                            else "of order %d" % self.options["order"]), NORMAL)
        self.fallback_return_time = tcur
        self._switch_table(tcur, y, tf, self.options["table"])
        return 1

    cdef int _retry_after_rhs_failure(self, int flag, double t, np.ndarray y, double tf) except -1:
        """A recoverable rhs failure at a directly evaluated stage (ARK_UNREC_RHSFUNC_ERR with
        the problem's failure counter advanced): restart at (t, y), the caller's last point --
        the last accepted step in one-step mode, the last output point otherwise (ARKODE's own
        current point may lie beyond a pending output time, which a reset there would make
        unreachable) -- with the step reduced by fail_factor. Returns 1 to go on, else 0."""
        cdef double tcur = t, hcur = 0.0
        cdef N_Vector ycur
        cdef int rflag
        if flag != ARK_UNREC_RHSFUNC_ERR or self.pData.nrhsfails == self._rhsfails_seen:
            return 0
        self._rhsfails_seen = self.pData.nrhsfails
        if self._consecutive_rhsfails >= int(self.options["fail_max"]):
            self.log_message("ARKODE: %d consecutive refused stage points at t = %g; giving up"
                             % (self._consecutive_rhsfails, t), NORMAL)
            return 0
        self._consecutive_rhsfails += 1
        ARK.ARKodeGetCurrentStep(self.ark_mem, &hcur)
        ycur = arr2nv(y, <void*>self.ctx)
        rflag = ARK.ARKodeReset(self.ark_mem, tcur, ycur)
        N_VDestroy(ycur)
        if rflag < 0: raise ARKODEError(rflag, tcur)
        if hcur <= 0.0:
            hcur = float(self.options["inith"])
        rflag = ARK.ARKodeSetInitStep(self.ark_mem, float(self.options["fail_factor"]) * hcur)
        if rflag < 0: raise ARKODEError(rflag, tcur)
        rflag = ARK.ARKodeSetStopTime(self.ark_mem, tf)
        if rflag < 0: raise ARKODEError(rflag, tcur)
        self.log_message("ARKODE: the rhs refused a stage point at t = %g (step %g); retrying from t = %g with the step %g"
                         % (t, hcur, tcur, float(self.options["fail_factor"]) * hcur), LOUD)
        return 1

    cpdef initialize_options(self):
        """Applies the options to the ARKODE memory (called on every (re)initialization)."""
        cdef int flag
        method = self.options["method"]

        # linear solver and Jacobian (implicit and imex only)
        if method in ("implicit", "imex"):
            if self.options["linear_solver"] in ("DENSE", "BLOCK"):
                if self.sun_matrix == NULL:
                    self.sun_matrix = SUNDIALS.SUNDenseMatrix(self.pData.dim, self.pData.dim, self.ctx)
                    if self.options["linear_solver"] == "DENSE":
                        self.sun_linearsolver = SUNDIALS.SUNLinSol_Dense(self.yTemp, self.sun_matrix, self.ctx)
                    else:
                        self._block_lu = BlockLU(self.problem.implicit_blocks, self.pData.dim,
                                                 1 if method == "implicit" else 0,
                                                 1 if self.options["block_check"] else 0,
                                                 float(self.options["block_check_tol"]))
                        self.sun_linearsolver = block_linear_solver(self._block_lu, self.ctx)
                    if self.sun_linearsolver == NULL:
                        raise ARKODEError(ARK_MEM_FAIL, self.t)
                    flag = ARK.ARKodeSetLinearSolver(self.ark_mem, self.sun_linearsolver, self.sun_matrix)
                    if flag < 0:
                        raise ARKODEError(flag, self.t)
                if self.pData.JAC != NULL and self.options["usejac"]:
                    flag = ARK.ARKodeSetJacFn(self.ark_mem, cv_jac)
                else:
                    flag = ARK.ARKodeSetJacFn(self.ark_mem, NULL)
                if flag < 0:
                    raise ARKODEError(flag, self.t)
            elif self.options["linear_solver"] == "SPGMR":
                if self.sun_linearsolver == NULL:
                    self.sun_linearsolver = SUNDIALS.SUNLinSol_SPGMR(self.yTemp, PREC_NONE, self.options["maxkrylov"], self.ctx)
                    flag = ARK.ARKodeSetLinearSolver(self.ark_mem, self.sun_linearsolver, NULL)
                    if flag < 0:
                        raise ARKODEError(flag, self.t)
                if self.pData.JACV != NULL and self.options["usejac"]:
                    flag = ARK.ARKodeSetJacTimes(self.ark_mem, SUNDIALS.cv_spils_jtsetup_dummy, cv_jacv)
                else:
                    flag = ARK.ARKodeSetJacTimes(self.ark_mem, NULL, NULL)
                if flag < 0:
                    raise ARKODEError(flag, self.t)
            else:
                raise AssimuloException("ARKODE: 'linear_solver' must be DENSE, BLOCK or SPGMR.")

            flag = ARK.ARKodeSetMaxNonlinIters(self.ark_mem, int(self.options["max_nonlin_iters"]))
            if flag < 0: raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetNonlinConvCoef(self.ark_mem, float(self.options["nonlin_conv_coef"]))
            if flag < 0: raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetPredictorMethod(self.ark_mem, int(self.options["predictor"]))
            if flag < 0: raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetDeduceImplicitRhs(self.ark_mem, 1 if self.options["deduce_implicit_rhs"] else 0)
            if flag < 0: raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetMaxConvFails(self.ark_mem, int(self.options["maxncf"]))
            if flag < 0: raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetLSetupFrequency(self.ark_mem, int(self.options["lsetup_frequency"]))
            if flag < 0: raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetJacEvalFrequency(self.ark_mem, int(self.options["jac_eval_frequency"]))
            if flag < 0: raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetDeltaGammaMax(self.ark_mem, float(self.options["delta_gamma_max"]))
            if flag < 0: raise ARKODEError(flag, self.t)

        # method: a named table, or ARKODE's default table of the requested order. Set once per
        # memory block: ARKStep frees its Butcher tables on these calls and re-creates them on the
        # first initialization only, not on ARKodeReset (the option setters recreate the memory).
        if self._fresh_memory:
            table = self._active_table
            if table is not None:
                if method == "imex":
                    # an ARK pair: (implicit, explicit) names, or the implicit name with the explicit
                    # partner by ARKODE's naming (ARKODE_ARK436L2SA_DIRK_6_3_4 / ..._ERK_6_3_4)
                    if isinstance(table, (tuple, list)):
                        itable, etable = str(table[0]), str(table[1])
                    else:
                        itable = str(table)
                        etable = itable.replace("_DIRK_", "_ERK_")
                    flag = ARK.ARKStepSetTableName(self.ark_mem, itable.encode("ascii"), etable.encode("ascii"))
                else:
                    tname = str(table).encode("ascii")
                    if method == "implicit":
                        flag = ARK.ARKStepSetTableName(self.ark_mem, tname, b"ARKODE_ERK_NONE")
                    elif method == "explicit":
                        flag = ARK.ARKStepSetTableName(self.ark_mem, b"ARKODE_DIRK_NONE", tname)
                    else:
                        flag = ARK.SPRKStepSetMethodName(self.ark_mem, tname)
                if flag < 0:
                    raise ARKODEError(flag, self.t, "unknown Butcher table '%s'" % (table,))
            else:
                flag = ARK.ARKodeSetOrder(self.ark_mem, int(self.options["order"]))
                if flag < 0:
                    raise ARKODEError(flag, self.t, "'order' %s is not available" % self.options["order"])
            if int(self.options["interpolant_degree"]) >= 0:
                flag = ARK.ARKodeSetInterpolantDegree(self.ark_mem, int(self.options["interpolant_degree"]))
                if flag < 0: raise ARKODEError(flag, self.t)
            self._fresh_memory = 0

        if method == "symplectic":
            # SPRKStep: no adaptivity, a fixed step, no tolerances (the rootfinding needs none)
            if float(self.options["fixed_h"]) <= 0.0:
                raise AssimuloException("ARKODE symplectic: 'fixed_h' (the fixed step size) must be positive.")
            flag = ARK.ARKodeSetFixedStep(self.ark_mem, float(self.options["fixed_h"]))
            if flag < 0: raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetUseCompensatedSums(self.ark_mem, 1 if self.options["compensated_sums"] else 0)
            if flag < 0: raise ARKODEError(flag, self.t)
            flag = ARK.ARKodeSetMaxNumSteps(self.ark_mem, int(self.options["maxsteps"]))
            if flag < 0: raise ARKODEError(flag, self.t)
            return

        # step-size adaptivity bounds
        flag = ARK.ARKodeSetSafetyFactor(self.ark_mem, float(self.options["safety"]))
        if flag < 0: raise ARKODEError(flag, self.t)
        flag = ARK.ARKodeSetMaxGrowth(self.ark_mem, float(self.options["max_growth"]))
        if flag < 0: raise ARKODEError(flag, self.t)
        flag = ARK.ARKodeSetMaxFirstGrowth(self.ark_mem, float(self.options["max_first_growth"]))
        if flag < 0: raise ARKODEError(flag, self.t)
        flag = ARK.ARKodeSetMaxEFailGrowth(self.ark_mem, float(self.options["max_efail_growth"]))
        if flag < 0: raise ARKODEError(flag, self.t)
        if method in ("implicit", "imex"):
            flag = ARK.ARKodeSetMaxCFailGrowth(self.ark_mem, float(self.options["max_cfail_growth"]))
            if flag < 0: raise ARKODEError(flag, self.t)
        # step limits
        flag = ARK.ARKodeSetMaxErrTestFails(self.ark_mem, int(self.options["maxnef"]))
        if flag < 0: raise ARKODEError(flag, self.t)
        flag = ARK.ARKodeSetMaxNumSteps(self.ark_mem, int(self.options["maxsteps"]))
        if flag < 0: raise ARKODEError(flag, self.t)
        flag = ARK.ARKodeSetMaxStep(self.ark_mem, float(self.options["maxh"]))
        if flag < 0: raise ARKODEError(flag, self.t)
        flag = ARK.ARKodeSetMinStep(self.ark_mem, float(self.options["minh"]))
        if flag < 0: raise ARKODEError(flag, self.t)
        # first step: the user's inith, else ARKODE's estimate; after an event the last step is
        # kept only with restart_h = "keep" (a full-size step over a fresh transient hides it
        # from the error estimator; measured on the native TRBDF2)
        if self.options["restart_h"] == "keep" and self.statistics["nsteps"] > 0:
            pass
        else:
            flag = ARK.ARKodeSetInitStep(self.ark_mem, float(self.options["inith"]))
            if flag < 0: raise ARKODEError(flag, self.t)

        # tolerances (atol may have changed with the nominals in handle_event)
        if self.nv_atol != NULL:
            N_VDestroy(self.nv_atol)
        self.nv_atol = arr2nv(self.options["atol"], <void*>self.ctx)
        flag = ARK.ARKodeSVtolerances(self.ark_mem, float(self.options["rtol"]), self.nv_atol)
        if flag < 0: raise ARKODEError(flag, self.t)

    def initialize_event_detection(self):
        if self.options["external_event_detection"]:
            def event_func(t, y):
                return 0, self.problem.state_events(t, y, self.sw)
            self.event_func = event_func
            _, self.g_old = self.event_func(self.t, self.y)
        else:
            def event_func(t, y):
                return self.problem.state_events(t, y, self.sw)
            self.event_func = event_func
            self.g_old = self.event_func(self.t, self.y)

    # ------------------------------------------------------------------ integration
    cpdef integrate(self, double t, np.ndarray[ndim=1, dtype=realtype] y, double tf, dict opts):
        cdef int flag, output_index
        cdef N_Vector yout
        cdef double tret = self.t, tout
        cdef list tr = [], yr = []
        cdef np.ndarray output_list

        yout = arr2nv(y, <void*>self.ctx)

        if opts["initialize"]:
            self.initialize_arkode()
            self.initialize_options()
            if self.options["external_event_detection"]:
                self.initialize_event_detection()

        flag = ARK.ARKodeSetStopTime(self.ark_mem, tf)
        if flag < 0:
            N_VDestroy(yout)
            raise ARKODEError(flag, t)

        if opts["report_continuously"] or opts["output_list"] is None:
            nsteps_call = 0
            while True:
                flag = ARK.ARKodeEvolve(self.ark_mem, tf, yout, &tret, ARK_ONE_STEP)
                if flag < 0:
                    if self._retry_after_failure(flag, tret, tf) or self._retry_after_rhs_failure(flag, t, y, tf):
                        continue
                    self.store_statistics(ARK_TSTOP_RETURN)
                    N_VDestroy(yout)
                    raise ARKODEError(flag, tret, self._block_message())
                self._consecutive_rhsfails = 0
                # ARKODE's maxsteps counts per ARKodeEvolve call, i.e. per step here -- as CVode's
                # mxstep in the same one-step mode, so no limit applies between output points
                nsteps_call += 1

                t = tret
                y = nv2arr(yout)
                if self.options["external_event_detection"] and self.problem_info["state_events"]:
                    told = self.t if (opts["report_continuously"] or not tr) else tr[-1]
                    event_flag, t, y = self.event_locator(told, t, y)
                    if event_flag == ID_PY_EVENT:
                        flag = ARK_ROOT_RETURN

                if opts["report_continuously"]:
                    try:
                        flag_initialize = self.report_solution(t, y, opts)
                    except Exception:
                        self.store_statistics(ARK_TSTOP_RETURN)
                        raise
                    if flag_initialize:
                        flag = CV_STEP_RETURN    # a step event: the integration is reinitialized
                else:
                    tr.append(t)
                    yr.append(y)

                if flag == ARK_ROOT_RETURN or flag == CV_STEP_RETURN:
                    self.store_statistics(flag)
                    flag = ID_EVENT
                    break
                if flag == ARK_TSTOP_RETURN:
                    flag = ID_COMPLETE
                    self.store_statistics(ARK_TSTOP_RETURN)
                    break
                self._fallback_check(t, y, tf)
        else:
            output_index = opts["output_index"]
            output_list = opts["output_list"][output_index:]
            for tout in output_list:
                flag = ARK.ARKodeEvolve(self.ark_mem, tout, yout, &tret, ARK_NORMAL)
                while flag < 0 and (self._retry_after_failure(flag, tret, tf) or
                                    self._retry_after_rhs_failure(flag, tr[-1] if tr else t, yr[-1] if yr else y, tf)):
                    flag = ARK.ARKodeEvolve(self.ark_mem, tout, yout, &tret, ARK_NORMAL)
                if flag < 0:
                    self.store_statistics(ARK_TSTOP_RETURN)
                    N_VDestroy(yout)
                    raise ARKODEError(flag, tret, self._block_message())
                tr.append(tret)
                yr.append(nv2arr(yout))
                if flag == ARK_ROOT_RETURN:
                    self.store_statistics(ARK_ROOT_RETURN)
                    flag = ID_EVENT
                    if tret == tout:
                        output_index += 1
                    break
                if flag == ARK_TSTOP_RETURN:
                    self.store_statistics(ARK_TSTOP_RETURN)
                    flag = ID_COMPLETE
                    if tret == tout:
                        output_index += 1
                    break
                output_index += 1
                self._fallback_check(tret, yr[-1], tf)
            else:
                flag = ID_COMPLETE
                self.store_statistics(ARK_TSTOP_RETURN)
            opts["output_index"] = output_index

        N_VDestroy(yout)
        return flag, tr, yr

    cpdef step(self, double t, np.ndarray y, double tf, dict opts):
        cdef int flag
        cdef N_Vector yout
        cdef double tret = t
        yout = arr2nv(y, <void*>self.ctx)
        if opts["initialize"]:
            self.initialize_arkode()
            self.initialize_options()
        flag = ARK.ARKodeSetStopTime(self.ark_mem, tf)
        if flag < 0:
            raise ARKODEError(flag, t)
        flag = ARK.ARKodeEvolve(self.ark_mem, tf, yout, &tret, ARK_ONE_STEP)
        if flag < 0:
            raise ARKODEError(flag, tret, self._block_message())
        tr = tret
        yr = nv2arr(yout)
        if flag == ARK_ROOT_RETURN:
            flag = ID_EVENT
            self.store_statistics(ARK_ROOT_RETURN)
        if flag == ARK_TSTOP_RETURN:
            flag = ID_COMPLETE
            self.store_statistics(ARK_TSTOP_RETURN)
        N_VDestroy(yout)
        return flag, tr, yr

    cpdef state_event_info(self):
        if self.options["external_event_detection"]:
            return self._event_info
        if self.pData.dimRoot == 0:
            return []           # a step event of a problem without state events: no rootfinding to ask
        cdef int* c_info
        cdef int flag
        c_info = <int*>malloc(self.pData.dimRoot * sizeof(int))
        for k in range(self.pData.dimRoot):
            c_info[k] = 0
        flag = ARK.ARKodeGetRootInfo(self.ark_mem, c_info)
        if flag < 0:
            free(c_info)
            raise ARKODEError(flag)
        event_info = [c_info[k] for k in range(self.pData.dimRoot)]
        free(c_info)
        return event_info

    def set_event_info(self, event_info):
        self._event_info = event_info

    cpdef np.ndarray interpolate(self, double t, int k = 0):
        """The dense output (ARKodeGetDky) at t within the last step; k = 0 for y, 1 for its derivative."""
        cdef int flag
        cdef N_Vector dky = N_VNew_Serial(self.pData.dim, self.ctx)
        flag = ARK.ARKodeGetDky(self.ark_mem, t, k, dky)
        if flag < 0:
            N_VDestroy(dky)
            raise ARKODEError(flag, t)
        res = nv2arr(dky)
        N_VDestroy(dky)
        return res

    cpdef get_local_errors(self):
        cdef int flag
        cdef N_Vector ele = N_VNew_Serial(self.pData.dim, self.ctx)
        flag = ARK.ARKodeGetEstLocalErrors(self.ark_mem, ele)
        if flag < 0:
            N_VDestroy(ele)
            raise ARKODEError(flag, self.t)
        res = nv2arr(ele)
        N_VDestroy(ele)
        return res

    cpdef get_error_weights(self):
        cdef int flag
        cdef N_Vector ew = N_VNew_Serial(self.pData.dim, self.ctx)
        flag = ARK.ARKodeGetErrWeights(self.ark_mem, ew)
        if flag < 0:
            N_VDestroy(ew)
            raise ARKODEError(flag, self.t)
        res = nv2arr(ew)
        N_VDestroy(ew)
        return res

    def get_weighted_local_errors(self):
        return np.abs(self.get_local_errors() * self.get_error_weights())

    cdef void store_statistics(self, int return_flag):
        """Adds ARKODE's counters (cumulative since the memory was created) as increments."""
        cdef long int nsteps = 0, nattempts = 0, netfails = 0, nfe = 0, nfi = 0, njevals = 0, nfevalsLS = 0
        cdef long int nniters = 0, nncfails = 0, nlinsetups = 0, ngevals = 0, njvevals = 0, nsolvefails = 0
        if self.ark_mem == NULL:
            raise ARKODEError(ARK_MEM_NULL)
        ARK.ARKodeGetNumSteps(self.ark_mem, &nsteps)
        ARK.ARKodeGetNumStepAttempts(self.ark_mem, &nattempts)
        ARK.ARKodeGetNumErrTestFails(self.ark_mem, &netfails)
        ARK.ARKodeGetNumRhsEvals(self.ark_mem, 0, &nfe)
        ARK.ARKodeGetNumRhsEvals(self.ark_mem, 1, &nfi)
        if self.options["method"] in ("implicit", "imex"):
            ARK.ARKodeGetNumNonlinSolvIters(self.ark_mem, &nniters)
            ARK.ARKodeGetNumNonlinSolvConvFails(self.ark_mem, &nncfails)
            ARK.ARKodeGetNumStepSolveFails(self.ark_mem, &nsolvefails)
            ARK.ARKodeGetNumLinSolvSetups(self.ark_mem, &nlinsetups)
            if self.options["linear_solver"] == "SPGMR":
                ARK.ARKodeGetNumJtimesEvals(self.ark_mem, &njvevals)
                ARK.ARKodeGetNumLinRhsEvals(self.ark_mem, &nfevalsLS)
            else:
                ARK.ARKodeGetNumJacEvals(self.ark_mem, &njevals)
                ARK.ARKodeGetNumLinRhsEvals(self.ark_mem, &nfevalsLS)
        if self.pData.ROOT != NULL:
            ARK.ARKodeGetNumGEvals(self.ark_mem, &ngevals)

        # the counters are cumulative over the life of the memory block (ARKodeReset keeps them):
        # store the difference to the last stored values
        cur = {"nsteps": nsteps, "nstepattempts": nattempts, "nerrfails": netfails, "nfcns": nfe + nfi,
               "nfcns_implicit": nfi,
               "nniters": nniters, "nnfails": nncfails, "nconvfails": nsolvefails, "nlus": nlinsetups,
               "njacs": njevals, "njacvecs": njvevals, "nfcnjacs": nfevalsLS, "nstatefcns": ngevals}
        cur["nrhsfails"] = self.pData.nrhsfails
        last = self._last_counters
        for key, val in cur.items():
            prev = last.get(key, 0)
            if val < prev:
                prev = 0                 # ARKodeReset zeroes some counters (the root evaluations)
            self.statistics[key] += val - prev
        self._last_counters = cur
        if return_flag == ARK_ROOT_RETURN and not self.options["external_event_detection"]:
            self.statistics["nstateevents"] += 1

    def print_statistics(self, verbose=NORMAL):
        Explicit_ODE.print_statistics(self, verbose)
        log = lambda msg: self.log_message(msg, verbose)
        log("\nSolver options:\n")
        log(" Solver                  : ARKODE (%s, %s)" % (self.options["method"],
            self.options["table"] if self.options["table"] else "order %d" % self.options["order"]))
        if self.fallback_time is not None:
            log(" Fallback                : %s from t = %g%s" % (self.options["fallback_table"], self.fallback_time,
                "" if self.fallback_return_time is None else ", back from t = %g" % self.fallback_return_time))
        if self.options["method"] in ("implicit", "imex"):
            log(" Linear solver           : " + self.options["linear_solver"]
                + ("" if self._block_lu is None else " (%d blocks of sizes %s, %d states outside)"
                   % (self._block_lu.nblocks, self._block_lu.sizes(), self._block_lu.nout)))
        log(" Tolerances (absolute)   : " + str(self._compact_tol(self.options["atol"])))
        log(" Tolerances (relative)   : " + str(self.options["rtol"]))
        log("")

    # ------------------------------------------------------------------ options as properties
    def _set_atol(self, atol):
        atol = set_type_shape_array(atol)
        if len(atol) == 1:
            atol = atol * np.ones(self.pData.dim)
        elif len(atol) != self.pData.dim:
            raise AssimuloException("atol must be of length one or same as the dimension of the problem.")
        if (atol < 0.0).any():
            raise AssimuloException("The absolute tolerance must be positive.")
        self.options["atol"] = atol
    def _get_atol(self):
        return self.options["atol"]
    atol = property(_get_atol, _set_atol)

    def _set_rtol(self, rtol):
        rtol = float(rtol)
        if rtol < 0.0:
            raise AssimuloException("The relative tolerance must be positive.")
        self.options["rtol"] = rtol
    def _get_rtol(self):
        return self.options["rtol"]
    rtol = property(_get_rtol, _set_rtol)

    def _set_maxh(self, maxh):
        self.options["maxh"] = float(maxh) if maxh else 0.0
    def _get_maxh(self):
        return self.options["maxh"]
    maxh = property(_get_maxh, _set_maxh)

    def _set_minh(self, minh):
        self.options["minh"] = float(minh)
    def _get_minh(self):
        return self.options["minh"]
    minh = property(_get_minh, _set_minh)

    def _set_inith(self, inith):
        self.options["inith"] = float(inith)
    def _get_inith(self):
        return self.options["inith"]
    inith = property(_get_inith, _set_inith)

    def _set_maxsteps(self, n):
        self.options["maxsteps"] = int(n)
    def _get_maxsteps(self):
        return self.options["maxsteps"]
    maxsteps = property(_get_maxsteps, _set_maxsteps)

    def _set_q_states(self, q):
        self.options["q_states"] = None if q is None else np.asarray(q)
    def _get_q_states(self):
        return self.options["q_states"]
    q_states = property(_get_q_states, _set_q_states)

    def _set_fixed_h(self, h):
        self.options["fixed_h"] = float(h)
    def _get_fixed_h(self):
        return self.options["fixed_h"]
    fixed_h = property(_get_fixed_h, _set_fixed_h)

    def _set_compensated_sums(self, b):
        self.options["compensated_sums"] = bool(b)
    def _get_compensated_sums(self):
        return self.options["compensated_sums"]
    compensated_sums = property(_get_compensated_sums, _set_compensated_sums)

    def _set_fallback_table(self, table):
        self.options["fallback_table"] = str(table) if table else None
    def _get_fallback_table(self):
        return self.options["fallback_table"]
    fallback_table = property(_get_fallback_table, _set_fallback_table)

    def _set_fallback_conv_fail_rate(self, rate):
        rate = float(rate)
        if not 0.0 < rate <= 1.0:
            raise AssimuloException("fallback_conv_fail_rate must be in (0, 1].")
        self.options["fallback_conv_fail_rate"] = rate
    def _get_fallback_conv_fail_rate(self):
        return self.options["fallback_conv_fail_rate"]
    fallback_conv_fail_rate = property(_get_fallback_conv_fail_rate, _set_fallback_conv_fail_rate)

    def _set_fallback_window(self, n):
        n = int(n)
        if n < 1:
            raise AssimuloException("fallback_window must be at least 1.")
        self.options["fallback_window"] = n
    def _get_fallback_window(self):
        return self.options["fallback_window"]
    fallback_window = property(_get_fallback_window, _set_fallback_window)

    def _set_usejac(self, jac):
        self.options["usejac"] = bool(jac)
    def _get_usejac(self):
        return self.options["usejac"]
    usejac = property(_get_usejac, _set_usejac)

    def _set_linear_solver(self, ls):
        ls = str(ls).upper()
        if ls not in ("DENSE", "BLOCK", "SPGMR"):
            raise AssimuloException("'linear_solver' must be DENSE, BLOCK or SPGMR.")
        if ls != self.options["linear_solver"]:
            self._free_memory()
        self.options["linear_solver"] = ls
    def _get_linear_solver(self):
        return self.options["linear_solver"]
    linear_solver = property(_get_linear_solver, _set_linear_solver)

    def _set_fail_factor(self, v):
        v = float(v)
        if not 0.0 < v < 1.0:
            raise AssimuloException("fail_factor must be in (0, 1).")
        self.options["fail_factor"] = v
    def _get_fail_factor(self):
        return self.options["fail_factor"]
    fail_factor = property(_get_fail_factor, _set_fail_factor)

    def _set_fail_max(self, n):
        n = int(n)
        if n < 0:
            raise AssimuloException("fail_max must be non-negative.")
        self.options["fail_max"] = n
    def _get_fail_max(self):
        return self.options["fail_max"]
    fail_max = property(_get_fail_max, _set_fail_max)

    def _set_block_check(self, b):
        self.options["block_check"] = bool(b)
    def _get_block_check(self):
        return self.options["block_check"]
    block_check = property(_get_block_check, _set_block_check)

    def _set_block_check_tol(self, tol):
        tol = float(tol)
        if tol <= 0.0:
            raise AssimuloException("block_check_tol must be positive.")
        self.options["block_check_tol"] = tol
    def _get_block_check_tol(self):
        return self.options["block_check_tol"]
    block_check_tol = property(_get_block_check_tol, _set_block_check_tol)

    def _set_maxkrylov(self, n):
        self.options["maxkrylov"] = int(n)
    def _get_maxkrylov(self):
        return self.options["maxkrylov"]
    maxkrylov = property(_get_maxkrylov, _set_maxkrylov)

    def _set_external_event_detection(self, v):
        v = bool(v)
        if v != self.options["external_event_detection"]:
            self._free_memory()                 # the root function set-up differs
        self.options["external_event_detection"] = v
    def _get_external_event_detection(self):
        return self.options["external_event_detection"]
    external_event_detection = property(_get_external_event_detection, _set_external_event_detection)

    def _set_method(self, m):
        m = str(m).lower()
        if m not in ("implicit", "explicit", "imex", "symplectic"):
            raise AssimuloException("'method' must be 'implicit', 'explicit', 'imex' or 'symplectic'.")
        if m != self.options["method"]:
            self._free_memory()
        self.options["method"] = m
    def _get_method(self):
        return self.options["method"]
    method = property(_get_method, _set_method)

    def _set_order(self, q):
        if int(q) != self.options["order"]:
            self._free_memory()
        self.options["order"] = int(q)
    def _get_order(self):
        return self.options["order"]
    order = property(_get_order, _set_order)

    def _set_table(self, name):
        if isinstance(name, (tuple, list)):
            name = tuple(str(x) for x in name)      # imex: (implicit, explicit) pair
        else:
            name = None if name is None else str(name)
        if name != self.options["table"]:
            self._free_memory()
        self.options["table"] = name
    def _get_table(self):
        return self.options["table"]
    table = property(_get_table, _set_table)

    def _set_predictor(self, p):
        self.options["predictor"] = int(p)
    def _get_predictor(self):
        return self.options["predictor"]
    predictor = property(_get_predictor, _set_predictor)

    def _set_nonlin_conv_coef(self, c):
        self.options["nonlin_conv_coef"] = float(c)
    def _get_nonlin_conv_coef(self):
        return self.options["nonlin_conv_coef"]
    nonlin_conv_coef = property(_get_nonlin_conv_coef, _set_nonlin_conv_coef)

    def _set_max_nonlin_iters(self, n):
        self.options["max_nonlin_iters"] = int(n)
    def _get_max_nonlin_iters(self):
        return self.options["max_nonlin_iters"]
    max_nonlin_iters = property(_get_max_nonlin_iters, _set_max_nonlin_iters)

    def _set_deduce_implicit_rhs(self, v):
        self.options["deduce_implicit_rhs"] = bool(v)
    def _get_deduce_implicit_rhs(self):
        return self.options["deduce_implicit_rhs"]
    deduce_implicit_rhs = property(_get_deduce_implicit_rhs, _set_deduce_implicit_rhs)

    def _set_interpolant_degree(self, d):
        if int(d) != self.options["interpolant_degree"]:
            self._free_memory()
        self.options["interpolant_degree"] = int(d)
    def _get_interpolant_degree(self):
        return self.options["interpolant_degree"]
    interpolant_degree = property(_get_interpolant_degree, _set_interpolant_degree)

    def _set_maxncf(self, n):
        self.options["maxncf"] = int(n)
    def _get_maxncf(self):
        return self.options["maxncf"]
    maxncf = property(_get_maxncf, _set_maxncf)

    def _set_maxnef(self, n):
        self.options["maxnef"] = int(n)
    def _get_maxnef(self):
        return self.options["maxnef"]
    maxnef = property(_get_maxnef, _set_maxnef)

    def _set_lsetup_frequency(self, n):
        self.options["lsetup_frequency"] = int(n)
    def _get_lsetup_frequency(self):
        return self.options["lsetup_frequency"]
    lsetup_frequency = property(_get_lsetup_frequency, _set_lsetup_frequency)

    def _set_jac_eval_frequency(self, n):
        self.options["jac_eval_frequency"] = int(n)
    def _get_jac_eval_frequency(self):
        return self.options["jac_eval_frequency"]
    jac_eval_frequency = property(_get_jac_eval_frequency, _set_jac_eval_frequency)

    def _set_delta_gamma_max(self, v):
        self.options["delta_gamma_max"] = float(v)
    def _get_delta_gamma_max(self):
        return self.options["delta_gamma_max"]
    delta_gamma_max = property(_get_delta_gamma_max, _set_delta_gamma_max)

    def _set_safety(self, v):
        self.options["safety"] = float(v)
    def _get_safety(self):
        return self.options["safety"]
    safety = property(_get_safety, _set_safety)

    def _set_max_growth(self, v):
        self.options["max_growth"] = float(v)
    def _get_max_growth(self):
        return self.options["max_growth"]
    max_growth = property(_get_max_growth, _set_max_growth)

    def _set_max_first_growth(self, v):
        self.options["max_first_growth"] = float(v)
    def _get_max_first_growth(self):
        return self.options["max_first_growth"]
    max_first_growth = property(_get_max_first_growth, _set_max_first_growth)

    def _set_max_efail_growth(self, v):
        self.options["max_efail_growth"] = float(v)
    def _get_max_efail_growth(self):
        return self.options["max_efail_growth"]
    max_efail_growth = property(_get_max_efail_growth, _set_max_efail_growth)

    def _set_max_cfail_growth(self, v):
        self.options["max_cfail_growth"] = float(v)
    def _get_max_cfail_growth(self):
        return self.options["max_cfail_growth"]
    max_cfail_growth = property(_get_max_cfail_growth, _set_max_cfail_growth)

    def _set_restart_h(self, v):
        v = str(v).lower()
        if v not in ("estimate", "keep"):
            raise AssimuloException("'restart_h' must be 'estimate' or 'keep'.")
        self.options["restart_h"] = v
    def _get_restart_h(self):
        return self.options["restart_h"]
    restart_h = property(_get_restart_h, _set_restart_h)


def _same_blocks(a, b):
    """Whether two block lists name the same index sets in the same order."""
    try:
        if len(a) != len(b):
            return False
        for x, y in zip(a, b):
            y = np.asarray(y)
            if y.dtype == bool:
                y = np.flatnonzero(y)
            if not np.array_equal(np.asarray(x), np.unique(y.astype(np.int64))):
                return False
        return True
    except Exception:
        return False


class ARKODEError(Exception):
    """An ARKODE return flag with its textual meaning."""
    msg = {ARK_TOO_MUCH_WORK: "The solver took the maximum number of internal steps but could not reach tout.",
           ARK_TOO_MUCH_ACC: "The solver could not satisfy the accuracy demanded for some internal step.",
           ARK_ERR_FAILURE: "Error test failures occurred too many times during one internal step or the minimum step size was reached.",
           ARK_CONV_FAILURE: "Convergence test failures occurred too many times during one internal step or the minimum step size was reached.",
           ARK_LINIT_FAIL: "The linear solver's initialization function failed.",
           ARK_LSETUP_FAIL: "The linear solver's setup function failed in an unrecoverable manner.",
           ARK_LSOLVE_FAIL: "The linear solver's solve function failed in an unrecoverable manner.",
           ARK_RHSFUNC_FAIL: "The right-hand side function failed in an unrecoverable manner.",
           ARK_FIRST_RHSFUNC_ERR: "The right-hand side function failed at the first call.",
           ARK_REPTD_RHSFUNC_ERR: "The right-hand side function had repeated recoverable errors.",
           ARK_UNREC_RHSFUNC_ERR: "The right-hand side function had a recoverable error, but no recovery is possible.",
           ARK_RTFUNC_FAIL: "The rootfinding function failed in an unrecoverable manner.",
           ARK_MEM_FAIL: "A memory allocation failed.",
           ARK_MEM_NULL: "The ARKODE memory was NULL.",
           ARK_ILL_INPUT: "One of the function inputs is illegal.",
           ARK_TOO_CLOSE: "The output and initial times are too close to each other.",
           ARK_INTERP_FAIL: "The interpolation module failed.",
           ARK_STEPPER_UNSUPPORTED: "The stepper does not support this operation."}

    def __init__(self, value, t=0.0, extra=None):
        self.value = value
        self.t = t
        self.extra = extra

    def __str__(self):
        text = self.msg.get(self.value, "Unknown error (%d)." % self.value)
        if self.extra:
            text += " " + str(self.extra)
        return "ARKODE failed with flag %d at time %g: %s" % (self.value, self.t, text)
