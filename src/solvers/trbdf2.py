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

import numpy as np
import scipy.sparse as sps

from assimulo.exception import AssimuloException, AssimuloRecoverableError, TimeLimitExceeded
from assimulo.ode import NORMAL, ID_PY_EVENT, ID_PY_COMPLETE
from assimulo.explicit_ode import Explicit_ODE


class TRBDF2Error(AssimuloException):
    """Error raised by the TR-BDF2 solver."""
    def __init__(self, value, t=None, err_msg=""):
        self.value = value
        self.t = t
        self.err_msg = err_msg

    def __str__(self):
        return "TR-BDF2 failed with flag %d%s: %s" % (self.value, "" if self.t is None else " at time %g" % self.t, self.err_msg)


class TRBDF2(Explicit_ODE):
    """
    TR-BDF2: a second-order, L-stable, stiffly accurate one-step method (a 3-stage ESDIRK)
    with a third-order embedded error estimate filtered for stiffness, a piecewise cubic
    Hermite dense output, and a CVode-like Jacobian policy (the Jacobian is kept until Newton
    fails or 'maxsteps_jac' steps have passed; the LU is reused while the step stays within
    'lu_band' of the step it was formed with). The C implementation is assimulo/lib/trbdf2.c.

    Reference: M. E. Hosea and L. F. Shampine, "Analysis and implementation of TR-BDF2",
    Applied Numerical Mathematics 20 (1996) 21-37.

    Events are located by Assimulo's event locator on the dense output after every step, so
    the solver supports state, time and step events (one-step mode).
    """

    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem)

        # Default values
        self.options["inith"] = 0.0             # 0: estimate the first step
        self.options["maxh"] = None             # None: no maximum step
        self.options["maxsteps"] = 100000
        self.options["rtol"] = 1.0e-6
        self.options["atol"] = 1.0e-6 * np.ones(self.problem_info["dim"])
        self.options["usejac"] = True if self.problem_info["jac_fcn"] else False
        self.options["newt"] = 6                # maximum Newton iterations per stage
        self.options["newton_tol"] = 0.1        # Newton stopping tolerance in the weighted norm (CVode's value)
        self.options["safe"] = 0.9              # step-size safety factor
        self.options["fac1"] = 0.2              # smallest step-size decrease factor
        self.options["fac2"] = 5.0              # largest step-size increase factor
        self.options["lu_band"] = (0.8, 1.25)   # reuse the LU while h/h_lu is within
        self.options["maxsteps_jac"] = 50       # steps between Jacobian evaluations without Newton failure
        self.options["fail_factor"] = 0.25      # step reduction on a recoverable rhs/jac failure
        self.options["fail_max"] = 40           # budget of consecutive recoverable failures
        self.options["report_continuously"] = False

        # Statistics ("nlus" is a default key of the base class)
        self.statistics.add_key("nrhsfails", "Number of recoverable rhs failures")
        self.statistics.add_key("njacfails", "Number of refused user Jacobians (differences used)")
        self.statistics.add_key("njacstale", "Number of refused Jacobians (previous one kept)")
        self.statistics.add_key("nnfails_div", "Number of Newton failures by divergence")
        self.statistics.add_key("nnfails_iter", "Number of Newton failures by the iteration limit")

        # Solver support
        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = True
        self.supports["state_events"] = True

        self._leny = len(self.y)
        self._memory = None
        self._py_err = None
        self._event_info = [0] * self.problem_info["dimRoot"] if self.problem_info["state_events"] else []
        self._werr = np.zeros(self._leny)

    # ------------------------------------------------------------------ set-up
    def initialize(self):
        self.statistics.reset()
        try:
            from assimulo.lib import trbdf2ode
            self._impl = trbdf2ode
        except Exception:
            raise TRBDF2Error(-1, err_msg="Failed to import the TR-BDF2 solver (assimulo.lib.trbdf2ode).") from None
        if self.usejac and not hasattr(self.problem, "jac"):
            raise TRBDF2Error(-1, err_msg="Use of an analytical Jacobian is enabled, but the problem has no 'jac' function.")
        self._memory = self._impl.TRBDF2Memory(self._leny)
        self._apply_options()

    def _apply_options(self):
        m = self._memory
        rtol = self.options["rtol"]
        rtol_vec = np.asarray(rtol, dtype=float) * np.ones(self._leny) if np.ndim(rtol) == 0 else np.asarray(rtol, dtype=float)
        atol_vec = np.asarray(self.options["atol"], dtype=float) * np.ones(self._leny)
        for name, ret in [("tolerances", m.set_tolerances(rtol_vec, atol_vec)),
                          ("maxh", m.set_hmax(self.options["maxh"] if self.options["maxh"] else 0.0)),
                          ("maxsteps", m.set_max_steps(int(self.options["maxsteps"]))),
                          ("newt/newton_tol", m.set_newton(int(self.options["newt"]), float(self.options["newton_tol"]))),
                          ("safe/fac1/fac2/lu_band", m.set_step_control(float(self.options["safe"]), float(self.options["fac1"]), float(self.options["fac2"]),
                                                                        float(self.options["lu_band"][0]), float(self.options["lu_band"][1]))),
                          ("maxsteps_jac/usejac", m.set_jac_policy(int(self.options["maxsteps_jac"]), int(bool(self.usejac)))),
                          ("fail_factor/fail_max", m.set_failure_policy(float(self.options["fail_factor"]), int(self.options["fail_max"])))]:
            if ret != 0:
                raise TRBDF2Error(ret, err_msg="Invalid option(s): %s" % name)

    def set_problem_data(self):
        if self.problem_info["state_events"]:
            def event_func(t, y):
                try:
                    res = self.problem.state_events(t, y, self.sw)
                except BaseException as E:
                    self._py_err = E
                    return -1, None
                return 0, res

            def f(t, y):
                try:
                    return self.problem.rhs(t, y, self.sw), [0]
                except BaseException as E:
                    return self._callback_failure(E, y)
            self.event_func = event_func
            self._event_info = [0] * self.problem_info["dimRoot"]
            ret, g0 = self.event_func(self.t, self.y)
            if ret < 0:
                raise self._py_err
            self.g_old = np.array(g0)
            self.statistics["nstatefcns"] += 1
        else:
            def f(t, y):
                try:
                    return self.problem.rhs(t, y), [0]
                except BaseException as E:
                    return self._callback_failure(E, y)
        self.f = f

    def _callback_failure(self, E, y):
        if isinstance(E, (np.linalg.LinAlgError, ZeroDivisionError, AssimuloRecoverableError)):
            return y.copy(), [1]        # recoverable: the trial point is refused
        self._py_err = E
        return y.copy(), [-1]           # unrecoverable

    def _jacobian(self, t, y):
        try:
            jac = self.problem.jac(t, y)
            if sps.issparse(jac):
                jac = jac.toarray()
            return np.asarray(jac, dtype=float), [0]
        except BaseException as E:
            return np.eye(self._leny), self._callback_failure(E, y)[1]

    # ------------------------------------------------------------------ integration
    def interpolate(self, time):
        y = np.empty(self._leny)
        self._memory.interpolate(time, y)
        return y

    def get_weighted_local_errors(self):
        """The weighted local error estimate of the last accepted step."""
        return np.abs(self._werr)

    def _solout(self, naccpt, told, t, y, werr):
        """Called by the C stepper after every accepted step."""
        try:
            self._werr = werr
            ret, flag = 0, 0
            if self.problem_info["state_events"]:
                flag, t, y = self.event_locator(told, t, y)
                if flag == ID_PY_EVENT:
                    ret = 1
                if flag < 0:
                    ret = -1
            if self._opts["report_continuously"]:
                try:
                    if self.report_solution(t, y.copy(), self._opts):
                        ret = 1
                except TimeLimitExceeded as e:
                    self._py_err = e
                    ret = -2
            else:
                if self._opts["output_list"] is None:
                    self._tlist.append(t)
                    self._ylist.append(y.copy())
                else:
                    output_list = self._opts["output_list"]
                    output_index = self._opts["output_index"]
                    try:
                        while output_list[output_index] <= t:
                            self._tlist.append(output_list[output_index])
                            self._ylist.append(self.interpolate(output_list[output_index]))
                            output_index += 1
                    except IndexError:
                        pass
                    self._opts["output_index"] = output_index
                    if self.problem_info["state_events"] and flag == ID_PY_EVENT and len(self._tlist) > 0 and self._tlist[-1] != t:
                        self._tlist.append(t)
                        self._ylist.append(y)
        except BaseException as E:
            self._py_err = E
            ret = -1
        return ret

    def integrate(self, t, y, tf, opts):
        if opts["initialize"]:
            self.set_problem_data()
            self._tlist = []
            self._ylist = []
            self._memory.reinit()
            self._apply_options()         # tolerances may have changed in handle_event (nominals)
        self._py_err = None
        self._opts = opts

        jac = self._jacobian if self.usejac else None
        flag, t, y = self._impl.trbdf2_py_solve(self._memory, self.f, jac, self._solout, t, np.array(y, dtype=float), tf,
                                                float(self.options["inith"]))

        st = self._memory.get_stats()
        self.statistics["nsteps"] += st["naccpt"]
        self.statistics["nerrfails"] += st["nreject"]
        self.statistics["nfcns"] += st["nfcn"]
        self.statistics["nfcnjacs"] += st["nfcnjac"]
        self.statistics["njacs"] += st["njac"]
        self.statistics["nlus"] += st["nlu"]
        self.statistics["nniters"] += st["nnewton"]
        self.statistics["nnfails"] += st["nnfail"]
        self.statistics["nrhsfails"] += st["nrhsfail"]
        self.statistics["njacfails"] += st["njacfail"]
        self.statistics["njacstale"] += st["njacstale"]
        self.statistics["nnfails_div"] += st["nnfail_div"]
        self.statistics["nnfails_iter"] += st["nnfail_iter"]
        self._memory.reset_stats()    # the counters were added to self.statistics

        if flag >= 0 and opts["output_list"] is not None and (not self._tlist or self._tlist[-1] != t):
            self._tlist.append(t)        # end point, so that a time event at tf is reached
            self._ylist.append(y)

        if flag == 0:
            flag = ID_PY_COMPLETE
        elif flag == 1:
            flag = ID_PY_EVENT
        else:
            msg = self._memory.get_err_msg()
            if isinstance(self._py_err, BaseException):
                raise self._py_err from None
            raise TRBDF2Error(flag, t, msg) from None
        return flag, self._tlist, self._ylist

    def state_event_info(self):
        return self._event_info

    def set_event_info(self, event_info):
        self._event_info = event_info

    def print_statistics(self, verbose=NORMAL):
        Explicit_ODE.print_statistics(self, verbose)
        log = lambda msg: self.log_message(msg, verbose)
        log("\nSolver options:\n")
        log(" Solver                  : TR-BDF2")
        log(" Tolerances (absolute)   : " + str(self._compact_tol(self.options["atol"])))
        log(" Tolerances (relative)   : " + str(self.options["rtol"]))
        log("")

    # ------------------------------------------------------------------ options as properties
    def _set_atol(self, atol):
        self.options["atol"] = np.array(atol, dtype=float) * np.ones(self._leny)
        if self.options["atol"].min() < 0:
            raise TRBDF2Error(-1, err_msg="The absolute tolerance must be non-negative.")

    def _get_atol(self):
        return self.options["atol"]
    atol = property(_get_atol, _set_atol)

    def _set_rtol(self, rtol):
        if np.ndim(rtol) == 0:
            if float(rtol) <= 0.0:
                raise TRBDF2Error(-1, err_msg="The relative tolerance must be positive.")
            self.options["rtol"] = float(rtol)
        else:
            self.options["rtol"] = np.array(rtol, dtype=float)
    def _get_rtol(self):
        return self.options["rtol"]
    rtol = property(_get_rtol, _set_rtol)

    def _set_maxh(self, maxh):
        self.options["maxh"] = None if maxh is None or float(maxh) <= 0.0 else float(maxh)
    def _get_maxh(self):
        return self.options["maxh"]
    maxh = property(_get_maxh, _set_maxh)

    def _set_inith(self, inith):
        self.options["inith"] = float(inith)
    def _get_inith(self):
        return self.options["inith"]
    inith = property(_get_inith, _set_inith)

    def _set_usejac(self, jac):
        self.options["usejac"] = bool(jac)
    def _get_usejac(self):
        return self.options["usejac"]
    usejac = property(_get_usejac, _set_usejac)

    def _set_maxsteps(self, n):
        self.options["maxsteps"] = int(n)
    def _get_maxsteps(self):
        return self.options["maxsteps"]
    maxsteps = property(_get_maxsteps, _set_maxsteps)

    def _set_newt(self, n):
        self.options["newt"] = int(n)
    def _get_newt(self):
        return self.options["newt"]
    newt = property(_get_newt, _set_newt)

    def _set_maxsteps_jac(self, n):
        self.options["maxsteps_jac"] = int(n)
    def _get_maxsteps_jac(self):
        return self.options["maxsteps_jac"]
    maxsteps_jac = property(_get_maxsteps_jac, _set_maxsteps_jac)

    def _set_safe(self, v):
        self.options["safe"] = float(v)
    def _get_safe(self):
        return self.options["safe"]
    safe = property(_get_safe, _set_safe)

    def _set_newton_tol(self, v):
        self.options["newton_tol"] = float(v)
    def _get_newton_tol(self):
        return self.options["newton_tol"]
    newton_tol = property(_get_newton_tol, _set_newton_tol)

    def _set_lu_band(self, band):
        lo, hi = band
        self.options["lu_band"] = (float(lo), float(hi))
    def _get_lu_band(self):
        return self.options["lu_band"]
    lu_band = property(_get_lu_band, _set_lu_band)

    def _set_fac1(self, v):
        self.options["fac1"] = float(v)
    def _get_fac1(self):
        return self.options["fac1"]
    fac1 = property(_get_fac1, _set_fac1)

    def _set_fac2(self, v):
        self.options["fac2"] = float(v)
    def _get_fac2(self):
        return self.options["fac2"]
    fac2 = property(_get_fac2, _set_fac2)

    def _set_fail_factor(self, v):
        self.options["fail_factor"] = float(v)
    def _get_fail_factor(self):
        return self.options["fail_factor"]
    fail_factor = property(_get_fail_factor, _set_fail_factor)

    def _set_fail_max(self, n):
        self.options["fail_max"] = int(n)
    def _get_fail_max(self):
        return self.options["fail_max"]
    fail_max = property(_get_fail_max, _set_fail_max)
