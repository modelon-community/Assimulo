#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2010-2024 Modelon AB
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
import scipy.linalg as spl

from assimulo.explicit_ode import Explicit_ODE
from assimulo.implicit_ode import Implicit_ODE
from assimulo.lib.radau_core import Radau_Common, Radau_Exception
from assimulo.exception import (
    AssimuloException,
    Explicit_ODE_Exception,
    Implicit_ODE_Exception,
    AssimuloRecoverableError,
)
from assimulo.ode import NORMAL, LOUD, SCREAM, ID_PY_OK, ID_PY_COMPLETE, ID_PY_EVENT


class Radau5Error(AssimuloException):
    """
    Defines the Radau5Error and provides the textual error message.
    """
    msg = { -1    : 'The input is not consistent.',
            -2    : 'The solver took max internal steps but could not reach the next output time.',
            -3    : 'The step size became too small.',
            -4    : 'The matrix is repeatedly singular.',
            -5    : 'Repeated unexpected step rejections.',
            -10   : 'Unrecoverable exception encountered during callback to problem (right-hand side/jacobian).'
            }

    def __init__(self, value=None, t=0.0, err_msg=None):
        self.value = value
        self.t = t
        self.err_msg = err_msg

    def __str__(self):
        if self.err_msg:
            return repr('Radau5 failed with flag %s. At time %f. Message: %s' % (self.value, self.t, self.err_msg))
        else:
            try:
                return repr(self.msg[self.value] + ' At time %f.' % self.t)
            except KeyError:
                return repr('Radau failed with flag %s. At time %f.' % (self.value, self.t))


class Radau5DAE(Radau_Common, Implicit_ODE):
    """
    Radau IIA fifth-order three-stages with step-size control and
    continuous output. Based on the FORTRAN code RADAU5 by E.Hairer and
    G.Wanner, which can be found here:
    http://www.unige.ch/~hairer/software.html

    Details about the implementation (FORTRAN) can be found in the book::

        Solving Ordinary Differential Equations II,
        Stiff and Differential-Algebraic Problems

        Authors: E. Hairer and G. Wanner
        Springer-Verlag, ISBN: 3-540-60452-9
    """

    def __init__(self, problem):
        Implicit_ODE.__init__(self, problem)

        self.options["inith"]    = 0.01
        self.options["newt"]     = 7
        self.options["thet"]     = 1.e-3
        self.options["fnewt"]    = 0.0
        self.options["quot1"]    = 1.0
        self.options["quot2"]    = 1.2
        self.options["fac1"]     = 0.2
        self.options["fac2"]     = 8.0
        self.options["maxh"]     = np.inf
        self.options["safe"]     = 0.9
        self.options["atol"]     = 1.0e-6 * np.ones(self.problem_info["dim"])
        self.options["rtol"]     = 1.0e-6
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 100000

        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = True
        self.supports["state_events"]        = True

        self._leny = len(self.y)
        self._type = '(implicit)'
        self._event_info = None

    def _get_linear_solver(self):
        return 'DENSE'

    def _set_linear_solver(self, linear_solver):
        raise Radau_Exception(
            "Radau5DAE does not support setting the 'linear_solver' attribute, "
            "since it only supports the DENSE linear solver in Fortran implementation of Radau5.")

    linear_solver = property(_get_linear_solver, _set_linear_solver)

    def initialize(self):
        self.statistics.reset()
        try:
            from assimulo.lib import radau5 as radau5_f
            self.radau5 = radau5_f
        except Exception:
            raise Radau_Exception("Failed to import the Fortran based Radau5 solver implementation.")

    def set_problem_data(self):
        if self.problem_info["state_events"]:
            if self.problem_info["type"] == 1:
                def event_func(t, y, yd):
                    return self.problem.state_events(t, y, yd, self.sw)
            else:
                def event_func(t, y, yd):
                    return self.problem.state_events(t, y, self.sw)

            def f(t, y):
                ret = 0
                try:
                    leny = self._leny
                    res = self.problem.res(t, y[:leny], y[leny:2*leny], self.sw)
                except BaseException as err:
                    res = y[:self._leny].copy()
                    if isinstance(err, (np.linalg.LinAlgError, ZeroDivisionError, AssimuloRecoverableError)):
                        ret = -1
                    else:
                        ret = -2
                return np.append(y[self._leny:2*self._leny], res), [ret]

            self._f = f
            self.event_func = event_func
            self._event_info = [0] * self.problem_info["dimRoot"]
            self.g_old = self.event_func(self.t, self.y, self.yd)
            self.statistics["nstatefcns"] += 1
        else:
            def f(t, y):
                ret = 0
                try:
                    leny = self._leny
                    res = self.problem.res(t, y[:leny], y[leny:2*leny])
                except BaseException as err:
                    res = y[:self._leny].copy()
                    if isinstance(err, (np.linalg.LinAlgError, ZeroDivisionError, AssimuloRecoverableError)):
                        ret = -1
                    else:
                        ret = -2
                return np.append(y[self._leny:2*self._leny], res), [ret]
            self._f = f

    def interpolate(self, time, k=0):
        y = np.empty(self._leny * 2)
        for i in range(self._leny * 2):
            y[i] = self.radau5.contr5(i + 1, time, self.cont)
        if k == 0:
            return y[:self._leny]
        elif k == 1:
            return y[self._leny:2*self._leny]

    def _solout(self, nrsol, told, t, y, cont, werr, lrc, irtrn):
        """Called after every successful step taken by Radau5."""
        self.cont = cont

        yd = y[self._leny:2*self._leny].copy()
        y  = y[:self._leny].copy()

        if self.problem_info["state_events"]:
            flag, t, y, yd = self.event_locator(told, t, y, yd)
            if flag == ID_PY_EVENT:
                irtrn = -1

        if self._opts["report_continuously"]:
            initialize_flag = self.report_solution(t, y, yd, self._opts)
            if initialize_flag:
                irtrn = -1
        else:
            if self._opts["output_list"] is None:
                self._tlist.append(t)
                self._ylist.append(y)
                self._ydlist.append(yd)
            else:
                output_list  = self._opts["output_list"]
                output_index = self._opts["output_index"]
                try:
                    while output_list[output_index] <= t:
                        self._tlist.append(output_list[output_index])
                        self._ylist.append(self.interpolate(output_list[output_index]))
                        self._ydlist.append(self.interpolate(output_list[output_index], 1))
                        output_index += 1
                except IndexError:
                    pass
                self._opts["output_index"] = output_index

                if (self.problem_info["state_events"] and flag == ID_PY_EVENT
                        and len(self._tlist) > 0 and self._tlist[-1] != t):
                    self._tlist.append(t)
                    self._ylist.append(y)
                    self._ydlist.append(yd)

        return irtrn

    def _mas_f(self, am):
        return self._mass_matrix

    def integrate(self, t, y, yd, tf, opts):
        if self.usejac:
            self.usejac = False
            self.log_message("Jacobians are not currently supported, disabling.", NORMAL)

        ITOL  = 1
        IJAC  = 1 if self.usejac else 0
        MLJAC = self.problem_info["dim"] * 2
        MUJAC = 0
        IMAS  = 1
        MLMAS = 0
        MUMAS = 0
        IOUT  = 1
        WORK  = np.array([0.0] * (5 * ((self.problem_info["dim"] * 2)**2 + 12) + 20))
        IWORK = np.array([0]   * (3 *  (self.problem_info["dim"] * 2)           + 20), dtype=np.intc)

        WORK[1] = self.safe
        WORK[2] = self.thet
        WORK[3] = self.fnewt
        WORK[4] = self.quot1
        WORK[5] = self.quot2
        WORK[6] = self.maxh
        WORK[7] = self.fac1
        WORK[8] = self.fac2

        IWORK[1] = self.maxsteps
        IWORK[2] = self.newt
        IWORK[4] = self._leny
        IWORK[5] = self._leny
        IWORK[8] = self._leny
        IWORK[9] = self._leny

        mas_dummy = lambda t: t
        jac_dummy = (lambda t: t) if not self.usejac else self.problem.jac

        if opts["initialize"]:
            self.set_problem_data()
            self._tlist  = []
            self._ylist  = []
            self._ydlist = []

        self._opts = opts
        y = np.append(y, yd)
        self._mass_matrix = np.array([[0] * self._leny])
        atol = np.append(self.atol, self.atol)

        t, y, h, iwork, flag = self.radau5.radau5(
            self._f, t, y.copy(), tf, self.inith,
            self.rtol * np.ones(self.problem_info["dim"] * 2), atol,
            ITOL, jac_dummy, IJAC, MLJAC, MUJAC,
            self._mas_f, IMAS, MLMAS, MUMAS,
            self._solout, IOUT, WORK, IWORK)

        if (flag >= 0 and opts["output_list"] is not None
                and (not self._tlist or self._tlist[-1] != t)):
            self._tlist.append(t)
            self._ylist.append(y[:self._leny].copy())
            self._ydlist.append(y[self._leny:2*self._leny].copy())

        if flag == 1:
            flag = ID_PY_COMPLETE
        elif flag == 2:
            flag = ID_PY_EVENT
        else:
            raise Radau5Error(flag, t)

        self.statistics["nsteps"]    += iwork[16]
        self.statistics["nfcns"]     += iwork[13]
        self.statistics["njacs"]     += iwork[14]
        self.statistics["nfcnjacs"]  += (iwork[14] * self.problem_info["dim"] if not self.usejac else 0)
        self.statistics["nerrfails"] += iwork[17]
        self.statistics["nlus"]      += iwork[18]

        return flag, self._tlist, self._ylist, self._ydlist

    def state_event_info(self):
        return self._event_info

    def set_event_info(self, event_info):
        self._event_info = event_info

    def print_statistics(self, verbose=NORMAL):
        """Prints the run-time statistics for the problem."""
        Implicit_ODE.print_statistics(self, verbose)

        log_message_verbose = lambda msg: self.log_message(msg, verbose)
        log_message_verbose('\nSolver options:\n')
        log_message_verbose(' Solver                  : Radau5' + self._type)
        log_message_verbose(' Tolerances (absolute)   : ' + str(self._compact_tol(self.options["atol"])))
        log_message_verbose(' Tolerances (relative)   : ' + str(self.options["rtol"]))
        log_message_verbose('')


class _Radau5ODE(Radau_Common, Explicit_ODE):
    """
    Radau IIA fifth-order three-stages with step-size control and continuous output.
    Based on the FORTRAN code by E.Hairer and G.Wanner, which can be found here:
    http://www.unige.ch/~hairer/software.html

    Details about the implementation (FORTRAN) can be found in the book::

        Solving Ordinary Differential Equations II,
        Stiff and Differential-Algebraic Problems

        Authors: E. Hairer and G. Wanner
        Springer-Verlag, ISBN: 3-540-60452-9

    This code is aimed at providing a Python implementation of the original code.
    """

    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem)

        self.options["inith"]    = 0.01
        self.options["newt"]     = 7
        self.options["thet"]     = 1.e-3
        self.options["fnewt"]    = 0
        self.options["quot1"]    = 1.0
        self.options["quot2"]    = 1.2
        self.options["fac1"]     = 0.2
        self.options["fac2"]     = 8.0
        self.options["maxh"]     = np.inf
        self.options["safe"]     = 0.9
        self.options["atol"]     = 1.0e-6
        self.options["rtol"]     = 1.0e-6
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 10000

        self._curjac   = False
        self._itfail   = False
        self._needjac  = True
        self._needLU   = True
        self._first    = True
        self._rejected = True
        self._leny     = len(self.y)
        self._oldh     = 0.0
        self._olderr   = 1.0
        self._eps      = np.finfo('double').eps
        self._col_poly = np.zeros(self._leny * 3)
        self._type     = '(explicit)'
        self._curiter  = 0

        self.f   = problem.rhs_internal
        self.Y1  = np.array([0.0] * len(self.y0))
        self.Y2  = np.array([0.0] * len(self.y0))
        self.Y3  = np.array([0.0] * len(self.y0))
        self._f0 = np.array([0.0] * len(self.y0))

        self.supports["one_step_mode"]      = True
        self.supports["interpolated_output"] = True

        self._load_parameters()

    def initialize(self):
        self.statistics.reset()

    def step_generator(self, t, y, tf, opts):

        if opts["initialize"]:
            self._oldh    = self.inith
            self.h        = self.inith
            self._fac_con = 1.0

        if self.fnewt == 0:
            self.fnewt = max(10. * self._eps / self.rtol, min(0.03, self.rtol**0.5))

        self.f(self._f0, t, y)
        self.statistics["nfcns"] += 1
        self._tc = t
        self._yc = y

        for i in range(self.maxsteps):
            if t < tf:
                t, y = self._step(t, y)
                self._tc = t
                self._yc = y

                if self.h > np.abs(tf - t):
                    self.h = np.abs(tf - t)

                if t < tf:
                    yield ID_PY_OK, t, y
                else:
                    yield ID_PY_COMPLETE, t, y
                    break

                self._first = False
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')

    def step(self, t, y, tf, opts):
        if opts["initialize"]:
            self._next_step = self.step_generator(t, y, tf, opts)
        return next(self._next_step)

    def integrate(self, t, y, tf, opts):
        if opts["output_list"] is not None:
            output_list  = opts["output_list"]
            output_index = opts["output_index"]
            next_step    = self.step_generator(t, y, tf, opts)
            tlist, ylist = [], []
            res = [ID_PY_OK]
            while res[0] != ID_PY_COMPLETE:
                res = next(next_step)
                try:
                    while output_list[output_index] <= res[1]:
                        tlist.append(output_list[output_index])
                        ylist.append(self.interpolate(output_list[output_index]))
                        output_index += 1
                except IndexError:
                    pass
            return res[0], tlist, ylist
        else:
            [flags, tlist, ylist] = list(zip(*list(self.step_generator(t, y, tf, opts))))
            return flags[-1], tlist, ylist

    def _step(self, t, y):
        self._scaling = np.array(abs(y) * self.rtol + self.atol)
        while True:
            self.newton(t, y)
            self._err = self.estimate_error()
            if self._err > 1.0:
                self._rejected = True
                self.statistics["nerrfails"] += 1
                ho     = self.h
                self.h = self.adjust_stepsize(self._err)
                self.log_message('Rejecting step at ' + str(t) + 'with old stepsize' + str(ho) + 'and new ' + str(self.h), SCREAM)
                if self._curjac or self._curiter == 1:
                    self._needjac = False
                    self._needLU  = True
                else:
                    self._needjac = True
                    self._needLU  = True
            else:
                self.log_message('Accepting step at ' + str(t) + 'with stepsize ' + str(self.h), SCREAM)
                self.statistics["nsteps"] += 1
                tn = t + self.h
                yn = y + self._Z[2*self._leny:3*self._leny]
                self.f(self._f0, tn, yn)
                self.statistics["nfcns"] += 1
                self._oldoldh = self._oldh
                self._oldh    = self.h
                self._oldt    = t
                self._newt    = tn
                ht     = self.adjust_stepsize(self._err, predict=True)
                self.h = min(self.h, ht) if self._rejected else ht
                self._rejected = False
                self._curjac   = False
                if self._oldoldh == self.h and (self._theta <= self.thet):
                    self._needjac = False
                    self._needLU  = False
                else:
                    if self._theta <= self.thet:
                        self._needjac = False
                        self._needLU  = True
                    else:
                        self._needjac = True
                        self._needLU  = True
                if self.thet < 0:
                    self._needjac = True
                    self._needLU  = True
                self._olderr = max(self._err, 1.e-2)
                break
        self._col_poly = self._collocation_pol(self._Z, self._col_poly, self._leny)
        return tn, yn

    def _collocation_pol(self, Z, col_poly, leny):
        col_poly[2*leny:3*leny] = Z[:leny] / self.C[0]
        col_poly[leny:2*leny]   = (Z[:leny] - Z[leny:2*leny]) / (self.C[0] - self.C[1])
        col_poly[:leny]         = (Z[leny:2*leny] - Z[2*leny:3*leny]) / (self.C[1] - 1.)
        col_poly[2*leny:3*leny] = (col_poly[leny:2*leny] - col_poly[2*leny:3*leny]) / self.C[1]
        col_poly[leny:2*leny]   = (col_poly[leny:2*leny] - col_poly[:leny]) / (self.C[0] - 1.)
        col_poly[2*leny:3*leny] = col_poly[leny:2*leny] - col_poly[2*leny:3*leny]
        return col_poly

    def _radau_F(self, Z, t, y):
        Z1 = Z[:self._leny]
        Z2 = Z[self._leny:2*self._leny]
        Z3 = Z[2*self._leny:3*self._leny]
        self.f(self.Y1, t + self.C[0] * self.h, y + Z1)
        self.f(self.Y2, t + self.C[1] * self.h, y + Z2)
        self.f(self.Y3, t + self.C[2] * self.h, y + Z3)
        self.statistics["nfcns"] += 3
        return np.hstack((np.hstack((self.Y1, self.Y2)), self.Y3))

    def calc_start_values(self):
        if self._first:
            Z = np.zeros(self._leny * 3)
            W = np.zeros(self._leny * 3)
        else:
            Z      = self._Z
            cq     = self.C * self.h / self._oldh
            newtval = self._col_poly
            leny   = self._leny
            Z[:leny]         = cq[0]*(newtval[:leny]+(cq[0]-self.C[1]+1.)*(newtval[leny:2*leny]+(cq[0]-self.C[0]+1.)*newtval[2*leny:3*leny]))
            Z[leny:2*leny]   = cq[1]*(newtval[:leny]+(cq[1]-self.C[1]+1.)*(newtval[leny:2*leny]+(cq[1]-self.C[0]+1.)*newtval[2*leny:3*leny]))
            Z[2*leny:3*leny] = cq[2]*(newtval[:leny]+(cq[2]-self.C[1]+1.)*(newtval[leny:2*leny]+(cq[2]-self.C[0]+1.)*newtval[2*leny:3*leny]))
            W = np.dot(self.T2, Z)
        return Z, W

    def newton(self, t, y):
        for k in range(20):
            self._curiter = 0
            self._fac_con = max(self._fac_con, self._eps)**0.8
            self._theta   = abs(self.thet)
            if self._needjac:
                self._jac = self.jacobian(t, y)
            if self._needLU:
                self.statistics["nlus"] += 1
                self._a = self._alpha / self.h
                self._b = self._beta  / self.h
                self._g = self._gamma / self.h
                self._B = self._g * self.I - self._jac
                self._P1, self._L1, self._U1 = spl.lu(self._B)
                self._P2, self._L2, self._U2 = spl.lu(self._a * self.I - self._jac)
                self._P3, self._L3, self._U3 = spl.lu(self._b * self.I - self._jac)
                self._needLU = False
                if min(abs(np.diag(self._U1))) < self._eps:
                    raise Explicit_ODE_Exception('Error, gI-J is singular.')
            Z, W = self.calc_start_values()
            for i in range(self.newt):
                self._curiter += 1
                self.statistics["nniters"] += 1
                Z = np.dot(self.T2, self._radau_F(Z.real, t, y))
                Z[:self._leny]               = Z[:self._leny]               - self._g * np.dot(self.I, W[:self._leny])
                Z[self._leny:2*self._leny]   = Z[self._leny:2*self._leny]   - self._a * np.dot(self.I, W[self._leny:2*self._leny])
                Z[2*self._leny:3*self._leny] = Z[2*self._leny:3*self._leny] - self._b * np.dot(self.I, W[2*self._leny:3*self._leny])
                Z[:self._leny]               = np.linalg.solve(self._U1, np.linalg.solve(self._L1, np.linalg.solve(self._P1, Z[:self._leny])))
                Z[self._leny:2*self._leny]   = np.linalg.solve(self._U2, np.linalg.solve(self._L2, np.linalg.solve(self._P2, Z[self._leny:2*self._leny])))
                Z[2*self._leny:3*self._leny] = np.linalg.solve(self._U3, np.linalg.solve(self._L3, np.linalg.solve(self._P3, Z[2*self._leny:3*self._leny])))
                newnrm = np.linalg.norm(Z.reshape(-1, self._leny) / self._scaling, 'fro') / np.sqrt(3. * self._leny)
                if i > 0:
                    thq = newnrm / oldnrm
                    if i == 1:
                        self._theta = thq
                    else:
                        self._theta = np.sqrt(thq * thqold)
                    thqold = thq
                    if self._theta < 0.99:
                        self._fac_con = self._theta / (1. - self._theta)
                        dyth = self._fac_con * newnrm * self._theta**(self.newt - (i+1) - 1) / self.fnewt
                        if dyth >= 1.0:
                            qnewt = max(1.e-4, min(20., dyth))
                            self.h = 0.8 * qnewt**(-1.0 / (4.0 + self.newt - (i+1) - 1)) * self.h
                            self._itfail   = True
                            self._rejected = True
                            break
                    else:
                        self._itfail = True
                        break
                oldnrm = max(newnrm, self._eps)
                W = W + Z
                Z = np.dot(self.T3, W)
                if self._fac_con * newnrm <= self.fnewt:
                    self._itfail = False
                    break
            else:
                self._itfail = True
            if not self._itfail:
                self._Z = Z.real
                break
            else:
                self.log_message('Iteration failed at time %e with step-size %e' % (t, self.h), SCREAM)
                self.statistics["nnfails"] += 1
                self._rejected = True
                if self._theta >= 0.99:
                    self.h = self.h / 2.0
                if self._curjac:
                    self._needjac = False
                    self._needLU  = True
                else:
                    self._needjac = True
                    self._needLU  = True
        else:
            raise Explicit_ODE_Exception('Newton iteration failed at time %e with step-size %e' % (t, self.h))

    def adjust_stepsize(self, err, predict=False):
        fac   = min(self.safe, self.safe * (2.*self.newt + 1.) / (2.*self.newt + self._curiter))
        quot  = max(1./self.fac2, min(1./self.fac1, (err**0.25) / fac))
        hnormal = self.h / quot
        if predict:
            if not self._first:
                facgus = (self._hacc / self.h) * (err**2 / self._olderr)**0.25 / self.safe
                facgus = max(1./self.fac2, min(1./self.fac1, facgus))
                quot   = max(quot, facgus)
                h      = self.h / quot
            else:
                h = hnormal
            self._hacc = self.h
        else:
            h = hnormal
        qt = h / self.h
        if (qt >= self.quot1) and (qt <= self.quot2):
            h = self.h
        if self._first and err >= 1.0:
            h = self.h / 10.
        if h < self._eps:
            raise Explicit_ODE_Exception('Step-size to small at %e with h = %e' % (self._tc, self.h))
        if h > self.maxh:
            h = self.maxh
        return h

    def estimate_error(self):
        temp  = 1./self.h * (self.E[0]*self._Z[:self._leny] + self.E[1]*self._Z[self._leny:2*self._leny] + self.E[2]*self._Z[2*self._leny:3*self._leny])
        scal  = self._scaling
        err_v = np.linalg.solve(self._U1, np.linalg.solve(self._L1, np.linalg.solve(self._P1, self._f0 + temp)))
        err   = np.linalg.norm(err_v / scal)
        err   = max(err / np.sqrt(self._leny), 1.e-10)
        if (self._rejected or self._first) and err >= 1.:
            self.statistics["nfcns"] += 1
            err_new = np.array([0.0] * self._leny)
            self.f(err_new, self._tc, self._yc + err_v)
            err_v = np.linalg.solve(self._U1, np.linalg.solve(self._L1, np.linalg.solve(self._P1, err_new + temp)))
            err   = np.linalg.norm(err_v / scal)
            err   = max(err / np.sqrt(self._leny), 1.e-10)
        return err

    def jacobian(self, t, y):
        self._curjac  = True
        self._needLU  = True
        self._needjac = False
        if self.usejac:
            cjac = self.problem.jac(t, y)
        else:
            delt  = np.array([(self._eps * max(abs(yi), 1.e-5))**0.5 for yi in y]) * np.identity(self._leny)
            Fdelt = np.array([self.problem.rhs(t, y + e) for e in delt])
            grad  = ((Fdelt - self.problem.rhs(t, y)).T / delt.diagonal()).T
            cjac  = np.array(grad).T
            self.statistics["nfcnjacs"] += 1 + self._leny
        self.statistics["njacs"] += 1
        return cjac

    def interpolate(self, t, k=0):
        leny = self._leny
        s    = (t - self._newt) / self._oldh
        Z    = self._col_poly
        return self._yc + s * (Z[:leny] + (s - self.C[1] + 1.) * (Z[leny:2*leny] + (s - self.C[0] + 1.) * Z[2*leny:3*leny]))

    def _load_parameters(self):
        A = np.zeros([3, 3])
        A[0, 0] = (88. - 7.*np.sqrt(6.)) / 360.0
        A[0, 1] = (296. - 169.*np.sqrt(6.)) / 1800.0
        A[0, 2] = (-2.0 + 3.0*np.sqrt(6.)) / 225.0
        A[1, 0] = (296.0 + 169.0*np.sqrt(6.)) / 1800.0
        A[1, 1] = (88. + 7.*np.sqrt(6.)) / 360.0
        A[1, 2] = (-2. - 3.*np.sqrt(6.)) / 225.0
        A[2, 0] = (16.0 - np.sqrt(6.)) / 36.0
        A[2, 1] = (16.0 + np.sqrt(6.)) / 36.0
        A[2, 2] = 1.0 / 9.0
        C = np.zeros([3])
        C[0] = (4.0 - np.sqrt(6.0)) / 10.0
        C[1] = (4.0 + np.sqrt(6.0)) / 10.0
        C[2] = 1.0
        B = np.zeros([3])
        B[0] = (16.0 - np.sqrt(6.0)) / 36.0
        B[1] = (16.0 + np.sqrt(6.0)) / 36.0
        B[2] = 1.0 / 9.0
        E    = np.zeros(3)
        E[0] = -13.0 - 7.*np.sqrt(6.)
        E[1] = -13.0 + 7.0*np.sqrt(6.)
        E[2] = -1.0
        E    = 1.0 / 3.0 * E
        Ainv       = np.linalg.inv(A)
        eig, T     = np.linalg.eig(Ainv)
        eig        = np.array([eig[2], eig[0], eig[1]])
        self._alpha = eig[1]
        self._beta  = eig[2]
        self._gamma = eig[0].real
        temp0 = T[:, 0].copy(); temp1 = T[:, 1].copy(); temp2 = T[:, 2].copy()
        T[:, 0] = temp2; T[:, 1] = temp0; T[:, 2] = temp1
        Tinv = np.linalg.inv(T)
        I  = np.eye(self._leny)
        T1 = np.kron(np.diag(eig), I)
        T2 = np.kron(Tinv, I)
        T3 = np.kron(T, I)
        self.A  = A; self.B = B; self.C = C; self.I = I; self.E = E
        self.T1 = T1; self.T2 = T2; self.T3 = T3
        self.I3 = np.eye(3); self.EIG = eig


class _Radau5DAE(Radau_Common, Implicit_ODE):
    """
    Radau IIA fifth-order three-stages with step-size control and continuous output.
    Based on the FORTRAN code by E.Hairer and G.Wanner, which can be found here:
    http://www.unige.ch/~hairer/software.html

    Details about the implementation (FORTRAN) can be found in the book::

        Solving Ordinary Differential Equations II,
        Stiff and Differential-Algebraic Problems

        Authors: E. Hairer and G. Wanner
        Springer-Verlag, ISBN: 3-540-60452-9

    This code is aimed at providing a Python implementation of the original code.
    """

    def __init__(self, problem):
        Implicit_ODE.__init__(self, problem)

        self._leny  = len(self.y)
        self._2leny = 2 * self._leny

        self.options["inith"]    = 0.01
        self.options["newt"]     = 7
        self.options["thet"]     = 1.e-3
        self.options["fnewt"]    = 0
        self.options["quot1"]    = 1.0
        self.options["quot2"]    = 1.2
        self.options["fac1"]     = 0.2
        self.options["fac2"]     = 8.0
        self.options["maxh"]     = np.inf
        self.options["safe"]     = 0.9
        self.options["atol"]     = np.array([1.0e-6] * self._leny)
        self.options["rtol"]     = 1.0e-6
        self.options["index"]    = np.array([1] * self._leny + [2] * self._leny)
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 10000

        self._curjac   = False
        self._itfail   = False
        self._needjac  = True
        self._needLU   = True
        self._first    = True
        self._rejected = True
        self._oldh     = 0.0
        self._olderr   = 1.0
        self._eps      = np.finfo('double').eps
        self._col_poly = np.zeros(self._2leny * 3)
        self._type     = '(implicit)'
        self._curiter  = 0

        self.f   = problem.res_internal
        self.RES = np.array([0.0] * len(self.y0))
        self.Y1  = np.array([0.0] * len(self.y0))
        self.Y2  = np.array([0.0] * len(self.y0))
        self.Y3  = np.array([0.0] * len(self.y0))
        self._f0 = np.array([0.0] * len(self.y0))

        self.supports["one_step_mode"]       = True
        self.supports["interpolated_output"] = True

        self._load_parameters()

    def _set_index(self, index):
        if len(index) == self._2leny:
            ind = np.array(index)
        elif len(index) == self._leny:
            ind = np.array(index + (np.array(index) + 1).tolist())
        else:
            raise Implicit_ODE_Exception('Wrong number of variables in the index vector.')
        self.options["index"] = ind

    def _get_index(self):
        return self.options["index"]

    index = property(_get_index, _set_index)

    def initialize(self):
        self.statistics.reset()

    def step_generator(self, t, y, yd, tf, opts):
        if opts["initialize"]:
            self._oldh    = self.inith
            self.h        = self.inith
            self._fac_con = 1.0
        if self.fnewt == 0:
            self.fnewt = max(10.*self._eps / self.rtol, min(0.03, self.rtol**0.5))
        self._f0 = self._ode_f(t, np.append(y, yd))
        self.statistics["nfcns"] += 1
        self._tc  = t
        self._yc  = y
        self._ydc = yd
        for i in range(self.maxsteps):
            if t < tf:
                t, y, yd = self._step(t, y, yd)
                self._tc  = t
                self._yc  = y
                self._ydc = yd
                if self.h > np.abs(tf - t):
                    self.h = np.abs(tf - t)
                if t < tf:
                    yield ID_PY_OK, t, y, yd
                else:
                    yield ID_PY_COMPLETE, t, y, yd
                    break
                self._first = False
        else:
            raise Implicit_ODE_Exception('Final time not reached within maximum number of steps')

    def step(self, t, y, yd, tf, opts):
        if opts["initialize"]:
            self._next_step = self.step_generator(t, y, yd, tf, opts)
        return next(self._next_step)

    def integrate(self, t, y, yd, tf, opts):
        if opts["output_list"] is not None:
            output_list  = opts["output_list"]
            output_index = opts["output_index"]
            next_step    = self.step_generator(t, y, yd, tf, opts)
            tlist, ylist, ydlist = [], [], []
            res = [ID_PY_OK]
            while res[0] != ID_PY_COMPLETE:
                res = next(next_step)
                try:
                    while output_list[output_index] <= res[1]:
                        tlist.append(output_list[output_index])
                        ylist.append(self.interpolate(output_list[output_index]))
                        ydlist.append(self.interpolate(output_list[output_index], k=1))
                        output_index += 1
                except IndexError:
                    pass
            return res[0], tlist, ylist, ydlist
        else:
            [flags, tlist, ylist, ydlist] = list(zip(*list(self.step_generator(t, y, yd, tf, opts))))
            return flags[-1], tlist, ylist, ydlist

    def _ode_f(self, t, y):
        self.f(self.RES, t, y[:self._leny], y[self._leny:])
        return np.hstack((y[self._leny:], self.RES))

    def _radau_F(self, Z, t, y, yd):
        Z1  = Z[:self._2leny]
        Z2  = Z[self._2leny:2*self._2leny]
        Z3  = Z[2*self._2leny:3*self._2leny]
        q   = np.append(y, yd)
        sol1 = self._ode_f(t + self.C[0]*self.h, q + Z1)
        sol2 = self._ode_f(t + self.C[1]*self.h, q + Z2)
        sol3 = self._ode_f(t + self.C[2]*self.h, q + Z3)
        self.statistics["nfcns"] += 3
        return np.hstack((np.hstack((sol1, sol2)), sol3))

    def _step(self, t, y, yd):
        self._scaling = np.array(abs(np.append(y, yd)) * self.rtol + self.atol.tolist() * 2)
        while True:
            self.newton(t, y, yd)
            self._err = self.estimate_error()
            if self._err > 1.0:
                self._rejected = True
                self.statistics["nerrfails"] += 1
                ho     = self.h
                self.h = self.adjust_stepsize(self._err)
                self.log_message('Rejecting step at ' + str(t) + 'with old stepsize' + str(ho) + 'and new ' + str(self.h) + '. Error: ' + str(self._err), SCREAM)
                if self._curjac or self._curiter == 1:
                    self._needjac = False; self._needLU = True
                else:
                    self._needjac = True;  self._needLU = True
            else:
                self.log_message("Accepting step at " + str(t) + ' with stepsize ' + str(self.h) + '. Error: ' + str(self._err), SCREAM)
                self.statistics["nsteps"] += 1
                tn  = t + self.h
                yn  = y  + self._Z[2*self._2leny:3*self._2leny][:self._leny]
                ydn = yd + self._Z[2*self._2leny:3*self._2leny][self._leny:]
                self._f0      = self._ode_f(t, np.append(yn, ydn))
                self.statistics["nfcns"] += 1
                self._oldoldh = self._oldh
                self._oldh    = self.h
                self._oldt    = t
                self._newt    = tn
                ht     = self.adjust_stepsize(self._err, predict=True)
                self.h = min(self.h, ht) if self._rejected else ht
                self._rejected = False
                self._curjac   = False
                if self._oldoldh == self.h and (self._theta <= self.thet or self._curiter == 1):
                    self._needjac = False; self._needLU = False
                else:
                    if self._theta <= self.thet or self._curiter == 1:
                        self._needjac = False; self._needLU = True
                    else:
                        self._needjac = True;  self._needLU = True
                if self.thet < 0:
                    self._needjac = True; self._needLU = True
                self._olderr = max(self._err, 1.e-2)
                break
        self._col_poly = self._collocation_pol(self._Z, self._col_poly, self._2leny)
        return tn, yn, ydn

    def newton(self, t, y, yd):
        for k in range(20):
            self._curiter = 0
            self._fac_con = max(self._fac_con, self._eps)**0.8
            self._theta   = abs(self.thet)
            if self._needjac:
                self._jac = self.jacobian(t, y, yd)
            if self._needLU:
                self.statistics["nlus"] += 1
                self._a = self._alpha / self.h
                self._b = self._beta  / self.h
                self._g = self._gamma / self.h
                self._B = self._g * self.M - self._jac
                self._P1, self._L1, self._U1 = spl.lu(self._B)
                self._P2, self._L2, self._U2 = spl.lu(self._a * self.M - self._jac)
                self._P3, self._L3, self._U3 = spl.lu(self._b * self.M - self._jac)
                self._needLU = False
                if min(abs(np.diag(self._U1))) < self._eps:
                    raise Implicit_ODE_Exception('Error, gM-J is singular at ', self._tc)
            Z, W = self.calc_start_values()
            for i in range(self.newt):
                self._curiter += 1
                self.statistics["nniters"] += 1
                Z = np.dot(self.T2, self._radau_F(Z.real, t, y, yd))
                Z[:self._2leny]               = Z[:self._2leny]               - self._g * np.dot(self.M, W[:self._2leny])
                Z[self._2leny:2*self._2leny]  = Z[self._2leny:2*self._2leny]  - self._a * np.dot(self.M, W[self._2leny:2*self._2leny])
                Z[2*self._2leny:3*self._2leny] = Z[2*self._2leny:3*self._2leny] - self._b * np.dot(self.M, W[2*self._2leny:3*self._2leny])
                Z[:self._2leny]               = np.linalg.solve(self._U1, np.linalg.solve(self._L1, np.linalg.solve(self._P1, Z[:self._2leny])))
                Z[self._2leny:2*self._2leny]  = np.linalg.solve(self._U2, np.linalg.solve(self._L2, np.linalg.solve(self._P2, Z[self._2leny:2*self._2leny])))
                Z[2*self._2leny:3*self._2leny] = np.linalg.solve(self._U3, np.linalg.solve(self._L3, np.linalg.solve(self._P3, Z[2*self._2leny:3*self._2leny])))
                self._scaling = self._scaling / self.h**(self.index - 1)
                newnrm = np.linalg.norm(Z.reshape(-1, self._2leny) / self._scaling, 'fro') / np.sqrt(3. * self._2leny)
                if i > 0:
                    thq = newnrm / oldnrm
                    if i == 1:
                        self._theta = thq
                    else:
                        self._theta = np.sqrt(thq * thqold)
                    thqold = thq
                    if self._theta < 0.99:
                        self._fac_con = self._theta / (1. - self._theta)
                        dyth = self._fac_con * newnrm * self._theta**(self.newt - (i+1) - 1) / self.fnewt
                        if dyth >= 1.0:
                            qnewt = max(1.e-4, min(20., dyth))
                            self._hhfac = 0.8 * qnewt**(-1.0 / (4.0 + self.newt - (i+1) - 1))
                            self.h      = self._hhfac * self.h
                            self._itfail   = True
                            self._rejected = True
                            break
                    else:
                        self._itfail = True
                        break
                oldnrm = max(newnrm, self._eps)
                W = W + Z
                Z = np.dot(self.T3, W)
                if self._fac_con * newnrm <= self.fnewt:
                    self._itfail = False
                    break
            else:
                self._itfail = True
            if not self._itfail:
                self._Z = Z.real
                break
            else:
                self.log_message("Iteration failed at time %e with step-size %e" % (t, self.h), SCREAM)
                self.statistics["nnfails"] += 1
                self._rejected = True
                if self._theta >= 0.99:
                    self._hhfac = 0.5
                    self.h      = self.h * self._hhfac
                if self._curjac:
                    self._needjac = False; self._needLU = True
                else:
                    self._needjac = True;  self._needLU = True
        else:
            raise Implicit_ODE_Exception('Newton iteration failed at time %e with step-size %e' % (t, self.h))

    def estimate_error(self):
        temp  = 1./self.h * (self.E[0]*self._Z[:self._2leny] + self.E[1]*self._Z[self._2leny:2*self._2leny] + self.E[2]*self._Z[2*self._2leny:3*self._2leny])
        temp  = np.dot(self.M, temp)
        self._scaling = self._scaling / self.h**(self.index - 1)
        scal  = self._scaling
        err_v = np.linalg.solve(self._U1, np.linalg.solve(self._L1, np.linalg.solve(self._P1, self._f0 + temp)))
        err   = np.linalg.norm(err_v / scal)
        err   = max(err / np.sqrt(self._2leny), 1.e-10)
        if (self._rejected or self._first) and err >= 1.:
            self.statistics["nfcns"] += 1
            err_v = self._ode_f(self._tc, np.append(self._yc, self._ydc) + err_v)
            err_v = np.linalg.solve(self._U1, np.linalg.solve(self._L1, np.linalg.solve(self._P1, err_v + temp)))
            err   = np.linalg.norm(err_v / scal)
            err   = max(err / np.sqrt(self._2leny), 1.e-10)
        return err

    def interpolate(self, t, k=0):
        leny = self._2leny
        s    = (t - self._newt) / self._oldh
        Z    = self._col_poly
        diff = s * (Z[:leny] + (s - self.C[1] + 1.) * (Z[leny:2*leny] + (s - self.C[0] + 1.) * Z[2*leny:3*leny]))
        yout  = self._yc  + diff[:self._leny]
        ydout = self._ydc + diff[self._leny:]
        if k == 0:
            return yout
        elif k == 1:
            return ydout
        else:
            raise Implicit_ODE_Exception('Unknown value of k. Should be either 0 or 1')

    def jacobian(self, t, y, yd):
        self._curjac  = True
        self._needLU  = True
        self._needjac = False
        q = np.append(y, yd)
        if self.usejac:
            cjac = self.problem.jac(t, y, yd)
        else:
            delt  = np.array([(self._eps * max(abs(yi), 1.e-5))**0.5 for yi in q]) * np.identity(self._2leny)
            Fdelt = np.array([self._ode_f(t, q + e) for e in delt])
            grad  = ((Fdelt - self._ode_f(t, q)).T / delt.diagonal()).T
            cjac  = np.array(grad).T
            self.statistics["nfcnjacs"] += 1 + self._2leny
        self.statistics["njacs"] += 1
        return cjac

    def adjust_stepsize(self, err, predict=False):
        fac   = min(self.safe, self.safe * (2.*self.newt + 1.) / (2.*self.newt + self._curiter))
        quot  = max(1./self.fac2, min(1./self.fac1, (err**0.25) / fac))
        hnormal = self.h / quot
        if predict:
            if not self._first:
                facgus = (self._hacc / self.h) * (err**2 / self._olderr)**0.25 / self.safe
                facgus = max(1./self.fac2, min(1./self.fac1, facgus))
                quot   = max(quot, facgus)
                h      = self.h / quot
            else:
                h = hnormal
            self._hacc = self.h
        else:
            h = hnormal
        qt = h / self.h
        if (qt >= self.quot1) and (qt <= self.quot2):
            h = self.h
        if h > self.maxh:
            h = self.maxh
        if self._first and err >= 1.0:
            self._hhfac = 0.1
            h = self.h * self._hhfac
        else:
            self._hhfac = h / self.h
        if h < self._eps:
            raise Implicit_ODE_Exception('Step-size to small at %e with h = %e' % (self._tc, self.h))
        return h

    def _collocation_pol(self, Z, col_poly, leny):
        col_poly[2*leny:3*leny] = Z[:leny] / self.C[0]
        col_poly[leny:2*leny]   = (Z[:leny] - Z[leny:2*leny]) / (self.C[0] - self.C[1])
        col_poly[:leny]         = (Z[leny:2*leny] - Z[2*leny:3*leny]) / (self.C[1] - 1.)
        col_poly[2*leny:3*leny] = (col_poly[leny:2*leny] - col_poly[2*leny:3*leny]) / self.C[1]
        col_poly[leny:2*leny]   = (col_poly[leny:2*leny] - col_poly[:leny]) / (self.C[0] - 1.)
        col_poly[2*leny:3*leny] = col_poly[leny:2*leny] - col_poly[2*leny:3*leny]
        return col_poly

    def calc_start_values(self):
        if self._first:
            Z = np.zeros(self._2leny * 3)
            W = np.zeros(self._2leny * 3)
        else:
            Z      = self._Z
            cq     = self.C * self.h / self._oldh
            newtval = self._col_poly
            leny   = self._2leny
            Z[:leny]         = cq[0]*(newtval[:leny]+(cq[0]-self.C[1]+1.)*(newtval[leny:2*leny]+(cq[0]-self.C[0]+1.)*newtval[2*leny:3*leny]))
            Z[leny:2*leny]   = cq[1]*(newtval[:leny]+(cq[1]-self.C[1]+1.)*(newtval[leny:2*leny]+(cq[1]-self.C[0]+1.)*newtval[2*leny:3*leny]))
            Z[2*leny:3*leny] = cq[2]*(newtval[:leny]+(cq[2]-self.C[1]+1.)*(newtval[leny:2*leny]+(cq[2]-self.C[0]+1.)*newtval[2*leny:3*leny]))
            W = np.dot(self.T2, Z)
        return Z, W

    def _load_parameters(self):
        A = np.zeros([3, 3])
        A[0, 0] = (88. - 7.*np.sqrt(6.)) / 360.0
        A[0, 1] = (296. - 169.*np.sqrt(6.)) / 1800.0
        A[0, 2] = (-2.0 + 3.0*np.sqrt(6.)) / 225.0
        A[1, 0] = (296.0 + 169.0*np.sqrt(6.)) / 1800.0
        A[1, 1] = (88. + 7.*np.sqrt(6.)) / 360.0
        A[1, 2] = (-2. - 3.*np.sqrt(6.)) / 225.0
        A[2, 0] = (16.0 - np.sqrt(6.)) / 36.0
        A[2, 1] = (16.0 + np.sqrt(6.)) / 36.0
        A[2, 2] = 1.0 / 9.0
        C = np.zeros([3])
        C[0] = (4.0 - np.sqrt(6.0)) / 10.0
        C[1] = (4.0 + np.sqrt(6.0)) / 10.0
        C[2] = 1.0
        B = np.zeros([3])
        B[0] = (16.0 - np.sqrt(6.0)) / 36.0
        B[1] = (16.0 + np.sqrt(6.0)) / 36.0
        B[2] = 1.0 / 9.0
        E    = np.zeros(3)
        E[0] = -13.0 - 7.*np.sqrt(6.)
        E[1] = -13.0 + 7.0*np.sqrt(6.)
        E[2] = -1.0
        E    = 1.0 / 3.0 * E
        M    = np.array([[1., 0.], [0., 0.]])
        Ainv       = np.linalg.inv(A)
        eig, T     = np.linalg.eig(Ainv)
        eig        = np.array([eig[2], eig[0], eig[1]])
        self._alpha = eig[1]
        self._beta  = eig[2]
        self._gamma = eig[0].real
        temp0 = T[:, 0].copy(); temp1 = T[:, 1].copy(); temp2 = T[:, 2].copy()
        T[:, 0] = temp2; T[:, 1] = temp0; T[:, 2] = temp1
        Tinv = np.linalg.inv(T)
        I  = np.eye(self._2leny)
        M  = np.kron(M, np.eye(self._leny))
        T1 = np.kron(np.diag(eig), M)
        T2 = np.kron(Tinv, I)
        T3 = np.kron(T, I)
        self.A  = A; self.B = B; self.C = C; self.I = I; self.E = E; self.M = M
        self.T1 = T1; self.T2 = T2; self.T3 = T3
        self.I3 = np.eye(3); self.EIG = eig
