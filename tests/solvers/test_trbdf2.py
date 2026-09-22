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
import pytest

from assimulo.problem import Explicit_Problem
from assimulo.solvers import Radau5ODE
from assimulo.solvers.trbdf2 import TRBDF2, TRBDF2Error
from assimulo.exception import AssimuloRecoverableError


def vanderpol(mu=1e6):
    def f(t, y):
        return np.array([y[1], mu * ((1. - y[0]**2) * y[1] - y[0])])

    def jac(t, y):
        return np.array([[0., 1.], [mu * (-2. * y[0] * y[1] - 1.), mu * (1. - y[0]**2)]])
    mod = Explicit_Problem(f, [2.0, -0.6])
    mod.jac = jac
    return mod


class Test_TRBDF2:

    def test_prothero_robinson(self):
        """Stiff linear problem with a known solution: the error follows the tolerance."""
        lam = -1000.0
        mod = Explicit_Problem(lambda t, y: np.array([lam * (y[0] - np.cos(t)) - np.sin(t)]), [1.0])
        mod.jac = lambda t, y: np.array([[lam]])
        errs = {}
        for tol in (1e-4, 1e-6, 1e-8):
            sim = TRBDF2(mod)
            sim.rtol = tol
            sim.atol = 1e-2 * tol
            sim.verbosity = 50
            t, y = sim.simulate(10.0)
            errs[tol] = abs(y[-1][0] - np.cos(10.0))
            assert sim.statistics["nnfails"] <= 2           # at most the first step's 1e4 growth overshooting the LU band
            assert sim.statistics["njacs"] <= 40          # Jacobian reuse: not one per step
        assert errs[1e-4] < 1e-4 and errs[1e-6] < 1e-6 and errs[1e-8] < 1e-8
        assert errs[1e-8] < errs[1e-6] < errs[1e-4]

    def test_vanderpol_against_radau5(self):
        mod = vanderpol()
        ref = Radau5ODE(mod)
        ref.rtol = 1e-10
        ref.atol = 1e-12
        ref.verbosity = 50
        _, yref = ref.simulate(2.0)
        sim = TRBDF2(mod)
        sim.rtol = 1e-6
        sim.atol = 1e-8
        sim.verbosity = 50
        _, y = sim.simulate(2.0)
        assert y[-1][0] == pytest.approx(yref[-1][0], abs=1e-4)
        assert sim.statistics["nfcnjacs"] == 0            # the user Jacobian was used
        assert sim.statistics["nsteps"] < 5000
        # stage derivatives come from the stage equations: no rhs call beyond the Newton
        # iterations except the start of the integration and the initial-step estimate
        assert sim.statistics["nfcns"] <= sim.statistics["nniters"] + 2

    def test_no_jac(self):
        mod = vanderpol(mu=100.0)
        del mod.jac
        sim = TRBDF2(mod)
        sim.verbosity = 50
        assert not sim.usejac
        sim.simulate(1.0)
        assert sim.statistics["nfcnjacs"] == 3 * sim.statistics["njacs"]   # forward differences: base point + one per column

    def test_interpolate_and_ncp(self):
        """Dense output: communication points are interpolated, not stepped to."""
        mod = Explicit_Problem(lambda t, y: np.array([-y[0], -100. * (y[1] - np.cos(t))]), [1.0, 1.0])
        sim = TRBDF2(mod)
        sim.rtol = 1e-8
        sim.atol = 1e-10
        sim.verbosity = 50
        t, y = sim.simulate(3.0, 300)
        assert len(t) == 301
        assert np.abs(np.array(y)[:, 0] - np.exp(-np.array(t))).max() < 1e-6
        sim.reset()
        sim.report_continuously = True
        t2, y2 = sim.simulate(3.0, 300)
        assert len(t2) == 301
        assert np.abs(np.array(y2)[:, 0] - np.exp(-np.array(t2))).max() < 1e-6

    def test_event_localizer(self, extended_problem):
        sim = TRBDF2(extended_problem)
        sim.verbosity = 50
        sim.report_continuously = True
        t, y = sim.simulate(10.0, 1000)
        assert y[-1][0] == pytest.approx(8.0)
        assert y[-1][1] == pytest.approx(3.0)
        assert y[-1][2] == pytest.approx(2.0)

    def test_time_event(self):
        events = [1.0, 2.0, 2.5, 3.0]
        seen = []

        def time_events(t, y, sw):
            for ev in events:
                if t < ev:
                    return ev
            return None

        def handle_event(solver, event_info):
            solver.y += 1.0
            seen.append(solver.t)
            assert event_info[0] == []
            assert event_info[1]
        mod = Explicit_Problem(lambda t, y: np.array([1.0]), [0.0])
        mod.time_events = time_events
        mod.handle_event = handle_event
        sim = TRBDF2(mod)
        sim.verbosity = 50
        t, y = sim.simulate(5.0, 100)
        assert [pytest.approx(e) for e in events] == seen
        assert y[-1][0] == pytest.approx(5.0 + 4.0)

    def test_switches(self):
        state_events = lambda t, x, sw: np.array([x[0] - 1.])

        def handle_event(solver, event_info):
            solver.sw = [False]
        mod = Explicit_Problem(lambda t, x, sw: np.array([1.0]), [0.0])
        mod.sw0 = [True]
        mod.state_events = state_events
        mod.handle_event = handle_event
        sim = TRBDF2(mod)
        sim.verbosity = 50
        assert sim.sw[0]
        sim.simulate(3)
        assert not sim.sw[0]
        assert sim.statistics["nstateevents"] == 1

    def test_recoverable_rhs_failure(self):
        """A trial point the model refuses is retried with a smaller step, like CVode."""
        def f(t, y):
            if y[0] > 1.05:          # only trial points beyond the solution's range
                raise AssimuloRecoverableError("out of range")
            return np.array([1.0 - y[0]])
        mod = Explicit_Problem(f, [0.0])
        sim = TRBDF2(mod)
        sim.verbosity = 50
        sim.inith = 5.0             # a first step that overshoots
        t, y = sim.simulate(3.0)
        assert y[-1][0] == pytest.approx(1.0 - np.exp(-3.0), abs=1e-4)
        assert sim.statistics["nrhsfails"] >= 1

    def test_maxsteps(self):
        sim = TRBDF2(vanderpol())
        sim.verbosity = 50
        sim.maxsteps = 5
        with pytest.raises(TRBDF2Error, match="maximum number of steps"):
            sim.simulate(2.0)

    def test_options_validation(self):
        sim = TRBDF2(vanderpol())
        with pytest.raises(TRBDF2Error):
            sim.rtol = -1.0
        with pytest.raises(TRBDF2Error):
            sim.atol = [-1.0, 1.0]
        sim.maxh = 0.0
        assert sim.maxh is None
        sim.maxh = 0.1
        assert sim.maxh == 0.1
        sim.atol = 1e-8
        sim.newton_tol = 0.05
        sim.lu_band = [0.9, 1.1]
        sim.fac2 = 10.0
        assert sim.newton_tol == 0.05 and sim.lu_band == (0.9, 1.1) and sim.fac2 == 10.0
        sim.verbosity = 50
        sim.simulate(0.1)                # the properties reach the C stepper without an error
        with pytest.raises(TRBDF2Error, match="Invalid option"):
            sim.lu_band = (1.5, 2.0)     # lower bound above 1
            sim.simulate(0.1)
