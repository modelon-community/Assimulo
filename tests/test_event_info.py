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

from assimulo import solvers
from assimulo.ode import QUIET
from assimulo.problem import Explicit_Problem, Implicit_Problem

TIME_EVENT = 1.0
STATE_EVENT_Y = 0.5
FINAL_TIME = 1.5
STATE_EVENT_ATOL = 1e-9


def never(solver):
    return False


def on_first_step(solver):
    return True


def at_time_event(solver):
    return np.isclose(solver.t, TIME_EVENT)


def at_state_event(solver):
    return np.isclose(solver.y[0], STATE_EVENT_Y, atol=STATE_EVENT_ATOL)


class Decay:
    """y = exp(-t) with a time event and one step event, when `signal` first holds."""

    def __init__(self, signal, with_state_event, **initial_values):
        super().__init__(y0=[1.0], t0=0.0, **initial_values)
        self.signal = signal
        self.signalled = False
        self.received = []
        if with_state_event:
            self.state_events = self.indicator

    def indicator(self, t, y, *args):
        return np.array([y[0] - STATE_EVENT_Y])

    def time_events(self, t, y, *args):
        return TIME_EVENT if t < TIME_EVENT else None

    def step_events(self, solver):
        if self.signalled or not self.signal(solver):
            return False
        self.signalled = True
        return True

    def handle_event(self, solver, event_info):
        state_info, time_event, step_event = event_info
        self.received.append([[int(i) for i in state_info], time_event, step_event])


class ExplicitDecay(Decay, Explicit_Problem):
    def rhs(self, t, y, *sw):
        return -y


class ImplicitDecay(Decay, Implicit_Problem):
    def __init__(self, signal, with_state_event):
        super().__init__(signal, with_state_event, yd0=[-1.0])

    def res(self, t, y, yd, *sw):
        return yd + y


SOLVER_CASES = [
    (ExplicitDecay, solvers.CVode),
    (ExplicitDecay, solvers.Radau5ODE),
    (ExplicitDecay, solvers.RungeKutta34),
    (ExplicitDecay, solvers.Dopri5),
    (ExplicitDecay, solvers.RodasODE),
    (ExplicitDecay, solvers.LSODAR),
    (ExplicitDecay, solvers.ExplicitEuler),
    (ExplicitDecay, solvers.ImplicitEuler),
    (ImplicitDecay, solvers.IDA),
    (ImplicitDecay, solvers.Radau5DAE),
]


def simulate(solver_class, problem):
    solver = solver_class(problem)
    solver.verbosity = QUIET
    solver.report_continuously = True
    solver.simulate(FINAL_TIME)
    return solver.get_statistics()["ntimeevents"]


def step_event_info(received):
    [event_info] = [event_info for event_info in received if event_info[2]]
    return event_info


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_time_event_without_step_event_reports_time_event_only(problem_class, solver_class):
    # Given
    problem = problem_class(never, with_state_event=False)

    # When
    ntimeevents = simulate(solver_class, problem)

    # Then
    assert problem.received == [[[], True, False]]
    assert ntimeevents == 1


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_step_event_alone_reports_step_event_only(problem_class, solver_class):
    # Given
    problem = problem_class(on_first_step, with_state_event=False)

    # When
    simulate(solver_class, problem)

    # Then
    assert step_event_info(problem.received) == [[], False, True]


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_step_event_alone_with_state_events_reports_no_located_root(problem_class, solver_class):
    # Given
    problem = problem_class(on_first_step, with_state_event=True)

    # When
    simulate(solver_class, problem)

    # Then
    assert step_event_info(problem.received) == [[0], False, True]


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_step_event_reaching_time_event_reports_both(problem_class, solver_class):
    # Given
    problem = problem_class(at_time_event, with_state_event=False)

    # When
    ntimeevents = simulate(solver_class, problem)

    # Then
    assert problem.received == [[[], True, True]]
    assert ntimeevents == 1


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_step_event_ending_at_state_event_reports_both(problem_class, solver_class):
    # Given
    problem = problem_class(at_state_event, with_state_event=True)

    # When
    simulate(solver_class, problem)

    # Then
    [state_info], time_event, step_event = step_event_info(problem.received)
    assert state_info != 0
    assert (time_event, step_event) == (False, True)
