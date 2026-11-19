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

from dataclasses import dataclass, field

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
        self.handled_events_info = []
        if with_state_event:
            self.state_events = self.indicator

    def indicator(self, t, y, *args):
        return np.array([y[0] - STATE_EVENT_Y])

    def time_events(self, t, y, *args):
        return TIME_EVENT if t < TIME_EVENT else None

    def step_events(self, solver):
        if not self.signalled and self.signal(solver):
            self.signalled = True
            return True

        return False

    def handle_event(self, solver, event_info):
        self.handled_events_info.append(list(event_info))


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


def simulate_and_get_nevents(solver_class, problem):
    solver = solver_class(problem)
    solver.verbosity = QUIET
    solver.report_continuously = True
    solver.simulate(FINAL_TIME)
    return solver.get_statistics()["ntimeevents"]


@dataclass
class Event:
    """Expected event_info passed to handle_event. Omitting `state` means no state event info."""

    time: bool = False
    step: bool = False
    state: list = field(default_factory=list)


def assert_events(problem, expects):
    received = [
        Event(time=time_event, step=step_event, state=[int(s) for s in state_info])
        for state_info, time_event, step_event in problem.handled_events_info
    ]
    assert received == expects


def downward_crossing(solver_class):
    # LSODAR reports the opposite sign for a root crossing compared to the other solvers.
    return 1 if solver_class is solvers.LSODAR else -1


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_time_event_without_step_event_reports_time_event_only(problem_class, solver_class):
    # Given
    problem = problem_class(never, with_state_event=False)

    # When
    ntimeevents = simulate_and_get_nevents(solver_class, problem)

    # Then
    assert_events(problem, expects=[Event(time=True)])
    assert ntimeevents == 1


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_step_event_alone_reports_step_event_only(problem_class, solver_class):
    # Given
    problem = problem_class(on_first_step, with_state_event=False)

    # When
    simulate_and_get_nevents(solver_class, problem)

    # Then
    assert_events(problem, expects=[
        Event(step=True),
        Event(time=True),
    ])


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_step_event_alone_with_state_events_reports_no_located_root(problem_class, solver_class):
    # Given
    problem = problem_class(on_first_step, with_state_event=True)

    # When
    simulate_and_get_nevents(solver_class, problem)

    # Then
    assert_events(problem, expects=[
        Event(step=True, state=[0]),
        Event(state=[downward_crossing(solver_class)]),
        Event(time=True),
    ])


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_step_event_reaching_time_event_reports_both(problem_class, solver_class):
    # Given
    problem = problem_class(at_time_event, with_state_event=False)

    # When
    ntimeevents = simulate_and_get_nevents(solver_class, problem)

    # Then
    assert_events(problem, expects=[Event(time=True, step=True)])
    assert ntimeevents == 1


@pytest.mark.parametrize("problem_class, solver_class", SOLVER_CASES)
def test_step_event_ending_at_state_event_reports_both(problem_class, solver_class):
    # Given
    problem = problem_class(at_state_event, with_state_event=True)

    # When
    simulate_and_get_nevents(solver_class, problem)

    # Then
    assert_events(problem, expects=[
        Event(step=True, state=[downward_crossing(solver_class)]),
        Event(time=True),
    ])
