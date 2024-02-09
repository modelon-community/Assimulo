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

import pytest
from assimulo import testattr
from assimulo.ode import ODE, NORMAL
from assimulo.problem import Explicit_Problem
from assimulo.exception import AssimuloException

class Test_ODE:
    @classmethod
    def setup_class(cls):
        cls.problem = Explicit_Problem(y0=4.0)
        cls.simulator = ODE(cls.problem)
    
    @testattr(stddist = True)
    def test_init(self):
        """
        This tests the functionality of the method __init__.
        """
        assert self.simulator.verbosity == NORMAL
        assert not self.simulator.report_continuously
    
    @testattr(stddist = True)
    def test_verbosity(self):
        """
        This tests the functionality of the property verbosity.
        """

        with pytest.raises(AssimuloException):
            self.simulator._set_verbosity('Test')
        with pytest.raises(AssimuloException):
            self.simulator._set_verbosity([1, 31])
        with pytest.raises(AssimuloException):
            self.simulator._set_verbosity([1])

        self.simulator.verbosity=1
        assert self.simulator.verbosity == 1
        assert self.simulator.options["verbosity"] == 1
        self.simulator.verbosity=4
        assert self.simulator.verbosity == 4
        assert self.simulator.options["verbosity"] == 4
        
    @testattr(stddist = True)    
    def test_report_continuously(self):
        """
        This tests the functionality of the property report_continuously.
        """
        assert not self.simulator.report_continuously #Test the default value
        
        self.simulator.report_continuously = True
        assert self.simulator.report_continuously
        assert self.simulator.options["report_continuously"]
    def test_step_events_report_continuously(self):
        """
        This test tests if report_continuously is set correctly, when step_events are present.
        """
        self.simulator.supports["report_continuously"] = True
        self.simulator.supports["interpolated_output"] = True
        self.simulator.problem_info["step_events"] = True
        self.simulator.problem=self.problem
        self.simulator(10.,ncp=10) # output points and step events should set report_continuously to True 
        assert self.simulator.report_continuously
