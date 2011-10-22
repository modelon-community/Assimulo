#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2011 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as N

from assimulo.ode import *
from assimulo.explicit_ode import Explicit_ODE

class ExplicitEuler(Explicit_ODE):
    """
    Explicit Euler.
    """
    def __init__(self, problem):
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Solver options
        self.solver_options["h"] = 0.01
        
        #Internal temporary result vector
        self.yd = N.array([0.0]*len(problem.y0))
        
        #RHS-Function
        self.f = problem.rhs_internal
        
        #Solver support
        self.solver_support["one_step_mode"] = True
        
    
    def one_step_mode(self, t, y, tf, initialize, *args):
        if initialize:
            self.solver_iterator = self.integrator(t,y,tf,*args)

        return self.solver_iterator.next()
    
    def integrator(self, t, y, tf, *args):
        """
        Integrates (t,y) values until t > tf
        """
        h = self.solver_options["h"]
        h = min(h, abs(tf-t))
 
        while t+h < tf:
            t, y = self.step(t, y, h)
            yield ID_OK, t,y
            h=min(h, abs(tf-t))
        else:
            t, y = self.step(t, y, h)
            yield ID_COMPLETE, t, y
    
    def step(self, t, y, h):
        """
        This calculates the next step in the integration.
        """
        self.f(self.yd,t,y) #The output is stored in yd
        return t + h, y + h*self.yd
        
    def print_statistics(self, verbose=NORMAL):
        """
        Should print the statistics.
        """
        self.logg_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        self.logg_message(' Step-length          : %s '%(self.solver_options["h"]), verbose)
        self.logg_message('\nSolver options:\n',                                    verbose)
        self.logg_message(' Solver            : ExplicitEuler',                     verbose)
        self.logg_message(' Solver type       : Fixed step\n',                      verbose)
