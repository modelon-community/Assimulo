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
Created on Tue Feb 21 14:03:12 2012
@author: tony
"""

import numpy as np

from assimulo.problem import Delay_Explicit_Problem
from assimulo.solvers import radar5

class Simple(Delay_Explicit_Problem):
    def __init__(self):
        Delay_Explicit_Problem.__init__(self)
        self.lagcompmap = [[0]]     # Delay 0 affects component 0
            
    def time_lags(self, t, y):
        return [t - 1.0]    # Only one constant time lag, -1
        
    def rhs(self, t, y, ydelay):
        ytm1 = ydelay[0][0] # y(t-1)

        return -y + ytm1
       
#    def jac(self, t, y, ydelay):
#        return J

    def phi(self, i, t):    # History function for $t \in [-1,0]$
        return np.sin(np.pi*t)
    
def run_example():
    """ Basic example of running Radar solver."""
    RTOL = 1.e-6
    ATOL = 1e-6
    H = 1.e-3
    mxst=1000


	# First solve
    p = Simple()

    
    p.y0 = 0    # sin(pi*0) = 0, so consistent with history function

    t0 = 0
    tf = 5.0
    #tf = 1.1
    
    s = radar5.Radar5ODE(p)
    s.grid = np.array([1.0])
    s.grid = np.array([0.5, 2.1])
    s.grid = np.array([0.7, 1.0, 1.1])
    s.grid = np.array([0.7, 1.0, 1.1, 2.5, 3.7, 4.3])
    #s.grid = array([0.7, 1.0, 1.1, 1.5, 3.7, 4.4])
    #s.grid = array([0.2, .7, 2.0])
    #s.grid = array([0.2, .7, 1.0, 2.0, 3.0, 4.0, 5.0])
    s.rtol = RTOL
    s.atol = ATOL
    s.mxst = mxst
    s.inith = H
    s.maxsteps = 1000
    t,y = s.simulate(tf)

    return t, y

if __name__ == '__main__':
    t, y = run_example()
