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
import pylab as P
import nose
from assimulo.solvers import DASP3ODE
from assimulo.problem import SingPerturbed_Problem

def run_example(with_plots=True):
    
    #Define the slow rhs
    def dydt(t,y,z):
        eps=(1./3)*1.e-3
        yp = N.array([-(.6*z[0]+.8*y[2])*y[0]+10.*y[1],
                      -10.*y[1]+ 1.6*z[0] *y[2],
                      -1.33*eps**2*y[2]*(y[0]+2.*z[0])])										 
        return yp
    
    #Define the fast rhs
    def dzdt(t,y,z):
        eps=(1./3)*1.e-3
        zp = N.array([1.6*z[0]*y[2]-.6*z[0]*y[0]
                      -45.*(eps*z[0])**2+.8*y[2]*y[0]])
        return zp
    
    #The initial values
    y0 = [ 3.0, 0.216, 1.0]
    z0 = [1.35]
    eps = N.array([.33333333e-3])
    
    #Define an Assimulo problem
    exp_mod = SingPerturbed_Problem(dydt, dzdt, yy0=y0, zz0=z0,eps=eps)
    exp_mod.name = 'Simple DASP3 Example'
    
    #Define an explicit solver
    exp_sim = DASP3ODE(exp_mod) #Create a CVode solver
    
    #Sets the parameters
    exp_sim.rtol = 1e-3 #Default 1e-6
    exp_sim.atol = 1e-3 #Default 1e-6

    #Simulate
    t, y = exp_sim.simulate(10) #Simulate 10 seconds

    #Basic test
    nose.tools.assert_almost_equal(y[-1,0], 10.87313065, 3)
    
    #Plot
    if with_plots:
        P.plot(t, y, color="b")
        P.show()

if __name__=='__main__':
    run_example()
