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

import numpy as np
import pytest
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem

def run_example(with_plots=True):
    r"""
    An example for IDA with scaled preconditioned GMRES method
    as a special linear solver.
    Note, how the operation Jacobian times vector is provided.
    
    ODE:
    
    .. math::
       
       \dot y_1 - y_2 &= 0\\
       \dot y_2 -9.82 &= 0
       
    
    on return:
    
       - :dfn:`imp_mod`    problem instance
    
       - :dfn:`imp_sim`    solver instance
       
    """
    
    #Defines the residual
    def res(t,y,yd):
        res_0 = yd[0] - y[1]
        res_1 = yd[1] + 9.82

        return np.array([res_0,res_1])
    
    #Defines the Jacobian*vector product
    def jacv(t,y,yd,res,v,c):
        jy = np.array([[0,-1.],[0,0]])
        jyd = np.array([[1,0.],[0,1]])
        j = jy+c*jyd
        return np.dot(j,v)
    
    #Initial conditions
    y0 = [1.0,0.0] 
    yd0 = [0.0, -9.82]
    
    #Defines an Assimulo implicit problem
    imp_mod = Implicit_Problem(res,y0,yd0,name = 'Example using the Jacobian Vector product')
    
    imp_mod.jacv = jacv #Sets the jacobian
    
    imp_sim = IDA(imp_mod) #Create an IDA solver instance
    
    #Set the parameters
    imp_sim.atol = 1e-5 #Default 1e-6
    imp_sim.rtol = 1e-5 #Default 1e-6
    imp_sim.linear_solver = 'SPGMR' #Change linear solver
    #imp_sim.options["usejac"] = False
    
    #Simulate
    t, y, yd = imp_sim.simulate(5, 1000) #Simulate 5 seconds with 1000 communication points
    
    #Basic tests
    assert y[-1][0] == pytest.approx(-121.75000000, abs = 1e-4)
    assert y[-1][1] == pytest.approx(-49.100000000)
    
    #Plot
    if with_plots:
        import pylab as pl
        pl.plot(t,y)
        pl.xlabel('Time')
        pl.ylabel('State')
        pl.title(imp_mod.name)
        pl.show()
    return imp_mod,imp_sim
    
if __name__=='__main__':
    mod,sim = run_example()
