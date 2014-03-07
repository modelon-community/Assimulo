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

from scipy import *
import pylab as P
import nose
from assimulo.problem import Explicit_Problem
from assimulo.solvers import CVode

def run_example(with_plots=True):
    """
    Simulations for the Gyro (Heavy Top) example in Celledoni/Safstrom: 
        Journal of Physics A, Vol 39, 5463-5478, 2006
        
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
    
    """

    def curl(v):
        return array([[0,v[2],-v[1]],[-v[2],0,v[0]],[v[1],-v[0],0]])

    #Defines the rhs
    def f(t,u):
        """
        Simulations for the Gyro (Heavy Top) example in Celledoni/Safstrom: 
        Journal of Physics A, Vol 39, 5463-5478, 2006
        """
        I1=1000.
        I2=5000.
        I3=6000.
        u0=[0,0,1.]
        pi=u[0:3]
        Q=(u[3:12]).reshape((3,3))
        Qu0=dot(Q,u0)
        f=array([Qu0[1],-Qu0[0],0.])
        f=0
        omega=array([pi[0]/I1,pi[1]/I2,pi[2]/I3])
        pid=dot(curl(omega),pi)+f
        Qd=dot(curl(omega),Q)
        return hstack([pid,Qd.reshape((9,))])

    def energi(state):
        energi=[]
        for st in state:
            Q=(st[3:12]).reshape((3,3))
            pi=st[0:3]
            u0=[0,0,1.]
            Qu0=dot(Q,u0)
            V=Qu0[2]  # potential energy
            T=0.5*(pi[0]**2/1000.+pi[1]**2/5000.+pi[2]**2/6000.)
            energi.append([T])
        return energi

    #Initial conditions
    y0=hstack([[1000.*10,5000.*10,6000*10],eye(3).reshape((9,))])
    
    #Create an Assimulo explicit problem
    exp_mod = Explicit_Problem(f,y0, name="Gyroscope Example")
    
    #Create an Assimulo explicit solver (CVode)
    exp_sim=CVode(exp_mod)
    
    #Sets the parameters
    exp_sim.discr='BDF'
    exp_sim.iter='Newton'
    exp_sim.maxord=2 #Sets the maxorder
    exp_sim.atol=1.e-10
    exp_sim.rtol=1.e-10
    
    #Simulate
    t, y = exp_sim.simulate(0.1)
    
    #Basic tests
    nose.tools.assert_almost_equal(y[-1][0],692.800241862)
    nose.tools.assert_almost_equal(y[-1][8],7.08468221e-1)
    
    #Plot
    if with_plots:
        P.plot(t,y/10000.)
        P.xlabel('Time')
        P.ylabel('States, scaled by $10^4$')
        P.title(exp_mod.name)
        P.show()
        
    return exp_mod, exp_sim    


if __name__=='__main__':
    mod,sim = run_example()
