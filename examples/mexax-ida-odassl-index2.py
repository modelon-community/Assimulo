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

from  __future__  import division
#import nose
import assimulo.problem as ap
import assimulo.special_systems as ass
import numpy as N
import scipy.linalg as sl
from assimulo.solvers import IDA, ODASSL, Mexax
from  matplotlib.pyplot import *


def pendulum():
    g=13.7503671
    n_p=2
    n_la=1
    def forces(t, p, v):
        return N.array([0.,-g])
    def GT(p):
        return N.array([p[0],p[1]]).reshape((2,1))
    def constr3(t,y):
        p=y[0:2]
        return N.array([p[0]**2+p[1]**2-1.])
    def constr2(t,y):
        p,v=y[0:2],y[2:4]
        return N.array([p[0]*v[0]+p[1]*v[1]])
    def constr1(t,y):
        p,v,la=y[0:2],y[2:4],y[4:5]
        return N.array([v[0]**2+v[1]**2 - la[0] * (p[0]**2 + p[1]**2) - p[1] * g])
    def mass_matrix(t,p):
        return N.eye(2)
    return ass.Mechanical_System(n_p, forces, n_la,  
                                   [1.,0.], [0.,0.],
                                   [0],
                                   [0.,0.], [0.,-g], GT=GT,
                                   constr3 = constr3,
                                   constr2 = constr2, 
                                   constr1 = constr1,
                                   mass_matrix=mass_matrix)

def run_example(eps,index, with_plots=False, with_test=False):
    my_pend_sys=pendulum()
    my_pend=my_pend_sys.generate_problem(index)
    my_pend.name='Index {} System'.format(index)
    if index in ('ind1','ind2','ind3','ggl2'):
        dae_pend = IDA(my_pend)  
    elif index in ('ovstab2','ovstab1'):
        dae_pend = ODASSL(my_pend) 
    elif index in ('oproj2') :
        dae_pend = Mexax(my_pend) 
    else:
        raise Exception('No method for index {}'.format(index))
    dae_pend.atol=eps
    dae_pend.rtol=1.e-10
    dae_pend.suppress_alg=True  
    t,y,yd=dae_pend.simulate(5.,100)   
    if index != 'oproj2':
        final_residual=my_pend.res(0.,dae_pend.y,dae_pend.yd)
        whichres='All'
    else:
        qflag=9*[0]
        qflag[6]=1
        parameterlist={'nl':my_pend.n_la,'ng':0,'nu':0,'t':dae_pend.t,
                                  'p':dae_pend.y[:2],'v':dae_pend.y[2:4],'u':N.array([0.]),
                                  'lam':dae_pend.y[4],'mass':N.empty((2,2)),'gp':N.empty((1,2)),
                                  'f':N.empty((2,)),'pdot':dae_pend.yd[:2],'udot':N.empty((1,)),
                                  'g':N.empty((1,)),'gi':N.empty((1,)),'fl':N.empty((1,)),
                                  'qflag':qflag}
        final_residual=my_pend.fprob(**parameterlist)[5]      # only position residual
        whichres='Position'
    print(my_pend.name+" {} residuals after the integration run\n".format(whichres))
    print final_residual, 'Norm:  ', sl.norm(final_residual) 
    if with_test:
        assert(sl.norm(final_residual) < 1.5e-2)
    if with_plots:
        dae_pend.plot(mask=[1,1]+(len(my_pend.y0)-2)*[0]) 
    nfv=dae_pend.statistics["nfevals"]
    return my_pend, dae_pend,nfv
        
if __name__=='__main__':
    eps=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13]
    index_values=['ind2','ggl2','ovstab2','oproj2']
    ind_color=['r','b','c','y']
    index_all=zip(index_values,ind_color)
    sim=[]
    mod=[]
    
    for ind,colour in index_all:
        nfv=[]
        for i in eps:
            m, s,nf=run_example(i,index=ind)
            mod.append(m)
            sim.append(s)
            nfv.append(nf)
        loglog(eps,nfv,'{}o'.format(colour))
        loglog(eps,nfv,'{}'.format(colour),label='Method={}'.format(ind))
        xlabel('Tolerance')
        ylabel('nfevals')
        legend()
    show()
        
        
                                  
                                       
    
        
