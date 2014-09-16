# -*- coding: utf-8 -*-
"""
@author: Najmeh
"""
from  __future__  import division
from scipy import *
from  matplotlib.pyplot import *

from assimulo.special_systems import Mechanical_System
from assimulo.exception import *
from assimulo.solvers import Mexax


def pendulum():
        g=13.7503671
        n_p=2
        n_la=1
        def forces(t, p, v):
            return array([0.,-g])
        def GT(p):
            return array([p[0],p[1]]).reshape((2,1))
        def constr3(t,y):
            p=y[0:2]
            return array([p[0]**2+p[1]**2-1.])
        def constr2(t,y):
            p,v=y[0:2],y[2:4]
            return array([p[0]*v[0]+p[1]*v[1]])
        def constr1(t,y):
            p,v,la=y[0:2],y[2:4],y[4:5]
            return array([v[0]**2+v[1]**2 - la[0] * (p[0]**2 + p[1]**2) - p[1] * g])
        def mass_matrix(t,p):
            return eye(2)
        return Mechanical_System(n_p, forces, n_la,
                                     [1.,0.], [0.,0.],
                                     [0],
                                     [0.,0.], [0.,-g], GT=GT,
                                     constr3 = constr3,
                                     constr2 = constr2,
                                     constr1 = constr1,
                                     mass_matrix=mass_matrix)
                                     
                                     
                                     
                                     
                                     

def run_example(index, with_plots=True, with_test=False):
    my_pend_sys=pendulum()
    my_pend=my_pend_sys.generate_problem(index)
    my_pend.name='Index = {}'.format(index)
    mexax_pend=Mexax(my_pend) 
    my_pend.atol=1.e-12
    my_pend.rtol=1.e-12  
    mexax_pend.simulate(10.,100)  
    print'y is {} and y0 {}'.format(mexax_pend.y,mexax_pend.y0)
    print 'final_residual={}'.format(my_pend.res(0.,mexax_pend.y,mexax_pend.yd))
    mexax_pend.plot(mask=[1,1]+(len(my_pend.y0)-2)*[0])  
    return my_pend, mexax_pend
        
if __name__=='__main__':
    index='oproj2'
    results=run_example(index)
    '''
    my_pend_sys=pendulum()
    my_pend=my_pend_sys.generate_problem(index)
    my_pend.name='Index = {}'.format(index)
    mexax_pend=Mexax(my_pend) 
    qflag=9*[False]
    parameterlist={'nl':1,'ng':0,'nu':0,'t':0.,
                                  'p':array([0.,1.]),'v':array([0.,0.]),'u':array([0.]),
                                  'lam':array([0.]),'mass':empty((2,2)),'gp':empty((1,2)),
                                  'f':empty((2,)),'pdot':empty((2,)),'udot':empty((1,)),
                                  'g':empty((1,)),'gi':empty((1,)),'fl':empty((1,)),
                                  'qflag':qflag}
                                  
    my_pend.fprob(**parameterlist)[0]                              
    '''





