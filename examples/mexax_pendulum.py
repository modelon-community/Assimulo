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
    mexax_pend.simulate(10.,100)   
    return my_pend, mexax_pend
        
if __name__=='__main__':
    index='oproj2'
    results=run_example(index)
    
