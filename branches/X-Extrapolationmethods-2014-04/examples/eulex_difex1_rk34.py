from  __future__  import division
from scipy import *
from  matplotlib.pyplot import *

from assimulo.solvers import Difex1
from assimulo.solvers import Eulex
from assimulo.solvers import RungeKutta34
from assimulo.problem import Explicit_Problem

def run_example(with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """

    
        
    #Defines the rhs
    def f(t,y):
        ydot = -y[0]
        return array([ydot])

    ###Eulex
    #Define an Assimulo problem
    
    eps=[1e-2,1e-4,1e-6,1e-8,1e-10,1e-12,1e-14]
    nfev_eulex=[]
    nfev_difex1=[]
    
    for j in eps:
        Explicit_Problem.reset
        exp_mod = Explicit_Problem(f, y0=array([4.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_eulex = Eulex(exp_mod) #Create an eulex solver
        exp_eulex.rtol = j
        exp_eulex.simulate(5,100) #Simulate 5 seconds
        #print exp_eulex.statistics["nfcn"]
        nfev_eulex.append(exp_eulex.statistics["nfcn"])
       
        
    '''    
    for j in eps:
        Explicit_Problem.reset
        exp_mod = Explicit_Problem(f, y0=array([4.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_difex1 = Difex1(exp_mod) #Create an difex1 solver
        exp_difex1.rtol = j
        exp_difex1.simulate(5,100) #Simulate 5 seconds
        nfev_difex1.append(exp_difex1.statistics["nfcn"])
        
    if with_plots:
        #print nfev_eulex,nfev_difex1,eps
        #loglog(eps,nfev_eulex)
        loglog(eps,nfev_difex1)
        
    '''
 
   
    
    
    
    
    
    
    
    
    
    
    
if __name__=='__main__':
    mod,sim = run_example()

