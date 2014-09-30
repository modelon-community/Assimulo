from  __future__  import division
from scipy import *
from  matplotlib.pyplot import *

from assimulo.solvers import Difex1
from assimulo.solvers import Eulex
from assimulo.solvers import RungeKutta34
from assimulo.solvers import Dopri5
from assimulo.problem import Explicit_Problem


def run_eulex(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """

    #Define an Assimulo problem   
    
    nfev_eulex=[]
    
    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_eulex = Eulex(exp_mod) #Create an eulex solver
        exp_eulex.rtol = j
        exp_eulex.simulate(5,100) #Simulate 5 seconds
        #Number of function evaluation 
        nfev_eulex.append(exp_eulex.statistics["nfcn"])
    
    return exp_mod, exp_eulex,nfev_eulex
    
    
def run_difex(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    nfev_difex1=[]
    
    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_difex1 = Difex1(exp_mod) #Create an difex1 solver
        exp_difex1.rtol = j
        exp_difex1.simulate(5,100) #Simulate 5 seconds
        nfev_difex1.append(exp_difex1.statistics["nfcn"])
   
   
    return exp_mod, exp_difex1,nfev_difex1
 
   
def run_rk34(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    nfev_rk34=[]

    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_rk34= RungeKutta34(exp_mod) #Create an difex1 solver
        exp_rk34.rtol = j 
        exp_rk34.simulate(5,100) #Simulate 5 seconds
        nfev_rk34.append(exp_rk34.statistics["nfcn"])
        
   
   
    return exp_mod, exp_rk34,nfev_rk34
    
def run_dopri5(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    nfev_Dopri5=[]
    
    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_Dopri5 = Dopri5(exp_mod) #Create an difex1 solver
        exp_Dopri5.rtol = j
        exp_Dopri5.simulate(5,100) #Simulate 5 seconds
        nfev_Dopri5.append(exp_Dopri5.statistics["nfcn"])
   
   
    return exp_mod, exp_Dopri5,nfev_Dopri5
    


    
    

     
    
if __name__=='__main__':
    '''
    def f(t,y):
        ydot = -y[0]
        return array([ydot])
    '''
    def f(t,y):
        l=1.0
        g=9.81
        yd_0=y[1]
        yd_1=-g/l*sin(y[0]) 
        
        return array([yd_0,yd_1])
    
    eps=linspace(1e-2,1e-14,100)
    mod1,sim1,nfev_eulex = run_eulex(f,eps)
    mod2,sim2 ,nfev_difex1= run_difex(f,eps)
    mod3,sim3 ,nfev_rk34= run_rk34(f,eps)
    mod4,sim4 ,nfev_Dopri5= run_dopri5(f,eps)
    fig=figure()
    loglog(eps,nfev_eulex,'r')
    loglog(eps,nfev_difex1,'b')
    loglog(eps,nfev_rk34,'c')
    loglog(eps,nfev_Dopri5,'g')
    title('Comparison of total number of function evaluations over whole test set(pendulum)',fontsize=12)
    xlabel('tolerance')
    ylabel('NFEV')
    show()
    #fig.saving('test.jpg')
    
  
