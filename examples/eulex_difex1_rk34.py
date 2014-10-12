from  __future__  import division
from scipy import *
from  matplotlib.pyplot import *
import time
from assimulo.solvers import Difex1
from assimulo.solvers import Eulex
from assimulo.solvers import RungeKutta34
from assimulo.solvers import Dopri5
from assimulo.problem import Explicit_Problem
from decimal import Decimal


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
    tim=[]
    nep=[]
    
    for j in eps:
        Explicit_Problem.reset
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        Eulex.reset
        exp_eulex = Eulex(exp_mod) #Create an eulex solver
        exp_eulex.rtol = j
        start=time.clock()
        exp_eulex.simulate(5,100) #Simulate 5 seconds
        stop=time.clock()
        #Number of function evaluation 
        nfev_eulex.append(exp_eulex.statistics["nfcn"])
        t=stop-start
        
        if t!=0:
            tim.append(t)
            nep.append(j)
    
    return exp_mod, exp_eulex,nfev_eulex,tim,nep
    
    
def run_difex(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    nfev_difex1=[]
    tim=[]
    nep=[]
    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_difex1 = Difex1(exp_mod) #Create an difex1 solver
        exp_difex1.rtol = j
        start=time.clock()
        exp_difex1.simulate(5,100) #Simulate 5 seconds
        stop=time.clock()
        nfev_difex1.append(exp_difex1.statistics["nfcn"])
        t=stop-start
        if t!=0:
            tim.append(t)
            nep.append(j)
   
    return exp_mod, exp_difex1,nfev_difex1,tim,nep
 
   
def run_rk34(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    nfev_rk34=[]
    tim=[]
    nep=[]

    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_rk34= RungeKutta34(exp_mod) #Create an difex1 solver
        exp_rk34.rtol = j 
        start=time.clock()
        exp_rk34.simulate(5,100) #Simulate 5 seconds
        stop=time.clock()
        nfev_rk34.append(exp_rk34.statistics["nfcn"])
        t=stop-start
        if t!=0:
            tim.append(t)
            nep.append(j)
   
   
    return exp_mod, exp_rk34,nfev_rk34,tim,nep
    
def run_dopri5(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    nfev_Dopri5=[]
    tim=[]
    nep=[]
    
    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_Dopri5 = Dopri5(exp_mod) #Create an difex1 solver
        exp_Dopri5.rtol = j
        start=time.clock()
        exp_Dopri5.simulate(5,100) #Simulate 5 seconds
        stop=time.clock()
        nfev_Dopri5.append(exp_Dopri5.statistics["nfcn"])
        t=stop-start
        if t!=0:
            tim.append(t)
            nep.append(j)
        
   
    return exp_mod, exp_Dopri5,nfev_Dopri5,tim,nep
    


    
    

     
    
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
    
    eps=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14]
    mod1,sim1,nfev_eulex ,t1,nep1 = run_eulex(f,eps)
    mod2,sim2 ,nfev_difex1,t2,nep2= run_difex(f,eps)
    mod3,sim3 ,nfev_rk34,t3,nep3= run_rk34(f,eps)
    mod4,sim4 ,nfev_Dopri5,t4,nep4= run_dopri5(f,eps)
    
    m=['eulex','difex1','RK34','Dopri5']
    fig=figure(1)
    loglog(eps,nfev_eulex,'r',label='method={}'.format(m[0]))
    loglog(eps,nfev_difex1,'b',label='method={}'.format(m[1]))
    loglog(eps,nfev_rk34,'c',label='method={}'.format(m[2]))
    loglog(eps,nfev_Dopri5,'g',label='method={}'.format(m[3]))
    axis([1e-14,1e-2,10,1e5])
    xlabel('Tolerance')
    ylabel('nfevals')
    legend()
    show()
    
