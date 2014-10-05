from  __future__  import division
from scipy import *
from  matplotlib.pyplot import *
import time
from assimulo.solvers import CVode
from assimulo.solvers import Metan1
from assimulo.solvers import Eulsim

from assimulo.problem import Explicit_Problem


def run_cvode(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """

    #Define an Assimulo problem   
    
    nfev_cvode=[]
    tim=[]
    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'CVode Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_cvode = CVode(exp_mod) #Create an eulex solver
        #exp_cvode.iter  = 'Newton' #Default 'FixedPoint'
        #exp_cvode.discr = 'BDF' #Default 'Adams'
        exp_cvode.rtol = j
        start=time.clock()
        for i in range(50):
            exp_cvode.reset()
            exp_cvode.simulate(5,100) #Simulate 5 seconds
        stop=time.clock()
        #Number of function evaluation 
        nfev_cvode.append(exp_cvode.statistics["nfevals"])
        t=stop-start
        tim.append(t)
    return exp_mod, exp_cvode,nfev_cvode,tim
    
    
def run_metan1(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    nfev_metan1=[]
    tim=[]
    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'metan1 Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_metan1 = Metan1(exp_mod) #Create an metan1 solver
        exp_metan1.rtol = j
        start=time.clock()
        exp_metan1.simulate(5,100) #Simulate 5 seconds
        stop=time.clock()
        nfev_metan1.append(exp_metan1.statistics["nfcn"])
        tim.append(stop-start)
   
    return exp_mod, exp_metan1,nfev_metan1,tim
 
   
def run_eulsim(f,eps,with_plots=True):
    r"""
    Demonstration of the use of eulex by solving the
    linear test equation :math:`\dot y = - y`
    
    on return:
    
       - :dfn:`exp_mod`    problem instance
    
       - :dfn:`exp_sim`    solver instance
       
    """
    nfev_eulsim=[]
    tim=[]
    for j in eps:
        exp_mod = Explicit_Problem(f, y0=array([1.0,1.0]), name = 'Eulex Test Example: $\dot y = - y$')
        #Define an explicit solver
        exp_eulsim= Eulsim(exp_mod) #Create an difex1 solver
        exp_eulsim.rtol = j 
        start=time.clock()
        exp_eulsim.simulate(5,100) #Simulate 5 seconds
        stop=time.clock()
        nfev_eulsim.append(exp_eulsim.statistics["nfcn"])
        tim.append(stop-start)
   
   
    return exp_mod, exp_eulsim,nfev_eulsim,tim

    
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
    eps1=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15]
    eps=linspace(1e-2,1e-14,100)
    mod1,sim1,nfev_cvode ,t1= run_cvode(f,eps1)
    mod2,sim2 ,nfev_metan1,t2= run_metan1(f,eps1)
    mod3,sim3 ,nfev_eulsim,t3= run_eulsim(f,eps1)
    
    fig=figure(1)
    loglog(eps1,nfev_cvode,'r')
    loglog(eps1,nfev_metan1,'b')
    loglog(eps1,nfev_eulsim,'c')
    fig=figure(2)
    plot(eps1,t1,'r')
    plot(eps1,t2,'b')
    plot(eps1,t3,'c')
    #axis([0.001,0.101,0.001,0.15])
    ##title('Comparison of total number of function evaluations over whole test set(pendulum)',fontsize=12)
    #xlabel('tolerance')
    ##ylabel('NFEV')
    show()
    #fig.saving('test.jpg')
    
  
