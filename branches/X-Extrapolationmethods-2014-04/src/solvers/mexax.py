__author__ = 'najmeh'

#ticket:346


import numpy as np

from assimulo.ode import *




from scipy import *
from assimulo.support import set_type_shape_array
from assimulo.exception import *
from assimulo.ode import *
from assimulo.implicit_ode import MexaxDAE


from assimulo.special_systems import Mechanical_System

from assimulo.exception import *



try:
    from assimulo.lib import mexax
except ImportError:
    print "Could not find extrapolation pack functions"

# cd -;sudo python setup.py install;cd -;ipython


class Mexax(MexaxDAE):
    '''
    add some documentation


    '''
    def __init__(self, problem):
        """
        Initiates the solver.

            Parameters::

                problem
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.
        """
        MexaxDAE.__init__(self, problem) #Calls the base class

         #Default values
        self.options["inith"]    = 1.e-4
        self.options["maxh"]     = 0.0
        self.options["safe"]     = 0.9 #Safety factor
        self.options["atol"]     = 1.0e-6*np.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 5000
        self.options["maxord"]   = 0
        


        # - Statistic values
        self.statistics["nsteps"]      = 0 #Number of steps
        self.statistics["nfcn"]        = 0 #Number of function evaluations
        self.statistics["njac"]        = 0 #Number of Jacobian evaluations
        self.statistics["njacfcn"]     = 0 #Number of function evaluations when evaluating the jacobian
        self.statistics["errfail"]     = 0 #Number of step rejections
        self.statistics["nlu"]         = 0 #Number of LU decompositions
        self.statistics["nstepstotal"] = 0 #Number of total computed steps (may NOT be equal to nsteps+nerrfail)
        self.statistics["nstateevents"]= 0 #Number of state events
        self.statistics["ngevals"]     = 0 #Root evaluations

        #Solver support
        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = False
        self.supports["state_events"] = False
        
        self.problem_info["dim"]
   

    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
            
   
    def _solout(self,t, p, v, u, a, rlam, infos, irtrn):   # <--- change pyf so that this becomes an out variable only
        self._tlist.append(t)
        y=vstack((p,v))
        yd=vstack((v,a))
        self._ylist.append(y.copy())   #  <---- return y and check if you want also a self.yd
        self._ydlist.append(yd.copy())
        return irtrn   #?? now istrn is output in pyf do we need it to be return?
    
    
              
    def _denout(self,t, p,v,u,a, lam,indp,indv,indu,inda,indl,info):
        
        self._tlist.append(t)
        self._ylist.append(self.p)  #  <---- return y and check if you want also a self.yd
        y=vstack((p,v))
        yd=vstack((v,a))
        self._ylist.append(y.copy())   #  <---- return y and check if you want also a self.yd
        self._ydlist.append(yd.copy())
    
        
    def fswit(self,t,p,v,u,a,rlam,g=9.81):
        raise Exception('fswit only in dummy mode')
        
    
   
    def integrate(self, t, y, yd, tfin):
                                                                            
        
        info = np.zeros((15,),np.int) 
        info[1] = 1  # Tolerances are vectors
        # set mxjob control parameters  
        mxjob=zeros((150,),np.int)
        normal_mode = 1 if opts["output_list"] == None or opts["report_continuously"] else 0  # intermediate output mode
        mxjob[30-1] = 1 if normal_mode == 0 else 0
        mxjob[31-1] = 1 if normal_mode else 0
        mxjob[32-1] = len(opts["output_list"])
                            
                               
        atol = self.options["atol"]
        rtol = self.options["rtol"]
        itol = 1
        
        
        self._tlist  = []
        self._ylist  = []
        self._ydlist = []
    
        #Store the opts
        self._opts = opts
    
        # Provide the workspace
        # a) Compute the size of the workspace
        np = nv = self.problem.n_p
        nu = 0
        nl = self.problem.n_la
        ny  = self.problem_info["dim"]  # should be np+nv+nu+nl
        if np+nv+nu+nl != ny:
           raise Exception('Dimension error: np+nv+nu+nl != ny')
        
        liwk=np + 4*nv + nl + nu + 60
        iwk=empty((liwk,),dtype=int)
        ngl=max(ng,nl)
        lo=(np+nv+nu+nv+nl)*156 if mxjob[31-1] else 0
        lrwk=(nv+nl)**2+np*(ngl+18)+nv*(nv+45)+28*nl+max(ng,1)+18*max(nu,1)+50+lo
        rwk=empty((lrwk,),dtype=float)
        p=y[:np].copy()
        v=y[np,np+nv].copy()
        a=yd[nv,2*nv].copy()   #  <---------------- take this from  ydot
        u=zeros((1,))
        lam=y[np+nv,:].copy()
        h=self.options["inith"]
        
        [t, p, v, u, a, lam, h, mxjob, ierr, iwk, rwk]= \
                                                 mexax.mexx(nl,ng,nu,self.fprob,t,tfin,                                            
                                                     p,v,u,a,lam,
                                                     itol,rtol,atol,
                                                     h,mxjob,iwk,rwk,
                                                     self.solout,self.denout,self.fswit,lrwk=lrwk,liwk=liwk)
        
        
        
        #Retrieving statistics

        self.statistics["nsteps"]        += mxjob[51-1]
        self.statistics["nfcn"]          += mxjob[55-1]
        self.statistics["errfail"]       += mxjob[73-1]
        self.statistics["nlu"]           += mxjob[58-1]


        return mxjob

    def state_event_info(self):
        return self._event_info

    def set_event_info(self, event_info):
        self._event_info = event_info

    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)

        self.log_message(' Number of steps                          : '+ str(self.statistics["nsteps"]),          verbose)
        self.log_message(' Number of function evaluations           : '+ str(self.statistics["nfcn"]),         verbose)
        self.log_message(' Number of Jacobian evaluations           : '+ str(self.statistics["njac"]),    verbose)
        self.log_message(' Number of error test failures            : '+ str(self.statistics["errfail"]),       verbose)
        self.log_message(' Number of LU decompositions              : '+ str(self.statistics["nlu"]),       verbose)
        if self.problem_info["state_events"]:
            self.log_message(' Number of event function evaluations     : '+ str(self.statistics["ngevals"]),        verbose)
            self.log_message(' Number of State-Events                   : '+ str(self.statistics["nstateevents"]),   verbose)

        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : mexax ' + self._type,          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self._compact_atol()),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)




