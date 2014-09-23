import numpy as N

from assimulo.exception import *
from assimulo.ode import *

from assimulo.explicit_ode import Explicit_ODE


try:
    from assimulo.lib import difex1
except ImportError:
    print "Could not find extrapolation pack functions"

# cd -;sudo python setup.py install;cd -;ipython


class Difex1(Explicit_ODE):
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
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Default values
        self.options["atol"]     = 1.e-6       #1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.e-6 #Relative tolerance
        self.options["usejac"]   = False
        self.options["maxsteps"] = 1.    #100000
        self.options["maxh"]     = 0.    #N.inf #Maximum step-size.
        self.options["maxordn"] = 12
        self.options["maxords"] =  5
        self.options["hmax"] = 0.
        self.options["inith"]    = 0.1
        
        
        
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
        self.supports["interpolated_output"] = True
        self.supports["state_events"] = True
        
        self._leny = len(self.y) #Dimension of the problem
        self._type = '(explicit)'
        self._event_info = None
        
    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0
            
    def set_problem_data(self):
        if self.problem_info["state_events"]:
            def event_func(t, y):
                return self.problem.state_events(t, y, self.sw)
            def f(t, y):
                return self.problem.rhs(t, y, self.sw)
            self.f = f
            self.event_func = event_func
            self._event_info = [0] * self.problem_info["dimRoot"]
            self.g_old = self.event_func(self.t, self.y)
        else:
            self.f = self.problem.rhs

    
    def _set_rtol(self, rtol):
        try:
            rtol = float(rtol)
        except (TypeError,ValueError):
            raise Explicit_ODE_Exception('Relative tolerance must be a float.')
        if rtol <= 0.0:
            raise Explicit_ODE_Exception('Relative tolerance must be a positive (scalar) float.')
        self.options["rtol"] = rtol
            
    def _get_rtol(self):
        """
        The relative tolerance to be used in the integration.
        
            Parameters::
            
                rtol    
                            - Default 1.0e-6
                            
                            - Should be a float.
        """
        return self.options["rtol"]
    
    rtol = property(_get_rtol, _set_rtol)
    
                                      
    
    def integrate(self, t, y, tf, opts):
       
        #Check for initialization
        if opts["initialize"]:
            self.set_problem_data()
            self._tlist = []
            self._ylist = []
            
        
####najmeh
        
        tresult=[]
        yresult=[]
        hresult=[]
        flag=[]
        output_index = opts["output_index"]
        output_list  = opts["output_list"][output_index:]   
            
        kflag = 0#ID_PY_COMPLETE

        for tout in output_list:
            output_index += 1
                  
            result=difex1.difex1(self.f,t,y.copy(),tout, self.options["rtol"] , self.options["maxh"] ,self.options["inith"],kflag)
            y=result[1]
            t=result[0]
            H=result[2]
            flag=result[3]
            tresult.append(t)
            hresult.append(H)
            yresult.append(y)
            
            
            self.statistics["nsteps"]        += difex1.statp.nstep
            self.statistics["nfcn"]          += difex1.statp.nfcn
            self.statistics["errfail"]       += difex1.statp.nrejct 
            self.statistics["nlu"]           += difex1.statp.ndec
             
            
            
        return flag, tresult, yresult
        
        
    ##############
        
        #Checking return
        if flag == 1:
            flag = ID_PY_COMPLETE
        elif flag == 2:
            flag = ID_PY_EVENT
        else:
            raise Exception("difex failed with flag %d"%flag)
      
        return flag, self._tlist, self._ylist
        
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
        self.log_message(' Solver                  : difex1 ' + self._type,          verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)
        

      
 
