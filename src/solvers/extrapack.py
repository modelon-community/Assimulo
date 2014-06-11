from  __future__  import division
from  scipy       import *
import numpy as N
from assimulo.ode import *
from assimulo.explicit_ode import Explicit_ODE
from assimulo.exception import *


try:
    from assimulo.lib import eulex
except ImportError:
    print "Could not find extrapolation pack functions"



class Eulex(Explicit_ODE):
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
        self.options["inith"]    = 0.01
        self.options["newt"]     = 7 #Maximum number of newton iterations
        self.options["thet"]     = 1.e-3 #Boundary for re-calculation of jac
        self.options["fnewt"]    = 0.0 #Stopping critera for Newtons Method
        self.options["quot1"]    = 1.0 #Parameters for changing step-size (lower bound)
        self.options["quot2"]    = 1.2 #Parameters for changing step-size (upper bound)
        self.options["fac1"]     = 0.2 #Parameters for step-size selection (lower bound)
        self.options["fac2"]     = 8.0 #Parameters for step-size selection (upper bound)
        self.options["maxh"]     = N.inf #Maximum step-size.
        self.options["safe"]     = 0.9 #Safety factor
        self.options["atol"]     = 1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["usejac"]   = True if self.problem_info["jac_fcn"] else False
        self.options["maxsteps"] = 10000
        
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
        
        
    '''    
        # Do we need this one?
    def _solout(self, nrsol, told, t, y, cont, lrc, irtrn):
        """
        This method is called after every successful step taken by Radau5
        """
        self.cont = cont #Saved to be used by the interpolation function.
        
        if self.problem_info["state_events"]:
            flag, t, y = self.event_locator(told, t, y)
            #Convert to Fortram indicator.
            if flag == ID_PY_EVENT: irtrn = -1
            
        if self._opts["report_continuously"]:
            initialize_flag = self.report_solution(t, y, self._opts)
            if initialize_flag: irtrn = -1
        else:
            if self._opts["output_list"] == None:
                self._tlist.append(t)
                self._ylist.append(y.copy())
            else:
                output_list = self._opts["output_list"]
                output_index = self._opts["output_index"]
                try:
                    while output_list[output_index] <= t:
                        self._tlist.append(output_list[output_index])
                        self._ylist.append(self.interpolate(output_list[output_index]))
                        
                        output_index += 1
                except IndexError:
                    pass
                self._opts["output_index"] = output_index
        
        return irtrn
    '''
   
   
   
                                     
    
    def integrate(self, t, y, tf, opts):
        ITOL  = 1 #Both atol and rtol are vectors
        #IJAC  = 1 if self.usejac else 0 #Switch for the jacobian, 0==NO JACOBIAN
        MLJAC = self.problem_info["dim"] #The jacobian is full
        MUJAC = self.problem_info["dim"] #See MLJAC
        IMAS  = 0 #The mass matrix is the identity
        MLMAS = self.problem_info["dim"] #The mass matrix is full
        MUMAS = self.problem_info["dim"] #See MLMAS
        IOUT  = 1 #solout is called after every step
        WORK  = N.array([0.0]*(4*self.problem_info["dim"]**2+12*self.problem_info["dim"]+20)) #Work (double) vector
        IWORK = N.array([0]*(3*self.problem_info["dim"]+20)) #Work (integer) vector
        '''
        #Setting work options
        WORK[1] = self.safe
        WORK[2] = self.thet
        WORK[3] = self.fnewt
        WORK[4] = self.quot1
        WORK[5] = self.quot2
        WORK[6] = self.maxh
        WORK[7] = self.fac1
        WORK[8] = self.fac2
        
        #Setting iwork options
        IWORK[1] = self.maxsteps
        IWORK[2] = self.newt
        
        #Dummy methods
        mas_dummy = lambda t:x
        jac_dummy = (lambda t:x) if not self.usejac else self.problem.jac
        '''
        #Check for initialization
        if opts["initialize"]:
            self.set_problem_data()
            self._tlist = []
            self._ylist = []
            
        
####najmeh
            
        output_list=linspace(0,10,100) 
        kflag=0
        h=1.e-15 #H
        hmax=1.
        y=array(1*[10.])
        told=0.  # T / starting point 
        tresult=[]
        yresult=[]
        hresult=[]
        
        for t in output_list:
            #output_index += 1
             
             
            '''            
            t,y,h, iwork, flag  = eulex.eulex(self.f, t, y.copy(),  tout, self.atol, self.options["maxh"] 
                                , self.h, kflag)
            #Store results
            tlist.append(t)
            ylist.append(y.copy())
            self._event_info = roots
            '''
            result=eulex.eulex(self.f,0.,y.copy(),t,1.e-4,1.,1.e-15,kflag)
            print eulex.statp.nfcn
            y=result[1]
            told=result[0]
            H=result[2]
            flag=result[3]
            tresult.append(t)
            hresult.append(H)
            yresult.append(y[0])
        
            #opts["output_index"] = output_index
       
        
        #Retrieving statistics
        #self.statistics["ngevals"]            += IWORK[9]
        #self.statistics["nsteps"]        += IWORK[10]
        self.statistics["nfcn"]          += eulex.statp.nfcn
        #self.statistics["njac"]          += IWORK[12]   
        #self.statistics["nevents"] += 1  if flag == ID_PY_EVENT else 0
        # save RWORK, IWORK for restarting feature
        
    
        return flag, tresult, yresult
        
        
    ##############
        
        
        
        
        
        #Checking return
        if flag == 1:
            flag = ID_PY_COMPLETE
        elif flag == 2:
            flag = ID_PY_EVENT
        else:
            raise Exception("eulex failed with flag %d"%flag)
        
        #Retrieving statistics
        self.statistics["nsteps"]      += iwork[16]
        self.statistics["nfcn"]        += iwork[13]
        self.statistics["njac"]        += iwork[14]
        self.statistics["nstepstotal"] += iwork[15]
        self.statistics["errfail"]     += iwork[17]
        self.statistics["nlu"]         += iwork[18]
        
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
        self.log_message(' Solver                  : eulex ' + self._type,          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self._compact_atol()),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)
        

      
            
        
        
    
