__author__ = 'najmeh'

#ticket:346

import numpy as np
from assimulo.exception import *
from assimulo.ode import *
from assimulo.implicit_ode import MexaxDAE


try:
    from assimulo.lib import mexax
except ImportError:
    print "Could not find extrapolation pack functions"

# cd -;sudo python setup.py install;cd -;ipython


class Mexax(ImplicitProblem):
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
        Implicit_ODE.__init__(self, problem) #Calls the base class

         #Default values
        self.options["inith"]    = 0.0
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
        
    def check_instance(self):
        if not isinstance(self.problem, MEXAX_Problem):
            raise Implicit_ODE_Exception('The problem needs to be a subclass of MEXAX_Problem.')


    def initialize(self):
        #Reset statistics
        for k in self.statistics.keys():
            self.statistics[k] = 0

    

    def integrate(self, t, y, yprime, tf, opts):
                                                                            
        
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
        tlist=[]
        ylist=[]
        ydlist=[]
        
        
        #Store the opts
        self._opts = opts
    
        # Provide the workspace
        # a) Compute the size of the workspace
        np = self.problem.n_p
        nv = np
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
        
        
        
        
        
        if opts["report_continuously"]:
            idid = 1
            while idid==1:




                mexax.mexx(self.problem.np,self.problem.nv,
                           self.n_la,
                           ng,nu,self.problem.fprob,
                           t,tfin,
                           y[:self.np],y[self.np:self.nv+self.np],  # p,v
                           u,a,   # ????
                           y[self.nv+self.np:],    # lambda
                           itol,rtol,atol,h,mxjob,ierr,liwk,iwk,lrwk,rwk,solout,denout,fswit)

                #?
                initialize_flag = self.report_solution(t, y, yprime, opts)
                if initialize_flag:
                    flag = ID_PY_EVENT
                    break
                if idid==2 or idid==3:
                    flag=ID_PY_COMPLETE
                elif idid < 0:
                    raise MEXAX_Exception("MEXAX failed with flag IDID {IDID}".format(IDID=idid))
        else:   
            if normal_mode == 1: # intermediate output mode
                idid = 1
                while idid==1:
                    t,y,yprime,tf,info,idid,rwork,iwork = \
                       mexax.mexax(callback_prob,neq,ny,t,y,yprime,
                             tf,info,rtol,atol,rwork,iwork,jac_dummy)
                    
                    tlist.append(t)
                    ylist.append(y.copy())
                    ydlist.append(yprime.copy())
                    if idid==2 or idid==3:
                        flag=ID_PY_COMPLETE
                    elif idid < 0:
                        raise MEXAX_Exception("MEXAX  failed with flag IDID {IDID}".format(IDID=idid))
                    
            else:   # mode with output_list          
                output_list  = opts["output_list"]
                for tout in output_list: 
                    t,y,yprime,tout,info,idid,rwork,iwork = \
                      mexax.mexax(callback_prob,neq,ny,t,y,yprime, \
                             tout,info,rtol,atol,rwork,iwork,jac_dummy)
                    tlist.append(t)
                    ylist.append(y.copy())
                    ydlist.append(yprime.copy())
                    if idid > 0 and t >= tf:
                        flag=ID_PY_COMPLETE
                    elif idid < 0:
                        raise MEXAX_Exception("MEXAX  failed with flag IDID {IDID}".format(IDID=idid))                                      
 
        
        
        #Retrieving statistics
        self.statistics["nsteps"]      += iwork[10]
        self.statistics["nfcn"]        += iwork[11]
        self.statistics["njac"]        += iwork[12]
        self.statistics["errfail"]     += iwork[13]
        self.statistics["convfail"]         += iwork[14]
        
        return flag, tlist, ylist, ydlist
        
    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        self.log_message('Final Run Statistics: %s \n' % self.problem.name,        verbose)
        
        self.log_message(' Number of steps                          : '+str(self.statistics["nsteps"]), verbose)               
        self.log_message(' Number of function evaluations           : '+str(self.statistics["nfcn"]), verbose)
        self.log_message(' Number of Jacobian evaluations           : '+ str(self.statistics["njac"]), verbose)
        self.log_message(' Number of error test failures            : '+ str(self.statistics["errfail"]), verbose)
        self.log_message(' Number of Convergence Test Failures      : '+ str(self.statistics["convfail"]), verbose)
        
        self.log_message('\nSolver options:\n', verbose)
        self.log_message(' Solver                  : mexax ',          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self._compact_atol()),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)



    
    
    
    
    
       
        

        tresult=[]
        yresult=[]
        hresult=[]
        flag=[]
        #opts["output_list"]=0
        print opts
        output_index = opts["output_index"]
        output_list  = opts["output_list"][output_index:]   #[0.,1.,2.,3.,4.]

        kflag = 0#ID_PY_COMPLETE

        for tout in output_list:
            output_index += 1
            print tout

            result=mexax.mexax(self.f,t,y.copy(),tout, self.atol , self.options["maxh"] ,self.options["inith"],kflag)
            y=result[1]
            t=result[0]
            H=result[2]
            flag=result[3]
            tresult.append(t)
            hresult.append(H)
            yresult.append(y)

            #opts["output_index"] = output_index


            #Retrieving statistics

            self.statistics["nsteps"]        += mexax.statp.nstep
            self.statistics["nfcn"]          += mexax.statp.nfcn
            self.statistics["errfail"]       += mexax.statp.nrejct
            self.statistics["nlu"]           += mexax.statp.ndec



        return flag, tresult, yresult


    ##############

        #Checking return
        if flag == 1:
            flag = ID_PY_COMPLETE
        elif flag == 2:
            flag = ID_PY_EVENT
        else:
            raise Exception("mexax failed with flag %d"%flag)

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
        self.log_message(' Solver                  : mexax ' + self._type,          verbose)
        self.log_message(' Tolerances (absolute)   : ' + str(self._compact_atol()),  verbose)
        self.log_message(' Tolerances (relative)   : ' + str(self.options["rtol"]),  verbose)
        self.log_message('',                                                         verbose)




