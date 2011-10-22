        
class IDA(Implicit_ODE, Sundials):
    """
    Sundials IDA.
    """
    
    def __init__(self, problem, y0=None, yd0=None, t0=None, switches0=None, p0=None):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Implicit_Problem' class.
                              
                y0          
                            - Default 'None'. The initial values for the states.
                              If 'None', the initial values are retrieved from
                              the problem.y0. If set they override problem.y0
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    y0 = [1.0, 0.0]
                                    
                yd0         
                            - Default 'None'. The initial values for the state
                              derivatives. If 'None', the initial values are
                              retrieved from the problem.yd0. If set they
                              override problem.yd0.
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    yd0 = [0.0, 0.0]
                                
                t0          
                            - Default 'None'. The initial time. If 'None'. the
                              initial time are retrieved from the problem.t0.
                              If set it override problem.t0. If NOT set and NOT
                              defined in problem.t0, t0 is set to 0.0.
                            
                            - Should be a float.
                            
                                Example:
                                    t0 = 1.0
                                    
                switches0   
                            - Default 'None'. Only used for hybrid (discontinuous)
                              systems. If 'None', the switches are retrieved from
                              the problem.switches0. If set, they override the
                              problem.switches0.
                            
                            - Should be a list of booleans.
                            
                                Example:
                                    switches0 = [True, False]
            
        """
        if y0 == None:
            if hasattr(problem, 'y0'):
                y0 = problem.y0
            else:
                raise Implicit_ODE_Exception('y0 must not be None.')
                
        Sundials.__init__(self, y0, 'IDA') #Creates an integrator
        Implicit_ODE.__init__(self, problem, y0, yd0, t0) #Calls the base class

        if hasattr(self._problem, 'algvar'): #Check if the algebraic components are defined in the Problem specifications
            self.algvar = self._problem.algvar
        else:
            self.algvar = [1.0]*len(self.y_cur) #No algebraic variables are set
        
        self.problem_data = {}
        
        if hasattr(self._problem, 'switches0') and switches0 == None:
            switches0 = self._problem.switches0
        
        if isinstance(switches0, list):
            for x in switches0:
                if not isinstance(x, bool):
                    raise Implicit_ODE_Exception('Switches must be a list of booleans.')
        elif switches0 is not None:
            raise Implicit_ODE_Exception('Switches must be a list of booleans.')
        
        self.switches = switches0
        self._problem.switches0 = switches0
        
        #Determine if we have an state event function and sets the integration data
        if hasattr(problem, 'state_events'):
            self.state_events = self._problem.state_events
            self.num_state_events=len(self.state_events(self._problem.t0,self._problem.y0,self._problem.yd0,self._problem.switches0))
            self.problem_data['ROOT'] = self.state_events
            self.problem_data['dimRoot'] = self.num_state_events
        else:
            self.num_state_events=0
        
        #Determine if we have a user supplied jacobian
        if hasattr(self._problem, 'jac'):
            if self.switches == None:
                trial = self._problem.jac(1,self._problem.t0,self.y_cur,self.yd_cur)
            else:
                trial = self._problem.jac(1,self._problem.t0,self.y_cur,self.yd_cur, self.switches)
            if trial.shape != (len(self.y_cur),len(self.y_cur)):
                raise Implicit_ODE_Exception('The Jacobian must be a numpy matrix of size len(f)*len(f).')
            
            self.jac = self._problem.jac
            self.Integrator.jacobian = True
            self.usejac = True
            self.problem_data['JAC']=self.jac
        else:
            self.Integrator.jacobian = False
            self.usejac = False
        
        self.problem_data['RHS']=self.res_fcn
        self.problem_data['dim']=len(self._problem.y0)
        
        #Check for sensitivites
        #----------------
        sens = False
        
        if p0 == None:
            if hasattr(problem, 'p0'):
                self.p = problem.p0
                sens = True
        else:
            self.p = p0
            sens = True

        if sens:
            self._problem.p0 = self.p
            #Set information to the solver IDAS
            self.Integrator.p = N.array(self.p)
            self.problem_data['dimSens'] = len(self.p)
            
        else:
            self.problem_data['dimSens'] = 0
        #-------------End Sensitivity initiation
        
        if hasattr(problem, 'completed_step'):
            self._completed_step = True
            self.Integrator.comp_step = True
        
        #Default values
        self.tout1 = 0.001
        self.lsoff = False #Use LineSearch
        self.suppress_alg = False #Don't suppres algebraic variables
        self.initstep = 0.0 #Setting the initial step to be estimated
        self.maxord = 5 #Maximal order is set to max
        self.maxh = 0.0 #Setting the maximum absolute step length to infinity
        self.atol = 1.0e-6 #Absolute tolerance
        self.rtol = 1.0e-6 #Relative tolerance
        if sens:
            if hasattr(self._problem, 'yS0'):
                self.yS0 = self._problem.yS0
            if hasattr(problem, 'pbar'):
                self.pbar = self._problem.pbar
            else:
                self.pbar = N.abs(self._problem.p0)
        
        # TEST METHODS
        try:
            jt = self.problem_data['JAC']
        except KeyError:
            jt = None
        try:
            rt = self.problem_data['ROOT']
        except KeyError:
            rt = None
        self._assert_return_types(rt,jt,self.switches, sens)
        
        #Sets the problem data to Sundials
        self.Integrator.set_problem_info(**self.problem_data)
    
    def _assert_return_types(self, root = None, jac = None, sw = None, sens = False):
        """
        Tests the provided user methods for correct return types. The
        methods should all return numpy arrays(matrix) of floats.
        """
        if sw == None and sens == False:
            testf = self.res_fcn(self._problem.t0, self.y_cur, self.yd_cur)
        elif sw==None and sens == True:
            testf = self.res_fcn(self._problem.t0, self.y_cur, self.yd_cur, self.p)
        elif sw and sens == False:
            testf = self.res_fcn(self._problem.t0, self.y_cur, self.yd_cur, sw)
        else:
            testf = self.res_fcn(self._problem.t0, self.y_cur, self.yd_cur, sw=sw, p=self.p)
            
        if not isinstance(testf, N.ndarray) or testf.dtype != float:
            raise Implicit_ODE_Exception('The residual function must return a numpy array of floats.')
            
        if root != None:
            testr = self.state_events(self._problem.t0, self.y_cur, self.yd_cur, sw)
            if not isinstance(testr, N.ndarray) or testr.dtype != float:
                raise Implicit_ODE_Exception('The state event function must return a numpy array of floats.')
            
        if jac != None:
            if sw == None:
                testj = self.jac(1.0, self._problem.t0, self.y_cur, self.yd_cur)
            else:
                testj = self.jac(1.0, self._problem.t0, self.y_cur, self.yd_cur, sw)
            if not isinstance(testj, N.ndarray) or testj.dtype != float:
                raise Implicit_ODE_Exception('The Jacobian function must return a numpy array of floats.')
    
    def _set_calcIC_tout1(self, tout1):
        """
        Sets the value used in the internal Sundials function
        for determine initial conditions. This value is needed
        in order to manipulate the existing machinery to be used
        in determining the initial conditions.
        
            Parameters::
            
                tout1       
                            - Default '0.001'.
                            
                            - Should be a float.
                            
                                Example:
                                    tout1 = 0.01
        """
        try:
            tout1 = float(tout1)
        except (ValueError,TypeError):
            raise Implicit_ODE_Exception('tout1 must be an integer or float.')
        self.__tout1 = tout1
        
    def _get_calcIC_tout1(self):
        """
        Sets the value used in the internal Sundials function
        for determine initial conditions. This value is needed
        in order to manipulate the existing machinery to be used
        in determining the initial conditions.
        
            Parameters::
            
                tout1       
                            - Default '0.001'.
                            
                            - Should be a float.
                            
                                Example:
                                    tout1 = 0.01
        """
        return self.__tout1

    tout1=property(_get_calcIC_tout1,_set_calcIC_tout1)
    
    def _set_usejac(self, jac):
        """
        This sets the option to use the user defined jacobian. If a
        user provided jacobian is implemented into the problem the
        default setting is to use that jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejac  
                        - True - use user defined jacobian
                          False - use an approximation
                    
                        - Should be a boolean.
                        
                            Example:
                                usejac = False
        """
        self.__usejac = bool(jac)
        self.Integrator.usejac = self.__usejac
    
    def _get_usejac(self):
        """
        This sets the option to use the user defined jacobian. If a
        user provided jacobian is implemented into the problem the
        default setting is to use that jacobian. If not, an
        approximation is used.
        
            Parameters::
            
                usejac  
                        - True - use user defined jacobian
                          False - use an approximation
                    
                        - Should be a boolean.
                        
                            Example:
                                usejac = False
        """
        return self.__usejac
    
    usejac = property(_get_usejac,_set_usejac)
    
    def _set_lsoff(self, lsoff):
        """
        Boolean value to turn OFF Sundials LineSearch when calculating
        initial conditions.
            
            Parameters::
            
                lsoff   
                        - Default 'False'. False indicates the use of
                          linesearch.

                        - Should be a boolean.
                        
                            Example:
                                lsoff = True
        """
        if not isinstance(lsoff, bool):
            raise Implicit_ODE_Exception('lsoff must be a boolean.')
        self.__lsoff = lsoff
        
    def _get_lsoff(self):
        """
        Boolean value to turn OFF Sundials LineSearch when calculating
        initial conditions.
            
            Parameters::
            
                lsoff   
                        - Default 'False'. False indicates the use of
                          linesearch.

                        - Should be a boolean.
                        
                            Example:
                                lsoff = True
        """
        return self.__lsoff
    
    lsoff = property(_get_lsoff, _set_lsoff)
    
    def make_consistent(self, method):
        """
        Directs IDA to try to calculate consistant initial conditions.
            
            Parameters::
            
                method  
                        - 'IDA_YA_YDP_INIT'
                                - This tries to calculate the
                                  algebraic components of y and the differential
                                  components of yd given the differential components
                                  of y. The algebraic components of y must have been
                                  specified with the property 'algvar'. The property
                                  'IDA.tout1' is  used in the calculations. It 
                                  represents the the next output point.
                                
                        - 'IDA_Y_INIT' 
                                - This tries to calculate all components
                                  of y given yd.
            
        See SUNDIALS IDA documentation 4.5.4 for more details.
        """
        self.Integrator.idinit(self.t_cur, self.y_cur, self.yd_cur, self.maxsteps,self.verbosity, self.switches)
        
        if method == 'IDA_YA_YDP_INIT':
            [flag, y, yd] = self.Integrator.calc_IC(method,self.tout1,self.lsoff)
        elif method == 'IDA_Y_INIT':
            [flag, y, yd] = self.Integrator.calc_IC(method,self.tout1,self.lsoff)
        else:
            raise Sundials_Exception('Make consistency must be one of the following:'\
                                            ' IDA_YA_YDP_INIT or IDA_Y_INIT')
        
        if flag < 0:
            raise Sundials_Exception('Calculation of initial conditions failed. IDA returned flag %d'%flag)
        else:
            self.y_cur = y
            self.yd_cur = yd
        
        return [self.y_cur, self.yd_cur]
        
    def make_consistency(self, method):
        """
        This method is being renamed to 'make_consistent', see docstring on that method.
        This method is being removed in version 1.4.
        """
        warnings.warn('Please use "make_consistent" instead of "make_consistency".')
        return self.make_consistent(method)
        
    
    def _integrator(self,t,y,yd,tfinal,nt=0):
        """
        Simulates the problem up until tfinal.
        """
        if self._flag_reset_statistics:
            self.Integrator.solver_stats = [0,0,0,0,0,0,0,0]
            self.Integrator.solver_sens_stats = [0,0,0,0,0,0]
            self._flag_reset_statistics = False
        
        self.Integrator.store_cont = self.store_cont
        if self._flag_init:
            self.Integrator.store_statistics()
            self.Integrator.idinit(t,y,yd, self.maxsteps, self.verbosity, self.switches)
        
        return self.Integrator.run(t,tfinal,nt)
    
    def _set_initial_step(self, initstep):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.0', which result in that the
                              the initial step is approximated.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        try:
            initstep = float(initstep)
        except (ValueError, TypeError):
            raise Implicit_ODE_Exception('The initial step must be an integer or float.')
        
        self.__initstep = initstep
        self.Integrator.inith = self.__initstep
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.0', which result in that the
                              the initial step is approximated.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        return self.__initstep
        
    initstep = property(_get_initial_step,_set_initial_step)
    
    def interpolate(self, t, k):
        """        
        Calls Sundials internal function IDAGetDky that computes the interpolated 
        values of the k-th derivative of y for any value of t in the last internal 
        step taken by IDA.
        
            Parameters::
            
                t
                    - Must be within tn − hu ≤ t ≤ tn  where tn denotes the current
                      internal time reached, and hu is the last internal step size used successfully.
                      
                    - Must be a float.
                      
                k
                    - Must be non-negative and samller than the last internal order used.
                    
                    - Must be an integer.
        """
        try:
            t = float(t)
        except (TypeError, ValueError):
            raise Implicit_ODE_Exception('t must be convertable to a float.')
        try:
            k = int(k)
        except (TypeError, ValueError):
            raise Implicit_ODE_Exception('k must be convertable to an integer.')
            
        return self.Integrator.interpolate(t,k)

    def _set_max_ord(self,maxord):
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - Default '5', which is the maximum.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
        """
        if not isinstance(maxord,int):
            raise Sundials_Exception('The maximal order must be an integer.')
        if maxord > 5:
            self.print_verbos(['The given maximal order is greater than the maximal order of the BDF. ', 
                      'Setting the maximal order to BDF maximum.'],self.WHISPER)
            self.__maxord=5
        elif maxord < 1:
            self.print_verbos(['The given maximal order is lower than the minimum of BDF. ', 
                      'Setting the maximal order to BDF minimum.'],self.WHISPER)
            self.__maxord=1
        else:
            self.__maxord=maxord
        
        self.Integrator.maxord = maxord #Sets the maximum order to the solver
            
    def _get_max_ord(self):
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - Default '5', which is the maximum.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
        """
        return self.__maxord

    maxord=property(_get_max_ord,_set_max_ord)
    
    def _set_suppress_alg(self,suppress_alg):
        """
        A boolean flag which indicates that the error-tests are 
        suppressed on algebraic variables. The algebraic variables
        are defined by setting the property 'algvar'.
        
            Parameters::
            
                suppress_alg    
                                - Default 'False'.
                
                                - Should be a boolean.
                                
                                    Example:
                                        suppress_alg = True
                                        
        See SUNDIALS IDA documentation 4.5.7 'IDASetSuppressAlg' 
        for more details.
        """
        if not isinstance(suppress_alg, bool):
            raise Sundials_Exception('Option suppress algebraic variables' \
                                              ' must be defined by a boolean.')
        self.Integrator.suppress_alg=suppress_alg
    def _get_suppress_alg(self):
        """
        A Boolean flag which indicates that the error-tests are 
        suppressed on algebraic variables. The algebraic variables
        are defined by setting the property 'algvar'.
        
            Parameters::
            
                suppress_alg    
                                - Default 'False'.
                
                                - Should be a boolean.
                                
                                    Example:
                                        suppress_alg = True
                                        
        See SUNDIALS IDA documentation 4.5.7 'IDASetSuppressAlg' 
        for more details.
        """
        return self.Integrator.suppress_alg    

    suppress_alg=property(_get_suppress_alg,_set_suppress_alg)
    
    def _set_algvar(self,algvar):
        """
        A list for defining which variables are differential and
        which are algebraic.
        This list is used when excluding algebraic variables from the error test
        by setting suppress_alg=True  and it is used, when computing consistent initial 
        values using the method make_consistency
        
            Parameters::
            
                algvar  
                        - The value True(1.0) indicates a differential
                          variable and the value False(0.0) indicates an
                          algebraic variable.
                          
                        - Should be a list or a numpy vector (ndarray)
                        
                            Example:
                                algvar = [1.0, 0.0, 1.0]
                                algvar = [True, False, True]
                                algvar = [1,0,1]
                                
        """
        try:
            algvar = N.array(algvar,dtype='float')
        except (ValueError, TypeError):
            raise Sundials_Exception('algvar needs to be a list or a numpy vector.')
            
        try:
            size = len(algvar)
        except TypeError:
            raise Sundials_Exception('algvar needs to be a list or a numpy vector.')
        
        if len(algvar) != self.Integrator.dim:
            raise Sundials_Exception('When setting the algebraic variables, the' \
                                ' vector must be of the same size as the problem dimension.')
        
        for alg in algvar:
            if alg != 1.0 and alg != 0.0:
                raise Sundials_Exception('The vector must consist of 1.0 or 0.0 .')
        self.Integrator.algvar=algvar

    def _get_algvar(self):
        """
        A vector for defining which variables are differential and
        which are algebraic.
        
            Parameters::
            
                algvar  
                        - The value True(1.0) indicates a differential
                          variable and the value False(0.0) indicates an
                          algebraic variable.
                          
                        - Should be a list or a numpy vector (ndarray)
                        
                            Example:
                                algvar = [1.0, 0.0, 1.0]
                                algvar = [True, False, True]
                                algvar = [1,0,1]
                                
        """
        return self.Integrator.algvar    

    algvar=property(_get_algvar,_set_algvar)
    
    def simulation_complete(self):
        """
        Method which returns a boolean value determining if the
        simulation completed to tfinal. Used for determining 
        time-events.
        
            Returns::
            
                sim_complete
                            - Boolean value
                                -True for success
                                -False for not complete
        """
        return self.Integrator.sim_complete
    
    def print_statistics(self,minimal_verbosity=0):
        """
        Prints the run-time statistics for the problem
        if self.verbosity >= minimal_verbosity
        
        """
        if self.verbosity < minimal_verbosity:
            return None
        print 'Final Run Statistics: %s \n' % self.problem_name
        statistics = self.stats
        
        if statistics!= None:
            print ' Number of Steps                          :', statistics[0]                         
            print ' Number of Function Evaluations           :', statistics[1]     
            print ' Number of Jacobian Evaluations           :', statistics[2]        
            print ' Number of F-Eval During Jac-Eval         :', statistics[3]    
            print ' Number of Root Evaluations               :', statistics[4]       
            print ' Number of Error Test Failures            :', statistics[5]       
            print ' Number of Nonlinear Iterations           :', statistics[6]     
            print ' Number of Nonlinear Convergence Failures :', statistics[7]
                
            if self.problem_data['dimSens'] > 0: #Senstivity calculations is on
                sens_stats = self.get_sensitivity_statistics()
                
                print '\nSensitivity Statistics:\n'
                print ' Number of Sensitivity Calculations             :', sens_stats[0]
                print ' Number of F-Evals Due to Finite Approximation  :', sens_stats[1]
                print ' Number of Local Error Test Failures            :', sens_stats[2]
                print ' Number of Linear Setups                        :', sens_stats[3]
                print ' Number of Nonlinear iterations                 :', sens_stats[4]
                print ' Number of Nonlinear Convergance Failures       :', sens_stats[5]
                
                print '\nSensitivity options:\n'
                print ' Method                   : ' ,self.sensmethod
                print ' Difference quotient type : ' ,self.dqtype
                print ' Suppress Sens            : ' ,self.suppress_sens
            
            print '\nSolver options:\n'
            print ' Solver                :  IDA (BDF)'
            print ' Maxord                : ' ,self.maxord
            print ' Suppress Alg          : ' ,self.suppress_alg
            print ' Tolerances (absolute) : ' ,self.atol
            print ' Tolerances (relative) : ' ,self.rtol
            print ''
        else:
            print 'No statistics available.'
    
    @property
    def is_disc(self):
        """Method to test if we are at an event."""
        return self.t_cur==self.Integrator.event_time
    
    def echo_options(self):
        """
        Echo the solver options.
        """
        print 'Solver options:\n'
        print ' Solver                        :  IDA'
        print ' Linear Multistep Method       :  BDF'
        print ' Maxord                        : ' ,self.maxord
        print ' Maximum step-size             : ' ,self.maxh
        print ' Initial step-size             : ' ,self.initstep
        print ' Maximum number of steps       : ' ,self.maxsteps
        print ' Use user-defined Jacobian     : ' ,self.usejac
        print ' Suppress algebraic components : ' ,self.suppress_alg
        print ' Differential components       : ' ,self.algvar
        print ' Tolerances (relative)         : ' ,self.rtol
        print ' Tolerances (absolute)         : ' ,self.atol
