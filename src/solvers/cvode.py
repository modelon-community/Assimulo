#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2011 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as N

from assimulo.explicit_ode import Explicit_ODE
from assimulo.sundials import Sundials, Sundials_Exception
    
class CVode(Explicit_ODE, Sundials):
    """
    Sundials CVode.
    """
    def __init__(self, problem, y0=None, t0=None, switches0=None, p0=None):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Explicit_Problem' class.
                              
                y0          
                            - Default 'None'. The initial values for the states.
                              If 'None', the initial values are retrieved from
                              the problem.y0. If set they override problem.y0
                            
                            - Should be a list of floats or a float.
                            
                                Example:
                                    y0 = [1.0, 0.0]
                                    
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
                raise Explicit_ODE_Exception('y0 must not be None.')
            
        Sundials.__init__(self, y0, 'CVode') #Creates an integrator
        Explicit_ODE.__init__(self, problem, y0, t0, switches0) #Calls the base class
        
        #Default values
        self.discr = 'Adams' #Setting default discretization to Adams
        self.iter = 'FixedPoint' #Setting default iteration to FixedPoint
        self.maxord = 12 #Setting default maxord to maximum
        self.initstep = 0.0 #Setting the initial step to be estimated
        self.maxh = 0.0 #Setting the maximum absolute step length to infinity
        self.atol = 1.0e-6 #Absolute tolerance
        self.rtol = 1.0e-6 #Relative tolerance
        self.pretype = 'PREC_NONE' #The preconditioner used (if any)
        self.maxkrylov = 5 #Max number of krylov subspace (if any)
        self.linearsolver = 'DENSE' #The linear solver to be used.
        self.problem_data = {}
        
        
        #Determine if we have an event function and sets the integration data
        if hasattr(problem, 'state_events'):
            self.state_events = self._problem.state_events #problem.state_events
            self.num_state_events=len(self.state_events(self._problem.t0,self._problem.y0,self._problem.switches0))
            self.problem_data['ROOT'] = self.state_events
            self.problem_data['dimRoot'] = self.num_state_events
        else:
            self.num_state_events=0
        
        #Determine if we have a user supplied jacobian
        if hasattr(self._problem, 'jac'):
            if self.switches == None:
                trial = self._problem.jac(self._problem.t0,self._problem.y0)
            else:
                trial = self._problem.jac(self._problem.t0,self._problem.y0, self.switches)
            if trial.shape != (len(self._problem.y0),len(self._problem.y0)):
                raise Explicit_ODE_Exception('The Jacobian must be a numpy matrix of size len(f)*len(f).')
            self.jac = self._problem.jac    
            self.Integrator.jacobian = True
            self.usejac = True
            self.problem_data['JAC']=self.jac
        else:
            self.Integrator.jacobian = False
            self.usejac = False
        
        #Look for the Jacobian times vector function
        if hasattr(problem, 'jacv'):
            self.usejac = True
            self.jacv = self._problem.jacv
            self.problem_data['JACV'] = self.jacv
        
        self.problem_data['RHS']=self.f
        self.problem_data['dim']=len(self._problem.y0)
        
        if hasattr(problem, 'completed_step'):
            self._completed_step = True
            self.Integrator.comp_step = True
        
        #Sensitivity
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
        
        #Defaul values
        if sens:
            if hasattr(problem, 'pbar'):
                self.pbar = self._problem.pbar
            else:
                self.pbar = N.abs(self._problem.p0)
            if hasattr(self._problem, 'yS0'):
                self.yS0 = self._problem.yS0
        
        # 
        # TEST METHODS
        #
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

    def _assert_return_types(self, root = None, jac = None, sw = None,sens = False):
        """
        Tests the provided user methods for correct return types. The
        methods should all return numpy arrays(matrix) of floats.
        """
        if sw == None and sens == False:
            testf = self.f(self._problem.t0, self.y_cur)
        elif sw == None and sens == True:
            testf = self.f(self._problem.t0, self.y_cur, self.p)
        elif sw and sens == False:
            testf = self.f(self._problem.t0, self.y_cur, sw)
        else:
            testf = self.f(self._problem.t0 , self.y_cur, sw=sw, p=self.p)
            
        if not isinstance(testf, N.ndarray) or testf.dtype != float:
            raise Explicit_ODE_Exception('The right-hand-side function must return a numpy array of floats.')
            
        if root != None:
            testr = self.state_events(self._problem.t0, self.y_cur, sw)
            if not isinstance(testr, N.ndarray) or testr.dtype != float:
                raise Explicit_ODE_Exception('The state event function must return a numpy array of floats.')
            
        if jac != None:
            if sw == None:
                testj = self.jac(self._problem.t0, self.y_cur)
            else:
                testj = self.jac(self._problem.t0, self.y_cur, sw)
            if not isinstance(testj, N.ndarray) or testj.dtype != float:
                raise Explicit_ODE_Exception('The Jacobian function must return a numpy array of floats.')
    
    def _integrator(self,t,y,tfinal,dt):
        """
        Simulates the problem up until tfinal.
        """
        if self._flag_reset_statistics:
            self.Integrator.solver_stats = [0,0,0,0,0,0,0,0]
            self.Integrator.solver_sens_stats = [0,0,0,0,0,0]
            for x in self.Integrator.statistics.keys():
                self.Integrator.statistics[x] = 0
            self._flag_reset_statistics = False
        
        self.Integrator.store_cont = self.store_cont
        if self._flag_init:
            self.Integrator.store_statistics()
            self.Integrator.cvinit(t,y,self.maxsteps,self.verbosity, self.switches)
            
        return self.Integrator.run(t,tfinal,dt)
    
    def _set_discr_method(self,discr='Adams'):
        """
        This determines the discretization method.
        
            Parameters::
            
                discr   
                        - Default 'Adams', which indicates the use
                          of the Adams method. Can also be set to
                          'BDF' which indicates the use of the BDF
                          method.
                
            Example::
                
                discr = 'BDF'
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        if discr=='BDF':
            self.Integrator.discr=2
            if self.maxord > 5:
                self.maxord = 5
        elif discr=='Adams':
            if self.Integrator.discr != 1 and self.maxord == 5:
                self.Integrator.discr=1
                self.maxord = 12
            self.Integrator.discr=1
            #self.maxord = 12
        else:
            raise Sundials_Exception('Discretization method must be either Adams or BDF')
            
    def _get_discr_method(self):
        """
        This determines the discretization method.
        
            Parameters::
            
                discr   
                        - Default 'Adams', which indicates the use
                          of the Adams method. Can also be set to
                          'BDF' which indicates the use of the BDF
                          method.
                
                    Example:
                        discr = 'BDF'
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        discr='Adams'
        if self.Integrator.discr==2:
            discr='BDF'
        return discr

    discr= property(_get_discr_method,_set_discr_method)
    
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
            raise Explicit_ODE_Exception('The initial step must be an integer or float.')
        
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
    
    def _set_iter_method(self,iter='FixedPoint'):
        """
        This determines the iteration method that is be used by the
        solver.
        
            Parameters::
            
                iter    
                        - Default 'FixedPoint', which indicates the
                          use of a fixedpoint iteration method. Can
                          also be set to 'Newton' which indicates
                          the use of a Newton method.
                          
                            Example:
                                iter = 'Newton'
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        if iter=='Newton':
            self.Integrator.iter=2
        elif iter=='FixedPoint':
            self.Integrator.iter=1
        else:
            raise Sundials_Exception('Iteration method must be either FixedPoint or Newton')
    
    def _get_iter_method(self):
        """
        This determines the iteration method that is be used by the
        solver.
        
            Parameters::
            
                iter    
                        - Default 'FixedPoint', which indicates the
                          use of a fixedpoint iteration method. Can
                          also be set to 'Newton' which indicates
                          the use of a Newton method.
                          
                            Example:
                                iter = 'Newton'
        
        See SUNDIALS CVODE documentation 2.1 for more details.
        """
        iter='FixedPoint'
        if self.Integrator.iter==2:
            iter='Newton'
        return iter
        
    iter = property(_get_iter_method,_set_iter_method)        
    
    def _set_max_ord(self,maxord):
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - Default '12', which is the maximum for the
                          Adams method, which is also default. For the
                          BDF method the maximum order is 5. 'maxord'
                          can be set in an interval from 1 to the
                          maximum order allowed.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
        """
        if not isinstance(maxord,int):
            raise Sundials_Exception('The maximal order must be an integer.')
        if self.Integrator.discr == 1: #Adams
            if maxord > 12:
                if self.verbosity > self.QUIET:
                    print 'The set maximal order is greater than that of the Adams. ' \
                           'Setting the maximal order to Adams maximum.'
                self.__maxord=12
            elif maxord < 1:
                if self.verbosity > self.QUIET:
                    print 'The set maximal order is lower than the minimum of Adams. ' \
                          'Setting the maximal order to Adams minimum.'
                self.__maxord=1
            else:
                self.__maxord=maxord
        if self.Integrator.discr == 2: #BDF
            if maxord > 5:
                if self.verbosity > self.QUIET:
                    print 'The set maximal order is greater than that of the BDF. ' \
                           'Setting the maximal order to BDF maximum.'
                self.__maxord=5
            elif maxord < 1:
                if self.verbosity > self.QUIET:
                    print 'The set maximal order is lower than the minimum of BDF. ' \
                          'Setting the maximal order to BDF minimum.'
                self.__maxord=1
            else:
                self.__maxord=maxord
        
        self.Integrator.maxord = maxord #Sets the maximum order to the solver
    
    def _get_max_ord(self):
        """
        This determines the maximal order that is be used by the solver.
        
            Parameters::
            
                maxord  
                        - Default '12', which is the maximum for the
                          Adams method, which is also default. For the
                          BDF method the maximum order is 5. 'maxord'
                          can be set in an interval from 1 to the
                          maximum order allowed.
                
                        - Should be an integer.
                        
                            Example:
                                maxord = 3
    
        
        An input value greater than the maximal order will result in the 
        maximum value.
        """
        return self.__maxord

    maxord=property(_get_max_ord,_set_max_ord)
    
    def _get_pre_type(self):
        """
        Specifies the preconditioning type. 
        
            Parameters::
            
                ptype
                        - Default 'PREC_NONE' which is the currently
                          only supported.
                          
        """
        return self.__pretype
        
    def _set_pre_type(self, ptype):
        """
        Specifies the preconditioning type. 
        
            Parameters::
            
                ptype
                        - Default 'PREC_NONE' which is the currently
                          only supported.
                          
        """
        if isinstance(ptype, str):
            if ptype.upper() == 'PREC_NONE':
                self.__pretype = ptype.upper()
                self.Integrator.pretype = 0
            else:
                raise Explicit_ODE_Exception('"PREC_NONE" is the only supported option.')
        else:
            raise Explicit_ODE_Exception('Must be a string, "PREC_NONE".')
    
    pretype=property(_get_pre_type,_set_pre_type)
    
    def _get_max_krylov(self):
        """
        Maximum dimension of the Krylov subspace to be used.
        
            Parameters::
            
                maxkrylov
                        - Default 5
                        
                        - Should be an integer.
        """
        return self.__maxkrylov
    
    def _set_max_krylov(self, mkrylov):
        """
        Maximum dimension of the Krylov subspace to be used.
        
            Parameters::
            
                maxkrylov
                        - Default 5
                        
                        - Should be an integer.
        """
        try:
            mkrylov = int(mkrylov)
            self.__maxkrylov = mkrylov
            self.Integrator.max_krylov = mkrylov
        except ValueError:
            raise Explicit_ODE_Exception('maxkrylov must be convertable to an integer.')
    
    maxkrylov = property(_get_max_krylov, _set_max_krylov)
    
    def _get_linear_solver(self):
        """
        Specifies the linear solver to be used.
        
            Parameters::
            
                linearsolver
                        - Default 'DENSE'. Can also be 'SPGMR'.
        """
        return self.__linearsolver
        
    def _set_linear_solver(self, lsolver):
        """
        Specifies the linear solver to be used.
        
            Parameters::
            
                linearsolver
                        - Default 'DENSE'. Can also be 'SPGMR'.
        """
        if isinstance(lsolver, str):
            if lsolver.upper() == 'DENSE':
                self.__linearsolver = lsolver.upper()
                self.Integrator.linear_solver = lsolver.upper()
            elif lsolver.upper() == 'SPGMR':
                self.__linearsolver = lsolver.upper()
                self.Integrator.linear_solver = lsolver.upper()
            else:
                raise Explicit_ODE_Exception('The linearsolver must be either "DENSE" or "SPGMR".')
        else:
            raise Explicit_ODE_Exception('The linearsolver must be a string.')
    
    linearsolver = property(_get_linear_solver, _set_linear_solver)
    
    def interpolate(self, t, k):
        """            
        Calls Sundials internal function CVodeGetDky that computes the interpolated 
        values of the k-th derivative of y for any value of t in the last internal 
        step taken by CVode.
        
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
            raise Explicit_ODE_Exception('t must be convertable to a float.')
        try:
            k = int(k)
        except (TypeError, ValueError):
            raise Explicit_ODE_Exception('k must be convertable to an integer.')
            
        return self.Integrator.interpolate(t,k)
            
    
    def print_statistics(self):
        """
        Prints the run-time statistics for the problem.
        """
        print 'Final Run Statistics: %s \n' % self.problem_name
        
        statistics = self.stats
        if statistics!= None:
            print ' Number of Steps                          :', statistics[0]                         
            print ' Number of Function Evaluations           :', statistics[1]
            if self.linearsolver == "SPGMR":
                print ' Number of Jacobian*Vector Evaluations    :', self.Integrator.statistics["JVEC"]
                print ' Number of F-Evals During Jac*Vec-Evals   :', self.Integrator.statistics["RHSJVEC"]
            else:     
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
            print ' Solver                  :  CVode'
            print ' Linear Multistep Method : ', self.discr
            print ' Nonlinear Solver        : ' ,self.iter
            print ' Maxord                  : ' ,self.maxord
            print ' Tolerances (absolute)   : ' ,self.atol
            print ' Tolerances (relative)   : ' ,self.rtol
            print ''
        else:
            print 'No statistics available.'
    
    @property
    def is_disc(self):
        """Method to test if we are at an event."""
        return self.t_cur==self.Integrator.event_time
        
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
        
    def echo_options(self):
        """
        Echo the solver options.
        """
        print 'Solver options:\n'
        print ' Solver                    :  CVode'
        print ' Linear Multistep Method   : ' ,self.discr
        print ' Nonlinear Solver          : ' ,self.iter
        print ' Maxord                    : ' ,self.maxord
        print ' Maximum step-size         : ' ,self.maxh
        print ' Initial step-size         : ' ,self.initstep
        print ' Maximum number of steps   : ' ,self.maxsteps
        print ' Use user-defined Jacobian : ' ,self.usejac
        print ' Tolerances (relative)     : ' ,self.rtol
        print ' Tolerances (absolute)     : ' ,self.atol
        
        

