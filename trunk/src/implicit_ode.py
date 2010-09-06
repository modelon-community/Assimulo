#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
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

from ode import *
from problem import Implicit_Problem
from sundials import Sundials, Sundials_Exception
from assimulo.lib.radau_core import Radau_Common
import warnings

class Implicit_ODE_Exception(Exception):
    """An integrator exception"""
    pass 

class Implicit_ODE(ODE):
    """
    Baseclass for our implicit ODE integrators.
    """
    
    def __init__(self, problem, y0=None, yd0=None, t0=None):
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
                                    
        """
        ODE.__init__(self) #Sets general attributes
        
        if problem == None:
            raise Implicit_ODE_Exception('Problem cannot be None. It has be a subclass of a Implicit_Problem')
        
        if isinstance(problem, Implicit_Problem):
            self._problem = problem
            self.res_fcn = problem.f
            self.problem_name = problem.problem_name
        else:
            raise Implicit_ODE_Exception('The problem needs to be a subclass of a Implicit_Problem.')
        
        if y0 == None:
            if hasattr(problem, 'y0'):
                y0 = problem.y0[:]
            else:
                raise Implicit_ODE_Exception('y0 must be specified. Either in the problem or in the initialization')
        
        if yd0 == None:
            if hasattr(problem, 'yd0'):
                yd0 = problem.yd0[:]
            else:
                raise Implicit_ODE_Exception('yd0 must be specified. Either in the problem or in the initialization')
        
        try:
            if isinstance(y0, int) or isinstance(y0, float):
                y0 = [y0]
            if isinstance(yd0, int) or isinstance(yd0, float):
                yd0 = [yd0]
            self._problem.y0 = y0[:]
            self._problem.yd0 = yd0[:]
            self.y_cur = N.array(y0, dtype=float)
            self.yd_cur = N.array(yd0, dtype=float)
            self.y = [N.array(y0, dtype=float)]
            self.yd = [N.array(yd0, dtype=float)]
        except ValueError:
            raise Implicit_ODE_Exception('Initial values must be a scalar/list/array of type int or float.')
        
        if t0 == None:
            if hasattr(problem, 't0'):
                t0 = problem.t0
            else:
                t0 = 0.0
        
        try:
            if isinstance(t0, list):
                self.t_cur = float(t0[0])
            else:
                self.t_cur = float(t0)
            self._problem.t0 = self.t_cur
        except ValueError:
            raise Implicit_ODE_Exception('Initial time must be an integer or float.')
        
        if hasattr(self._problem, 'time_events'):
            self._time_function  = True
        
        self.t_cur  = N.array(self.t_cur)
        self.y_cur  = N.array(self.y_cur)
        self.yd_cur = N.array(self.yd_cur)
        
        self.t  = []
        self.y  = []
        self.yd = []
    
    
    def integrate(self, t, y, yd, tf,nt):
        pass
        
    def reset(self):
        """
        
        Resets the problem. If the problem is defined with a reset method, its called
        and then the method re_init. The re_init method is called with the initial
        values set in the problem, problem.t0, problem.y0 and problem.yd0
        
        """
        self._problem.reset()
        self.re_init(self._problem.t0, self._problem.y0, self._problem.yd0)
        
    def re_init(self,t0, y0, yd0):
        """
        Reinitiates the solver.
        
            Parameters::
                
                t0  - The initial time.
                y0  - The initial values for the states
                yd0 - The initial values for the state derivatives.
                
        See information in the __init__ method.
        """
        if len(self.y_cur) != len(y0) or len(self.yd_cur) != len(yd0):
            raise Implicit_ODE_Exception('y0/yd0 must be of the same length as the original problem.')
        
        Implicit_ODE.__init__(self, self._problem,y0,yd0,t0)

    def __call__(self, tfinal, ncp=0):
        """
        Calls the integrator to perform the simulation over the given time-interval.
        If a second call to simulate is performed, the simulation starts from the last
        given final time.
        
            Parameters::
            
                tfinal  
                        - Final time for the simulation
                
                        - Should be a float or integer greater than the initial time.
                        
                ncp     
                        - Default '0'. Number of communication points where the 
                          solution is returned. If '0', the integrator will return 
                          at its internal steps.
                          
                        - Should be an integer.
                          
                    Example:
                    
                        __call__(10.0, 100), 10.0 is the final time and 100 is the number
                                             communication points.
        """
        try:
            tfinal = float(tfinal)
            tfinal = N.array(tfinal)
        except ValueError:
            raise Implicit_ODE_Exception('Final time must be an integer or float.')
            
        if self.t_cur > tfinal:
            raise Implicit_ODE_Exception('Final time must be greater than start time.')
        if not isinstance(ncp, int):
            raise Implicit_ODE_Exception('Number of communication points must be an integer.')
        if ncp < 0:
            ncp = 0
            self.print_verbos(['Number of communication points must be a positive integer, setting' \
                        ' nt = 0.'], self.WHISPER)
        
        t0  = self.t_cur
        y0  = self.y_cur
        yd0 = self.yd_cur
        
        if ncp != 0 and self._completed_step:
            mode = 'SPECIAL'
            dist_space = [(x+1)*(tfinal-self.t_cur)/ncp for x in range(ncp+1)]
            dt = 0.0
        elif ncp != 0:
            dt = (tfinal-t0)/ncp
            mode = 'NORMAL'
        else:
            dt = 0.0
            mode = 'ONE_STEP'
        
        ncp_ori = ncp
        tfinal_ori = tfinal
        time_start = time.clock()
        
        self._problem.handle_result(self,t0,y0,yd0) #Logg the first point
        self._flag_init = True #Reinitiate the solver
        
        while self.t_cur < tfinal_ori:
            
            #Time event function is specified.
            if self._time_function:
                tevent = self._problem.time_events(self.t_cur, self.y_cur, self.yd_cur, self.switches)
                if tevent == None:
                    tfinal = tfinal_ori
                else:
                    tfinal = tevent if tevent < tfinal_ori else tfinal_ori
            
            solution = list(self.integrate(self.t_cur, self.y_cur, self.yd_cur, tfinal,dt))

            temp_t, temp_y, temp_yd = solution[-1]
            
            self.t_cur  = temp_t.copy()
            self.y_cur  = temp_y.copy()
            self.yd_cur = temp_y.copy()
            
            if mode == 'SPECIAL':
                while dist_space[0] <= self.t_cur:
                    self._problem.handle_result(self, dist_space[0], self.interpolate(dist_space[0],0),self.interpolate(dist_space[0],1))
                    last_logg = dist_space[0].copy()
                    dist_space.pop(0)
            else:
                for q in solution:
                    self._problem.handle_result(self,q[0],q[1],q[2])
                last_logg = self.t_cur

            #Check if there is a time event
            if tfinal == tfinal_ori:
                time_event = False
            else:
                time_event = self.simulation_complete()
            
            if self.is_disc or time_event: #Is discontinious?
                
                event_info = [[],time_event]
                if self.is_disc:
                    event_info[0] = self.disc_info[1]
                
                #Log the information
                self._log_event_info.append([self.t_cur, event_info])
                
                self.print_verbos(['A discontinuity occured at t = %e.'%self.t_cur,'\n',
                              'Current switches: ', self.switches,'\n',
                              'Event info: ', event_info],self.LOUD)

                self.print_statistics(self.SCREAM) #Prints statistics
                    
                self.print_verbos(['Calling problem specified event handling...'],self.LOUD)
                
                self._problem.handle_event(self, event_info) #self corresponds to the solver
                #self.event_iteration(event_info) #Handles the event iteration
                self._flag_init = True
            else:
                self._flag_init = False
                
            if self._completed_step: #If the option completed is set.
                self._flag_init = self._problem.completed_step(self) or self._flag_init
            
            if self._flag_init and last_logg == self.t_cur: #Logg after the event handling if there was a communication point there.
                self._problem.handle_result(self, self.t_cur, self.y_cur, self.yd_cur)
                
        #Simulation complete, call finalize
        self._problem.finalize(self)            
        
        time_stop = time.clock()
        
        self.print_statistics(self.NORMAL)
        self.print_verbos(['Elapsed simulation time:', time_stop-time_start, 'seconds.'],self.NORMAL)
        
        
        #return [self.t, self.y, self.yd]
        
        
    def plot(self, mask=None, der=False, **kwargs):
        """
        Plot the computed solution.
        
            Parameters::
            
                mask    
                        - Default 'None'. Used to determine which variables that is to be plotted.
                          Used as a list of integers, ones represents the variable that is to be
                          plotted and zeros that is not. 
                        
                        - Should be a list of integers.
                        
                            Example:
                                mask = [1,0] , plots the first variable.
                        
                der     
                        - Default 'False'. When 'True' plots the derivative variables also.
                        
                        - Should be a boolean.
                        
                            Example:
                                der = True
                **kwargs
                        - See http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
                          for information about the available options for **kwargs.
        """
        P.figure(1)
        if not mask:
            P.plot(self.t, self.y, **kwargs)
        else:
            if not isinstance(mask, list):
                raise Implicit_ODE_Exception('Mask must be a list of integers')
            if not len(mask)==len(self.y[-1]):
                raise Implicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                             'the number of variables.')
            for i in range(len(mask)):
                if mask[i]:
                    P.plot(self.t, N.array(self.y)[:,i], **kwargs)

        P.xlabel('time')
        P.ylabel('state')
        P.title(self.problem_name)

        
        if der and not mask:
            P.figure(2)
            P.plot(self.t, self.yd, **kwargs)
            P.xlabel('time')
            P.ylabel('state derivatives')
            P.title(self.problem_name)
        elif mask and der:
            P.figure(2)
            if not isinstance(mask, list):
                raise Implicit_ODE_Exception('Mask must be a list of integers')
            if not len(mask)==len(self.yd[-1]):
                raise Implicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                             'the number of variables.')
            for i in range(len(mask)):
                if mask[i]:
                    P.plot(self.t, N.array(self.yd)[:,i], **kwargs)
                    
            P.xlabel('time')
            P.ylabel('state derivatives')
            P.title(self.problem_name)
        
        
        P.show()

        
        
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
            self.Integrator.num_state_events=len(self.state_events(self._problem.t0,self._problem.y0,self._problem.yd0,self._problem.switches0))
            self._ROOT = [self.state_events, self.switches]
            self.problem_data['ROOT'] = self.state_events
            self.problem_data['dimRoot'] = self.Integrator.num_state_events
        else:
            self.Integrator.num_state_events=0
        
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
            self._RES = [self.res_fcn, self._problem.jac]
            self.problem_data['JAC']=self.jac
        else:
            self.Integrator.jacobian = False
            self.usejac = False
            self._RES = [self.res_fcn]
        
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
            self.Integrator.nbr_params = len(self.p)
            self.Integrator.p = N.array(self.p)
            self.problem_data['dimSens'] = len(self.p)
            self._RES = [len(self.p)]+self._RES #Indicator for Cython
        else:
            self._RES = [0]+self._RES #Indicator for Cython
        #-------------End Sensitivity initiation

        if hasattr(self, '_ROOT'):
            self.problem_spec = [self._RES, self._ROOT]
        else:
            self.problem_spec = [self._RES]
        
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
        
        try:
            sens = [self.problem_spec[0][0]]
        except AttributeError:
            sens = [0]
        
        if not bool(jac):
            self.Integrator.jacobian = False
            self._RES = sens+[self.res_fcn]
        else:
            self.Integrator.jacobian = True
            if not hasattr(self, 'jac'):
                raise Implicit_ODE_Exception('No jacobian defined.')
            self._RES = sens+[self.res_fcn, self.jac]
            
        if hasattr(self, '_ROOT'):
            self.problem_spec = [self._RES, self._ROOT]
        else:
            self.problem_spec = [self._RES]
    
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
        self.Integrator.idinit(self.t_cur, self.problem_spec, self.y_cur, self.yd_cur, self.maxord, self.maxsteps,self.initstep,self.maxh,self.verbosity, self.switches)
        
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
        
    
    def integrate(self,t,y,yd,tfinal,nt=0):
        """
        Simulates the problem up until tfinal.
        """
        self.Integrator.store_cont = self.store_cont
        if self._flag_init:
            self.Integrator.store_statistics()
            self.Integrator.idinit(t,self.problem_spec,y,yd,self.maxord, self.maxsteps, self.initstep, self.maxh, self.verbosity, self.switches)
        
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
            keys = statistics.keys()
            keys.sort()
            
            for x in keys:
                print ' %s = %s'%(x, statistics[x])
                
            if self.problem_spec[0][0]: #Senstivity calculations is on
                sens_stats = self.get_sensitivity_statistics()
                
                print '\nSensitivity Statistics:\n'
                print ' Number of Sensitivity Calculations             :', sens_stats[0]
                print ' Number of F-Evals Due to Finite Approximation  :', sens_stats[1]
                print ' Number of Local Error Test Failures            :', sens_stats[2]
                print ' Number of Linear Setups                        :', sens_stats[3]
                print ' Number of Nonlinear iterations                 :', sens_stats[4]
                print ' Number of Nonlinear Convergance Failures       :', sens_stats[5]
            
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
        
    def get_sensitivity_statistics(self):
        """
        Returns the sensitivity statistics.
        
            Returns::
            
                A list of statistics.
                
            Example::
            
                stats = IDA.get_sensitivity_statistics()
                
                stats[0] #Number of Sensitivity Calculations
                stats[1] #Number of F-Evals Due to Finite Approximation
                stats[2] #Number of Local Error Test Failures
                stats[3] #Number of Linear Setups            
                stats[4] #Number of Nonlinear iterations     
                stats[5] #Number of Nonlinear Convergance Failures
        """
        return self.Integrator.sens_stats

    
    #SENSITIVITY METHODS
    def interpolate_sensitivity(self, t, k, p=-1):
        """        
        Calls Sundials internal function IDAGetSensDky1 that computes the interpolated 
        values of the k-th derivative of the sensitivity p for any value of t in the 
        last internal step taken by IDA.
        
            Parameters::
            
                t
                    - Must be within tn − hu ≤ t ≤ tn  where tn denotes the current
                      internal time reached, and hu is the last internal step size used successfully.
                      
                    - Must be a float.
                      
                k
                    - Must be non-negative and samller than the last internal order used.
                    
                    - Must be an integer.
                    
                p
                    - An integer to determine of which parameter the solution is to be computed.
                    
                    - Default value -1 indicates that all is to be calculated.
                    
                    - Must be an integer.
                    
            Returns::
            
                A numpy array of the calculated sensitivites.
        """
        try:
            t = float(t)
        except (TypeError, ValueError):
            raise Implicit_ODE_Exception('t must be convertable to a float.')
        try:
            k = int(k)
        except (TypeError, ValueError):
            raise Implicit_ODE_Exception('k must be convertable to an integer.')
        try:
            p = int(p)
        except (TypeError, ValueError):
            raise Implicit_ODE_Exception('p must be convertable to an integer.')
            
        return self.Integrator.interpolate_sensitivity(t,k,p)
        
    def _set_DQtype(self, dqtype):
        """
        Specifies the difference quotient type in the sensitivity calculations
        and can be either 'IDA_CENTERED' or 'IDA_FORWARD'.
        
            Parameters::
            
                DQtype 
                        - A string of either 'IDA_CENTERED' or 'IDA_FORWARD'
                        - Default 'IDA_CENTERED'
                        
            Returns::
            
                The current value of DQtype.
        
        See SUNDIALS IDAS documentation 5.2.6 'IDASetSensDQMethod' 
        """
        if not isinstance(dqtype, str):
            raise Implicit_ODE_Exception('DQtype must be string.')
        
        if dqtype.upper() == 'IDA_CENTERED':
            self.Integrator.DQtype = 1
        elif dqtype.upper() == 'IDA_FORWARD':
            self.Integrator.DQtype = 2
        else:
            raise Implicit_ODE_Exception('DQtype must be either "IDA_CENTERED" or "IDA_FORWARD".')
            
    def _get_DQtype(self):
        """
        Specifies the difference quotient type in the sensitivity calculations
        and can be either 'IDA_CENTERED' or 'IDA_FORWARD'.
        
            Parameters::
            
                DQtype 
                        - A string of either 'IDA_CENTERED' or 'IDA_FORWARD'
                        - Default 'IDA_CENTERED'
                        
            Returns::
            
                The current value of DQtype.
        
        See SUNDIALS IDAS documentation 5.2.6 'IDASetSensDQMethod' 
        """
        if self.Integrator.DQtype == 1:
            return 'IDA_CENTERED'
        elif self.Integrator.DQtype == 2:
            return 'IDA_FORWARD'
        else:
            raise Implicit_ODE_Exception('Unknown value of DQtype.')
    
    DQtype = property(_get_DQtype, _set_DQtype)
    
    
    def _set_DQrhomax(self, dqrhomax):
        """
        Specifies the selection parameters used in deciding switching between a simultaneous
        or separate approximation of the two terms in the sensitivity residual.
        
            Parameters::
            
                DQrhomax
                        - A postive float.
                        - Default 0.0
                        
            Returns::
            
                The current value of DQrhomax (float)
        
        See SUNDIALS IDAS documentation 5.2.6 'IDASetSensDQMethod' 
        """
        try:
            dqrhomax = float(dqrhomax)
            if dqrhomax < 0.0:
                raise Implicit_ODE_Exception('DQrhomax must be a positive float.')
        except (TypeError, ValueError):
            raise Implicit_ODE_Exception('DQrhomax must be convertable to a float.')
            
        self.Integrator.DQrhomax = dqrhomax
        
    def _get_DQrhomax(self):
        """
        Specifies the selection parameters used in deciding switching between a simultaneous
        or separate approximation of the two terms in the sensitivity residual.
        
            Parameters::
            
                DQrhomax
                        - A postive float.
                        - Default 0.0
                        
            Returns::
            
                The current value of DQrhomax (float)
        
        See SUNDIALS IDAS documentation 5.2.6 'IDASetSensDQMethod' 
        """
        return self.Integrator.DQrhomax
    
    DQrhomax = property(_get_DQrhomax, _set_DQrhomax)
    
    def _set_usesens(self, usesens):
        """
        Specifies if the sensitivity calculations should be used or turned off.
        
            Parameters::
            
                usesens
                            - A boolean type.
                            - Default True
                            
            Returns::
            
                The current value of usesens (boolean)
        
        See SUNDIALS IDAS documentation 5.2.1 'IDASensToggleOff' 
        """
        self.Integrator.sensToggleOff = not bool(usesens)
        
    def _get_usesens(self):
        """
        Specifies if the sensitivity calculations should be used or turned off.
        
            Parameters::
            
                usesens
                            - A boolean type.
                            - Default True
                            
            Returns::
            
                The current value of usesens (boolean)
        
        See SUNDIALS IDAS documentation 5.2.1 'IDASensToggleOff' 
        """
        return not self.Integrator.sensToggleOff
        
    usesens = property(_get_usesens, _set_usesens)
    
    def _set_sensitivity_method(self, ism):
        """
        Specifies the sensitivity solution method. Can be either
        'IDA_SIMULTANEOUS' or 'IDA_STAGGERED'.
        
            Parameters::
            
                ism
                        - A string of either 'IDA_SIMULTANEOUS' or 'IDA_STAGGERED'
                        - Default 'IDA_STAGGERED'
                        
            Returns::
            
                The current value of sensmethod (string)
        
        See SUNDIALS IDAS documentation 5.2.1 'IDASensInit'
        """
        if not isinstance(ism, str):
            raise Implicit_ODE_Exception('sensmethod must be string.')
        
        if ism.upper() == 'IDA_SIMULTANEOUS':
            self.Integrator.ism = 1
        elif ism.upper() == 'IDA_STAGGERED':
            self.Integrator.ism = 2
        else:
            raise Implicit_ODE_Exception('sensmethod must be either "IDA_SIMULTANEOUS" or "IDA_STAGGERED".')
        
    def _get_sensitivity_method(self):
        """
        Specifies the sensitivity solution method. Can be either
        'IDA_SIMULTANEOUS' or 'IDA_STAGGERED'.
        
            Parameters::
            
                ism
                        - A string of either 'IDA_SIMULTANEOUS' or 'IDA_STAGGERED'
                        - Default 'IDA_STAGGERED'
                        
            Returns::
            
                The current value of sensmethod (string)
        
        See SUNDIALS IDAS documentation 5.2.1 'IDASensInit'
        """
        if self.Integrator.ism == 1:
            return 'IDA_SIMULTANEOUS'
        elif self.Integrator.ism == 2:
            return 'IDA_STAGGERED'
        else:
            raise Implicit_ODE_Exception('Unknown value of DQtype.')
    
    sensmethod = property(_get_sensitivity_method, _set_sensitivity_method)
    
    def _set_suppress_sens(self, suppress_sens):
        """
        Specifies whether sensitivity variables are included in the error test
        or not. True means that the variables are suppressed and not included 
        in the error test.
        
            Parameters::
            
                suppress_sens
                        - A boolean
                        - Default True
                        
            Returns::
                
                The current value of suppress_sens (boolean)
                
        See SUNDIALS IDAS documentation 5.2.6 'IDASetSensErrCon'.
        
        NOTE:: 
        
            That this method does the opposite of IDASetSensERRCon to have 
            the same meaning as suppress_alg.
        """
        self.Integrator.errconS = not bool(suppress_sens)
        
    def _get_suppress_sens(self):
        """
        Specifies whether sensitivity variables are included in the error test
        or not. True means that the variables are suppressed and not included 
        in the error test.
        
            Parameters::
            
                suppress_sens
                        - A boolean
                        - Default True
                        
            Returns::
                
                The current value of suppress_sens (boolean)
                
        See SUNDIALS IDAS documentation 5.2.6 'IDASetSensErrCon'.
        
        NOTE:: 
        
            That this method does the opposite of IDASetSensERRCon to have 
            the same meaning as suppress_alg.
        """
        return not self.Integrator.errconS
    
    suppress_sens = property(_get_suppress_sens, _set_suppress_sens)
    
    def _set_max_nonlin(self, maxsensiter):
        """
        Specifies the maximum number of nonlinear solver iterations for
        sensitivity variables per step. (>0)
        
            Parameters::
            
                maxsensiter 
                            - An integer
                            - Default 3
                            
            Returns::
            
                The current value of maxsensiter.
        
        See SUNDIALS IDAS documentation 5.2.6 'IDASetSensMaxNonlinIters'.
        """
        try:
            maxsensiter = int(maxsensiter)
            if maxsensiter < 1:
                raise Implicit_ODE_Exception('maxsensiter must be greater than zero.')
        except (TypeError, ValueError):
            raise Implicit_ODE_Exception('maxsensiter must be convertable to an integer.')
        
        self.Integrator.maxcorS = maxsensiter
        
    def _get_max_nonlin(self):
        """
        Specifies the maximum number of nonlinear solver iterations for
        sensitivity variables per step. (>0)
        
            Parameters::
            
                maxsensiter 
                            - An integer
                            - Default 3
                            
            Returns::
            
                The current value of maxsensiter.
        
        See SUNDIALS IDAS documentation 5.2.6 'IDASetSensMaxNonlinIters'.
        """
        return self.Integrator.maxcorS
    
    maxsensiter = property(_get_max_nonlin, _set_max_nonlin)
    
    def _set_pbar(self, pbar):
        """
        Specifies the order of magnitude for the parameters. This is useful if IDAS is
        to estimate tolerances for the sensitivity solution vectors.
        
            Parameters::
            
                pbar
                        - An array of positive floats equal to the number of parameters.
                        - Default array of ones.
                        
            Returns::
            
                The current value of pbar.
                
        See SUNDIALS IDA documentation 5.2.6 'IDASetSensParams'
        """
        if len(pbar) != self.Integrator.nbr_params:
            raise Implicit_ODE_Exception('pbar must be of equal length as the parameters.')
        
        self.Integrator.pbar = pbar
    
    def _get_pbar(self):
        """
        Specifies the order of magnitude for the parameters. This is useful if IDAS is
        to estimate tolerances for the sensitivity solution vectors.
        
            Parameters::
            
                pbar
                        - An array of positive floats equal to the number of parameters.
                        - Default array of ones.
                        
            Returns::
            
                The current value of pbar.
                
        See SUNDIALS IDA documentation 5.2.6 'IDASetSensParams'
        """
        return self.Integrator.pbar
    
    pbar = property(_get_pbar, _set_pbar)

class Radau5(Radau_Common,Implicit_ODE):
    """
    Radau IIA fifth-order three-stages with step-size control and continuous output.
    Based on the FORTRAN code by E.Hairer and G.Wanner, which can be found here: 
    http://www.unige.ch/~hairer/software.html
    
    Details about the implementation (FORTRAN) can be found in the book,
    
    Solving Ordinary Differential Equations II,
    Stiff and Differential-Algebraic Problems
    
    Authors: E. Hairer and G. Wanner
    Springer-Verlag, ISBN: 3-540-60452-9
    
    This code is aimed at providing a Python implementation of the original code.
    
    Input and bug reports are very welcome.
    
    HOMEPAGE:  http://www.jmodelica.org/assimulo
    FORUM:     http://www.jmodelica.org/forums/jmodelicaorg-users/assimulo
    """
    def __init__(self, problem, y0=None, yd0=None, t0=None):
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
                                    
        """
        Implicit_ODE.__init__(self, problem, y0, yd0, t0) #Calls the base class
        
        #Determine if we have a user supplied jacobian
        if hasattr(self._problem, 'jac'):
            self.usejac = True
        else:
            self.usejac = False
        
        #Internal values
        self._leny = len(y0) #Dimension of the problem
        self._2leny = 2*self._leny
        
        #Default values
        self.initstep = 0.01
        self.newt = 7 #Maximum number of newton iterations
        self.thet = 1.e-3 #Boundary for re-calculation of jac
        self.fnewt = 0 #Stopping critera for Newtons Method
        self.quot1 = 1.0 #Parameters for changing step-size (lower bound)
        self.quot2 = 1.2 #Parameters for changing step-size (upper bound)
        self.fac1 = 0.2 #Parameters for step-size selection (lower bound)
        self.fac2 = 8.0 #Parameters for step-size selection (upper bound)
        self.maxh = N.inf #Maximum step-size.
        self.safe = 0.9 #Safety factor
        self.index = [2]*len(y0)
        
        #Internal values
        self._curjac = False #Current jacobian?
        self._itfail = False #Iteration failed?
        self._needjac = True #Need to update the jacobian?
        self._needLU = True #Need new LU-factorisation?
        self._first = True #First step?
        self._rejected = True #Is the last step rejected?
        self._oldh = 0.0 #Old stepsize
        self._olderr = 1.0 #Old error
        self._eps = N.finfo('double').eps
        self._col_poly = N.zeros(self._2leny*3)
        
        # - Statistic values
        self._nsteps = 0 #Number of steps
        self._nfcn = 0 #Number of function evaluations
        self._njac = 0 #Number of jacobian evaluations
        self._njacfcn = 0 #Number of function evaluations when evaluating the jacobian
        self._nniter = 0 #Number of nonlinear iterations
        self._nniterfail = 0 #Number of nonlinear failures
        self._errfail = 0 #Number of step rejections
        self._nlu = 0 #Number of LU decompositions
        self._curiter = 0 #Number of current iterations
        
        # - Retrieve the Radau5 parameters
        self._load_parameters() #Set the Radau5 parameters
    
    def _set_index(self, index):
        """
        Sets the index of the variables in the problem which in turn
        determine the error estimations.
        
            Parameters::
            
                    index - A list of integers, indicating the index
                            (1,2,3) of the variable.
                            
                            Example:
                                Radau5.index = [2,1]
                            
        """
        if len(index) == self._2leny:
            self._index = N.array(index)
        elif len(index) == self._leny:
            self._index = N.array([1]*self._leny+index)
        else:
            raise Implicit_ODE_Exception('Wrong number of variables in the index vector.')
            
    def _get_index(self):
        """
        Sets the index of the variables in the problem which in turn
        determine the error estimations.
        
            Parameters::
            
                    index - A list of integers, indicating the index
                            (1,2,3) of the variable.
                            
                            Example:
                                Radau5.index = [2,1]
                            
        """
        return self._index
        
    index = property(_get_index,_set_index)
    
    def integrate(self, t, y, yd, tf,dt):
        """
        Integrates (t,y,yd) values until t > tf
        """
        self._oldh = self.initstep
        self.h = self.initstep
        self._hhfac = self.h
        
        self._fac_con = 1.0
        
        if self.fnewt == 0:
            self.fnewt = max(10.*self._eps/self.rtol,min(0.03,self.rtol**0.5))
            
        self._f0 = self._ode_f(t,N.append(y,yd))
        self._nfcn +=1
        self._tc = t
        self._yc = y
        self._ydc = yd 
        
        if dt > 0.0:
            dist_space = [(x+1)*dt for x in range(int((tf-t)/dt)+1)]
        
        for i in range(self.maxsteps):
            if t >= tf:
                break
            t, y, yd = self.step(t, y, yd)
            self._tc = t
            self._yc = y
            self._ydc = yd
            
            if dt > 0.0:
                while dist_space[0] <= t:
                    yy,yyd=self.interpolate(dist_space[0],y)
                    yield dist_space[0], yy, yyd
                    dist_space.pop(0)
            else:
                yield t,y,yd
            
            if self.h > N.abs(tf-t):
                self.h = N.abs(tf-t)
            self._hhfac = self.h

            self._first = False
        else:
            raise Implicit_ODE_Exception('Final time not reached within maximum number of steps')
    
    def _ode_f(self, t, y):
        return N.hstack((y[self._leny:],self.res_fcn(t,y[:self._leny],y[self._leny:])))
    
    def _radau_F(self, Z, t, y, yd):
        
        Z1 = Z[:self._2leny]
        Z2 = Z[self._2leny:2*self._2leny]
        Z3 = Z[2*self._2leny:3*self._2leny]
        
        q = N.append(y,yd)
        
        sol1 = self._ode_f(t+self.C[0]*self.h, q+Z1)
        sol2 = self._ode_f(t+self.C[1]*self.h, q+Z2)
        sol3 = self._ode_f(t+self.C[2]*self.h, q+Z3)
        
        self._nfcn += 3
        
        return N.hstack((N.hstack((sol1,sol2)),sol3))
    
    def step(self, t, y, yd):
        """
        This calculates the next step in the integration.
        """
        self._scaling = N.array(abs(N.append(y,yd))*self.rtol + self.atol) #The scaling used.
        
        while True: #Loop for integrating one step.
            
            self.newton(t,y,yd)
            self._err = self.estimate_error()
            
            if self._err > 1.0: #Step was rejected.
                self._rejected = True
                self._errfail += 1
                ho = self.h
                self.h = self.adjust_stepsize(self._err)
                
                self.print_verbos(['Rejecting step at ', t, 'with old stepsize', ho, 'and new ',
                                   self.h, '. Error: ', self._err],self.SCREAM)
                
                if self._curjac or self._curiter == 1:
                    self._needjac = False
                    self._needLU = True
                else:
                    self._needjac = True
                    self._needLU = True
            else:
                
                self.print_verbos(['Accepting step at ', t,'with stepsize ', self.h, '. Error: ', self._err],self.SCREAM)
                self._nsteps += 1
                
                tn = t+self.h #Preform the step
                yn = y+self._Z[2*self._2leny:3*self._2leny][:self._leny]
                ydn = yd+self._Z[2*self._2leny:3*self._2leny][self._leny:]
                self._f0 = self._ode_f(t,N.append(yn,ydn))
                self._nfcn += 1
                
                self._oldoldh = self._oldh #Store the old(old) step-size for use in the test below.
                self._oldh = self.h #Store the old step-size
                self._oldt = t #Store the old time-point
                self._newt = tn #Store the new time-point
                
                #Adjust the new step-size
                ht = self.adjust_stepsize(self._err, predict=True)
                self.h = min(self.h,ht) if self._rejected else ht
                
                self._rejected = False
                self._curjac = False
                
                if self._oldoldh == self.h and (self._theta <= self.thet or self._curiter==1):
                    self._needjac = False
                    self._needLU = False
                else:
                    if self._theta <= self.thet or self._curiter == 1:
                        self._needjac = False
                        self._needLU = True
                    else:
                        self._needjac = True
                        self._needLU = True
                        
                self._olderr = max(self._err,1.e-2) #Store the old error
                break
                
        self._col_poly = self._collocation_pol(self._Z, self._col_poly, self._2leny) #Calculate the new collocation polynomial
        
        return tn, yn, ydn #Return the step
    
    def newton(self,t,y,yd):
        """
        The newton iteration. 
        """
        
        for k in xrange(20):
            
            self._curiter = 0 #Reset the iteration
            self._fac_con = max(self._fac_con, self._eps)**0.8;
            self._theta = abs(self.thet);
            
            if self._needjac:
                self._jac = self.jacobian(t,y,yd)
            
            if self._needLU:
                self._nlu += 1
                self._a = self._alpha/self.h
                self._b = self._beta/self.h
                self._g = self._gamma/self.h
                self._B = self._g*self.M - self._jac
                
                self._P1,self._L1,self._U1 = S.linalg.lu(self._B) #LU decomposition
                self._P2,self._L2,self._U2 = S.linalg.lu(self._a*self.M-self._jac)
                self._P3,self._L3,self._U3 = S.linalg.lu(self._b*self.M-self._jac)
                
                self._needLU = False
                
                if min(abs(N.diag(self._U1)))<self._eps:
                    raise Implicit_ODE_Exception('Error, gM-J is singular at ',self._tc)
                    
            Z, W = self.calc_start_values()

            for i in xrange(self.newt):
                self._curiter += 1 #The current iteration
                self._nniter += 1 #Adding one iteration

                #Solve the system
                Z = N.dot(self.T2,self._radau_F(Z.real,t,y,yd))

                Z[:self._2leny]               =Z[:self._2leny]               -self._g*N.dot(self.M,W[:self._2leny])
                Z[self._2leny:2*self._2leny]  =Z[self._2leny:2*self._2leny]  -self._a*N.dot(self.M,W[self._2leny:2*self._2leny])   #+self._b*N.dot(self.I,W[2*self._leny:3*self._leny])
                Z[2*self._2leny:3*self._2leny]=Z[2*self._2leny:3*self._2leny]-self._b*N.dot(self.M,W[2*self._2leny:3*self._2leny]) #-self._a*N.dot(self.I,W[2*self._leny:3*self._leny])
                
                Z[:self._2leny]               =N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,Z[:self._2leny])))
                Z[self._2leny:2*self._2leny]  =N.linalg.solve(self._U2,N.linalg.solve(self._L2,N.linalg.solve(self._P2,Z[self._2leny:2*self._2leny])))
                Z[2*self._2leny:3*self._2leny]=N.linalg.solve(self._U3,N.linalg.solve(self._L3,N.linalg.solve(self._P3,Z[2*self._2leny:3*self._2leny])))
                #----
                
                self._scaling = self._scaling/self.h**(self.index-1)#hfac
                
                newnrm = N.linalg.norm(Z.reshape(-1,self._2leny)/self._scaling,'fro')/N.sqrt(3.*self._2leny)
                
                if i > 0:
                    thq = newnrm/oldnrm
                    if i == 1:
                        self._theta = thq
                    else:
                        self._theta = N.sqrt(thq*thqold)
                    thqold = thq
                    
                    if self._theta < 0.99: #Convergence
                        self._fac_con = self._theta/(1.-self._theta)
                        dyth = self._fac_con*newnrm*self._theta**(self.newt-(i+1)-1)/self.fnewt
                        
                        if dyth >= 1.0: #Too slow convergence
                            qnewt = max(1.e-4,min(20.,dyth))
                            self._hhfac = 0.8*qnewt**(-1.0/(4.0+self.newt-(i+1)-1))
                            self.h = self._hhfac*self.h
                            self._itfail = True
                            self._rejected = True
                            break
                    else: #Not convergence, abort
                        self._itfail = True
                        break
                
                oldnrm = max(newnrm,self._eps) #Store oldnorm
                W = W+Z #Perform the iteration
                
                Z = N.dot(self.T3,W) #Calculate the new Z values
                
                if self._fac_con*newnrm <= self.fnewt: #Convergence?
                    self._itfail = False;
                    break
                
            else: #Iteration failed
                self._itfail = True
                
            if not self._itfail: #Newton iteration converged
                self._Z = Z.real
                break
            else: #Iteration failed
                self.print_verbos(['Iteration failed at time %e with step-size %e'%(t,self.h)],self.SCREAM)
                self._nniterfail += 1
                self._rejected = True #The step is rejected
                
                if self._theta >= 0.99:
                    self._hhfac = 0.5
                    self.h = self.h*self._hhfac
                if self._curjac:
                    self._needjac = False
                    self._needLU = True
                else:
                    self._needjac = True
                    self._needLU = True
        else:
            raise Implicit_ODE_Exception('Newton iteration failed at time %e with step-size %e'%(t,self.h))
    
    def estimate_error(self):
        
        temp = 1./self.h*(self.E[0]*self._Z[:self._2leny]+self.E[1]*self._Z[self._2leny:2*self._2leny]+self.E[2]*self._Z[2*self._2leny:3*self._2leny])
        temp = N.dot(self.M,temp)
        
        self._scaling = self._scaling/self.h**(self.index-1)#hfac
        
        scal = self._scaling#/self.h
        err_v = N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,self._f0+temp)))
        err = N.linalg.norm(err_v/scal)
        err = max(err/N.sqrt(self._2leny),1.e-10)

        if (self._rejected or self._first) and err >= 1.: #If the step was rejected, use the more expensive error estimation
            self._nfcn += 1
            err_v = self._ode_f(self._tc,N.append(self._yc,self._ydc)+err_v)
            err_v = N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,err_v+temp)))
            err = N.linalg.norm(err_v/scal)
            err = max(err/N.sqrt(self._2leny),1.e-10)
            
        return err
    
    def interpolate(self, t, k):
        """
        Calculates the continuous output from Radau5.
        """
        leny = self._2leny
        s = (t-self._newt)/self._oldh
        Z = self._col_poly
        
        diff = s*(Z[:leny]+(s-self.C[1,0]+1.)*(Z[leny:2*leny]+(s-self.C[0,0]+1.)*Z[2*leny:3*leny]))
        
        yout  = self._yc + diff[:self._leny]
        ydout = self._ydc+ diff[self._leny:]

        return yout, ydout
    
    def jacobian(self, t, y, yd):
        """
        Calculates the Jacobian, either by an approximation or by the user
        defined (jac specified in the problem class).
        """
        self._curjac = True #The jacobian is up to date
        self._needLU = True #A new LU-decomposition is needed
        self._needjac = False #A new jacobian is not needed
        
        q = N.append(y,yd)
        
        if self.usejac: #Retrieve the user-defined jacobian
            cjac = self._problem.jac(t,y,yd)
        else:           #Calculate a numeric jacobian
            delt = N.array([(self._eps*max(abs(yi),1.e-5))**0.5 for yi in q])*N.identity(self._2leny) #Calculate a disturbance
            Fdelt = N.array([self._ode_f(t,q+e) for e in delt]) #Add the disturbance (row by row) 
            grad = ((Fdelt-self._ode_f(t,q)).T/delt.diagonal()).T
            cjac = N.array(grad).T
            self._njacfcn += 1+self._2leny #Add the number of function evaluations

        self._njac += 1 #add the number of jacobian evaluation
        return cjac
    
    def adjust_stepsize(self, err, predict=False):
        
        fac = min(self.safe, self.safe*(2.*self.newt+1.)/(2.*self.newt+self._curiter))
        quot = max(1./self.fac2,min(1./self.fac1,(err**0.25)/fac))        
        hnormal = self.h/quot
        
        if predict:
            if not self._first:
                facgus = (self._hacc/self.h)*(err**2/self._olderr)**0.25/self.safe
                facgus = max(1./self.fac2,min(1./self.fac1,facgus))
                quot = max(quot,facgus)
                h = self.h/quot
            else:
                h = hnormal
            self._hacc = self.h
        else:
            h = hnormal
        
        qt = h/self.h
        
        if (qt >= self.quot1) and (qt <= self.quot2):
            h = self.h
        
        if h > self.maxh:
            h = self.maxh
        
        if self._first and err>=1.0:
            self._hhfac = 0.1
            h = self.h*self._hhfac
        else:
            self._hhfac = h/self.h
        
        if h < self._eps:
            raise Implicit_ODE_Exception('Step-size to small at %e with h = %e'%(self._tc,self.h))
    
        return h
    
    def _collocation_pol(self, Z, col_poly, leny):

        col_poly[2*leny:3*leny] = Z[:leny] / self.C[0,0]
        col_poly[leny:2*leny]   = ( Z[:leny] - Z[leny:2*leny] ) / (self.C[0,0]-self.C[1,0])
        col_poly[:leny]         = ( Z[leny:2*leny] -Z[2*leny:3*leny] ) / (self.C[1,0]-1.)
        col_poly[2*leny:3*leny] = ( col_poly[leny:2*leny] - col_poly[2*leny:3*leny] ) / self.C[1,0]
        col_poly[leny:2*leny]   = ( col_poly[leny:2*leny] - col_poly[:leny] ) / (self.C[0,0]-1.)
        col_poly[2*leny:3*leny] =   col_poly[leny:2*leny]-col_poly[2*leny:3*leny]
        
        return col_poly
    
    def calc_start_values(self):
        """
        Calculate newton starting values.
        """
        if self._first:
            Z = N.zeros(self._2leny*3)
            W = N.zeros(self._2leny*3)
        else:
            Z = self._Z
            cq = self.C*self.h/self._oldh#self._oldoldh#self._oldh
            newtval = self._col_poly
            leny = self._2leny
            
            Z[:leny]        = cq[0,0]*(newtval[:leny]+(cq[0,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[0,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            Z[leny:2*leny]  = cq[1,0]*(newtval[:leny]+(cq[1,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[1,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            Z[2*leny:3*leny]= cq[2,0]*(newtval[:leny]+(cq[2,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[2,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            
            W = N.dot(self.T2,Z)
            
        return Z, W
    
    def _load_parameters(self):
        
        #Parameters
        A = N.zeros([3,3])
        A[0,0] = (88.-7.*N.sqrt(6.))/360.0
        A[0,1] = (296.-169.*N.sqrt(6.))/1800.0
        A[0,2] = (-2.0+3.0*N.sqrt(6.))/225.0
        A[1,0] = (296.0+169.0*N.sqrt(6.))/1800.0
        A[1,1] = (88.+7.*N.sqrt(6.))/360.0
        A[1,2] = (-2.-3.*N.sqrt(6.))/225.0
        A[2,0] = (16.0-N.sqrt(6.))/36.0
        A[2,1] = (16.0+N.sqrt(6.))/36.0
        A[2,2] = (1.0/9.0)
        
        C = N.zeros([3,1])
        C[0,0]=(4.0-N.sqrt(6.0))/10.0
        C[1,0]=(4.0+N.sqrt(6.0))/10.0
        C[2,0]=1.0
        
        B = N.zeros([1,3])
        B[0,0]=(16.0-N.sqrt(6.0))/36.0
        B[0,1]=(16.0+N.sqrt(6.0))/36.0
        B[0,2]=1.0/9.0
        
        E = N.zeros(3)
        E[0] = -13.0-7.*N.sqrt(6.)
        E[1] = -13.0+7.0*N.sqrt(6.)
        E[2] = -1.0
        E = 1.0/3.0*E
        
        M = N.array([[1.,0.],[0.,0.]])
        
        Ainv = N.linalg.inv(A)
        [eig, T] = N.linalg.eig(Ainv)
        eig = N.array([eig[2],eig[0],eig[1]])
        J = N.diag(eig)

        self._alpha = eig[1]
        self._beta  = eig[2]
        self._gamma = eig[0].real
        
        temp0 = T[:,0].copy()
        temp1 = T[:,1].copy()
        temp2 = T[:,2].copy()
        T[:,0] = temp2
        T[:,1] = temp0
        T[:,2] = temp1
        Tinv = N.linalg.inv(T)
        
        I = N.eye(self._2leny)
        M = N.kron(M,N.eye(self._leny))
        I3 = N.eye(3)
        T1 = N.kron(J,M)
        T2 = N.kron(Tinv,I)
        T3 = N.kron(T,I)
        
        self.A = A
        self.B = B
        self.C = C
        self.I = I
        self.E = E
        self.M = M
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.I3 = I3
        self.EIG = eig
