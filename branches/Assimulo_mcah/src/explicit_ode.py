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
from problem import Explicit_Problem
from sundials import Sundials, Sundials_Exception
from assimulo.lib.radau_core import Radau_Common

class Explicit_ODE_Exception(Exception):
    """ An integrator exception. """
    pass

class Explicit_ODE(ODE):
    """
    Baseclass for our explicit ODE integrators.
    """
    
    def __init__(self, problem, y0=None, t0=None, switches0=None):
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
                                    
        """
        ODE.__init__(self) #Sets general attributes
        
        if problem == None:
            raise Explicit_ODE_Exception('The problem needs to be a subclass of a Explicit_Problem')
        
        if isinstance(problem, Explicit_Problem):
            self._problem = problem
            self.f = problem.f
            self.problem_name = problem.problem_name
        else:
            raise Explicit_ODE_Exception('The problem needs to be a subclass of a Explicit_Problem.')
        
        if y0 == None:
            if hasattr(problem, 'y0'):
                y0 = problem.y0
            else:
                raise Explicit_ODE_Exception('y0 must be specified. Either in the problem or in the initialization')

        try:
            if isinstance(y0, int) or isinstance(y0, float):
                y0 = [y0]
            self._problem.y0 = y0[:]
            self.y_cur = N.array(y0, dtype=float)
        except ValueError:
            raise Explicit_ODE_Exception('Initial values must be a scalar/list/array of type int or float.')
        
        if len(self.y_cur)==0:
            raise Explicit_ODE_Exception('Initial values must be provided.')
        
        if hasattr(self._problem, 'switches0') and switches0 == None:
            switches0 = self._problem.switches0
        
        if isinstance(switches0, list):
            for x in switches0:
                if not isinstance(x, bool):
                    raise Explicit_ODE_Exception('Switches must be a list of booleans.')
        elif switches0 is not None:
            raise Explicit_ODE_Exception('Switches must be a list of booleans.')
        
        self._problem.switches0 = switches0
        self.switches = switches0
        
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
            raise Explicit_ODE_Exception('Initial time must be an integer or float.')
        
        if hasattr(self._problem, 'time_events'):
            self._time_function  = True
        
        self.t_cur = N.array(self.t_cur)
        self.y_cur = N.array(self.y_cur)
        
        self.y = []
        self.t = []
            
    def reset(self):
        """
        
        Resets the problem. If the problem is defined with a reset method, its called
        and then the method re_init. The re_init method is called with the initial
        values set in the problem, problem.t0 and problem.y0.
        
        """
        self._problem.reset()
        self.re_init(self._problem.t0, self._problem.y0)
    
    def _integrator(self, t, y, tf, nt):
        pass 

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
        except ValueError:
            raise Explicit_ODE_Exception('Final time must be an integer or float.')
            
        if self.t_cur > tfinal:
            raise Explicit_ODE_Exception('Final time must be greater than start time.')
        
        if not isinstance(ncp, int):
            raise Explicit_ODE_Exception('Number of communication points must be an integer')
        if ncp < 0:
            ncp = 0
            if self.verbosity > self.QUIET:
                print 'Number of communication points must be a positive integer, setting' \
                      ' nt = 0.'

        ncp_ori = ncp
        tfinal_ori = tfinal
        time_start = time.clock()
        
        t0 = self.t_cur
        y0 = self.y_cur
        last_logg = t0
        
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
        
        self._problem.handle_result(self,t0,y0) #Logg the first point
        self._flag_init = True #Reinitiate the solver
        
        while self.t_cur < tfinal_ori:

            #Time event function is specified.
            if self._time_function:
                tevent = self._problem.time_events(self.t_cur, self.y_cur, self.switches)
                if tevent == None:
                    tfinal = tfinal_ori
                else:
                    tfinal = tevent if tevent < tfinal_ori else tfinal_ori

            solution = list(self._integrator(self.t_cur, self.y_cur, tfinal,dt))
            tt, yy = solution[-1]
            self.t_cur = tt.copy()
            self.y_cur = yy.copy()
            
            if mode == 'SPECIAL':
                while dist_space[0] <= self.t_cur:
                    self._problem.handle_result(self, dist_space[0], self.interpolate(dist_space[0],0))
                    last_logg = dist_space[0].copy()
                    dist_space.pop(0)
            else:
                for q in solution:
                    self._problem.handle_result(self,q[0],q[1])
                last_logg = self.t_cur
            
            if self._completed_step: #If the option completed is set.
                self._flag_init = self._problem.completed_step(self)
            else:
                self._flag_init = False
            
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
                
                if self.verbosity > self.NORMAL:
                    print 'A discontinuity occured at t = %e.'%self.t_cur
                if self.verbosity >= self.LOUD:
                    print 'Current switches: ', self.switches
                    print 'Event info: ', event_info
                    
                if self.verbosity >= self.SCREAM:
                    self.print_statistics() #Prints statistics
                    
                if self.verbosity > self.NORMAL:
                    print 'Calling problem specified event handling...'
                
                self._problem.handle_event(self, event_info) #self corresponds to the solver
                self._flag_init = True
            
            if self._flag_init and last_logg == self.t_cur: #Logg after the event handling if there was a communication point there.
                self._problem.handle_result(self, self.t_cur, self.y_cur)
        
        #Simulation complete, call finalize
        self._problem.finalize(self)
        
        time_stop = time.clock()
        
        if self.verbosity >= self.NORMAL:
            self.print_statistics()
            print 'Elapsed simulation time:', time_stop-time_start, 'seconds.'
        
        #return [self.t, self.y]
    
    def re_init(self,t0, y0):
        """
        Reinitiates the solver.
        
            Parameters::
                
                t0  
                    - The initial time.
                y0  
                    - The initial values for the states
                
        See information in the __init__ method.
        """
        if len(self.y_cur) != len(y0):
            raise Explicit_ODE_Exception('y0 must be of the same length as the original problem.')
        Explicit_ODE.__init__(self, self._problem,y0,t0)
    
    def plot(self, mask=None, **kwargs):
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
                
                **kwargs
                        - See http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
                          for information about the available options for **kwargs.
        """
        P.xlabel('time')
        P.ylabel('state')
        P.title(self.problem_name)
        
        if not mask:
            P.plot(self.t, self.y, **kwargs)
        else:
            if not isinstance(mask, list):
                raise Explicit_ODE_Exception('Mask must be a list of integers')
            if not len(mask)==len(self.y[-1]):
                raise Explicit_ODE_Exception('Mask must be a list of integers of equal length as '\
                                             'the number of variables.')
            for i in range(len(mask)):
                if mask[i]:
                    P.plot(self.t, N.array(self.y)[:,i],**kwargs)
        
        
        P.show()
            
            
    
class Explicit_Euler(Explicit_ODE):
    """
    Explicit Euler.
    """
    def _integrator(self, t, y, tf, dt):
        """
        _integrates (t,y) values until t > tf
        """
        if dt <= 0.0:
            raise Explicit_ODE_Exception('Explicit Euler is a fixed step-size method. Provide' \
                                         ' the number of communication points.')
        self.h = dt
        self._hlength = dt
        
        for i in range(self.maxsteps):
            if t >= tf:
                break
            t, y = self.step(t, y)
            yield t,y
            self.h=min(self.h,N.abs(tf-t))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
    
    def step(self, t, y):
        """
        This calculates the next step in the integration.
        """
        f = self.f
        h = self.h
        return t + h, y + h*f(t, y)
        
    def print_statistics(self):
        """
        Should print the statistics.
        """
        print 'Final Run Statistics: %s \n' % self.problem_name
        print 'Step-length          : %s'%(self._hlength)
        
        print '\nSolver options:\n'
        print ' Solver            : Explicit_Euler'
        print ' Solver type       : Fixed step'
        print ''
    
    
class RungeKutta34(Explicit_ODE):
    """
    Adaptive Runge-Kutta of order four.
    Obs. Step rejection not implemented.
    """
    def __init__(self, problem, y0=None, t0=None, switches0=None):
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
                                    
        """
        Explicit_ODE.__init__(self, problem, y0, t0, switches0) #Calls the base class
        
        #Default values
        self.initstep = 0.01
        self.atol = 1.e-6
        self.rtol = 1.e-6
        
        #Internal values
        # - Statistic values
        self._nsteps = 0 #Number of steps
        self._nfcn = 0 #Number of function evaluations
        
    def _set_initial_step(self, initstep):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.01'.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        try:
            initstep = float(initstep)
        except (ValueError, TypeError):
            raise Explicit_ODE_Exception('The initial step must be an integer or float.')
        
        self.__initstep = initstep
        
    def _get_initial_step(self):
        """
        This determines the initial step-size to be used in the integration.
        
            Parameters::
            
                initstep    
                            - Default '0.01'.
                            
                            - Should be float.
                            
                                Example:
                                    initstep = 0.01
        """
        return self.__initstep
        
    initstep = property(_get_initial_step,_set_initial_step)
    
    def _set_atol(self,atol):
        """
        Sets the absolute tolerance to be used in the integration.
        
            Parameters::
            
                atol    
                            - Default 1.0e-6.
                            
                            - Should be float or an array/list of len(y)
                            
                                Example:
                                    atol=1.e5
                                    atol=[1.e-5,1.e-4]
        """
        try:
            atol_arr = N.array(atol, dtype=float)
            if (atol_arr <= 0.0).any():
                raise Explicit_ODE_Exception('Absolute tolerance must be a positive float or a float vector.')
        except (ValueError,TypeError):
            raise Explicit_ODE_Exception('Absolute tolerance must be a positive float or a float vector.')
        if atol_arr.size == 1:
            self.__atol = float(atol)
        elif atol_arr.size == len(self.y_cur):
            self.__atol = [float(x) for x in atol]
        else:
            raise Explicit_ODE_Exception('Absolute tolerance must be a float vector of same dimension as the problem or a scalar.')

    def _get_atol(self):
        """
        Sets the absolute tolerance to be used in the integration.
        
            Parameters::
            
                atol    
                            - Default 1.0e-6.
                            
                            - Should be float or an array/list of len(y)
                            
                                Example:
                                    atol=1.e5
                                    atol=[1.e-5,1.e-4]
        """
        return self.__atol
    
    atol = property(_get_atol,_set_atol)
    
    def _set_rtol(self, rtol):
        """
        The relative tolerance to be used in the integration.
        
            Parameters::
            
                rtol    
                            - Default 1.0e-6
                            
                            - Should be a float.
        """
        try:
            rtol = float(rtol)
        except (TypeError,ValueError):
            raise Explicit_ODE_Exception('Relative tolerance must be a float.')
        if rtol <= 0.0:
            raise Explicit_ODE_Exception('Relative tolerance must be a positive (scalar) float.')
        self.__rtol = rtol
            
    def _get_rtol(self):
        """
        The relative tolerance to be used in the integration.
        
            Parameters::
            
                rtol    
                            - Default 1.0e-6
                            
                            - Should be a float.
        """
        return self.__rtol
    
    rtol = property(_get_rtol, _set_rtol)
    
    def _integrator(self, t, y, tf, dt):
        """
        Integrates (t,y) values until t > tf
        """
        self.h = self.initstep
        
        for i in range(self.maxsteps):
            if t >= tf:
                break
            t, y = self.step(t, y)
            self._nsteps += 1
            yield t,y
            self.adjust_stepsize()
            self.h=min(self.h,N.abs(tf-t))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
    
    def adjust_stepsize(self):
        """
        Adjusts the stepsize.
        """
        fac=min((1./self.error)**(1.0/4.0),2.)
        self.h *= fac
    
    def step(self, t, y):
        """
        This calculates the next step in the integration.
        """
        self._nfcn += 5
        self._scaling = N.array(abs(y)*self.rtol + self.atol) # to normalize the error 
        f = self.f
        h = self.h
        Y1 = f(t, y)
        Y2 = f(t + h/2., y + h*Y1/2.)
        Y3 = f(t + h/2, y + h*Y2/2)
        Z3 = f(t + h, y - h*Y1 + 2*h*Y2)
        Y4 = f(t + h, y + h*Y3)
        self.error = N.linalg.norm(h/6*(2*Y2 + Z3 - 2*Y3 - Y4)/self._scaling) #normalized 
        return t+h, y + h/6*(Y1 + 2*Y2 + 2*Y3 + Y4)
    
    def print_statistics(self):
        """
        Should print the statistics.
        """
        print 'Final Run Statistics: %s \n' % self.problem_name
        print 'Number of Steps                 :', self._nsteps
        print 'Number of Function Evaluations  :', self._nfcn
        
        print '\nSolver options:\n'
        print ' Solver               : RungeKutta34'
        print ' Solver type          : Adaptive'
        print ' Relative tolerance   : ', self.rtol
        print ' Absolute tolerance   : ', self.atol
        print ''
    
    
class RungeKutta4(Explicit_ODE):
    """
    Runge-Kutta of order 4.
    """
    def _integrator(self, t, y, tf, dt):
        """
        Integrates (t,y) values until t > tf
        """
        if dt <= 0.0:
            raise Explicit_ODE_Exception('RungeKutta4 is a fixed step-size method. Provide' \
                                         ' the number of communication points.')
        
        self.h = dt
        self._hlength = dt

        for i in range(self.maxsteps):
            if t >= tf:
                break
            t, y = self.step(t, y)
            yield t,y
            self.h=min(self.h,N.abs(tf-t))
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
    
    def step(self, t, y):
        """
        This calculates the next step in the integration.
        """
        f = self.f
        h = self.h
        Y1 = f(t, y)
        Y2 = f(t + h/2., y + h*Y1/2.)
        Y3 = f(t + h/2., y + h*Y2/2.)
        Y4 = f(t + h, y + h*Y3)
        return t+h, y + h/6.*(Y1 + 2.*Y2 + 2.*Y3 + Y4)
        
    def print_statistics(self):
        """
        Should print the statistics.
        """
        print 'Final Run Statistics: %s \n' % self.problem_name
        print 'Step-length        :', self._hlength
        
        print '\nSolver options:\n'
        print ' Solver            : RungeKutta4'
        print ' Solver type       : Fixed step'
        print ''
    
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
            self.pbar = N.abs(self._problem.p0)
        
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
        
        

class Radau5(Radau_Common,Explicit_ODE):
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
    
    def __init__(self, problem, y0=None, t0=None):
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
                                    
        """
        Explicit_ODE.__init__(self, problem, y0, t0) #Calls the base class
        
        #Determine if we have a user supplied jacobian
        if hasattr(self._problem, 'jac'):
            self.usejac = True
        else:
            self.usejac = False

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
        self.atol = 1.0e-6 #Absolute tolerance
        self.rtol = 1.0e-6 #Relative tolerance
        
        #Internal values
        self._curjac = False #Current jacobian?
        self._itfail = False #Iteration failed?
        self._needjac = True #Need to update the jacobian?
        self._needLU = True #Need new LU-factorisation?
        self._first = True #First step?
        self._rejected = True #Is the last step rejected?
        self._leny = len(self.y_cur) #Dimension of the problem
        self._oldh = 0.0 #Old stepsize
        self._olderr = 1.0 #Old error
        self._eps = N.finfo('double').eps
        self._col_poly = N.zeros(self._leny*3)
        self._type = '(explicit)'
        
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
        
    def _integrator(self, t, y, tf, dt):
        """
        Integrates (t,y) values until t > tf
        """
        if self._flag_reset_statistics:
            self._nsteps = 0 #Number of steps
            self._nfcn = 0 #Number of function evaluations
            self._njac = 0 #Number of jacobian evaluations
            self._njacfcn = 0 #Number of function evaluations when evaluating the jacobian
            self._nniter = 0 #Number of nonlinear iterations
            self._nniterfail = 0 #Number of nonlinear failures
            self._errfail = 0 #Number of step rejections
            self._nlu = 0 #Number of LU decompositions
            self._curiter = 0 #Number of current iterations
            self._flag_reset_statistics = False
        
        self._oldh = self.initstep
        self.h = self.initstep
        
        self._fac_con = 1.0
        
        if self.fnewt == 0:
            self.fnewt = max(10.*self._eps/self.rtol,min(0.03,self.rtol**0.5))
            
        self._f0 = self.f(t,y)
        self._nfcn +=1
        self._tc = t
        self._yc = y
        
        if dt > 0.0:
            ncp = (tf-t)/dt
            dist_space = [(x+1)*(tf-t)/ncp for x in range(int(ncp)+1)]
        
        for i in range(self.maxsteps):
            if t >= tf:
                break
            t, y = self.step(t, y)
            self._tc = t
            self._yc = y
            
            if dt > 0.0:
                while dist_space[0] <= t:
                    yy=self.interpolate(dist_space[0],y)
                    yield dist_space[0], yy
                    dist_space.pop(0)
            else:
                yield t,y
            if self.h > N.abs(tf-t):
                self.h = N.abs(tf-t)

            self._first = False
        else:
            raise Explicit_ODE_Exception('Final time not reached within maximum number of steps')
        
    def step(self, t, y):
        """
        This calculates the next step in the integration.
        """
        self._scaling = N.array(abs(y)*self.rtol + self.atol) #The scaling used.
        
        while True: #Loop for integrating one step.
            
            self.newton(t,y)
            self._err = self.estimate_error()
            
            if self._err > 1.0: #Step was rejected.
                self._rejected = True
                self._errfail += 1
                ho = self.h
                self.h = self.adjust_stepsize(self._err)
                
                if self.verbosity >= self.SCREAM:
                    print 'Rejecting step at ', t, 'with old stepsize', ho, 'and new ', self.h
                
                if self._curjac or self._curiter == 1:
                    self._needjac = False
                    self._needLU = True
                else:
                    self._needjac = True
                    self._needLU = True
            else:
                if self.verbosity >= self.SCREAM:
                    print 'Accepting step at ', t,'with stepsize ', self.h
                self._nsteps += 1
                
                tn = t+self.h #Preform the step
                yn = y+self._Z[2*self._leny:3*self._leny]
                self._f0 = self.f(tn,yn)
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
                
                if self._oldoldh == self.h and (self._theta <= self.thet):# or self._curiter==1):
                    self._needjac = False
                    self._needLU = False
                else:
                    if self._theta <= self.thet: #or self._curiter == 1:
                        self._needjac = False
                        self._needLU = True
                    else:
                        self._needjac = True
                        self._needLU = True
                if self.thet < 0:
                    self._needjac = True
                    self._needLU = True
                        
                self._olderr = max(self._err,1.e-2) #Store the old error
                break
                
        self._col_poly = self._collocation_pol(self._Z, self._col_poly, self._leny) #Calculate the new collocation polynomial
        
        return tn, yn #Return the step
    
    def _collocation_pol(self, Z, col_poly, leny):
        
        col_poly[2*leny:3*leny] = Z[:leny] / self.C[0,0]
        col_poly[leny:2*leny]   = ( Z[:leny] - Z[leny:2*leny] ) / (self.C[0,0]-self.C[1,0])
        col_poly[:leny]         = ( Z[leny:2*leny] -Z[2*leny:3*leny] ) / (self.C[1,0]-1.)
        col_poly[2*leny:3*leny] = ( col_poly[leny:2*leny] - col_poly[2*leny:3*leny] ) / self.C[1,0]
        col_poly[leny:2*leny]   = ( col_poly[leny:2*leny] - col_poly[:leny] ) / (self.C[0,0]-1.)
        col_poly[2*leny:3*leny] =   col_poly[leny:2*leny]-col_poly[2*leny:3*leny]
        
        return col_poly
    
    def _radau_F(self, Z, t, y):
        
        Z1 = Z[:self._leny]
        Z2 = Z[self._leny:2*self._leny]
        Z3 = Z[2*self._leny:3*self._leny]

        sol1 = self.f(t+self.C[0]*self.h, y+Z1)
        sol2 = self.f(t+self.C[1]*self.h, y+Z2)
        sol3 = self.f(t+self.C[2]*self.h, y+Z3)
        
        self._nfcn += 3
        
        return N.hstack((N.hstack((sol1,sol2)),sol3))
    
    def calc_start_values(self):
        """
        Calculate newton starting values.
        """
        if self._first:
            Z = N.zeros(self._leny*3)
            W = N.zeros(self._leny*3)
        else:
            Z = self._Z
            cq = self.C*self.h/self._oldh#self._oldoldh#self._oldh
            newtval = self._col_poly
            leny = self._leny
            
            Z[:leny]        = cq[0,0]*(newtval[:leny]+(cq[0,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[0,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            Z[leny:2*leny]  = cq[1,0]*(newtval[:leny]+(cq[1,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[1,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            Z[2*leny:3*leny]= cq[2,0]*(newtval[:leny]+(cq[2,0]-self.C[1,0]+1.)*(newtval[leny:2*leny]+(cq[2,0]-self.C[0,0]+1.)*newtval[2*leny:3*leny]))
            
            W = N.dot(self.T2,Z)
            
        return Z, W
    
    def newton(self,t,y):
        """
        The newton iteration. 
        """
        
        for k in xrange(20):
            
            self._curiter = 0 #Reset the iteration
            self._fac_con = max(self._fac_con, self._eps)**0.8;
            self._theta = abs(self.thet);
            
            if self._needjac:
                self._jac = self.jacobian(t,y)
            
            if self._needLU:
                self._nlu += 1
                self._a = self._alpha/self.h
                self._b = self._beta/self.h
                self._g = self._gamma/self.h
                self._B = self._g*self.I - self._jac
                
                self._P1,self._L1,self._U1 = S.linalg.lu(self._B) #LU decomposition
                self._P2,self._L2,self._U2 = S.linalg.lu(self._a*self.I-self._jac)
                self._P3,self._L3,self._U3 = S.linalg.lu(self._b*self.I-self._jac)
                
                self._needLU = False
                
                if min(abs(N.diag(self._U1)))<self._eps:
                    raise Explicit_ODE_Exception('Error, gI-J is singular.')
                    
            Z, W = self.calc_start_values()
        
            for i in xrange(self.newt):
                self._curiter += 1 #The current iteration
                self._nniter += 1 #Adding one iteration
                
                #Solve the system
                Z = N.dot(self.T2,self._radau_F(Z.real,t,y))

                Z[:self._leny]              =Z[:self._leny]              -self._g*N.dot(self.I,W[:self._leny])
                Z[self._leny:2*self._leny]  =Z[self._leny:2*self._leny]  -self._a*N.dot(self.I,W[self._leny:2*self._leny])   #+self._b*N.dot(self.I,W[2*self._leny:3*self._leny])
                Z[2*self._leny:3*self._leny]=Z[2*self._leny:3*self._leny]-self._b*N.dot(self.I,W[2*self._leny:3*self._leny]) #-self._a*N.dot(self.I,W[2*self._leny:3*self._leny])
                
                Z[:self._leny]              =N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,Z[:self._leny])))
                Z[self._leny:2*self._leny]  =N.linalg.solve(self._U2,N.linalg.solve(self._L2,N.linalg.solve(self._P2,Z[self._leny:2*self._leny])))
                Z[2*self._leny:3*self._leny]=N.linalg.solve(self._U3,N.linalg.solve(self._L3,N.linalg.solve(self._P3,Z[2*self._leny:3*self._leny])))
                #----
                newnrm = N.linalg.norm(Z.reshape(-1,self._leny)/self._scaling,'fro')/N.sqrt(3.*self._leny)
                      
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
                            self.h = 0.8*qnewt**(-1.0/(4.0+self.newt-(i+1)-1))*self.h
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
                if self.verbosity >= self.SCREAM:
                    print 'Iteration failed at time %e with step-size %e'%(t,self.h)
                self._nniterfail += 1
                self._rejected = True #The step is rejected
                
                if self._theta >= 0.99:
                    self.h = self.h/2.0
                if self._curjac:
                    self._needjac = False
                    self._needLU = True
                else:
                    self._needjac = True
                    self._needLU = True
        else:
            raise Explicit_ODE_Exception('Newton iteration failed at time %e with step-size %e'%(t,self.h))
        
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
            
        if self._first and err>=1.0:
            h = self.h/10.
        
        if h < self._eps:
            raise Explicit_ODE_Exception('Step-size to small at %e with h = %e'%(self._tc,self.h))
        
        if h > self.maxh:
            h = self.maxh
        
        return h
        
    def estimate_error(self):
        
        temp = 1./self.h*(self.E[0]*self._Z[:self._leny]+self.E[1]*self._Z[self._leny:2*self._leny]+self.E[2]*self._Z[2*self._leny:3*self._leny])

        scal = self._scaling#/self.h
        err_v = N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,self._f0+temp)))
        err = N.linalg.norm(err_v/scal)
        err = max(err/N.sqrt(self._leny),1.e-10)

        if (self._rejected or self._first) and err >= 1.: #If the step was rejected, use the more expensive error estimation
            self._nfcn += 1
            err_v = self.f(self._tc,self._yc+err_v)
            err_v =  N.linalg.solve(self._U1,N.linalg.solve(self._L1,N.linalg.solve(self._P1,err_v+temp)))
            err = N.linalg.norm(err_v/scal)
            err = max(err/N.sqrt(self._leny),1.e-10)
            
        return err
    
    def jacobian(self, t, y):
        """
        Calculates the Jacobian, either by an approximation or by the user
        defined (jac specified in the problem class).
        """
        self._curjac = True #The jacobian is up to date
        self._needLU = True #A new LU-decomposition is needed
        self._needjac = False #A new jacobian is not needed
        
        if self.usejac: #Retrieve the user-defined jacobian
            cjac = self._problem.jac(t,y)
        else:           #Calculate a numeric jacobian
            delt = N.array([(self._eps*max(abs(yi),1.e-5))**0.5 for yi in y])*N.identity(self._leny) #Calculate a disturbance
            Fdelt = N.array([self.f(t,y+e) for e in delt]) #Add the disturbance (row by row) 
            grad = ((Fdelt-self.f(t,y)).T/delt.diagonal()).T
            cjac = N.array(grad).T

            self._njacfcn += 1+self._leny #Add the number of function evaluations
        
        self._njac += 1 #add the number of jacobian evaluation
        return cjac
    
    def interpolate(self, t, k):
        """
        Calculates the continuous output from Radau5.
        """
        leny = self._leny
        s = (t-self._newt)/self._oldh
        Z = self._col_poly
        
        yout = self._yc+s*(Z[:leny]+(s-self.C[1,0]+1.)*(Z[leny:2*leny]+(s-self.C[0,0]+1.)*Z[2*leny:3*leny]))
        return yout
    
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
        
        I = N.eye(self._leny)
        I3 = N.eye(3)
        T1 = N.kron(J,I)
        T2 = N.kron(Tinv,I)
        T3 = N.kron(T,I)
        
        self.A = A
        self.B = B
        self.C = C
        self.I = I
        self.E = E
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3
        self.I3 = I3
        self.EIG = eig
