#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import numpy as N
import scipy.linalg as Sc
import scipy.sparse as sp
import sys
from assimulo.exception import *
from assimulo.ode import *
import logging

from assimulo.explicit_ode import Explicit_ODE

try:
    from assimulo.lib.odepack import dlsodar, dcfode, dintdy
    from assimulo.lib.odepack import set_lsod_common, get_lsod_common
except ImportError:
    sys.stderr.write("Could not find ODEPACK functions.\n")

# Real work array  
class common_like(object):
    def __call__(self):
         return self.__dict__

def g_dummy(t,y):
    return y

def jac_dummy(t,y):
    return N.zeros((len(y),len(y)))

class LSODAR(Explicit_ODE):
    """
        LOSDAR is a multistep method for solving explicit ordinary 
        differential equations on the form,
        
        .. math::
    
            \dot{y} = f(t,y), \quad y(t_0) = y_0.
            
        LSODAR automatically switches between using an ADAMS method
        or an BDF method and is also able to monitor events.
        
        LSODAR is part of ODEPACK, http://www.netlib.org/odepack/opkd-sum
    """

    def __init__(self, problem):
        """
        Initiates the solver.
        
            Parameters::
            
                problem     
                            - The problem to be solved. Should be an instance
                              of the 'Implicit_Problem' class.
        """
        Explicit_ODE.__init__(self, problem) #Calls the base class
        
        #Default values
        self.options["atol"]     = 1.0e-6*N.ones(self.problem_info["dim"]) #Absolute tolerance
        self.options["rtol"]     = 1.0e-6 #Relative tolerance
        self.options["usejac"]   = False
        self.options["maxsteps"] = 100000
        self.options["rkstarter"] = 1
        self.options["maxordn"] = 12
        self.options["maxords"] =  5
        self.options["maxh"] = 0.
        
        self._leny = len(self.y) #Dimension of the problem
        self._nordsieck_array = []
        self._nordsieck_order = 0
        self._nordsieck_time  = 0.0
        self._nordsieck_h  = 0.0
        self._update_nordsieck = False
        
        # Solver support
        self.supports["state_events"] = True
        self.supports["report_continuously"] = True
        self.supports["interpolated_output"] = True
        
        self._RWORK = N.array([0.0]*(22 + self.problem_info["dim"] * 
                               max(16,self.problem_info["dim"]+9) + 
                               3*self.problem_info["dimRoot"]))
        self._IWORK = N.array([0]*(20 + self.problem_info["dim"]))
        
        
    def initialize(self):
        """
        Initializes the overall simulation process
        (called before _simulate) 
        """ 
        #Reset statistics
        self.statistics.reset()
            
        #self._tlist = []
        #self._ylist = []
        
        #starts simulation with classical multistep starting procedure
        # Runge-Kutta starter will be started if desired (see options) 
        # only after an event occured.
        self._rkstarter_active = False
        
    def interpolate(self, t):
        """
        Helper method to interpolate the solution at time t using the Nordsieck history
        array. Wrapper to ODEPACK's subroutine DINTDY.
        """
        #print 'interpolate at t={} and nyh={}'.format(t,self._nyh)
        if self._update_nordsieck:
            #Nordsieck start index
            nordsieck_start_index = 21+3*self.problem_info["dimRoot"] - 1
            
            hu, nqu ,nq ,nyh, nqnyh = get_lsod_common()
            self._nordsieck_array = \
                     self._RWORK[nordsieck_start_index:nordsieck_start_index+(nq+1)*nyh].reshape((nyh,-1),order='F') 
            self._nyh = nyh
            self._update_nordsieck = False
                    
        dky, iflag = dintdy(t, 0, self._nordsieck_array, self._nyh)
        if iflag!= 0 and iflag!=-2:
            raise ODEPACK_Exception("DINTDY returned with iflag={} (see ODEPACK documentation).".format(iflag))   
        elif iflag==-2:
            dky=self.y.copy()
        return dky
     
    def autostart(self,t,y,sw0=[]):
        """
        autostart determines the initial stepsize for Runge--Kutta solvers 
        """
        RWORK=self._RWORK.copy()
        IWORK=self._IWORK.copy()
        pr=self.rkstarter+1 
        tol=self.options["atol"]  
        tolscale=tol[0]**(1./pr)
        normscale=1.
        f=self.problem.rhs
        
        t0=t
        tf=RWORK[0]
        T=abs(tf-t0)
        direction=N.sign(tf-t0)
        
        #Perturb initial condition and compute rough Lipschitz constant
        cent=Sc.norm(y)/normscale/100.
        v0=y+cent*N.random.rand(len(y),1)
        u0prime=f(t,y,sw0)
        v0prime=f(t,v0,sw0)
        Lip=Sc.norm(u0prime-v0prime)/Sc.norm(y-v0)
        h=direction*min(1e-3*T,max(1e-8*T,0.05/Lip))
        #step 1: fwd Euler step
        u1=y+h*u0prime
        t1=t+h
        #step 2: fwd Euler step in reverse
        u1prime=f(t1,u1,sw0)
        u0comp=u1-h*u1prime
        #step 3: estimate of local error
        du=u0comp-y
        dunorm=Sc.norm(du)
        errnorm=dunorm/normscale
        #step 4: new estimate of Lipschitz constant
        u0comprime=f(t0,u0comp,sw0)
        L=Sc.norm(u0comprime-u0prime)/dunorm
        M=N.dot(du,u0comprime-u0prime)/dunorm**2
        #step 5: construct a refined starting stepsize
        theta1=tolscale/N.sqrt(errnorm)
        theta2=tolscale/abs(h*(L+M/2))
        h=h*(theta1+theta2)/2
        h=direction*min(3e-3*T,abs(h))
        return h
        
        
    def integrate_start(self, t, y):
        """
        Helper program for the initialization of LSODAR
        """
        #print ' We have rkstarter {} and rkstarter_active {}'.format(self.rkstarter, self._rkstarter_active)
        if not(self.rkstarter>1 and self._rkstarter_active):
            # first call or classical restart after a discontinuity
            ISTATE=1
            # reset work arrays
            RWORK = 0.0*self._RWORK
            IWORK = 0.0*self._IWORK
        else: #self.rkstarter and self._rkstarter_active
            # RK restart
            RWORK=self._RWORK.copy()
            IWORK=self._IWORK.copy()
            ISTATE=2   #  should be 2  
            dls001=common_like()
            dlsr01=common_like()
            # invoke rkstarter
            # a) get previous stepsize if any
            hu, nqu ,nq ,nyh, nqnyh = get_lsod_common()
            #H = hu if hu != 0. else 1.e-4  # this needs some reflections 
            #H =(abs(RWORK[0]-t)*((self.options["rtol"])**(1/(self.rkstarter+1))))/(100*Sc.norm(self.problem.rhs(t,y,self.sw))+10)#if hu != 0. else 1.e-4
            H=1e-2
            #H=self.autostart(t,y)
            #H=3*H
            # b) compute the Nordsieck array and put it into RWORK
            rkNordsieck = RKStarterNordsieck(self.problem.rhs,H,number_of_steps=self.rkstarter)
            t,nordsieck = rkNordsieck(t,y,self.sw)
            nordsieck=nordsieck.T
            nordsieck_start_index = 21+3*self.problem_info["dimRoot"] - 1
            RWORK[nordsieck_start_index:nordsieck_start_index+nordsieck.size] = \
                                       nordsieck.flatten(order='F')
                        
            # c) compute method coefficients and update the common blocks
            dls001.init = 1
            #dls001.jstart = -1.0    #take the next step with a new value of H,n,meth,..
            mf = 20
            nq = self.rkstarter
            dls001.meth = meth = mf // 10
            dls001.miter =mf % 10
            elco,tesco =dcfode(meth)  #  
            dls001.maxord= 5      #max order 
            dls001.nq= self.rkstarter          #Next step order 
            dls001.nqu=self.rkstarter           #Method order last used  (check if this is needed)
            dls001.meo= meth      #meth
            dls001.nqnyh= nq*self.problem_info["dim"]    #nqnyh
            dls001.conit= 0.5/(nq+2)                     #conit   
            dls001.el= elco[:,nq-1]  # el0 is set internally
            dls001.hu=H  # this sets also hold and h internally
            dls001.jstart=1
            
            # IWORK[...] =  
            #IWORK[13]=dls001.nqu
            IWORK[14]=dls001.nq
            #IWORK[18]=dls001.meth
            #IWORK[7]=dlsa01.mxordn    #max allowed order for Adams methods
            #IWORK[8]=dlsa01.mxords    #max allowed order for BDF
            IWORK[19]=meth         #the current method indicator
            #RWORK[...]
            dls001.tn=t
            RWORK[12]=t
            RWORK[10]=hu         #step-size used successfully
            RWORK[11]=H         #step-size to be attempted for the next step 
            #RWORK[6]=dls001.hmin
            #RWORK[5]=dls001.hmxi
            
            number_of_fevals=N.array([1,2,4,7,11])
            # d) Reset statistics
            IWORK[9:13]=[0]*4
            dls001.nst=1
            dls001.nfe=number_of_fevals[self.rkstarter-1]   # from the starter
            dls001.nje=0
            dlsr01.nge=0
            # set common block
            commonblocks={}
            commonblocks.update(dls001())
            commonblocks.update(dlsr01())
            set_lsod_common(**commonblocks)
            
 
        return ISTATE, RWORK, IWORK
                                    
    def _jacobian(self, t, y):
        """
        Calculates the Jacobian, either by an approximation or by the user
        defined (jac specified in the problem class).
        """
        jac = self.problem.jac(t,y)
        
        if isinstance(jac, sp.csc_matrix):
            jac = jac.toarray()
        
        return jac
    
    def integrate(self, t, y, tf, opts):
        ITOL  = 2 #Only  atol is a  vector
        ITASK = 5 #For one step mode and hitting exactly tcrit, normally tf
        IOPT = 1 #optional inputs are used
        
        # provide work arrays and set common blocks (if needed)
        ISTATE, RWORK, IWORK = self.integrate_start( t, y)
        
        JT = 1 if self.usejac else 2#Jacobian type indicator
        JROOT = N.array([0]*self.problem_info["dimRoot"])
        
        #Setting work options
        RWORK[0] = tf #Do not integrate past tf
        RWORK[5] = self.options["maxh"]
        
        #Setting iwork options
        IWORK[5] = self.maxsteps
        
        #Setting maxord to IWORK
        IWORK[7] = self.maxordn
        IWORK[8] = self.maxords

        
        #Dummy methods
        #g_dummy = (lambda t:x) if not self.problem_info["state_events"] else self.problem.state_events
        g_fcn = g_dummy if not self.problem_info["state_events"] else self.problem.state_events
        #jac_dummy = (lambda t,y:N.zeros((len(y),len(y)))) if not self.usejac else self.problem.jac
        jac_fcn = jac_dummy if not self.usejac else self._jacobian
        
        #Extra args to rhs and state_events
        rhs_extra_args = (self.sw,) if self.problem_info["switches"] else ()
        g_extra_args = (self.sw,) if self.problem_info["switches"] else ()
        
        #Store the opts
        self._opts = opts
        
        #Outputs
        tlist = []
        ylist = []
        
        
        
        #Run in normal mode?
        normal_mode = 1 if opts["output_list"] is not None else 0
        
        #Tolerances:
        atol = self.atol
        rtol = self.rtol*N.ones(self.problem_info["dim"])
        rhs = self.problem.rhs
        
        #if normal_mode == 0:
        if opts["report_continuously"] or opts["output_list"] is None:
            
            while (ISTATE == 2 or ISTATE == 1) and t < tf:
            
                y, t, ISTATE, RWORK, IWORK, roots = dlsodar(rhs, y.copy(), t, tf, ITOL, 
                        rtol, atol,
                        ITASK, ISTATE, IOPT, RWORK, IWORK, jac_fcn, JT, g_fcn, JROOT,
                        f_extra_args = rhs_extra_args, g_extra_args = g_extra_args)
                
                self._update_nordsieck = True
                self._IWORK = IWORK
                self._RWORK = RWORK
                #hu, nqu ,nq ,nyh, nqnyh = get_lsod_common()
                #self._nordsieck_array = \
                #     RWORK[nordsieck_start_index:nordsieck_start_index+(nq+1)*nyh].reshape((nyh,-1),order='F') 
                #self._nyh = nyh              
                self._event_info = roots
                
                if opts["report_continuously"]:
                    flag_initialize = self.report_solution(t, y, opts)
                    if flag_initialize:
                        #If a step event has occured the integration has to be reinitialized
                        ISTATE = 3
                else:
                    #Store results
                    tlist.append(t)
                    ylist.append(y.copy())
            
                #Checking return
                if ISTATE == 2:
                    flag = ID_PY_COMPLETE
                elif ISTATE == 3:
                    flag = ID_PY_EVENT
                else:
                    raise ODEPACK_Exception("LSODAR failed with flag %d"%ISTATE)
            
        else:
            
            #Change the ITASK
            ITASK = 4 #For computation of yout
            
            output_index = opts["output_index"]
            output_list  = opts["output_list"][output_index:]
            
            flag = ID_PY_COMPLETE

            for tout in output_list:
                output_index += 1

                y, t, ISTATE, RWORK, IWORK, roots = dlsodar(rhs, y.copy(), t, tout, ITOL, 
                    rtol, atol,
                    ITASK, ISTATE, IOPT, RWORK, IWORK, jac_fcn, JT, g_fcn, JROOT,
                    f_extra_args = rhs_extra_args, g_extra_args = g_extra_args)
                
                #Store results
                tlist.append(t)
                ylist.append(y.copy())
                self._event_info = roots
                
                #Checking return
                if ISTATE == 2 and t >= tf:
                    flag = ID_PY_COMPLETE
                    break
                elif ISTATE == 3:
                    flag = ID_PY_EVENT
                    break
                elif ISTATE < 0:
                    raise ODEPACK_Exception("LSODAR failed with flag %d"%ISTATE)
            
            opts["output_index"] = output_index
        # deciding on restarting options
        self._rkstarter_active = True if ISTATE == 3 and self.rkstarter > 1 else False
        #print 'rkstarter_active set to {} and ISTATE={}'.format(self._rkstarter_active, ISTATE)
        
        #Retrieving statistics
        self.statistics["nstatefcns"]            += IWORK[9]
        self.statistics["nsteps"]        += IWORK[10]
        self.statistics["nfcns"]          += IWORK[11]
        self.statistics["njacs"]          += IWORK[12]
        self.statistics["nstateevents"] += 1  if flag == ID_PY_EVENT else 0
        # save RWORK, IWORK for restarting feature
        if self.rkstarter>1:
            self._RWORK=RWORK
            self._IWORK=IWORK        
       
        return flag, tlist, ylist
        
    def get_algorithm_data(self):
        """
        Returns the order and step size used in the last successful step.
        """
        hu, nqu ,nq ,nyh, nqnyh = get_lsod_common()
            
        return hu, nqu
    
    def state_event_info(self):
        """
        Returns the state events indicator as a list where a one (1)
        indicates that the event function have been activated and a (0)
        if not.
        """
        return self._event_info

    def print_statistics(self, verbose=NORMAL):
        """
        Prints the run-time statistics for the problem.
        """
        Explicit_ODE.print_statistics(self, verbose) #Calls the base class

        self.log_message('\nSolver options:\n',                                      verbose)
        self.log_message(' Solver                  : LSODAR ',         verbose)
        self.log_message(' Absolute tolerances     : {}'.format(self.options["atol"]),  verbose)
        self.log_message(' Relative tolerances     : {}'.format(self.options["rtol"]),  verbose)
        if self.rkstarter==1:
           self.log_message(' Starter                 : {}'.format('classical'),  verbose)
        else:
           self.log_message(' Starter                 : Runge-Kutta order {}'.format(self.rkstarter),  verbose)
        if self.maxordn < 12 or self.maxords < 5:
            self.log_message(' Maximal order Adams     : {}'.format(self.maxordn),  verbose)
            self.log_message(' Maximal order BDF       : {}'.format(self.maxords),  verbose)
        if self.maxh > 0. :
            self.log_message(' Maximal stepsize maxh   : {}'.format(self.maxh),  verbose)
        self.log_message('',                                                         verbose)
    
    def _set_usejac(self, jac):
        self.options["usejac"] = bool(jac)
    
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
        return self.options["usejac"]
    
    usejac = property(_get_usejac,_set_usejac)
    
    def _set_atol(self,atol):
        
        self.options["atol"] = N.array(atol,dtype=N.float) if len(N.array(atol,dtype=N.float).shape)>0 else N.array([atol],dtype=N.float)
    
        if len(self.options["atol"]) == 1:
            self.options["atol"] = self.options["atol"]*N.ones(self._leny)
        elif len(self.options["atol"]) != self._leny:
            raise ODEPACK_Exception("atol must be of length one or same as the dimension of the problem.")

    def _get_atol(self):
        """
        Defines the absolute tolerance(s) that is to be used by the solver.
        Can be set differently for each variable.
        
            Parameters::
            
                atol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float or a numpy vector
                          of floats.
                        
                            Example:
                                atol = [1.0e-4, 1.0e-6]
        """
        return self.options["atol"]
    
    atol=property(_get_atol,_set_atol)
    
    def _set_rtol(self,rtol):
        try:
            self.options["rtol"] = float(rtol)
        except (ValueError, TypeError):
            raise ODEPACK_Exception('Relative tolerance must be a (scalar) float.')
        if self.options["rtol"] <= 0.0:
            raise ODEPACK_Exception('Relative tolerance must be a positive (scalar) float.')
    
    def _get_rtol(self):
        """
        Defines the relative tolerance that is to be used by the solver.
        
            Parameters::
            
                rtol    
                        - Default '1.0e-6'.
                
                        - Should be a positive float.
                        
                            Example:
                                rtol = 1.0e-4
        """
        return self.options["rtol"]
        
    rtol=property(_get_rtol,_set_rtol)

    def _get_maxsteps(self):
        """
        The maximum number of steps allowed to be taken to reach the
        final time.
        
            Parameters::
            
                maxsteps
                            - Default 10000
                            
                            - Should be a positive integer
        """
        return self.options["maxsteps"]
    
    def _set_maxsteps(self, max_steps):
        try:
            max_steps = int(max_steps)
        except (TypeError, ValueError):
            raise ODEPACK_Exception("Maximum number of steps must be a positive integer.")
        self.options["maxsteps"] = max_steps
    
    maxsteps = property(_get_maxsteps, _set_maxsteps)
    def _get_hmax(self):
        """
        The absolute value of the maximal stepsize for all methods
        
        Parameters::
        
               hmax
                          - Default:  0.  (no maximal step size)
                          
                          - Should be a positive float
        """
        logging.warning("The option 'hmax' is deprecated and will be removed, please use 'maxh' instead.")
        return self.options["maxh"]
    def _set_hmax(self,hmax):
        if not (isinstance(hmax,float) and hmax >= 0.):
           raise ODEPACK_Exception("Maximal step size hmax should be a positive float")
        logging.warning("The option 'hmax' is deprecated and will be removed, please use 'maxh' instead.")
        self.options["maxh"]=hmax
    hmax = property(_get_hmax, _set_hmax)
    def _get_maxh(self):
        """
        The absolute value of the maximal stepsize for all methods
        
        Parameters::
        
               maxh
                          - Default:  0.  (no maximal step size)
                          
                          - Should be a positive float
        """
        return self.options["maxh"]
    def _set_maxh(self,maxh):
        if not (isinstance(maxh,float) and maxh >= 0.):
           raise ODEPACK_Exception("Maximal step size maxh should be a positive float")
        self.options["maxh"]=maxh
    maxh = property(_get_maxh, _set_maxh)
    def _get_maxordn(self):
        """
        The maximum order used by the Adams-Moulton method (nonstiff case)
        
            Parameters::
            
                maxordn
                            - Default 12
                            
                            - Should be a positive integer
        """
        return self.options["maxordn"]
    
    def _set_maxordn(self, maxordn):
        try:
            maxordn = int(maxordn)
        except (TypeError, ValueError):
            raise ODEPACK_Exception("Maximum order must be a positive integer.")
        if maxordn > 12:
            raise ODEPACK_Exception("Maximum order should not exceed 12.")    
        self.options["maxordn"] = maxordn
    
    maxordn = property(_get_maxordn, _set_maxordn)
    def _get_maxords(self):
        """
        The maximum order used by the BDF method (stiff case)
        
            Parameters::
            
                maxords
                            - Default 5
                            
                            - Should be a positive integer
        """
        return self.options["maxords"]
    
    def _set_maxords(self, maxords):
        try:
            maxords = int(maxords)
        except (TypeError, ValueError):
            raise ODEPACK_Exception("Maximum order must be a positive integer.")
        if maxords > 5:
            raise ODEPACK_Exception("Maximum order should not exceed 5.")    
        self.options["maxords"] = maxords
    
    maxords = property(_get_maxords, _set_maxords)
    def _get_rkstarter(self):
        """
        This defines how LSODAR is started. 
        (classically or with a fourth order Runge-Kutta starter)
        
            Parameters::
            
                rkstarter
                            - Default False  starts LSODAR in the classical multistep way
                            
                            - Should be a Boolean
        """
        return self.options["rkstarter"]
    
    def _set_rkstarter(self, rkstarter):
        if not rkstarter in {1,2,3,4,5}:
            raise ODEPACK_Exception("Must be a positive integer less than 6")
        self.options["rkstarter"] = rkstarter
    
    rkstarter = property(_get_rkstarter, _set_rkstarter)

class RKStarterNordsieck(object):
    """
    A family of Runge-Kutta starters producing a 
    Nordsieck array to (re-)start a Nordsieck based multistep
    method with a given order.
    
    See: Mohammadi (2013): https://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=4196026&fileOId=4196027
    """
    # Gamma matrix of Gear's RK starter which produce Nordsieck vector at t0
    Gamma_0=[N.array([[1.,0.],                        # 1st order
                      [0.,1.]]),
             N.array([[1.,0.,0.],                     # 2nd order
                      [0.,1.,-1.],
                      [0.,0.,1.]]),
             N.array([[1.,0.,0.,0.],                  # 3rd order
                      [0.,1.,-5./3.,1.],         
                      [0.,0.,3.,-2.],
                      [0.,0.,0.,1.],
                      [0.,0.,-4./3.,0.]]),
             N.array([[1.,0.,0.,0.,0.],               # 4th order
                      [0.,1.,-5./6.,4./9.,-1./9.],         
                      [0.,0.,0.,0.,0.],
                      [0.,0.,1./2.,-4./9.,1./9.],
                      [0.,0.,7./3.,-19./9.,7./9.],
                      [0.,0.,-3.,10./3.,-4./3.],
                      [0.,0.,1.,-11./9.,5./9.]])]
                      
    # A matrices of RK starter with equidistanced states 
    
    A_s=[    N.array([1]),                                           # 1st order
             N.array([[0.,0],[1,0]]),                                # 2nd order
             N.array([[0.,0.,0.,0.,0.],[1./2,0.,0.,0.,0.],          # 3rd order
                     [0.,3./4,0.,0.,0.],[2./9,1./3,4./9,0.,0.],
                     [17./72,1./6,2./9,-1./8,0.]]),
             N.array([[0.,0.,0.,0.,0.,0.,0.],                       # 4th order
                      [1./6.,0.,0.,0.,0.,0.,0.],
                      [0.,1./6.,0.,0.,0.,0.,0.],
                      [0.,0.,1./3.,0.,0.,0.,0.],
                      [1./18.,1./9.,1./9.,1./18.,0.,0.,0.],
                      [2.5,-3.,-3.,2.25,2.25,0.,0.],
                      [10./45.,-8./45.,-8./45.,-4./45.,13./15.,1./45.,0.]]),
             N.array([[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],  # 5th order
                      [1./20,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                      [3./160,9./160,0.,0.,0.,0.,0., 0.,0.,0.,0.,0.,0.,0.],
                      [3./40,-9./40,6./20,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                      [-11./216,5/8,-70/108,35./108,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                      [1631/221184,175./2048,575./55296,44275./442368,253./16384,0.,0.,0.,0.,0.,0.,0.,0.,0.],
                      [37./1512,0.,250/2484,125./2376,0.,512/7084,0.,0.,0.,0.,0.,0.,0.,0.],
                      [32.2687682,-39.4956646,-20.7849921,22.2296308,45.4543383,51.9961815,-91.1682621,0.,0.,0.,0.,0.,0.,0.],
                      [-23.867154567764317,41+31195543/31462660,8+6080897/76520184,-23-7543649/15003372,-63-91662407/130833420,-55-27823937/32947321,117+23615874/27978401,0.,0.,0.,0.,0.,0.,0.],
                      [-0.5118403967750724,0.,2+379808/6250467,-2-866339/6430593,-1-151963/7810968,-18723895/33225736,8/3,0.,303385/143008499,0.,0.,0.,0.,0.],
                      [-118.50992840087292,4+1586290/42052551,7+76050433/143964928,11+61385198/111742353,16+54354964/63299173,15+84679766/214253769,15+71618817/71932238,24+14108046/149713813,
                       2+12784603/81260661,21+38484710/59808171,0.,0.,0.,0.],
                      [10.183365237189525,0.,-36-34947657/185923229,38+50251318/55929787,18+114017528/325586295,10+78529445/98887874,-43.5,0.,-6714159/175827524,2.25,0.,0.,0.,0.],
                      [-0.1491588850008383,0.,19145766/113939551,8434687/36458574,2012005/39716421,8409989/57254530,10739409/94504714,0.,-3321/86525909,-30352753/150092385,0.,70257074/109630355,
                       0.,0.],
                      [0.03877310906055409,0.,3021245/89251943,5956469/58978530,851373/32201684,11559106/149527791,11325471/112382620,0.,-12983/235976962,17692261/82454251,0.,
                       38892959/120069679,11804845./141497517]])]  

    co_ord_s=[[],[],[4],[4,6],[6,9,11]]                       
    b_s=[    N.array([1]),
             N.array([1./2,1./2]),
             N.array([1./6,0.,0.,1./6,2./3]),
             N.array([29./1062.,83./531.,83./531.,83./1062.,2./531.,56./531.,251./531.]),  # b-vectors of of 1st to 5th order RK starter method
             N.array([0.03877310906055409,0.,3021245/89251943,5956469/58978530,851373/32201684,11559106/149527791,11325471/112382620,0.,-12983/235976962,17692261/82454251,0.,
                       38892959/120069679,11804845./141497517,0.])] 
                               
    C_s=[    N.array([0.]),
             N.array([0.]),
             N.array([0.,1./2,3./4,1.,1./2,1.]),
             N.array([0.,1./6,1./6,1./3,1./3,1.,1.,2./3,1.]),
             N.array([0.,1./20,3./40,3./20,1./4,7./28,1./4,2./4,1.,2./4,3./4,3./4,1.,1.])]    # C-values in Butcher tableau of 8-stages Runge-Kutta                           
            
    # A matrices of RK starter with non-equidistanced stages        
    A_n=[    N.array([1]),
             N.array([[0.,0.],[0.,0.]]),
             N.array([[0.,0.,0.,0.],[1./2,0.,0.,0.],[0.,1./2,0.,0.],[0.,0.,1.,0.]]),           # 3rd order
             N.array([[0.,0.,0.,0.,0.,0.],[2./5,0.,0.,0.,0.,0.],[-3./20,3./4,0.,0.,0.,0.],    # 4th order
                      [19./44,-15./44,40./44,0.,0.,0.],[-31./64,185./192,5./64,-11./192,0.,0.],[11./72,25./72,25./72,11./72,0.,0.]])]
             
    b_n=[    N.array([0.]),
             N.array([1./2.,1./2.]),
             N.array([[5./24,1./6,1./6,-1/24],[1./6,1./3,1./3,1./6]]),
             N.array([[802./5625,68./225,-67./225,-143./5625,144./625,6./125],[699./5000,81./200,-39./200,99./5000,144./625,0.],[11./72,25./72,25./72,11./72,0.,0.]])]   
             
    C_n=[    N.array([0.]),
             N.array([0.]),
             N.array([0.,1./2,1./2,1.]),
             N.array([0.,2./5,3./5,1.,1./2,1.])]

    #co_ord_n=[[],[],[1./2,1.],[2./5,3./5,1.]]
    #A=N.array([[1.,0.,0.,0.],
    #          [1.,1./9.,1./27,1./81.],
    #           [1.,4./9.,8./27,16./81.],
    #           [1.,1.,1.,1.]]) 
    A=[N.array([0.]),                       # Convert the state values to Nordsieck vector
       N.array([0.]),
       N.array([1.]),
       N.array([[1./4,1./8],[1.,1.]]),
       N.array([[1./9,1./27,1./81.],
                [4./9.,8./27,16./81],
                [1.,1.,1.]]),
       N.array([[1./16,1./64,1./256,1./1024],
                [1./4,1./8,1./16,1./32],
                [9./16,27./64,81./256,243./1024],
                [1.,1.,1.,1.]])]
                      
    scale=N.array([1, 1, 1/2., 1./6., 1./24., 1./120.]).reshape(-1,1)
    
    
    def __init__(self,  rhs, H, method='RKs_f', eval_at=0., number_of_steps=4):
        """
        Initiates the Runge-Kutta Starter.
        
            Parameters::
                            
                rhs     
                            - The problem's rhs function
                             
                H
                            - starter step size
                            
                eval_at
                            - Where to evaluate the Nordiseck vector
                              evaluation point is (eval_at)*H
                              Possible choices 
                              eval_at=0.
                              eval_at=1.
             
                number_of_steps
                            - the number of steps :math:`(k)` to start the multistep method with
                              This will determine the initial order of the method, which is
                              :math:`k` for BDF methods and :math:`k+1` for Adams Moulton Methods   .
        """
        self.f = rhs
        self.H = H
        self.method=method
        
        if not 1 < number_of_steps < 6:
            raise RKStarter_Exception('Step number larget than 4 not yet implemented')            
        self.number_of_steps = number_of_steps 
        self.eval_at = float(eval_at)
        if not self.eval_at == 0.:
           raise RKStarter_Exception("Parameter eval_at different from 0 not yet implemented.")

    def RKs_f(self,t0,y0,sw0):

        s=self.number_of_steps
        A_s,C_s,b_s,co_ord_s=self.A_s,self.C_s,self.b_s,self.co_ord_s
        A_s=A_s[s-1]
        C_s=C_s[s-1]
        b_s=b_s[s-1]
        co_ord_s=co_ord_s[s-1]
        H=(s-1)*self.H
        K=N.zeros((N.size(A_s,0),len(y0)))
        for i in range(N.size(A_s,0)):
            K[i,:]=self.f(t0+C_s[i]*H,y0+H*N.dot(A_s[i,:],K),sw0)   
        y=N.zeros((s,len(y0)))
        y[0,:]=y0
        for i in range(1,s):
            if i==s-1:
                y[i,:]=y0+H*N.dot(b_s,K)
            else:
                y[i,:]=y0+H*N.dot(A_s[co_ord_s[i-1],:],K)
        return y    
    def RKn_f(self,t0,y0,sw0):
        s=self.number_of_steps
        H=(s-1)*self.H
        A_n,C_n,b_n=self.A_n,self.C_n,self.b_n
        A_n=A_n[s-1]
        C_n=C_n[s-1]
        b_n=b_n[s-1]
        
        K=N.zeros((N.size(A_n,0),len(y0)))
        for i in range(N.size(A_n,0)):
            K[i,:]=self.f(t0+C_n[i]*H,y0+H*N.dot(A_n[i,:],K),sw0)
        y=N.zeros((s,len(y0)))  
        y[0,:]=y0
        for i in range(1,s):
                y[i,:]=y0+H*N.dot(b_n[i-1],K)
        return y
  
            
    def rk_like4(self, t0, y0, sw0): 
        """
        rk_like computes Runge-Kutta stages
        Note, the currently implementation is **only** correct for
        autonomous systems.
        """
        f = lambda y: self.f(t0, y, sw0)
        h = self.H/4.
        k1 = h*f(y0)
        k2 = h*f(y0 + k1)
        k3 = h*f(y0 + 2. * k2)
        k4 = h*f(y0 + 3./4. * k1 + 9./4. * k3)
        k5 = h*f(y0 + k1/2. + k2 + k3/2. + 2. * k4)
        k6 = h*f(y0+k1/12.+2. * k2 + k3/4. + 2./3. * k4 + 2. * k5)
        return N.array([y0,k1,k2,k3,k4,k5,k6])
    def rk_like3(self, t0, y0, sw0): 
        """
        rk_like computes Runge-Kutta stages
        Note, the currently implementation is **only** correct for
        autonomous systems.

        """
        f = lambda y: self.f(t0, y, sw0)
        h = self.H/3.
        k1 = h*f(y0)
        k2 = h*f(y0 + k1)
        k3 = h*f(y0 + k1+ k2)
        k4 = h*f(y0 + 3./2. * k1)
        return N.array([y0,k1,k2,k3,k4])
    def rk_like2(self, t0, y0, sw0):
        """
        rk_like2 computes Runge-Kutta 2nd-stages
        Note, the currently implementation is **only** correct for
        autonomous systems.
        """
        f=lambda y: self.f(t0, y, sw0)
        h=self.H/2.
        k1=h*f(y0)
        k2=h*f(y0+k1)
        return N.array([y0,k1,k2])
    def rk_like13(self, t0, y0, sw0):
        """
        rk_like6 computes Runge-Kutta 8th-stages 
        """
        h = self.H
        self.Gamma_2=self.Gamma_0[3]
        f=lambda y: self.f(t0 , y , sw0)
        K=N.zeros((6,len(y0)))
        sol=N.zeros((3,len(y0)))
        b=N.zeros((2,len(y0)))          #remove the fifth stage value that is for error estimation
        nord = N.zeros((4,len(y0)))     #Nordsieck vector
        for i in range(5):
            K[i,:]= f(y0+h*N.dot(self.Gamma_2[i,:],K))
        c=0
        for i in range(3):
            sol[i,:]=y0+h*N.dot(self.Gamma_2[i+3,:],K)
            if i!=0:
                b[c,:]=sol[i,:]-y0-(c+1)*h/2*K[0,:]
                c+=1
        nord[0,:] = y0
        nord[1,:] = h*K[0,:]
        nord[2:,:] = Sc.solve(self.A[self.number_of_steps],b)
        return nord     
    def rk_like14(self, t0, y0, sw0):
        """
        rk_like6 computes Runge-Kutta 8th-stages 
        """
        h = self.H
        Gamma_2=self.Gamma_0[4]
        f=lambda y: self.f(t0 , y , sw0)
        K=N.zeros((8,len(y0)))
        sol=N.zeros((4,len(y0)))
        b=N.zeros((3,len(y0)))          #remove the fifth stage value that is for error estimation
        nord = N.zeros((5,len(y0)))     #Nordsieck vector
        for i in range(7):
            K[i,:]= f(y0+h*N.dot(Gamma_2[i,:],K))
        c=0
        for i in range(4):
            sol[i,:]=y0+h*N.dot(Gamma_2[i+4,:],K)
            if i!=1:
                b[c,:]=sol[i,:]-y0-(c+1)*h/3*K[0,:]
                c+=1
        nord[0,:] = y0
        nord[1,:] = h*K[0,:]
        nord[2:,:] = Sc.solve(self.A[self.number_of_steps],b)
        return nord       
    def rk_like15(self, t0, y0, sw0):
        """
        rk_like6 computes Runge-Kutta 5th-stages ****needs to be modified****
        """
        h = self.H
        Gamma_2=self.Gamma_0[5]
        f=lambda y: self.f(t0 , y , sw0)
        K=N.zeros((14,len(y0)))
        sol=N.zeros((8,len(y0)))
        b=N.zeros((4,len(y0)))          #remove the fifth stage value that is for error estimation
        nord = N.zeros((6,len(y0)))     #Nordsieck vector
        for i in range(13):
            K[i,:]= f(y0+h*N.dot(Gamma_2[i,:],K))
        c=0
        for i in range(8):
            sol[i,:]=y0+h*N.dot(Gamma_2[i+6,:],K)
            if (i!=1) and (i!=2) and (i!=4) and (i!=6):
                b[c,:]=sol[i,:]-y0-(c+1)*h/4*K[0,:]
                c+=1
        nord[0,:] = y0
        nord[1,:] = h*K[0,:]
        nord[2:,:] = Sc.solve(self.A[self.number_of_steps],b)        
        return nord     
    def nordsieck(self,k):
        """
        Nordsieck array computed at initial point
        """
        nord=self.scale[:self.number_of_steps+1]*N.dot(self.Gamma_0[self.number_of_steps-1].T,k)
 
        return nord  
    def Nordsieck_RKn(self,t0,y,sw0):
        s=self.number_of_steps
        H=(s-1)*self.H
        co_nord=[N.array([1./2,1.]),N.array([2./5,3./5,1.])]
        l=size(y,0)
        y0=y[0,:]
        yf=self.f(t0,y0,sw0)
        
        if l==3:
            co=N.array([co_nord[0]])
            nord_n=N.vander(co_nord[0],self.number_of_steps+1)
            b=y[1:]-y0-co.T*yf
            nord=Sc.solve(nord_n[0:2,0:2],b)
        elif l==4:
            co=N.array([co_nord[1]])
            nord_n=N.vander(co_nord[1],self.number_of_steps+1)
            b=y[1:]-y0-H*co.T*yf
            nord=Sc.solve(nord_n[0:3,0:3],b)
        nord=N.vstack((y0,H*yf,nord[::-1]))       
        return nord
    def Nordsieck_RKs(self,t0,y,sw0):
        s=self.number_of_steps
        H=(s-1)*self.H
        co_nord=[N.array([1]),N.array([1./2,1]),N.array([1./3,2./3,1]),
                 N.array([1./4,2./4,3./4,1.])]
        A=self.A
        y0=y[0,:]
        yf=self.f(t0,y0,sw0)
        co=co_nord[s-2]
        co=N.array([co])
        b=y[1:]-y0-H*co.T*yf
        nord=Sc.solve(A[s],b)
        nord=N.vstack((y0,H*yf,nord))
        return nord
        
        
        
    def __call__(self, t0 , y0, sw0=[]):
        """
        Evaluates the Runge-Kutta starting values
        
            Parameters::
            
                y0   
                    - starting value
        """
        if self.method=='RK_G':
        # We construct a call like: rk_like4(self, t0, y0, sw0) 
            k=self.__getattribute__('rk_like{}'.format(self.number_of_steps))(t0, y0, sw0)
            t = t0+self.eval_at*self.H
            #t= t0 + self.H
            k=self.nordsieck(k)
            
        elif self.method=='RKs_f':
            y=self.RKs_f(t0, y0, sw0)
            t = t0+self.eval_at*self.H
            k=self.Nordsieck_RKs(t0,y,sw0)

        elif self.method=='RKn_f':
            y=self.RKn_f(t0,y0,sw0)
            t=t0+self.eval_at*self.H
            k=self.Nordsieck_RKn(t0,y,sw0)
        return t,k 


                
