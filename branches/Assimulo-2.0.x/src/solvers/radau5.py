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
from assimulo.lib.radau_core import Radau_Common

class ExplicitRadau5(Radau_Common,Explicit_ODE):
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
            dist_space = [t+(x+1)*(tf-t)/ncp for x in range(int(ncp)+1)]
        
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



class ImplicitRadau5(Radau_Common,Implicit_ODE):
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
        self._leny = len(self.y_cur) #Dimension of the problem
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
        self.index = [1]*self._leny
        self.atol = 1.0e-6 #Absolute tolerance
        self.rtol = 1.0e-6 #Relative tolerance
        
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
        self._type = '(implicit)'
        
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
            self._index = N.array(index+(N.array(index)+1).tolist())
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
    
    def _integrator(self, t, y, yd, tf,dt):
        """
        Integrates (t,y,yd) values until t > tf
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
            dist_space = [t+(x+1)*dt for x in range(int((tf-t)/dt)+1)]
        
        for i in range(self.maxsteps):
            if t >= tf:
                break
            t, y, yd = self.step(t, y, yd)
            self._tc = t
            self._yc = y
            self._ydc = yd
            
            if dt > 0.0:
                while dist_space[0] <= t:
                    yy  = self.interpolate(dist_space[0],0)
                    yyd = self.interpolate(dist_space[0],1)
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
                if self.thet < 0:
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
        
        if k==0:
            return yout
        elif k==1:
            return ydout
        else:
            raise Implicit_ODE_Exception('Unknown value of k. Should be either 0 or 1')

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
