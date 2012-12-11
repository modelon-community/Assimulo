from scipy import *
import assimulo.lib.dasp3dp as dasp
import numpy as np

# Just calling DASP3 from Python (no Assimulo)

def dydt(t,y,z):
	eps=(1./3)*1.e-3
	yp = array([-(.6*z[0]+.8*y[2])*y[0]+10.*y[1],
				  -10.*y[1]+ 1.6*z[0] *y[2],
				  -1.33*eps**2*y[2]*(y[0]+2.*z[0])])										 
	return yp
def dzdt(t,y,z):
	eps=(1./3)*1.e-3
	zp = array([1.6*z[0]*y[2]-.6*z[0]*y[0]
				  -45.*(eps*z[0])**2+.8*y[2]*y[0]])
	return zp
def outpda(t,y,z,n,m,jstop):
	print 'outpda  :', t, y[0:3], z[0]

#	Main Program to call DASP3
n = 3
m = 1
t = 0.
tend = 10.

#...INITIAL VALUEs
wsy=empty((10*n,))
wsy[0:3]=[ 3., .216, 1.]
wsz=empty((9*m,))
wsz[0] = 1.35

#...ERROR TOLERANCE PARAMETERS
tol = 1.e-3
absrel=ones((4,))
wght=ones((4,))

#...SINGULAR PERTURBATION PARAMETER
eps = array([.33333333e-3])

#...SIMULATION OF BELOUSOV-ZHABOTINSKII REACTION

a = empty((m,m))
w = empty((m,m))
slu= empty((2*m,))
ips= empty((m,),'int32')
ind = empty((2*m,),'int32')
eq= empty((m,),'bool')
print tol
t,lflag=dasp.dasp3(dydt,dzdt,outpda,t,tend,wsy,wsz,n,m,tol,absrel,wght,eps,a,w,slu,ips,eq,ind)
print 't, lflag', t, lflag
stats= """ 
Statistics
==========
NSTEPS {nsteps}  NYDERIV {nyderiv}  NZDERIV {nzderiv}
Rejected steps {nrej}
"""
print stats.format(nsteps=dasp.COUNTS.NSTEP, 
                   nyderiv=dasp.COUNTS.NYDER,
                   nzderiv=dasp.COUNTS.NZDER,
                   nrej=dasp.COUNTS.NREJ)
