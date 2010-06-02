from Assimulo.Explicit_ODE import *
from Assimulo.Problem import Explicit_Problem
from numpy import array, pi, exp, sin, cos, abs


def bridge(t,u):
	"""
	Example problem:
		
	Tacoma bridge (Sauer, p. 328 ff)
	"""
	len=6 # half width
	a=0.2 # Hooke's nonlinearity coeff
	W=80  # wind speed km/h
	d=0.01 # damping ratio
	k=1000 # stiffness
	m=2500 # mass
	omega=2*pi*38/60
	kma=k/(m*a)
	y,yd,theta,thetad=u
	a1=exp(a*(y-len*sin(theta)))
	a2=exp(a*(y+len*sin(theta)))
	ydd=-d*yd-kma*(a1+a2-2)+0.2*W*sin(omega*t)
	thetadd=-d*thetad+3*cos(theta)/len*kma*(a1-a2)
	return array([yd,ydd,thetad,thetadd])
bridge.name='Tacoma Narrows Bridge'

problem = Explicit_Problem()
problem.f = bridge
problem.problem_name = bridge.name

def test_cvode_BN_tacoma():
	"""
	Tacoma bridge CVode BDF Newton
	See Exercise 1 at Sauer page 330
	"""
	cv=CVode(problem,[0,0,0.001,0])
	cv.discr='BDF'
	cv.iter='Newton'
	cv.atol=1.e-6
	cv.rtol=1.e-9
	print cv.discr,cv.iter,cv.atol,cv.rtol,'\n'
	cv.simulate(10,200)
	print cv.y[-1][0]
	assert abs(cv.y[-1][0]- (-0.3615)) < 2.e-3
	assert cv.t[-1]==10
def test_cvode_BF_tacoma():
	"""
	Tacoma bridge CVode BDF Fixed Point
	See Exercise 1 at Sauer page 330
	"""
	cv=CVode(problem,[0,0,0.001,0])
	cv.discr='BDF'
	cv.atol=1.e-4
	cv.rtol=1.e-4
	print cv.discr,cv.iter,cv.atol,cv.rtol,'\n'
	cv.simulate(10,200)
	#cv.plot()
	assert abs(cv.y[-1][0]- (-0.3615)) < 2.e-3
	assert cv.t[-1]==10
def test_cvode_AF_tacoma():
	"""
	Tacoma bridge CVode Adams Fixed Point
	See Exercise 1 at Sauer page 330
	"""
	cv=CVode(problem,[0,0,0.001,0])
	#cv.atol=1.e-4
	cv.rtol=1.e-4
	print cv.discr,cv.iter,cv.atol,cv.rtol,'\n'
	cv.simulate(10,200)
	cv.print_statistics()
	#cv.plot()
	assert abs(cv.y[-1][0]- (-0.3615)) < 2.e-3
	assert cv.t[-1]==10
def test_cvode_AN_tacoma():
	"""
	Tacoma bridge CVode Adams Newton
	See Exercise 1 at Sauer page 330
	"""
	cv=CVode(problem,[0,0,0.001,0])
	cv.discr='Adams'
	cv.iter='FixedPoint'
	# we use default atols, rtols
	print cv.discr,cv.iter,cv.atol,cv.rtol,'\n'
	cv.simulate(10,200)
	#cv.plot()
	assert abs(cv.y[-1][0]- (-0.3615)) < 2.e-3
	assert cv.t[-1]==10
def test_rk_tacoma():
	"""
	Tacoma bridge RK4
	See Exercise 1 at Sauer page 330
	"""
	rk4=RungeKutta4(problem,[0,0,0.001,0])
	rk4.simulate(10,1000)
	print rk4.y[-1][0],rk4.t
	assert abs(rk4.y[-1][0]- (-0.3615)) < 1.e-3
	assert rk4.t[-1]==10
def test_rk34_tacoma():
	"""
	Tacoma bridge RK34
	See Exercise 1 at Sauer page 330
	"""
	rk34=RungeKutta34(problem,[0,0,0.001,0])
	rk34.simulate(10)
	assert abs(rk34.y[-1][0]- (-0.3615)) < 1.e-3
	assert rk34.t[-1]==10
