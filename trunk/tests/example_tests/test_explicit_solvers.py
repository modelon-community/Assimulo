import nose
import numpy as N
from scipy import *
from Assimulo.Explicit_ODE import *
from Assimulo.Problem import Explicit_Problem

def run_basic_test():
    
    def f(t,y):
        ydot = -y[0]
        return N.array([ydot])
    global exp_mod
    global exp_sim
    exp_mod = Explicit_Problem()
    exp_mod.f = f
    exp_mod.y0 = 4.0
    exp_mod.problem_name = 'Simple Explicit Example'

    #Explicit Euler
    exp_sim = Explicit_Euler(exp_mod) #Create a explicit Euler solver
    exp_sim.simulate(3,100) #Simulate 3 seconds
    exp_sim.simulate(5,100) #Simulate 2 second more
    #exp_sim.plot() #Plot the solution
    
    assert exp_sim.problem_name == exp_mod.problem_name
    assert exp_sim.y[0] == 4.0
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.0252255, places=4)
    
    #RungeKutta4
    exp_sim = RungeKutta4(exp_mod)
    exp_sim.simulate(5, 100) #Simulate 5 seconds
    #exp_sim.plot()
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.02697622, places=4)
    
    #RungeKutta34
    exp_sim = RungeKutta34(exp_mod)
    exp_sim.simulate(5) #Simulate 5 seconds
    #exp_sim.plot()
    print exp_sim.t[-1]
    assert len(exp_sim.t) == 62
    assert exp_sim.t[-1] == 5.0
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.02697622, places=4)
    exp_sim.reset()
    exp_sim.initstep = 0.1
    exp_sim.simulate(5) #Simulate 5 seconds
    #exp_sim.plot()
    print exp_sim.t[-1]
    assert exp_sim.t[-1] == 5.0
    assert len(exp_sim.t) == 60
    nose.tools.assert_almost_equal(exp_sim.y[-1], 0.02697622, places=4)

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

modTacoma = Explicit_Problem()
modTacoma.f = bridge
modTacoma.problem_name = bridge.name

def test_cvode_BN_tacoma():
	"""
	Tacoma bridge CVode BDF Newton
	See Exercise 1 at Sauer page 330
	"""
	cv=CVode(modTacoma,[0,0,0.001,0])
	cv.discr='BDF'
	cv.iter='Newton'
	cv.atol=1.e-6
	cv.rtol=1.e-9
	cv.simulate(10,200)
	assert abs(cv.y[-1][0]- (-0.3615)) < 2.e-3
	assert cv.t[-1]==10
    
def test_cvode_BF_tacoma():
	"""
	Tacoma bridge CVode BDF Fixed Point
	See Exercise 1 at Sauer page 330
	"""
	cv=CVode(modTacoma,[0,0,0.001,0])
	cv.discr='BDF'
	cv.atol=1.e-4
	cv.rtol=1.e-4
	cv.simulate(10,200)
	assert abs(cv.y[-1][0]- (-0.3615)) < 2.e-3
	assert cv.t[-1]==10
    
def test_cvode_AF_tacoma():
	"""
	Tacoma bridge CVode Adams Fixed Point
	See Exercise 1 at Sauer page 330
	"""
	cv=CVode(modTacoma,[0,0,0.001,0])
	cv.rtol=1.e-4
	cv.simulate(10,200)
	assert abs(cv.y[-1][0]- (-0.3615)) < 2.e-3
	assert cv.t[-1]==10

def test_cvode_AN_tacoma():
	"""
	Tacoma bridge CVode Adams Newton
	See Exercise 1 at Sauer page 330
	"""
	cv=CVode(modTacoma,[0,0,0.001,0])
	cv.discr='Adams'
	cv.iter='FixedPoint'
	cv.simulate(10,200)
	assert abs(cv.y[-1][0]- (-0.3615)) < 2.e-3
	assert cv.t[-1]==10

def test_rk_tacoma():
	"""
	Tacoma bridge RK4
	See Exercise 1 at Sauer page 330
	"""
	rk4=RungeKutta4(modTacoma,[0,0,0.001,0])
	rk4.simulate(10,1000)
	assert abs(rk4.y[-1][0]- (-0.3615)) < 1.e-3
	assert rk4.t[-1]==10

def test_rk34_tacoma():
	"""
	Tacoma bridge RK34
	See Exercise 1 at Sauer page 330
	"""
	rk34=RungeKutta34(modTacoma,[0,0,0.001,0])
	rk34.simulate(10)
	assert abs(rk34.y[-1][0]- (-0.3615)) < 1.e-3
	assert rk34.t[-1]==10

def test_explicit():
    run_basic_test()

if __name__=='__main__':
    run_example()
