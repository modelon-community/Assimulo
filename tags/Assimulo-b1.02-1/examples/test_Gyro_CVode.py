from Assimulo.Problem import Explicit_Problem
from Assimulo.Explicit_ODE import *
from scipy import *

def curl(v):
		return array([[0,v[2],-v[1]],[-v[2],0,v[0]],[v[1],-v[0],0]])

def f(t,u):
	"""
	Simulations for the Gyro (Heavy Top) example in Celledoni/Safstrom: 
	Journal of Physics A, Vol 39, 5463-5478, 2006
	"""
	I1=1000.
	I2=5000.
	I3=6000.
	u0=[0,0,1.]
	pi=u[0:3]
	Q=(u[3:12]).reshape((3,3))
	Qu0=dot(Q,u0)
	f=array([Qu0[1],-Qu0[0],0.])
	f=0
	omega=array([pi[0]/I1,pi[1]/I2,pi[2]/I3])
	pid=dot(curl(omega),pi)+f
	Qd=dot(curl(omega),Q)
	return hstack([pid,Qd.reshape((9,))])

def energi(state):
	energi=[]
	for st in state:
		Q=(st[3:12]).reshape((3,3))
		pi=st[0:3]
		u0=[0,0,1.]
		Qu0=dot(Q,u0)
		V=Qu0[2]  # potential energy
		T=0.5*(pi[0]**2/1000.+pi[1]**2/5000.+pi[2]**2/6000.)
		energi.append([T])
	return energi


if __name__=='__main__':
    gyro_mod = Explicit_Problem()
    gyro_mod.f = f

    y0=hstack([[1000.*10,5000.*10,6000*10],eye(3).reshape((9,))])
    gyro_sim=CVode(gyro_mod,y0)
    
    gyro_sim.discr='BDF'
    gyro_sim.iter='Newton'
    gyro_sim.maxord=2 #Sets the maxorder
    gyro_sim.atol=1.e-10
    gyro_sim.rtol=1.e-10
    gyro_sim.simulate(0.1)
    en=energi(gyro_sim.y)
