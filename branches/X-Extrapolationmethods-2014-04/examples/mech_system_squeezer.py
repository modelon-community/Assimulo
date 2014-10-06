
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

from  __future__  import division
#import nose
import assimulo.problem as ap
import assimulo.special_systems as ass
import numpy as N
import scipy.linalg as sl
from assimulo.solvers import IDA, ODASSL, Mexax
from pylab import figure
from scipy import *

def squeezer():
    
    n_p=7
    n_la=6
    m1,m2,m3,m4,m5,m6,m7=.04325,.00365,.02373,.00706,.07050,.00706,.05498
    i1,i2,i3,i4,i5,i6,i7=2.194e-6,4.410e-7,5.255e-6,5.667e-7,1.169e-5,5.667e-7,1.912e-5
    # Geometry
    xa,ya=-.06934,-.00227
    xb,yb=-0.03635,.03273
    xc,yc=.014,.072
    d,da,e,ea=28.e-3,115.e-4,2.e-2,1421.e-5
    rr,ra=7.e-3,92.e-5
    ss,sa,sb,sc,sd=35.e-3,1874.e-5,1043.e-5,18.e-3,2.e-2
    ta,tb=2308.e-5,916.e-5
    u,ua,ub=4.e-2,1228.e-5,449.e-5
    zf,zt=2.e-2,4.e-2
    fa=1421.e-5
    # Driving torque
    mom=0.033
    # Spring data
    c0=4530.
    lo=0.07785
        
    def forces(t, p, v):
        beta,theta,gamma,phi,delta,omega,epsilon=p
        bep,thp,gap,php,dep,omp,epp=v
        sibe,sith,siga,siph,side,siom,siep=sin(p)
        cobe,coth,coga,coph,code,coom,coep=cos(p)
        
        xd = sd*coga + sc*siga + xb
        yd = sd*siga - sc*coga + yb
        lang  = sqrt ((xd-xc)**2 + (yd-yc)**2)
        force = - c0 * (lang - lo)/lang
        fx = force * (xd-xc)
        fy = force * (yd-yc)
        ff=array([
            mom - m2*da*rr*thp*(thp+2*bep)*sith,    
            m2*da*rr*bep**2*sith,
            fx*(sc*coga - sd*siga) + fy*(sd*coga + sc*siga),
            m4*zt*(e-ea)*dep**2*coph,
            - m4*zt*(e-ea)*php*(php+2*dep)*coph,
            - m6*u*(zf-fa)*epp**2*coom,
            m6*u*(zf-fa)*omp*(omp+2*epp)*coom])        
        return ff
        
    def GT(p):
        beta,theta,gamma,phi,delta,omega,epsilon=p
        sibe,sith,siga,siph,side,siom,siep=sin(p)
        cobe,coth,coga,coph,code,coom,coep=cos(p)
        gp=zeros((6,7))
        
        sibeth = sin(beta+theta);cobeth = cos(beta+theta)
        siphde = sin(phi+delta);cophde = cos(phi+delta)
        siomep = sin(omega+epsilon);coomep = cos(omega+epsilon)

    
        gp[0,0] = - rr*sibe + d*sibeth
        gp[0,1] = d*sibeth
        gp[0,2] = - ss*coga
        gp[1,0] = rr*cobe - d*cobeth
        gp[1,1] = - d*cobeth
        gp[1,2] = - ss*siga
        gp[2,0] = - rr*sibe + d*sibeth
        gp[2,1] = d*sibeth
        gp[2,3] = - e*cophde
        gp[2,4] = - e*cophde + zt*side
        gp[3,0] = rr*cobe - d*cobeth
        gp[3,1] = - d*cobeth
        gp[3,3] = - e*siphde
        gp[3,4] = - e*siphde - zt*code
        gp[4,0] = - rr*sibe + d*sibeth
        gp[4,1] = d*sibeth
        gp[4,5] = zf*siomep
        gp[4,6] = zf*siomep - u*coep
        gp[5,0] = rr*cobe - d*cobeth
        gp[5,1] = - d*cobeth
        gp[5,5] = - zf*coomep
        gp[5,6] = - zf*coomep - u*siep
        return N.array(gp).T
        
    def constr3(t,y):
        p=y[0:2]
        return N.array([p[0]**2+p[1]**2-1.])
    def constr2(t,y):
        g=dot(GT,y[7:14])
        return N.array(g)
    def constr1(t,y):
        #p,v,la=y[0:2],y[2:4],y[4:5]
        return N.array([v[0]**2+v[1]**2 - la[0] * (p[0]**2 + p[1]**2) - p[1] * g])
    def mass_matrix(t,p):
        sibe,sith,siga,siph,side,siom,siep=sin(p)
        cobe,coth,coga,coph,code,coom,coep=cos(p)
               
        m=zeros((7,7))
        m[0,0] = m1*ra**2 + m2*(rr**2-2*da*rr*coth+da**2) + i1 + i2
        m[1,0] = m[0,1] = m2*(da**2-da*rr*coth) + i2
        m[1,1] = m2*da**2 + i2
        m[2,2] = m3*(sa**2+sb**2) + i3
        m[3,3] = m4*(e-ea)**2 + i4
        m[4,3] = m[3,4] = m4*((e-ea)**2+zt*(e-ea)*siph) + i4
        m[4,4] = m4*(zt**2+2*zt*(e-ea)*siph+(e-ea)**2) + m5*(ta**2+tb**2)+ i4 + i5
        m[5,5] = m6*(zf-fa)**2 + i6
        m[6,5] = m[5,6] = m6*((zf-fa)**2-u*(zf-fa)*siom) + i6
        m[6,6] = m6*((zf-fa)**2-2*u*(zf-fa)*siom+u**2) + m7*(ua**2+ub**2)+ i6 + i7
        return N.array(m)
    return ass.Mechanical_System(n_p, forces, n_la,  
                                   [-0.0617138900142764496358948458001,  #  beta
                                    0.,                                 #  theta
                                    0.455279819163070380255912382449,   # gamma
                                    0.222668390165885884674473185609,   # phi
                                    0.487364979543842550225598953530,   # delta
                                    -0.222668390165885884674473185609,  # Omega
                                    1.23054744454982119249735015568],   #epsilon 
                                   zeros(7,),
                                   [98.5668703962410896057654982170,        # lambda[0]
                                    -6.12268834425566265503114393122]+4*[0.],       # lambda[1:7]
                                   zeros(7,),
                                   [14222.4439199541138705911625887,        #  betadotdot
                                    -10666.8329399655854029433719415,       #  Thetadotdot
                                    0.,0.,0.,0.,0.],
                                   GT=GT,
                                   constr3 = constr3,
                                   constr2 = constr2, 
                                   constr1 = constr1,
                                   mass_matrix=mass_matrix)
                                   
                                   
                                   
def run_example(index, with_plots=True, with_test=False):
    my_sque_sys=squeezer()
    my_sque=my_sque_sys.generate_problem(index)
    my_sque.name='Index {} System'.format(index)
    if index in ('ind1','ind2','ind3','ggl2'):
        dae_sque = IDA(my_sque)  
    elif index in ('ovstab2','ovstab1'):
        dae_sque = ODASSL(my_sque) 
    elif index in ('oproj2') :
        dae_sque = Mexax(my_sque) 
    else:
        raise Exception('No method for index {}'.format(index))
    dae_sque.atol=1.e-4
    dae_sque.rtol=1.e-4
    dae_sque.suppress_alg=True  
    t,y,yd=dae_sque.simulate(10.,1)   
    if index != 'oproj2':
        final_residual=my_sque.res(0.,dae_sque.y,dae_sque.yd)
        whichres='All'
    else:
        qflag=9*[0]
        qflag[6]=1
        parameterlist={'nl':my_sque.n_la,'ng':0,'nu':0,'t':dae_sque.t,
                                  'p':dae_sque.y[:2],'v':dae_sque.y[2:4],'u':N.array([0.]),
                                  'lam':dae_sque.y[4],'mass':N.empty((2,2)),'gp':N.empty((1,2)),
                                  'f':N.empty((2,)),'pdot':dae_sque.yd[:2],'udot':N.empty((1,)),
                                  'g':N.empty((1,)),'gi':N.empty((1,)),'fl':N.empty((1,)),
                                  'qflag':qflag}
        final_residual=my_sque.fprob(**parameterlist)[5]      # only position residual
        whichres='Position'
    print(my_sque.name+" {} residuals after the integration run\n".format(whichres))
    print final_residual, 'Norm:  ', sl.norm(final_residual) 
    if with_test:
        assert(sl.norm(final_residual) < 1.5e-2)
    if with_plots:
        dae_sque.plot(mask=[1,1]+(len(my_sque.y0)-2)*[0]) 
    return my_sque, dae_sque
        
if __name__=='__main__':
    index_values=[['ind1','ind2','ind3','ggl2','ovstab2','ovstab1','oproj2'][-1]]
    sim={}
    mod={}
    for ind in index_values:
        mod[ind], sim[ind]=run_example(index=ind)
    
    
                              
                                   
    
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
