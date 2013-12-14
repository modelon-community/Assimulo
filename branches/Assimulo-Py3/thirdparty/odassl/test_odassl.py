from scipy import *
import odassl as od

g=13.7503716373294544 # gravitation constant such that the period is 2

def tolerance_set(irun,tol):
    rtol=atol=tol*ones((5,))
    if irun == 1 or irun == 3:
       rtol[4]=atol[4]=1.e7
    elif irun == 0:                                                                
       rtol[2:]=atol[2:]=1.e7                                                    
    info=1
    return rtol,atol,info            
def res1(t,p,pp):
    ires=0                                                                                                                      
    delta=empty((5,))                                                      
    # kinematics                                                                               
    delta[0:2] = [pp[0] - p[2], pp[1] - p[3]]                                                                         
    # dynamics                                                  
    delta[2:4] =  [pp[2] + p[4]*p[0], pp[3]+p[4]*p[1]+g]                                             
    # lambda constraint                                                              
    delta[4] = p[2]**2 + p[3]**2 - (p[0]**2 + p[1]**2)*p[4] - p[1]*g                                          
    return delta, ires
def res2(t,p,pp):
    ires=0                                                                      
    delta=empty((5,))                                                      
    # kinematics                                                                               
    delta[0:2] = [pp[0] - p[2], pp[1] - p[3]]                                                                         
    # dynamics                                                  
    delta[2:4] =  [pp[2] + p[4]*p[0], pp[3]+p[4]*p[1]+g]                                             
    # velocity constraint                                                              
    delta[4] = p[0]*p[2] + p[1]*p[3]                                          
    return delta, ires
def res3(t,p,pp):
    ires=0                                                                      
    delta=empty((5,))   
    #print t,p,pp                                                   
    # kinematics                                                                               
    delta[0:2] = [pp[0] - p[2], pp[1] - p[3]]                                                                         
    # dynamics                                                  
    delta[2:4] =  [pp[2] + p[4]*p[0], pp[3]+p[4]*p[1]+g]                                             
    # position constraint                                                              
    delta[4] = p[0]**2 + p[1]**2 - 1.0                                           
    return delta, ires
def resStab2(t,p,pp):
    ires=0                                                                      
    delta=empty((neq,))   
    #print t,p,pp                                                   
    # kinematics                                                                               
    delta[0:2] = [pp[0] - p[2], pp[1] - p[3]]                                                                         
    # dynamics                                                  
    delta[2:4] =  [pp[2] + p[4]*p[0], pp[3]+p[4]*p[1]+g]
    # velocity constraint  
    delta[4] = p[0]*p[2] + p[1]*p[3]                                            
    # position constraint                                                              
    delta[5] = p[0]**2 + p[1]**2 - 1.0                                           
    return delta, ires    
def resStab1(t,p,pp):
    ires=0                                                                      
    delta=empty((neq,))   
    #print t,p,pp                                                   
    # kinematics                                                                               
    delta[0:2] = [pp[0] - p[2], pp[1] - p[3]]                                                                         
    # dynamics                                                  
    delta[2:4] =  [pp[2] + p[4]*p[0], pp[3]+p[4]*p[1]+g]
    # lambda constraint
    delta[4] = p[2]**2 + p[3]**2 - (p[0]**2 + p[1]**2)*p[4] - p[1]*g
    # velocity constraint  
    delta[6] = p[0]*p[2] + p[1]*p[3]                                            
    # position constraint                                                              
    delta[5] = p[0]**2 + p[1]**2 - 1.0                                           
    return delta, ires    

p=empty(5)
pp=empty(5)
rwork=empty((1000,))
iwork=empty((200,),dtype=int32)  

                                                  
ny = 5                                                                    

tol= 1.e-5       

res=[res3,res2,res1, resStab2, resStab1]                                                 

for irun, task in enumerate(['INDEX-3', 'INDEX-2',
             'INDEX-1', 
             'STAB. INDEX-2','STAB. INDEX-1']):
     
     info=zeros((15,),dtype=int32)                                                                                         
     # initial conditions                                                                                                                                               
     p[0] = 1.0                                                             
     p[1] = -sqrt(1.0 - p[0]**2)                                           
     p[2:] = 0.0                                                                                                                         
     pp[0:2] = p[2:4]                                                           
     pp[2:5] = [-p[4]*p[0], -g-p[4]*p[1],0.]                                                                                                                        
     # other initialisations                                                                                                                                                                        
     t,th = 0., 2.0                                                            
                                                             
     rtol,atol,info[1] = tolerance_set(irun,tol)                                                                             
     ny=len(p)
     if irun <= 2:                                                                    
        neq=ny                                 
     elif irun==3:  # stab index 2                                                                   
        neq=ny+1                                        
     elif irun==4:  # stab index 1
        neq=ny+2                                                                 

#          INTEGRATION LOOP                                                     
     for i in range(50):
         tout=t+th                        
         t,y,yprime,tout,info,idid,rwork,iwork = od.odassl(res[irun],neq,ny,t,p,pp,tout,info,rtol,atol,
                                                         rwork,iwork,res[irun])                                                                                
         if idid<0:                              
            print task+':  ', idid                                                  
            break                                                         
         else:                                                                
            print '{0}  : t= {1: >8.2f}  p={2[0]: >-12.8f} {2[1]: > 12.8f} {2[2]: >-12.8f} {2[3]: >-12.8f} {2[4]: >-12.8f}'.format(task, t, p)                            
     print iwork[10:15], irun  
