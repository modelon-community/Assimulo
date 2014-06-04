      SUBROUTINE DDATRP(X,XOUT,YOUT,YPOUT,NEQ,KOLD,PHI,PSI)                     
C                                                                               
C***BEGIN PROLOGUE  DDATRP                                                      
C***REFER TO  DDASSL                                                            
C***ROUTINES CALLED  (NONE)                                                     
C***DATE WRITTEN   830315   (YYMMDD)                                            
C***REVISION DATE  830315   (YYMMDD)                                            
C***END PROLOGUE  DDATRP                                                        
C                                                                               
C-----------------------------------------------------------------------        
C     THE METHODS IN SUBROUTINE DASTEP USE POLYNOMIALS                          
C     TO APPROXIMATE THE SOLUTION. DDATRP APPROXIMATES THE                      
C     SOLUTION AND ITS DERIVATIVE AT TIME XOUT BY EVALUATING                    
C     ONE OF THESE POLYNOMIALS,AND ITS DERIVATIVE,THERE.                        
C     INFORMATION DEFINING THIS POLYNOMIAL IS PASSED FROM                       
C     DASTEP, SO DDATRP CANNOT BE USED ALONE.                                   
C                                                                               
C     THE PARAMETERS ARE%                                                       
C     X     THE CURRENT TIME IN THE INTEGRATION.                                
C     XOUT  THE TIME AT WHICH THE SOLUTION IS DESIRED                           
C     YOUT  THE INTERPOLATED APPROXIMATION TO Y AT XOUT                         
C           (THIS IS OUTPUT)                                                    
C     YPOUT THE INTERPOLATED APPROXIMATION TO YPRIME AT XOUT                    
C           (THIS IS OUTPUT)                                                    
C     NEQ   NUMBER OF EQUATIONS                                                 
C     KOLD  ORDER USED ON LAST SUCCESSFUL STEP                                  
C     PHI   ARRAY OF SCALED DIVIDED DIFFERENCES OF Y                            
C     PSI   ARRAY OF PAST STEPSIZE HISTORY                                      
C-----------------------------------------------------------------------        
C                                                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      DIMENSION YOUT(1),YPOUT(1)                                                
      DIMENSION PHI(NEQ,1),PSI(1)                                               
      KOLDP1=KOLD+1                                                             
      TEMP1=XOUT-X                                                              
      DO 10 I=1,NEQ                                                             
         YOUT(I)=PHI(I,1)                                                       
10       YPOUT(I)=0.0D0                                                         
      C=1.0D0                                                                   
      D=0.0D0                                                                   
      GAMMA=TEMP1/PSI(1)                                                        
      DO 30 J=2,KOLDP1                                                          
         D=D*GAMMA+C/PSI(J-1)                                                   
         C=C*GAMMA                                                              
         GAMMA=(TEMP1+PSI(J-1))/PSI(J)                                          
         DO 20 I=1,NEQ                                                          
            YOUT(I)=YOUT(I)+C*PHI(I,J)                                          
20          YPOUT(I)=YPOUT(I)+D*PHI(I,J)                                        
30       CONTINUE                                                               
      RETURN                                                                    
C                                                                               
C------END OF SUBROUTINE DDATRP------                                           
      END                                                                       
