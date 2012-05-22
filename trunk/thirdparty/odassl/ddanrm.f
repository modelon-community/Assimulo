      DOUBLE PRECISION FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR)                      
C                                                                               
C***BEGIN PROLOGUE  DDANRM                                                      
C***REFER TO  DDASSL                                                            
C***ROUTINES CALLED  (NONE)                                                     
C***DATE WRITTEN   830315   (YYMMDD)                                            
C***REVISION DATE  830315   (YYMMDD)                                            
C***END PROLOGUE  DDANRM                                                        
C-----------------------------------------------------------------------        
C     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED                               
C     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH                             
C     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS                                 
C     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.                                  
C        DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)                                
C-----------------------------------------------------------------------        
C                                                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      DIMENSION V(NEQ),WT(NEQ)                                                  
      DIMENSION RPAR(1),IPAR(1)                                                 
      DDANRM = 0.0D0                                                            
      VMAX = 0.0D0                                                              
      DO 10 I = 1,NEQ                                                           
10      IF(DABS(V(I)/WT(I)) .GT. VMAX) VMAX = DABS(V(I)/WT(I))                  
      IF(VMAX .LE. 0.0D0) GO TO 30                                              
      SUM = 0.0D0                                                               
      DO 20 I = 1,NEQ                                                           
20      SUM = SUM + ((V(I)/WT(I))/VMAX)**2                                      
      DDANRM = VMAX*DSQRT(SUM/DFLOAT(NEQ))                                      
30    CONTINUE                                                                  
      RETURN                                                                    
C------END OF FUNCTION DDANRM------                                             
      END                                                                       
