      SUBROUTINE DDAWTS(NEQ,IWT,RTOL,ATOL,Y,WT,RPAR,IPAR)                       
C                                                                               
C***BEGIN PROLOGUE  DDAWTS                                                      
C***REFER TO  DDASSL                                                            
C***ROUTINES CALLED  (NONE)                                                     
C***DATE WRITTEN   830315   (YYMMDD)                                            
C***REVISION DATE  830315   (YYMMDD)                                            
C***END PROLOGUE  DDAWTS                                                        
C-----------------------------------------------------------------------        
C     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR                              
C     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),                          
C     I=1,-,N.                                                                  
C     RTOL AND ATOL ARE SCALARS IF IWT = 0,                                     
C     AND VECTORS IF IWT = 1.                                                   
C-----------------------------------------------------------------------        
C                                                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      DIMENSION RTOL(1),ATOL(1),Y(1),WT(1)                                      
      DIMENSION RPAR(1),IPAR(1)                                                 
      RTOLI=RTOL(1)                                                             
      ATOLI=ATOL(1)                                                             
      DO 20 I=1,NEQ                                                             
         IF (IWT .EQ.0) GO TO 10                                                
           RTOLI=RTOL(I)                                                        
           ATOLI=ATOL(I)                                                        
10         WT(I)=RTOLI*DABS(Y(I))+ATOLI                                         
20         CONTINUE                                                             
      RETURN                                                                    
C-----------END OF SUBROUTINE DDAWTS-------------------------------------       
      END                                                                       
