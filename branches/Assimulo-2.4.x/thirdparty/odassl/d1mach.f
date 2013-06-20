      DOUBLE PRECISION FUNCTION D1MACH (IDUM)                                   
      INTEGER IDUM                                                              
C-----------------------------------------------------------------------        
C THIS ROUTINE COMPUTES THE UNIT ROUNDOFF OF THE MACHINE IN DOUBLE              
C PRECISION.  THIS IS DEFINED AS THE SMALLEST POSITIVE MACHINE NUMBER           
C U SUCH THAT  1.0E0 + U .NE. 1.0E0 (IN SINGLE PRECISION).                      
C-----------------------------------------------------------------------        
      DOUBLE PRECISION U, COMP                                                  
C     U = 1.0E0                                                                 
C10   U = U*0.5E0                                                               
C     COMP = 1.0E0 + U                                                          
C     IF (COMP .NE. 1.0E0) GO TO 10                                             
C     D1MACH = U*2.0E0 
C              FUER APPOLLO
      D1MACH = 2.2E-16
      RETURN                                                                    
      END                                                                       
