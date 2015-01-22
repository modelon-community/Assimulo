      INTEGER FUNCTION IDAMAX(N,DX,INCX)                                        
C                                                                               
C     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.                    
C     JACK DONGARRA, LINPACK, 3/11/78.                                          
C                                                                               
      DOUBLE PRECISION DX(1),DMAX                                               
      INTEGER I,INCX,IX,N                                                       
C                                                                               
      IDAMAX = 0                                                                
      IF( N .LT. 1 ) RETURN                                                     
      IDAMAX = 1                                                                
      IF(N.EQ.1)RETURN                                                          
      IF(INCX.EQ.1)GO TO 20                                                     
C                                                                               
C        CODE FOR INCREMENT NOT EQUAL TO 1                                      
C                                                                               
      IX = 1                                                                    
      DMAX = DABS(DX(1))                                                        
      IX = IX + INCX                                                            
      DO 10 I = 2,N                                                             
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5                                       
         IDAMAX = I                                                             
         DMAX = DABS(DX(IX))                                                    
    5    IX = IX + INCX                                                         
   10 CONTINUE                                                                  
      RETURN                                                                    
C                                                                               
C        CODE FOR INCREMENT EQUAL TO 1                                          
C                                                                               
   20 DMAX = DABS(DX(1))                                                        
      DO 30 I = 2,N                                                             
         IF(DABS(DX(I)).LE.DMAX) GO TO 30                                       
         IDAMAX = I                                                             
         DMAX = DABS(DX(I))                                                     
   30 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
