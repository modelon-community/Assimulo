      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)                         
C                                                                               
C     FORMS THE DOT PRODUCT OF TWO VECTORS.                                     
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.                          
C     JACK DONGARRA, LINPACK, 3/11/78.                                          
C                                                                               
      DOUBLE PRECISION DX(1),DY(1),DTEMP                                        
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N                                         
C                                                                               
      DDOT = 0.0D0                                                              
      DTEMP = 0.0D0                                                             
      IF(N.LE.0)RETURN                                                          
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20                                       
C                                                                               
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS                        
C          NOT EQUAL TO 1                                                       
C                                                                               
      IX = 1                                                                    
      IY = 1                                                                    
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1                                         
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1                                         
      DO 10 I = 1,N                                                             
        DTEMP = DTEMP + DX(IX)*DY(IY)                                           
        IX = IX + INCX                                                          
        IY = IY + INCY                                                          
   10 CONTINUE                                                                  
      DDOT = DTEMP                                                              
      RETURN                                                                    
C                                                                               
C        CODE FOR BOTH INCREMENTS EQUAL TO 1                                    
C                                                                               
C                                                                               
C        CLEAN-UP LOOP                                                          
C                                                                               
   20 M = MOD(N,5)                                                              
      IF( M .EQ. 0 ) GO TO 40                                                   
      DO 30 I = 1,M                                                             
        DTEMP = DTEMP + DX(I)*DY(I)                                             
   30 CONTINUE                                                                  
      IF( N .LT. 5 ) GO TO 60                                                   
   40 MP1 = M + 1                                                               
      DO 50 I = MP1,N,5                                                         
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +                     
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)        
   50 CONTINUE                                                                  
   60 DDOT = DTEMP                                                              
      RETURN                                                                    
      END                                                                       
