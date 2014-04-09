C     ==================================================================
      SUBROUTINE SPAPAT (N,A,IND,IW,R)        
C     ==================================================================
      DIMENSION A(1),IND(1),R(1)      
      LOGICAL IW(1)   
C                     
C    THIS SUBROUTINE COMPUTES THE 'SPARSITY PATTERN" OF THE
C    MATRIX A. THE OUTPUT OF THE ROUTINE IS CONTAINED IN THE
C    VECTOR IND WHICH GIVES INFORMATION ON WHICH COLUMNS IN
C    A THAT DO NOT HAVE NONZEROS IN THE SAME ROW. THESE COLUMNS
C    OF THE JACOBIAN MATRIX CAN THEN BE ESTIMATED SIMULTANEOUSLY
C    BY NUMERICAL DIFFERENCING.      
C                     
C    AUTHOR J OPPELSTRUP 1977        
C    REVISED C SODERLIND 1980-05-29  
C                     
      DO 1  I=1, N  
      IW(1)=.FALSE.   
    1 CONTINUE        
      IPEK=0  
      K=0     
      DO 2 I=1,N     
      IF(IW(I)) GOTO 2 
      IW(I)=.TRUE.   
      K=K+1
      IC=(I-1)*N     
      DO 21 I1=1,N  
      R(I1)=ABS(A(IC+I1))     
   21 CONTINUE        
      IND(K)=I        
      IP1=I+1 
      IF(IP1.GT.N) GOTO 2  
      DO 3 J=IP1,N    
      IF(IW(J)) GOTO 3        
      IK=(J-1)*N    
      DO 31 J1=1,N  
      IF(R(J1).NE.0.0.AND.A(IK+J1).NE.0.0) GOTO 3  
   31 CONTINUE        
C                     
C    COLUMN #J BELONGS TO GROUP #K  
C                     
       K=K+1   
       IND(K)=J        
       IW(J)=.TRUE.    
       DO 32 L=1,N       
       R(L)=AMAX1(R(L),ABS(A(IK+L)))  
   32  CONTINUE        
    3  CONTINUE        
       IPEK=IPEK+1    
       IND(N+IPEK)=K+1 
    2  CONTINUE        
       RETURN  
       END     
C                     
C     * * *   END SPAPAT      * * *   
C 
