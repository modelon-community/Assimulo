      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)                                       
      INTEGER LDA,N,IPVT(1),INFO                                                
      DOUBLE PRECISION A(LDA,1)                                                 
C                                                                               
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.          
C                                                                               
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED                    
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.                  
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .                           
C                                                                               
C     ON ENTRY                                                                  
C                                                                               
C        A       DOUBLE PRECISION(LDA, N)                                       
C                THE MATRIX TO BE FACTORED.                                     
C                                                                               
C        LDA     INTEGER                                                        
C                THE LEADING DIMENSION OF THE ARRAY  A .                        
C                                                                               
C        N       INTEGER                                                        
C                THE ORDER OF THE MATRIX  A .                                   
C                                                                               
C     ON RETURN                                                                 
C                                                                               
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS                 
C                WHICH WERE USED TO OBTAIN IT.                                  
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE               
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER                  
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.               
C                                                                               
C        IPVT    INTEGER(N)                                                     
C                AN INTEGER VECTOR OF PIVOT INDICES.                            
C                                                                               
C        INFO    INTEGER                                                        
C                = 0  NORMAL VALUE.                                             
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR               
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES                
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO          
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE           
C                     INDICATION OF SINGULARITY.                                
C                                                                               
C     LINPACK. THIS VERSION DATED 08/14/78 .                                    
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.              
C                                                                               
C     SUBROUTINES AND FUNCTIONS                                                 
C                                                                               
C     BLAS DAXPY,DSCAL,IDAMAX                                                   
C                                                                               
C     INTERNAL VARIABLES                                                        
C                                                                               
      DOUBLE PRECISION T                                                        
      INTEGER IDAMAX,J,K,KP1,L,NM1                                              
C                                                                               
C                                                                               
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                                
C                                                                               
      INFO = 0                                                                  
      NM1 = N - 1                                                               
      IF (NM1 .LT. 1) GO TO 70                                                  
      DO 60 K = 1, NM1                                                          
         KP1 = K + 1                                                            
C                                                                               
C        FIND L = PIVOT INDEX                                                   
C                                                                               
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1                                     
         IPVT(K) = L                                                            
C                                                                               
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED                  
C                                                                               
         IF (A(L,K) .EQ. 0.0D0) GO TO 40                                        
C                                                                               
C           INTERCHANGE IF NECESSARY                                            
C                                                                               
            IF (L .EQ. K) GO TO 10                                              
               T = A(L,K)                                                       
               A(L,K) = A(K,K)                                                  
               A(K,K) = T                                                       
   10       CONTINUE                                                            
C                                                                               
C           COMPUTE MULTIPLIERS                                                 
C                                                                               
            T = -1.0D0/A(K,K)                                                   
            CALL DSCAL(N-K,T,A(K+1,K),1)                                        
C                                                                               
C           ROW ELIMINATION WITH COLUMN INDEXING                                
C                                                                               
            DO 30 J = KP1, N                                                    
               T = A(L,J)                                                       
               IF (L .EQ. K) GO TO 20                                           
                  A(L,J) = A(K,J)                                               
                  A(K,J) = T                                                    
   20          CONTINUE                                                         
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)                          
   30       CONTINUE                                                            
         GO TO 50                                                               
   40    CONTINUE                                                               
            INFO = K                                                            
   50    CONTINUE                                                               
   60 CONTINUE                                                                  
   70 CONTINUE                                                                  
      IPVT(N) = N                                                               
      IF (A(N,N) .EQ. 0.0D0) INFO = N                                           
      RETURN                                                                    
      END                                                                       
