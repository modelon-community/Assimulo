C                                                                               
C     DEMONSTRATION - PROGRAM FOR ODASSL  AND SOME SHORT COMMENTS                                      
C     ================================== =========================                                       
C                                                                               
C----------------------------------------------------------------------         
C      ODASSL IS A MODIFIED VERSION OF THE DAE-SOLVER DASSL. IT HAS             
C      BEEN DESIGNED TO INTEGRATE MECHANICAL SYSTEMS IN THEIR                   
C      DESCRIPTORFORM.                                                          
C      THE EQUATIONS OF MOTION OF MECHANICAL SYSTEMS CAN BE GIVEN IN            
C      VARIOUS WAYS. HERE FOR THE PENDULUM EXAMPLE THE EQUATIONS OF             
C      MOTION ARE INTEGRATED USING THE CONSTRAINT EQUATION ON                   
C       -- POSITION LEVEL  (P)  INDEX-3 CASE                                         
C       -- VELOCITY LEVEL  (V)  INDEX-2 CASE                                        
C       -- ACCELERATION (A)/  OR CONSTRAINT FORCE LEVEL (L)  INDEX-1 CASE                 
C      OR USING                                                                 
C       -- THE BAUMGARTE-APPROACH WITH THE STABILIZING               
C          COEFFICIENTS BEEING STORED IN RPAR(1),RPAR(2)                        
C       -- ALL THREE TYPES OF CONSTRAINT EQUATIONS AND SOLVING                  
C          THE RESULTING OVERDETERMINED SYSTEM BY                               
C          THE SSF-ITERATION (SEE LITERATURE). FOR THIS IT IS IMPORTANT
C           THAT THE LAST EQUATIONS ARE ORDERED IN THE FOLLOWING SENSE:
C           D,K,L,P,V OR K,D,L,P,V 
C           WHERE THE FOLLOWING ABBREVIATIONS ARE USED 
C                   D:   DYNAMICS    (RELATION BETWEEN MASS, ACCEL. AND FORCES)
C                   K:   KINEMTATICS (RELATION BETWEEN THE DERIVATIVE OF
C                                     THE POSITION AND THE VELOCITIES)
C                   P:   POSTITON CONSTRAINT (INDEX-3 CONSTRAINT)
C                   V:   VELOCITY CONSTRAINT
C                   L:   LAMBDA CONSTRAINT  
C           VARIOUS EXPERIMENTS HAVE
C           SHOWN THAT THE FIRST OF THESE TWO VARIANTS ARE MORE EFFECTIVE.   
C          THIS FORMULATION IS CALLED THE STABILIZED INDEX-1 APPROACH            
C       -- TWO TYPES OF CONSTRAINT EQUATIONS IN THE FOLLOWING ORDER
C             K,D,V,P  OR EQUIVALENTLY  D,K,V,P                                                                         
C          THIS FORMULATION IS CALLED THE STABILIZED INDEX-2 APPROACH.
C
C     IMPORTANT:
C     ==========  
C     IN THE CURRENT VERSION THE ERROR ESTIMATION IN ODASSL IS
C     BASED ON THE LOCAL TRUNCATION ERROR ONLY. THIS MEANS THAT FOR HIGHER INDEX
C     SYSTEMS SPECIAL ACTIONS MUST BE TAKEN TO EXCLUDE CERTAIN VARIABLES FROM
C     THE ERROR ESTIMATING PROCESS:
C            FOR INDEX-3 PROBLEMS THE VELOCITY AND LAMBDA-VARIABLES MUST
C                                                  BE EXCLUDED
C            FOR INDEX-2 PROBLEMS THE (ALGEBRAIC) LAMBDA VARIABLES MUST BE
C                                                  EXCLUDED
C     THIS IS DONE IN THE EXAMPLE BELOW BY JUST SETTING THE CORRESPONDING
C     COMPONENTS OF THE ERROR-VECTOR TO A HIGH NUMBER.
C                                                                               
      IMPLICIT REAL*8(A-H,O-Z) 
      EXTERNAL RES1,RES2,RES3                                                   
      DIMENSION P(5), PP(5), INFO(15), RWORK(1000), IWORK(200),                   
     *          RTOL(5), ATOL(5), IPAR(10), RPAR(10)
      DATA     INFO/15*0/                                                       
      NY = 5                                                                    
C                GRAVITY - CONSTANT (WITH THIS VALUE THE EXACT                  
C                SOLUTION HAS A PERIOD OF 2 SECS.)                              
      RPAR(4) = 13.7503716373294544D0                                           
      G       = RPAR(4) 
      TOL = 1.E-5                                                        
C                                                                               
      DO 100 IRUN = 1,6 
         IF (IRUN .EQ. 1) THEN 
            WRITE(6,*) 'INDEX-3  FORMULATION'
         ELSE IF (IRUN .EQ. 2) THEN
            WRITE(6,*) 'INDEX-2  FORMULATION'    
         ELSE IF (IRUN .EQ. 3) THEN
            WRITE(6,*) 'INDEX-1  FORMULATION'
         ELSE IF (IRUN .EQ. 4) THEN
            WRITE(6,*) 'BAUMGARTE FORMULATION'
         ELSE IF (IRUN .EQ. 5) THEN
            WRITE(6,*) 'STAB. INDEX-1  FORMULATION'
         ELSE IF (IRUN .EQ. 6) THEN
            WRITE(6,*) 'STAB. INDEX-2  FORMULATION'
         END IF
C
         INFO(1)  = 0                                                           
C                                                                               
C                INITIAL CONDITIONS                                             
C                                                                               
         P(1) = 1.0                                                             
         P(2) = -SQRT(1.D0 - P(1)**2)                                           
         P(3) = 0.0                                                             
         P(4) = 0.0                                                             
         P(5) = 0.0                                                             
         PP(1) = P(3)                                                           
         PP(2) = P(4)                                                           
         PP(3) = 0.                                                             
         PP(4) = -G - P(5)*P(2)                                                 
         PP(5) = 0.                                                             
C                                                                               
C          BAUMGARTE PARAMETER OR CONSTRAINT EQUATION SELECTION PARAM.          
C                                                                               
         IF (IRUN .EQ. 1) THEN                                                  
C        POSITION CONSTRAINT                                                    
            RPAR(1) = 1.0                                                       
            RPAR(2) = 0.                                                        
            RPAR(3) = 0.                                                        
         ELSE IF (IRUN .EQ. 2) THEN                                             
C        VELOCITY CONSTRAINT                                                    
            RPAR(1) = 0.                                                        
            RPAR(2) = 1.                                                        
            RPAR(3) = 0.                                                        
         ELSE IF (IRUN .EQ. 3) THEN                                             
C        ACCELERATION CONSTRAINT                                                
            RPAR(1) = 0.                                                        
            RPAR(2) = 0.                                                        
            RPAR(3) = 1.                                                        
         ELSE IF (IRUN .EQ. 4) THEN                                             
C        BAUMGARTE STABILIZATION                                                
            RPAR(1) = 6.25D+2                                                   
            RPAR(2) = 100.                                                      
            RPAR(3) = 2.0                                                       
         END IF                                                                 
C           SCALING PARAMETERS (NECESSARY FOR ODAES)                            
         IF (IRUN .EQ. 5) THEN                                                  
C           SCALE EQUATION 1 TO 4                                               
            INFO(12) = 0                                                        
            IWORK(1) = 1                                                        
            IWORK(2) = 4                                                        
         END IF                                                                 
C                                                                               
C           OTHER INITIALISATIONS                                               
C                                                                               
         T = 0.0D0                                                              
         TH = 2.0D0                                                             
C                                                                               
C                                                                               
C             SET ERROR TOLERANCES AND FILTER OUT APPROPRIATE                   
C             SOLUTION COMPONENTS FROM THE ERROR TEST BY SETTING                
C             THE CORRESPONDING TOLERANCES TO INFINITY                          
C                                                                               
         DO 20 I=1,5                                                            
           RTOL(I) = TOL                                                      
           ATOL(I) = TOL                                                      
 20      CONTINUE                                                               
         IF (IRUN.EQ.2 .OR. IRUN.EQ.6) THEN                                                    
           RTOL(5) = 1.D+7                                                      
           ATOL(5) = 1.D+7                                                      
           INFO(2) = 1                                                          
         ELSE IF (IRUN.EQ.1) THEN                                               
           DO 30 I=3,5                                                          
              RTOL(I) = 1.D+7                                                   
              ATOL(I) = 1.D+7                                                   
  30       CONTINUE                                                             
           INFO(2) = 1                                                          
         ELSE 
           INFO(2) = 0                                                          
         END IF                                                                 
C                                                                               
C        SET DIMENSION VARIABLES:                                               
C                                                                               
      IF (IRUN .LE. 4 ) THEN                                                
C        DAE - FORMULATION                                                     
         NEQ = 5                                                                
         NY  = 5                                                                
      ELSE IF (IRUN .EQ.5) THEN                                                                     
C        STABILIZED INDEX-1 ODAE - FORMULATION                                                      
         NEQ = 7                                                                
         NY  = 5                                           
      ELSE
C        STABILIZED INDEX-1 ODAE - FORMULATION                                                      
         NEQ = 6                                                                
         NY  = 5                                                                
      END IF                                                                    
C                                                                               
C          INTEGRATION LOOP                                                     
C                             
         DO 200 I=1,50                                                          
            TOUT = T + TH                                                       
            IF (IRUN.LE.4) THEN                                                 
                CALL ODASSL (RES1, NEQ, NY, T, P, PP, TOUT,                        
     *                INFO, RTOL, ATOL,                                         
     *                IDID, RWORK, 370, IWORK,  29, RPAR, IPAR, RES1)           
            ELSE IF (IRUN.EQ.5) THEN                                          
                CALL ODASSL (RES2, NEQ, NY, T, P, PP, TOUT,                        
     *                INFO, RTOL, ATOL,                                         
     *                IDID, RWORK, 370, IWORK,  29, RPAR, IPAR, RES2)      
            ELSE 
                CALL ODASSL (RES3, NEQ, NY, T, P, PP, TOUT,                        
     *                INFO, RTOL, ATOL,                                         
     *                IDID, RWORK, 370, IWORK,  29, RPAR, IPAR, RES3)
            END IF           
         IF (IDID .LT. 0) THEN                                                  
            WRITE(6,1000) IDID                                                  
            GOTO  100                                                           
         END IF                                                                 
C                                                                               
C           OUTPUT OF THE RESULTS                                               
C                                 
         WRITE (6, 1001) T, P                                                  
 200  CONTINUE  
      WRITE (6,1002)  (IWORK(JJ),JJ=11,15),IRUN                               
 100  CONTINUE                                                                  
 999  CONTINUE                                                                  
1000  FORMAT(1H ,35HODASSL TERMINATED WITH ERROR CODE: ,I3)                     
1001  FORMAT(1H ,6(F8.4,3H I ))                                                 
1002  FORMAT(1H ,4HNST= ,I6,5H NFE= ,I6,5H NJE= ,I6,                            
     *       5H NEF= ,I6,5H NCF= ,I6,6H IRUN=,I6)                               
      END                                                                       
      SUBROUTINE RES1(T,P,PP,DELTA,IRES,RPAR,IPAR)                              
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION P(*), PP(*), DELTA(*), RPAR(*), IPAR(*)                         
C                                                                               
      IRES     = 0                                                              
C             GRAVITATION CONSTANT                                            
      G = RPAR(4)                                                               
C                                                                               
C           KINEMATICS                                                          
C                                                                               
      DELTA(1) =  PP(1) - P(3)                                                  
      DELTA(2) =  PP(2) - P(4)                                                  
C                                                                               
C           DYNAMICS                                                            
C                                                                               
      DELTA(3) =  PP(3) + P(5)*P(1)                                             
      DELTA(4) =  PP(4) + P(5)*P(2) + G                                         
C                                                                               
C           POSITION LEVEL CONSTRAINTS                                          
C                                                                               
      PLC = P(1)**2 + P(2)**2 - 1.0D0                                           
C                                                                               
C           VELOCITY LEVEL CONSTRAINTS                                          
C                                                                               
      VLC = P(1)*P(3) + P(2)*P(4)                                               
C                                                                               
C           ACCELERATION LEVEL CONSTRAINTS                                      
C                                                                               
      ALC = P(3)**2 + P(4)**2 + P(1)*PP(3) + P(2)*PP(4)                         
C                                                                               
C           ALGEBRAIC CONDITION (BY RPAR THE BAUMGARTE - COEFFICIENTS           
C           OR THE TYPE OF CONSTRAINT IS SPECIFIED)                             
C                                                                               
      DELTA(5) = RPAR(3) * ALC + RPAR(2) * VLC + RPAR(1) * PLC                  
C                                                                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE RES2(T,P,PP,DELTA,IRES,RPAR,IPAR)                              
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION P(*), PP(*), DELTA(*), RPAR(*), IPAR(*)                         
C                                                                               
      IRES     = 0                                                              
C             GRAVITATIONAL CONSTANT                                            
      G = RPAR(4)                                                               
C                                                                               
C           KINEMATICS                   
C                                                                               
      DELTA(1) =  PP(1) - P(3)                                                  
      DELTA(2) =  PP(2) - P(4)                                                  
C                                                                               
C           DYNAMICS                                                            
C                                                                               
      DELTA(3) =  PP(3) + P(5)*P(1)                                             
      DELTA(4) =  PP(4) + P(5)*P(2) + G  
C                                                                               
C           CONSTRAINT FORCE LEVEL CONSTRAINTS                                      
C                                                                               
      DELTA(5) = P(3)**2 + P(4)**2 - P(1)**2*P(5) - P(2)**2*P(5) -              
     *           P(2)*G 
C                                                                 
C                                                                               
C           POSITION LEVEL CONSTRAINTS                                          
C                                                                               
      DELTA(6) = P(1)**2 + P(2)**2 - 1.0D0                                      
C                                                                               
C           VELOCITY LEVEL CONSTRAINTS                                          
C                                                                               
      DELTA(7) = P(1)*P(3) + P(2)*P(4)                                          
                                                     
      RETURN                                                                    
      END                       
      SUBROUTINE RES3(T,P,PP,DELTA,IRES,RPAR,IPAR)                              
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION P(*), PP(*), DELTA(*), RPAR(*), IPAR(*)                         
C                                                                               
      IRES     = 0                                                              
C             GRAVITATIONAL CONSTANT                                            
      G = RPAR(4)                                                               
C                                                                               
C           KINEMATICS                   
C                                                                               
      DELTA(1) =  PP(1) - P(3)                                                  
      DELTA(2) =  PP(2) - P(4)                                                  
C                                                                               
C           DYNAMICS                                                            
C                                                                               
      DELTA(3) =  PP(3) + P(5)*P(1)                                             
      DELTA(4) =  PP(4) + P(5)*P(2) + G 
C                                                                               
C           VELOCITY LEVEL CONSTRAINTS                                          
C                                                                               
      DELTA(5) = P(1)*P(3) + P(2)*P(4)       
C                                                                               
C           POSITION LEVEL CONSTRAINTS                                          
C                                                                               
      DELTA(6) = P(1)**2 + P(2)**2 - 1.0D0                                      
C                                                                               
      RETURN                                                                    
      END                        
                                          
