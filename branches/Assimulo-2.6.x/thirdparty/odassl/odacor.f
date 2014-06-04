      SUBROUTINE ODACOR(NEQ, NY, M1, X, Y, YPRIME, RPAR, IPAR,             
     *                  JCALC,H, CJ, CJOLD, UROUND, NRE , NJE,                  
     *                  MAXIT, RES, WT, MTYPE, INORM,  ML, MU,                          
     *                  E, DELTA, WM, IWM, IRETC, IRES, JAC)                         
C----------------------------------------------------------------------         
C      MODIFIED VERSION OF THE DASSL-CORRECTOR. THIS VERSION                    
C      HAS EXTENSIONS FOR THE TREATMENT OF OVERDETERMINED SYSTEMS               
C      OF (SINGULARLY) IMPLICIT ODES               
C      (FOR DETAILS SEE REFERENCE).                                  
C      DASSL-AUTHOR:  PETZOLD, LINDA                                            
C             APPLIED MATHEMATICS DIVISION 8331                                 
C             SANDIA NATIONAL LABORATORIES                                      
C             LIVERMORE, CA.    94550                                           
C      MODIFICATIONS BY: FUEHRER, CLAUS                                         
C             DEUTSCHE FORSCHUNGSANSTALT                          
C             FUER LUFT- UND RAUMFAHRT (DLR)                                  
C             INST. DYNAMIC DER FLUGSYSTEME                                     
C             D-8031 WESSLING  (F.R.G)                                          
C                                                                               
C      DATE OF LAST MODIFICATION 900103                                       
C                                                                               
C----------------------------------------------------------------------         
C                                                                               
C***REFER TO  ODACOR                                                            
C***ROUTINES CALLED  RES, ODAJAC, DGESL
C                                                                               
C                                                                               
C     THE PARAMETERS REPRESENT                                                  
C     NEQ --       NUMBER OF EQUATIONS TO BE INTEGRATED                         
C     NY  --       NUMBER OF COMPONENTS OF Y OR YPRIME                          
C     M1 --        SEE COMMENTS TO INFO(13) IN THE HEAD OF ODASSL               
C     X  --        INDEPENDENT VARIABLE                                         
C     Y  --        SOLUTION VECTOR AT X                                         
C     YPRIME --    DERIVATIVE OF SOLUTION VECTOR                                
C                  AFTER SUCCESSFUL STEP                                        
C     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS THAT                       
C                  ARE USED FOR COMMUNICATION BETWEEN THE                       
C                  CALLING PROGRAM AND EXTERNAL USER ROUTINES                   
C                  THEY ARE NOT ALTERED BY DASTEP                               
C     JCALC --     STATE OF THE JACOBIAN:                                       
C                  JCALC =  0  THE JACOBIAN HAS BEEN ACTUALLY COMPUTED          
C                        =  1  JACOBIAN CORRESPONDS TO A PAST VALUE OF          
C                              T, BUT WAS SUFFICIENT FOR THE ITERATION          
C                        = -1  JACOBIAN MUST BE RECOMPUTED                      
C     H --         STEP SIZE                                                    
C     CJ --        LEADING COEFFICIENT FOR THE CURRENT STEP                     
C     CJOLD --     LEADING COEFFICIENT USED FOR THE COMPUTATION                 
C                  OF THE ITERATION MATRIX                                      
C     UROUND --    MACINE PRECISION                                             
C     NRE --       NUMBER OF CALLS TO RES                                       
C     NJE --       NUMBER OF EVALUATIONS OF THE JACOBIAN                        
C     MAXIT --     MAXIMUM NUMBER OF ITERATION STEPS ALLOWED                    
C     RES --       EXTERNAL USER-SUPPLIED SUBROUTINE                            
C                  TO EVALUATE THE RESIDUAL. FOR THE CALL                       
C                  SEE DRIVER - PROGRAMM ODASSL.                                
C     WT --        VECTOR OF WEIGHTS FOR ERROR CRITERION. 
C     MTYPE --     = 1 IF JACOBIAN IS GIVEN BY SUBROUTINE JAC
C                  = 2 IS TO BE COMPUTED BY FINITE DIFFERENCES                        
C     INORM, ML,                                                                
C     MU --        IF INORM .NE. 0 THE ROWS ML UP TO MU OF THE                  
C                  ITERATION MATRIX ARE TO BE SCALED                            
C     IRETC --     RETURN-CODE: = 0  EVERYTHING IS O.K.                         
C                               = 1  NO CONVERGENCE, TRY                        
C                                    AGAIN WITH SMALLER STEPSIZE                
C                               = 2  ITERATION MATRIX IS SINGULAR               
C                               = 3  MORE THAN 3 TIMES IRETC=2  
C     JAC  --      IF MTYPE = 1: SUBROUTINE FOR GENERATING THE JACOBIAN
C                  DUMMY OTHERWISE                
C                                                                               
C     THE OTHER PARAMETERS ARE INFORMATION                                      
C     WHICH IS NEEDED INTERNALLY BY ODACOR                                      
C                                                                               
C-----------------------------------------------------------------------        
C                                                                               
C                                                                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      LOGICAL REDO                                                              
      EXTERNAL RES, JAC                                                              
      DIMENSION Y(NY), YPRIME(NY), DELTA(NEQ), E(NY), WM(*), IWM(*),            
     *          WT(NY), RPAR(*), IPAR(*)                                        
      SAVE  NSF                                                                 
C                                                                               
C           INITIALIZATIONS                                                     
C     
      DO 10 I=1,NY                                                              
            E(I) = 0.0D0                                                        
 10   CONTINUE                                                                  
      IRETC = 0                                                                 
      IRES  = 0                                                                 
      NITER = 0                                                                 
      S     = 100.D0                                                            
      PNORM = DDANRM(NY, Y, WT, RPAR, IPAR)    
C      POINTER OF WM  (TOTAL STORAGE REQUIRED FOR WM:(NY+M1)**2 + NEQ + NEQ
C                                                     =  NEQ**2 + 2*NEQ + NEQ
      LAUX  = 1
      LA1   = LAUX + NEQ 
      LA2   = LA1+NY**2
      LA2T  = LA2 + M1*NY
      LXNORM = LA2T + NY*M1                                                  
      LQ    = LXNORM + NEQ                                                      
C     LRSUM = LQ + M1**2    
C      POINTER OF IWM (TOTAL STORAGE REQUIRED FOR IWM: NEQ = NY+M1
      LIPVA1 = 1
      LIPVQ  = LIPVA1+NY
C     LISUM  = LIPVQ + M1                                                                 
C                                                                               
C     DO WHILE (REDO) ----- ITERATIONLOOP                                       
C                                                                               
 100  CONTINUE                                                                  
C                                                                               
         NRE = NRE + 1                                                          
         CALL RES(X, NEQ, NY, Y, YPRIME,  DELTA, IRES, RPAR, IPAR)                       
         IF (IRES.NE.0) GOTO 999                                                
C                                                                               
         IF (JCALC.EQ.-1 ) THEN                                                 
C                           COMPUTE THE ACTUAL JACOBIAN                         
             JCALC=0                                                            
             NJE = NJE + 1                                                      
             CALL ODAJAC(NEQ, NY,M1,X,Y, YPRIME, RPAR, IPAR,  RES,         
     *                  DELTA, UROUND, CJ, H,                                   
     *                  WT,  WM(LAUX), WM(LA1), WM(LA2), WM(LA2T),                  
     *                  WM(LQ), WM(LXNORM),MTYPE,INORM, ML, MU,                              
     *                  IWM(LIPVA1),IWM(LIPVQ), IRES, IERR, JAC)                        
             CJOLD = CJ                                                         
             IF (IERR.NE. 0 ) THEN                                              
C               MATRIX SINGULAR -- SET COUNTER AND RETURN                       
C                                                                               
                NSF = NSF+1                                                     
                IF (NSF .LT. 3) THEN                                            
                    IRETC = 2                                                   
                ELSE                                                            
                    IRETC = 3                                                   
                END IF                                                          
                GOTO 999                                                        
             END IF                                                             
             NSF = 0                                                            
             IF (IRES .NE. 0 ) GOTO 999                                         
         END IF                                                                 
         IF (CJ .NE. CJOLD) THEN                                                
C                           DAMPING OF THE NEWTON STEP                          
             TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)                                   
             DO 20 I=1,NEQ                                                      
                DELTA(I) = DELTA(I)*TEMP1                                       
   20        CONTINUE                                                           
         END IF                                                                 
         IF (INORM .NE. 1) GOTO 22                                              
         DO 21 I = ML, MU                                                       
            DELTA(I) = DELTA(I)/WM(LXNORM-1+I)                                  
   21    CONTINUE                                                               
   22    CONTINUE                                                               
C
C            BEGIN OF LINEAR ALGEBRA STEP (NOTATION SEE DFVLR FB 89-08)
C
C            1.) A1 INV * A1 
C
         DO 23 I=1,NY
              WM(LAUX+I-1) = DELTA(I)
   23    CONTINUE
C                                          
        IF (M1 .EQ. 0) THEN   
           CALL DGESL(WM(LA1), NY, NY, IWM(LIPVA1), DELTA, 0) 
        ELSE
           CALL DGESL(WM(LA1), NY, NY, IWM(LIPVA1), WM(LAUX), 0) 
C
C            2.) A2 - A2 A1 INV * A1   (NUR FALLS M1 > 0)               
C                                                                               
           DO 30 I=1,M1 
              DO 30 J=1,NY
                 DELTA(NY+I) = DELTA(NY+I) - WM(LA2+(J-1)*M1+I-1)*
     *                                     WM(LAUX + J - 1)
   30      CONTINUE
C
C            3.) Q INV (....)           (NUR FALLS M1 > 0)               

C
           CALL DGESL(WM(LQ),M1,M1,IWM(LIPVQ),DELTA(NY+1),0)
C
C            4.) A1 + A2T*(....)        (NUR FALLS M1 > 0)               

C   
           DO 40 I=1,NY  
              DO 40 J=1,M1
                DELTA(I) = DELTA(I) + WM(LA2+(I-1)*M1 + J-1)*DELTA(NY+J)
   40      CONTINUE
C
C             5.) A1 INV *(...)         (NUR FALLS M1 > 0)               

C               
            CALL DGESL(WM(LA1),NY,NY,IWM(LIPVA1),DELTA,0) 

        END IF
C
C             END OF LINEAR ALGEBRA STEP                                        
C  
C  
C           UPDATING OF THE SOLUTION, DERIVATIVE AND ERROR VECTOR               
         DO 60 I=1,NY                                                           
            Y(I)      = Y(I) -         DELTA(I)  
            E(I)      = E(I) -         DELTA(I) 
            YPRIME(I) = YPRIME(I) - CJ*DELTA(I)                                 
   60    CONTINUE                                                               
C        TEST FOR CONVERGENCE                                                   
         DELNRM = DDANRM(NY,DELTA,WT,RPAR,IPAR)   
         IF (DELNRM .GT. 100.0D0*PNORM*UROUND) THEN                             
            IF (NITER .EQ. 0) THEN                                              
C                       THEN FIRST STEP OF ITERATION                            
               OLDNRM = DELNRM                                                  
               RATE = 0.9D0                                                     
            ELSE                                                                
C                       THEN SUBSEQUENT ITERATION STEP                          
               RATE = (DELNRM/OLDNRM)**(1.0D0/DBLE(NITER))                      
            END IF                                                              
            IF (RATE .GT. 0.9D0 ) THEN                                          
C                        CONVERGENCE RATE TOO BIG                               
               IF (JCALC .EQ. 0) THEN                                           
C                        COMPUTATION HAS BEEN DONE WITH THE ACTUAL              
C                        JACOBIAN - ITERATION TERMINATED WITH                   
C                        RETURNCODE = 1                                         
                   IRETC = 1                                                    
                   REDO = .FALSE.                                               
               ELSE                                                             
C                        TRY AGAIN WITH A NEW JACOBIAN                          
                   JCALC = -1                                                   
                   REDO  = .TRUE.                                               
C                        RECONSTRUCTION OF THE PREDICTED VALUES                 
                   DO 70 I=1,NY                                                 
                      Y(I)      = Y(I)      -    E(I)                           
                      YPRIME(I) = YPRIME(I) - CJ*E(I)                           
                      E(I)      =  0.0D0                                        
   70              CONTINUE                                                     
                   NITER = 0                                                    
                   S = 100.D0                                                   
               END IF                                                           
            ELSE                                                                
               IF (NITER .GT. 0) S = RATE / (1.0D0 - RATE)                      
               IF (S*DELNRM .GT. 0.33D0 ) THEN                                  
C                             ITERATION HAS NOT YET CONVERGED                   
                  NITER = NITER+1                                               
                  IF (NITER .GE. MAXIT) THEN                                    
C                             TOO MANY ITERATION STEPS                          
                     IF (JCALC .EQ. 0) THEN                                     
C                          COMPUTATION HAS BEEN DONE WITH THE ACTUAL            
C                          JACOBIAN - ITERATION TERMINATED WITH                 
C                          RETURNCODE = 1                                       
                        IRETC = 1                                               
                        REDO = .FALSE.                                          
                     ELSE                                                       
C                          TRY AGAIN WITH A NEW JACOBIAN                        
                        JCALC = -1                                              
                        REDO  = .TRUE.                                          
C                          RECONSTRUCTION OF THE PREDICTED VALUES               
                        DO 80 I=1,NY                                            
                           Y(I)      = Y(I)      -    E(I)                      
                           YPRIME(I) = YPRIME(I) - CJ*E(I)                      
                           E(I)      =  0.0D0                                   
   80                   CONTINUE                                                
                        NITER = 0                                               
                        S = 100.D0                                              
                     END IF                                                     
                  ELSE                                                          
C                          PROCEED WITH THE ITERATION                           
                     REDO = .TRUE.                                              
                  END IF                                                        
               ELSE                                                             
C                         ITERATION HAS CONVERGED                               
                  REDO = .FALSE.                                                
                  IRETC = 0                                                     
                  JCALC = 1                                                     
               END IF                                                           
            END IF                                                              
         ELSE                                                                   
C                         ITERATION HAS CONVERGED 
            REDO = .FALSE.                                                      
            IRETC = 0                                                           
            JCALC = 1 
         END IF                                                                 
      IF (REDO) GOTO 100                                                        
C     END DO WHILE                                                              
 999  CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
