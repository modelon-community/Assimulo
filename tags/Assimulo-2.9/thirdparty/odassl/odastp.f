      SUBROUTINE ODASTP(X,Y,YPRIME,NEQ,NY,                                 
     *                  RES,JAC,H,WT,JSTART,IDID,RPAR,IPAR,                     
     *                  PHI,DELTA,E,WM,IWM,                                     
     *                  ALPHA,BETA,GAMMA,PSI,SIGMA,                             
     *                  CJ,CJOLD,HOLD,S,HMIN,UROUND,                            
     *                  IPHASE,JCALC,K,KOLD,NS)                                 
C----------------------------------------------------------------------         
C      MODIFIED VERSION OF DASTP, A SUBROUTINE OF THE DRIVER PROGRAM            
C      ODASSL FOR SOLVING OVERDETERMINED SYSTEMS                                
C      OF (SINGULARLY) IMPLICIT ODES.                                           
C      DASSL-AUTHOR:  PETZOLD, LINDA                                            
C             APPLIED MATHEMATICS DIVISION 8331                                 
C             SANDIA NATIONAL LABORATORIES                                      
C             LIVERMORE, CA.    94550                                           
C      MODIFICATIONS BY: FUEHRER, CLAUS                                         
C             DEUTSCHE FORSCHUNGS- UND VERSUCHSANSTALT                          
C             FUER LUFT- UND RAUMFAHRT (DFVLR)                                  
C             INST. DYNAMIC DER FLUGSYSTEME                                     
C             D-8031 WESSLING  (F.R.G)                                          
C                                                                               
C      DATE OF LAST MODIFICATION 900103                                         
C                                                                               
C----------------------------------------------------------------------         
C                                                                               
C***REFER TO  ODASSL                                                            
C***ROUTINES CALLED  DDANRM,ODACOR,DDATRP                                       
C***COMMON BLOCKS    ODA001                                                     
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C     ODASTEP SOLVES A SYSTEM OF DIFFERENTIAL/                                  
C     ALGEBRAIC EQUATIONS OF THE FORM                                           
C     G(X,Y,YPRIME) = 0,  FOR ONE STEP (NORMALLY                                
C     FROM X TO X+H).                                                           
C                                                                               
C     THE METHODS USED ARE MODIFIED DIVIDED                                     
C     DIFFERENCE,FIXED LEADING COEFFICIENT                                      
C     FORMS OF BACKWARD DIFFERENTIATION                                         
C     FORMULAS. THE CODE ADJUSTS THE STEPSIZE                                   
C     AND ORDER TO CONTROL THE LOCAL ERROR PER                                  
C     STEP.                                                                     
C                                                                               
C                                                                               
C     THE PARAMETERS REPRESENT                                                  
C     X  --        INDEPENDENT VARIABLE                                         
C     Y  --        SOLUTION VECTOR AT X                                         
C     YPRIME --    DERIVATIVE OF SOLUTION VECTOR                                
C                  AFTER SUCCESSFUL STEP                                        
C     NEQ --       NUMBER OF EQUATIONS TO BE INTEGRATED                         
C     NY  --       NUMBER OF COMPONENTS OF Y OR YPRIME                          
C     RES --       EXTERNAL USER-SUPPLIED SUBROUTINE                            
C                  TO EVALUATE THE RESIDUAL. FOR THE CALL                       
C                  SEE DRIVER - PROGRAMM ODASSL.                                
C     JAC --       EXTERNAL USER-SUPPLIED ROUTINE -                             
C                  NOT REQUIRED BY THIS VERSION OF ODASSL                       
C     H --         APPROPRIATE STEP SIZE FOR NEXT STEP.                         
C                  NORMALLY DETERMINED BY THE CODE                              
C     WT --        VECTOR OF WEIGHTS FOR ERROR CRITERION.                       
C     JSTART --    INTEGER VARIABLE SET 0 FOR                                   
C                  FIRST STEP, 1 OTHERWISE.                                     
C     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS                  
C                  IDID= 1 -- THE STEP WAS COMPLETED SUCCESSFULLY               
C                  IDID=-6 -- THE ERROR TEST FAILED REPEATEDLY                  
C                  IDID=-7 -- THE CORRECTOR COULD NOT CONVERGE                  
C                  IDID=-8 -- THE ITERATION MATRIX IS SINGULAR                  
C                  IDID=-9 -- THE CORRECTOR COULD NOT CONVERGE.                 
C                             THERE WERE REPEATED ERROR TEST                    
C                             FAILURES ON THIS STEP.                            
C                  IDID=-10-- THE CORRECTOR COULD NOT CONVERGE                  
C                             BECAUSE IRES WAS EQUAL TO MINUS ONE               
C                  IDID=-11-- IRES EQUAL TO -2 WAS ENCOUNTERED,                 
C                             AND CONTROL IS BEING RETURNED TO                  
C                             THE CALLING PROGRAM                               
C     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS THAT                       
C                  ARE USED FOR COMMUNICATION BETWEEN THE                       
C                  CALLING PROGRAM AND EXTERNAL USER ROUTINES                   
C                  THEY ARE NOT ALTERED BY DASTEP                               
C     PHI --       ARRAY OF DIVIDED DIFFERENCES USED BY                         
C                  DASTEP. THE LENGTH IS NY*(K+1),WHERE                         
C                  K IS THE MAXIMUM ORDER                                       
C     DELTA,E --   WORK VECTORS FOR ODASTEP OF LENGTH NEQ AND NY                
C     WM,IWM --    REAL AND INTEGER ARRAYS STORING                              
C                  MATRIX INFORMATION SUCH AS THE MATRIX                        
C                  OF PARTIAL DERIVATIVES,PERMUTATION                           
C                  VECTOR,AND VARIOUS OTHER INFORMATION.                        
C                                                                               
C     THE OTHER PARAMETERS ARE INFORMATION                                      
C     WHICH IS NEEDED INTERNALLY BY DASTEP TO                                   
C     CONTINUE FROM STEP TO STEP.                                               
C                                                                               
C-----------------------------------------------------------------------        
C                                                                               
C                                                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      LOGICAL CONVGD                                                            
      DIMENSION Y(NY),YPRIME(NY),WT(NY)                                         
      DIMENSION PHI(NY,*),DELTA(NEQ),E(NY)                                      
      DIMENSION WM(*),IWM(*)                                                    
      DIMENSION PSI(*),ALPHA(*),BETA(*),GAMMA(*),SIGMA(*)                       
      DIMENSION RPAR(*),IPAR(*)                                                 
      EXTERNAL RES,JAC 
      COMMON/ODA001/NPD,NTEMP,INORM,                                            
     *   LML,LMU,LMXORD,LNY,LM1,LMTYPE,                                         
     *   LNST,LNRE,LNJE,LETF,LCTF,LIPVT                                         
      DATA MAXIT/4/                                                             
      DATA XRATE/0.25D0/                                                        
C                                                                               
C-----------------------------------------------------------------------        
C     BLOCK 1.                                                                  
C     INITIALIZE. ON THE FIRST CALL,SET                                         
C     THE ORDER TO 1 AND INITIALIZE                                             
C     OTHER VARIABLES.                                                          
C-----------------------------------------------------------------------        
C                                                                               
C     INITIALIZATIONS FOR ALL CALLS                                             
      IDID=1                                                                    
      NCF=0                                                                     
      NEF=0                                                                     
      M1 = NEQ - NY                                                              
      IF(JSTART .NE. 0) GO TO 120                                               
C                                                                               
C     IF THIS IS THE FIRST STEP,PERFORM                                         
C     OTHER INITIALIZATIONS                                                     
      IWM(LETF) = 0                                                             
      IWM(LCTF) = 0                                                             
      K=1                                                                       
      KOLD=0                                                                    
      HOLD=0.0D0                                                                
      JSTART=1                                                                  
      PSI(1)=H                                                                  
      CJOLD = 1.0D0/H                                                           
      CJ = CJOLD                                                                
      JCALC = -1                                                                
      IPHASE = 0                                                                
      NS=0                                                                      
120   CONTINUE       
                                                             
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C     BLOCK 2                                                                   
C     COMPUTE COEFFICIENTS OF FORMULAS FOR                                      
C     THIS STEP.                                                                
C-----------------------------------------------------------------------        
200   CONTINUE                                                                  
      KP1=K+1                                                                   
      KP2=K+2                                                                   
      KM1=K-1                                                                   
      XOLD=X                                                                    
      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0                                       
      NS=MIN0(NS+1,KOLD+2)                                                      
      NSP1=NS+1                                                                 
      IF(KP1 .LT. NS)GO TO 230                                                  
C                                                                               
      BETA(1)=1.0D0                                                             
      ALPHA(1)=1.0D0                                                            
      TEMP1=H                                                                   
      GAMMA(1)=0.0D0                                                            
      SIGMA(1)=1.0D0                                                            
      DO 210 I=2,KP1                                                            
         TEMP2=PSI(I-1)                                                         
         PSI(I-1)=TEMP1                                                         
         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2                                       
         TEMP1=TEMP2+H                                                          
         ALPHA(I)=H/TEMP1                                                       
         SIGMA(I)=DFLOAT(I-1)*SIGMA(I-1)*ALPHA(I)                               
         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H                                       
210      CONTINUE                                                               
      PSI(KP1)=TEMP1                                                            
230   CONTINUE                                                                  
C                                                                               
C     COMPUTE ALPHAS, ALPHA0                                                    
      ALPHAS = 0.0D0                                                            
      ALPHA0 = 0.0D0                                                            
      DO 240 I = 1,K                                                            
        ALPHAS = ALPHAS - 1.0D0/DFLOAT(I)                                       
        ALPHA0 = ALPHA0 - ALPHA(I)                                              
240     CONTINUE                                                                
C                                                                               
C     COMPUTE LEADING COEFFICIENT CJ                                            
      CJLAST = CJ                                                               
      CJ = -ALPHAS/H                                                            
C                                                                               
C     COMPUTE VARIABLE STEPSIZE ERROR COEFFICIENT CK                            
      CK = DABS(ALPHA(KP1) + ALPHAS - ALPHA0)                                   
      CK = DMAX1(CK,ALPHA(KP1))                                                 
C                                                                               
C     DECIDE WHETHER NEW JACOBIAN IS NEEDED                                     
      TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)                                   
      TEMP2 = 1.0D0/TEMP1                                                       
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1              
      IF (CJ .NE. CJLAST) S = 100.D0                                            
C                                                                               
C     CHANGE PHI TO PHI STAR                                                    
      IF(KP1 .LT. NSP1) GO TO 280                                               
      DO 270 J=NSP1,KP1                                                         
         DO 260 I=1,NY                                                          
260         PHI(I,J)=BETA(J)*PHI(I,J)                                           
270      CONTINUE                                                               
280   CONTINUE                                                                  
C                                                                               
C     UPDATE TIME                                                               
      X=X+H                                                                     
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C     BLOCK 3                                                                   
C     PREDICT THE SOLUTION AND DERIVATIVE,                                      
C     AND SOLVE THE CORRECTOR EQUATION                                          
C-----------------------------------------------------------------------        
C                                                                               
C     FIRST,PREDICT THE SOLUTION AND DERIVATIVE                                 
300   CONTINUE                                                                  
      DO 310 I=1,NY                                                             
         Y(I)=PHI(I,1)                                                          
310      YPRIME(I)=0.0D0                                                        
      DO 330 J=2,KP1                                                            
         DO 320 I=1,NY                                                          
            Y(I)=Y(I)+PHI(I,J)                                                  
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)                               
330   CONTINUE                                                                  
C                                                                               
C     SOLVE THE CORRECTOR EQUATION                                              
C                                                                               
      CALL ODACOR(NEQ, NY, M1, X, Y, YPRIME, RPAR, IPAR, JCALC,           
     *            H, CJ, CJOLD, UROUND, IWM(LNRE), IWM(LNJE), MAXIT,            
     *            RES, WT, IWM(LMTYPE), INORM, IWM(LML), IWM(LMU),                           
     *            E, DELTA, WM, IWM(LIPVT), IRETC, IRES, JAC)                        
      IF (IRETC.EQ.0 .AND. IRES.EQ.0) THEN                                      
          CONVGD = .TRUE.                                                       
      ELSE                                                                      
          CONVGD = .FALSE.                                                      
      END IF                                                                    
      IF(.NOT.CONVGD) GO TO 600                                                 
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C     BLOCK 4                                                                   
C     ESTIMATE THE ERRORS AT ORDERS K,K-1,K-2                                   
C     AS IF CONSTANT STEPSIZE WAS USED. ESTIMATE                                
C     THE LOCAL ERROR AT ORDER K AND TEST                                       
C     WHETHER THE CURRENT STEP IS SUCCESSFUL.                                   
C-----------------------------------------------------------------------        
C                                                                               
C     ESTIMATE ERRORS AT ORDERS K,K-1,K-2                                       
      ENORM = DDANRM(NY,E,WT,RPAR,IPAR)                                         
      ERK = SIGMA(K+1)*ENORM                                                    
      TERK = FLOAT(K+1)*ERK                                                     
      EST = ERK                                                                 
      KNEW=K                                                                    
      IF(K .EQ. 1)GO TO 430                                                     
      DO 405 I = 1,NY                                                           
405     DELTA(I) = PHI(I,KP1) + E(I)                                            
      ERKM1=SIGMA(K)*DDANRM(NY,DELTA,WT,RPAR,IPAR)                              
      TERKM1 = FLOAT(K)*ERKM1                                                   
      IF(K .GT. 2)GO TO 410                                                     
      IF(TERKM1 .LE. 0.5*TERK)GO TO 420                                         
      GO TO 430                                                                 
410   CONTINUE                                                                  
      DO 415 I = 1,NY                                                           
415     DELTA(I) = PHI(I,K) + DELTA(I)                                          
      ERKM2=SIGMA(K-1)*DDANRM(NY,DELTA,WT,RPAR,IPAR)                            
      TERKM2 = FLOAT(K-1)*ERKM2                                                 
      IF(DMAX1(TERKM1,TERKM2).GT.TERK)GO TO 430                                 
C     LOWER THE ORDER                                                           
420   CONTINUE                                                                  
      KNEW=K-1                                                                  
      EST = ERKM1                                                               
C                                                                               
C                                                                               
C     CALCULATE THE LOCAL ERROR FOR THE CURRENT STEP                            
C     TO SEE IF THE STEP WAS SUCCESSFUL                                         
430   CONTINUE                                                                  
      ERR = CK * ENORM                                                          
      IF(ERR .GT. 1.0D0)GO TO 600                                               
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C     BLOCK 5                                                                   
C     THE STEP IS SUCCESSFUL. DETERMINE                                         
C     THE BEST ORDER AND STEPSIZE FOR                                           
C     THE NEXT STEP. UPDATE THE DIFFERENCES                                     
C     FOR THE NEXT STEP.                                                        
C-----------------------------------------------------------------------        
      IDID=1                                                                    
      IWM(LNST)=IWM(LNST)+1                                                     
      KDIFF=K-KOLD                                                              
      KOLD=K                                                                    
      HOLD=H                                                                    
C                                                                               
C                                                                               
C     ESTIMATE THE ERROR AT ORDER K+1 UNLESS                                    
C        ALREADY DECIDED TO LOWER ORDER, OR                                     
C        ALREADY USING MAXIMUM ORDER, OR                                        
C        STEPSIZE NOT CONSTANT, OR                                              
C        ORDER RAISED IN PREVIOUS STEP                                          
      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1                               
      IF(IPHASE .EQ. 0)GO TO 545                                                
      IF(KNEW.EQ.KM1)GO TO 540                                                  
      IF(K.EQ.IWM(LMXORD)) GO TO 550                                            
      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550                                      
      DO 510 I=1,NY                                                             
510      DELTA(I)=E(I)-PHI(I,KP2)                                               
      ERKP1 = (1.0D0/DFLOAT(K+2))*DDANRM(NY,DELTA,WT,RPAR,IPAR)                 
      TERKP1 = FLOAT(K+2)*ERKP1                                                 
      IF(K.GT.1)GO TO 520                                                       
      IF(TERKP1.GE.0.5D0*TERK)GO TO 550                                         
      GO TO 530                                                                 
520   IF(TERKM1.LE.DMIN1(TERK,TERKP1))GO TO 540                                 
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550                           
C                                                                               
C     RAISE ORDER                                                               
530   K=KP1                                                                     
      EST = ERKP1                                                               
      GO TO 550                                                                 
C                                                                               
C     LOWER ORDER                                                               
540   K=KM1                                                                     
      EST = ERKM1                                                               
      GO TO 550                                                                 
C                                                                               
C     IF IPHASE = 0, INCREASE ORDER BY ONE AND MULTIPLY STEPSIZE BY             
C     FACTOR TWO                                                                
545   K = KP1                                                                   
      HNEW = H*2.0D0                                                            
      H = HNEW                                                                  
      GO TO 575                                                                 
C                                                                               
C                                                                               
C     DETERMINE THE APPROPRIATE STEPSIZE FOR                                    
C     THE NEXT STEP.                                                            
550   HNEW=H                                                                    
      TEMP2=K+1                                                                 
      R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)                                    
      IF(R .LT. 2.0D0) GO TO 555                                                
      HNEW = 2.0D0*H                                                            
      GO TO 560                                                                 
555   IF(R .GT. 1.0D0) GO TO 560                                                
      R = DMAX1(0.5D0,DMIN1(0.9D0,R))                                           
      HNEW = H*R                                                                
560   H=HNEW                                                                    
C                                                                               
C                                                                               
C     UPDATE DIFFERENCES FOR NEXT STEP                                          
575   CONTINUE                                                                  
      IF(KOLD.EQ.IWM(LMXORD))GO TO 585                                          
      DO 580 I=1,NY                                                             
580      PHI(I,KP2)=E(I)                                                        
585   CONTINUE                                                                  
      DO 590 I=1,NY                                                             
590      PHI(I,KP1)=PHI(I,KP1)+E(I)                                             
      DO 595 J1=2,KP1                                                           
         J=KP1-J1+1                                                             
         DO 595 I=1,NY                                                          
595      PHI(I,J)=PHI(I,J)+PHI(I,J+1)                                           
      RETURN                                                                    
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C-----------------------------------------------------------------------        
C     BLOCK 6                                                                   
C     THE STEP IS UNSUCCESSFUL. RESTORE X,PSI,PHI                               
C     DETERMINE APPROPRIATE STEPSIZE FOR                                        
C     CONTINUING THE INTEGRATION, OR EXIT WITH                                  
C     AN ERROR FLAG IF THERE HAVE BEEN MANY                                     
C     FAILURES.                                                                 
C-----------------------------------------------------------------------        
600   IPHASE = 1                                                                
C                                                                               
C     RESTORE X,PHI,PSI                                                         
      X=XOLD                                                                    
      IF(KP1.LT.NSP1)GO TO 630                                                  
      DO 620 J=NSP1,KP1                                                         
         TEMP1=1.0D0/BETA(J)                                                    
         DO 610 I=1,NY                                                          
610         PHI(I,J)=TEMP1*PHI(I,J)                                             
620      CONTINUE                                                               
630   CONTINUE                                                                  
      DO 640 I=2,KP1                                                            
640      PSI(I-1)=PSI(I)-H                                                      
C                                                                               
C                                                                               
C     TEST WHETHER FAILURE IS DUE TO CORRECTOR ITERATION                        
C     OR ERROR TEST                                                             
      IF(CONVGD)GO TO 660                                                       
      IWM(LCTF)=IWM(LCTF)+1                                                     
C                                                                               
C                                                                               
C     THE NEWTON ITERATION FAILED TO CONVERGE WITH                              
C     A CURRENT ITERATION MATRIX.  DETERMINE THE CAUSE                          
C     OF THE FAILURE AND TAKE APPROPRIATE ACTION.                               
      IF(IRETC .LT. 2) GOTO 650                                                 
C                                                                               
C     THE ITERATION MATRIX IS SINGULAR. REDUCE                                  
C     THE STEPSIZE BY A FACTOR OF 4. IF                                         
C     THIS HAPPENS THREE TIMES IN A ROW ON                                      
C     THE SAME STEP, RETURN WITH AN ERROR FLAG                                  
      R = 0.25D0                                                                
      H=H*R                                                                     
      IF (IRETC .LT. 3 .AND. DABS(H) .GE. HMIN) GO TO 690                       
      IDID=-8                                                                   
      GO TO 675                                                                 
C                                                                               
C                                                                               
C     THE NEWTON ITERATION FAILED TO CONVERGE FOR A REASON                      
C     OTHER THAN A SINGULAR ITERATION MATRIX.  IF IRES = -2, THEN               
C     RETURN.  OTHERWISE, REDUCE THE STEPSIZE AND TRY AGAIN, UNLESS             
C     TOO MANY FAILURES HAVE OCCURED.                                           
650   CONTINUE                                                                  
      IF (IRES .EQ. 0) GO TO 655                                                
      IDID = -11                                                                
      GO TO 675                                                                 
655   NCF = NCF + 1                                                             
      R = 0.25D0                                                                
      H = H*R                                                                   
      IF (NCF .LT. 10 .AND. DABS(H) .GE. HMIN) GO TO 690                        
      IDID = -7                                                                 
      IF (IRES .NE. 0) IDID = -10                                               
      IF (NEF .GE. 3) IDID = -9                                                 
      GO TO 675                                                                 
C                                                                               
C                                                                               
C     THE NEWTON SCHEME CONVERGED,AND THE CAUSE                                 
C     OF THE FAILURE WAS THE ERROR ESTIMATE                                     
C     EXCEEDING THE TOLERANCE.                                                  
660   NEF=NEF+1                                                                 
      IWM(LETF)=IWM(LETF)+1                                                     
      IF (NEF .GT. 1) GO TO 665                                                 
C                                                                               
C     ON FIRST ERROR TEST FAILURE, KEEP CURRENT ORDER OR LOWER                  
C     ORDER BY ONE.  COMPUTE NEW STEPSIZE BASED ON DIFFERENCES                  
C     OF THE SOLUTION.                                                          
      K = KNEW                                                                  
      TEMP2 = K + 1                                                             
      R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)                           
      R = DMAX1(0.25D0,DMIN1(0.9D0,R))                                          
      H = H*R                                                                   
      IF (DABS(H) .GE. HMIN) GO TO 690                                          
      IDID = -6                                                                 
      GO TO 675                                                                 
C                                                                               
C     ON SECOND ERROR TEST FAILURE, USE THE CURRENT ORDER OR                    
C     DECREASE ORDER BY ONE.  REDUCE THE STEPSIZE BY A FACTOR OF                
C     ONE QUARTER.                                                              
665   IF (NEF .GT. 2) GO TO 670                                                 
      K = KNEW                                                                  
      H = 0.25D0*H                                                              
      IF (DABS(H) .GE. HMIN) GO TO 690                                          
      IDID = -6                                                                 
      GO TO 675                                                                 
C                                                                               
C     ON THIRD AND SUBSEQUENT ERROR TEST FAILURES, SET THE ORDER TO             
C     ONE AND REDUCE THE STEPSIZE BY A FACTOR OF ONE QUARTER                    
670   K = 1                                                                     
      H = 0.25D0*H                                                              
      IF (DABS(H) .GE. HMIN) GO TO 690                                          
      IDID = -6                                                                 
      GO TO 675                                                                 
C                                                                               
C                                                                               
C                                                                               
C                                                                               
C     FOR ALL CRASHES, RESTORE Y TO ITS LAST VALUE,                             
C     INTERPOLATE TO FIND YPRIME AT LAST X, AND RETURN                          
675   CONTINUE                                                                  
      CALL DDATRP(X,X,Y,YPRIME,NY,K,PHI,PSI)                                    
      RETURN                                                                    
C                                                                               
C                                                                               
C     GO BACK AND TRY THIS STEP AGAIN                                           
690   GO TO 200                                                                 
C                                                                               
C------END OF SUBROUTINE DASTEP------                                           
      END                                                                       
