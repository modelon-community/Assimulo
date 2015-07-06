      SUBROUTINE ODAJAC (NEQ, NY, M1, X, Y, YPRIME, RPAR, IPAR,            
     *                  RES, DELTA, UROUND, CJ, H,                              
     *                  WT, AUX, A1 , A2, A2T, Q,                                 
     *                  XNORM, MTYPE, INORM, ML, MU,                                   
     *                  IPVTA1, IPVTQ,                                    
     *                  IRES, IERR, JAC)                                  
C----------------------------------------------------------------------         
C      MODIFIED VERSION OF THE DASSL-SUBROUTINE DAJAC. 
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
C***REFER TO  ODAJAC                                                            
C***ROUTINES CALLED  RES, DGEFA, DGESL                           
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
C     RES --       EXTERNAL USER-SUPPLIED SUBROUTINE                            
C                  TO EVALUATE THE RESIDUAL. FOR THE CALL                       
C                  SEE DRIVER - PROGRAMM ODASSL.                                
C     DELTA --     RESIDUEVECTOR (RESULT OF THE LAST CALL TO RES)               
C     UROUND --    MACINE PRECISION                                             
C     CJ --        LEADING COEFFICIENT FOR THE CURRENT STEP                     
C     H --         STEP SIZE                                                    
C     WT --        VECTOR OF WEIGHTS FOR ERROR CRITERION.                       
C     A1, A2 --    LU-FACTORED ITERATION MATRICES          
C     A2T --       A2 TRANSPOSE  
C     Q----        LU FACTORED MATRIX A2*A1^INV*A2T  
C     IVPTA1, 
C     IPVTQ--      THE CORRESPONDING PIVOT VECTORS
C     MTYPE --     = 1 IF JACOBIAN IS GIVEN BY SUBROUTINE JAC
C                  = 2 IS TO BE COMPUTED BY FINITE DIFFERENCES    
C     XNORM --     VECTOR WITH THE SCALING FACTORS (SEE INORM BELOW)            
C     INORM, ML,                                                                
C     MU --        IF INORM .NE. 0 THE ROWS ML UP TO MU OF THE                  
C                  ITERATION MATRIX ARE TO BE SCALED                            
C     IERR --     RETURN-CODE: = 0  EVERYTHING IS O.K.                          
C                              = 1  SINGULARITY DETECTED                        
C                                                                               
C     JAC  --      IF MTYPE = 1: SUBROUTINE FOR GENERATING THE JACOBIAN
C                  DUMMY OTHERWISE   C                                                                               
C-----------------------------------------------------------------------        
C                                                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      EXTERNAL RES                                                             
      DIMENSION Y(NY),YPRIME(NY),DELTA(NEQ),WT(NY),AUX(NEQ),                  
     *          A1(NY, NY), A2(M1,NY), A2T(NY,M1),RPAR(*),IPAR(*),            
     *          Q(M1,M1),IPVTA1(NY),IPVTQ(M1),XNORM(NEQ)                          
C
      IERR = 0                                                                
      IRES = 0  
      SQUR = SQRT(UROUND) 
      IF (MTYPE .EQ. 1) THEN 
C                                                                               
C        1A. COMPUTATION OF THE ITERATION MATRIX                                 
C                  DG/DY + CJ*DG/DYPRIME                                        
C                  BY CALLING JAC
         CALL JAC(T,Y,YPRIME,A1,A2,CJ,RPAR,IPAR) 
         DO 5, I=1,M1
            DO 5,J=1,NY
               A2T(J,I) = A2(I,J)
   5     CONTINUE
      ELSE                                    
C                                                       
C                                                                               
C        1B. COMPUTATION OF THE ITERATION MATRIX                                 
C                  DG/DY + CJ*DG/DYPRIME                                        
C                  BY FINITE DIFFERENCES                                        
C                                                                               
         DO 30 I=1,NY                                                              
            DEL=SQUR * MAX( ABS(Y(I)), ABS(H*YPRIME(I)), ABS(WT(I)))               
            DEL=SIGN(DEL,H*YPRIME(I))                                              
            DEL=(Y(I)+DEL)-Y(I)                                                    
            YSAVE=Y(I)                                                             
            YPSAVE=YPRIME(I)                                                       
            Y(I)=Y(I)+DEL                                                          
            YPRIME(I)=YPRIME(I)+CJ*DEL                                             
            CALL RES(X,NEQ, NY, Y,YPRIME,AUX,IRES,RPAR,IPAR)                              
            IF (IRES .NE. 0) GOTO 999                                              
            DELINV=1.0D0/DEL                                                       
C                COMPUTATION OF THE ITERATION MATRIX A1                       
            DO 10 L=1,NY                                                       
               A1(L,I) = (AUX(L)-DELTA(L))*DELINV                           
 10         CONTINUE  
            IF (M1.GT.0) THEN
C                COMPUTATION OF THE ITERATION MATRIX A2,A2T                       
               DO 20  L=1,M1
                 A2(L,I) = (AUX(L+NY)-DELTA(L+NY))*DELINV                           
                 A2T(I,L) = A2(L,I)
  20           CONTINUE
            END IF                                    
            Y(I)=YSAVE                                                             
            YPRIME(I)=YPSAVE                                                       
 30      CONTINUE 
      END IF 
C        2. SCALING OF THE ITERATION IF REQUIRED -                                 
C           THE ROWS ML UP TO MU ARE SCALED IN SUCH A WAY THAT                     
C           THEY HAVE UNIT TWO-NORM AFTERWARDS                                     
      IF (INORM .NE. 1) GOTO 100 
         IF (ML.LE.NY) THEN
C           SCALING OF A1 IS REQUIRED
            MUA1 = MIN(MU,NY)
            DO 50 I=ML,MUA1                                                          
               XNORM(I) = 0.0D0                                                     
               DO 45 L=1,NY                                                         
                 XNORM(I) = XNORM(I) + A1(I,L)**2                                
 45            CONTINUE                                                             
               XNORM(I) = SQRT(XNORM(I))                                            
               IF (XNORM(I) .LT. UROUND) THEN                                       
C                  MATRIX SINGULAR                                                      
                  IERR = 1                                                          
                  GOTO 999                                                          
                END IF                                                               
                DO 48 L=1,NY                                                         
                   A1(I,L) = A1(I,L)/XNORM(I)                                    
 48             CONTINUE                                                             
 50          CONTINUE 
          END IF
          IF (ML .GT. NY .OR. MU .GT. NY) THEN  
C           SCALING OF A2 IS REQUIRED
            MLA2 = MAX(ML,NY) - NY 
            MUA2 = MU - NY 
            DO 70 I=MLA2,MUA2
               XNORM(I+NY) = 0.0D0                                                     
               DO 55 L=1,NY                                                         
                 XNORM(I+NY) = XNORM(I+NY) + A2(I,L)**2                                
 55            CONTINUE                                                             
               XNORM(I+NY) = SQRT(XNORM(I+NY))                                            
               IF (XNORM(I+NY) .LT. UROUND) THEN                                       
C                  MATRIX SINGULAR                                                      
                  IERR = 1                                                          
                  GOTO 999                                                          
               END IF                                                               
               DO 60 L=1,NY                                                         
                  A2(I,L) = A2(I,L)/XNORM(I+NY)  
                  A2T(L,I) = A2(I,L)                                  
 60            CONTINUE                                                             
 70         CONTINUE 
         END IF  
100   CONTINUE                                                                  
C                                                                               
C                     3. LU - DECOMPOSITION OF A1                                     
      CALL DGEFA(A1,NY,NY,IPVTA1,IERR) 
      IF (IERR .NE. 0) GOTO 999
      IF (M1.EQ. 0) GOTO 999
C
C                     4.(A1)-1 * (A2)T   
      DO 110 I=1,M1
         CALL DGESL(A1,NY,NY,IPVTA1,A2T(1,I),0)
 110  CONTINUE
C
C                     5.Q=A2*(A1)-1 * (A2)T)   
      DO 140 L=1,M1
         DO 130 I=1,M1
            Q(L,I) = 0.0D0
            DO 120 J=1,NY
               Q(L,I) = Q(L,I) + A2(L,J)*A2T(J,I)
 120        CONTINUE
 130     CONTINUE
 140  CONTINUE
C
C                     6. LU - DECOMPOSITION OF Q                                     
      CALL DGEFA(Q,M1,M1,IPVTQ,IERR)
      IF (IERR .NE. 0 ) GOTO 999
C  ===============================================================
 999  CONTINUE  
      RETURN                                                                    
      END                                                                       
