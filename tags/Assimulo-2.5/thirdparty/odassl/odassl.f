      SUBROUTINE ODASSL (RES,NEQ,NY,T,Y,YPRIME,TOUT,                               
     *  INFO,RTOL,ATOL,IDID,                                                    
     *  RWORK,LRW,IWORK,LIW,RPAR,IPAR,                                          
     *  JAC)                                                                    
C                                                                               
C----------------------------------------------------------------------         
C      MODIFIED VERSION OF DASSL FOR SOLVING OVERDETERMINED SYSTEMS             
C      OF (SINGULARLY) IMPLICIT ODES. THE MAIN DIFFERENCE TO DASSL IS           
C      IN THE CORRECTOR ITERATION PART.                                         
C                                                                               
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
C----------------------------------------------------------------------         
C                                                                               
C  THIS CODE SOLVES A SYSTEM OF OVERDETERMINED DIFFERENTIAL/                    
C  ALGEBRAIC EQUATIONS OF THE FORM                                              
C  G(T,Y,YPRIME) = 0. IT IS SOLVED IN SUCH A WAY, THAT THE
C  THE ADDITIONAL LAST NEQ-NY EQUATIONS (IF THERE ARE ANY)
C  ARE TAKEN TO DEFINE SOLUTION INVARIANTS TO THE "SQUARE"
C  NY X NY - SYSTEM OF THE FIRST NY -EQUATIONS. THE ALGEBRAIC SYSTEMS
C  ARE THEN SOLVED BY USING THE SSF-ITERATION IF NEQ > NY OTHERWISE BY
C  A MODIFIED NEWTON-SCHEME. IN THE LAST CASE ODASSL IS ALGORITHMICALLY IDENTICAL
C  TO DASSL.                                                           
C                                                                               
C  SUBROUTINE ODASSL USES THE BACKWARD                                          
C  DIFFERENTIATION FORMULAS OF ORDERS ONE                                       
C  THROUGH FIVE TO SOLVE A  SYSTEM OF THE ABOVE                                  
C  FORM FOR Y AND YPRIME. VALUES FOR Y                                          
C  AND YPRIME AT THE INITIAL TIME MUST                                          
C  BE GIVEN AS INPUT. THESE VALUES MUST                                         
C  BE CONSISTENT, (THAT IS. IF T,Y,YPRIME                                       
C  ARE THE GIVEN INITIAL VALUES, THEY MUST                                      
C  SATISFY G(T,Y,YPRIME) = 0.)                                                  
C  THE SUBROUTINE SOLVES THE SYSTEM FROM T TO TOUT. IT IS                       
C  EASY TO CONTINUE THE SOLUTION TO GET RESULTS                                 
C  AT ADDITIONAL TOUT. THIS IS THE INTERVAL                                     
C  MODE OF OPERATION. INTERMEDIATE RESULTS CAN                                  
C  ALSO BE OBTAINED EASILY BY USING THE INTERMEDIATE-                           
C  OUTPUT CAPABILITY.                                                           
C                                                                               
C                                                                               
CC ------------DESCRIPTION OF ARGUMENTS TO ODASSL-----------------------        
C  ------------(AN OVERVIEW)--------------------------------------------        
C                                                                               
C  THE PARAMETERS ARE                                                           
C                                                                               
C  RES -- THIS IS A SUBROUTINE WHICH YOU PROVIDE                                
C         TO DEFINE THE DIFFERENTIAL/ALGEBRAIC                                  
C         SYSTEM                                                                
C                                                                               
C  NEQ -- THIS IS THE NUMBER OF EQUATIONS                                       
C         TO BE SOLVED   
C
C  NY -- NUMBER OF COMPONENTS OF THE SOLUTION VECTOR Y                                                       
C                                                                               
C                                                                               
C  T -- THIS IS THE CURRENT VALUE OF THE                                        
C       INDEPENDENT VARIABLE.                                                   
C                                                                               
C  TOUT -- THIS IS A POINT AT WHICH A SOLUTION                                  
C      IS DESIRED.                                                              
C                                                                               
C  INFO(*) -- THE BASIC TASK OF THE CODE IS                                     
C             TO SOLVE THE SYSTEM FROM T TO                                     
C             TOUT AND RETURN AN ANSWER AT TOUT.                                
C             INFO(*) IS AN INTEGER ARRAY WHICH IS                              
C             USED TO COMMUNICATE EXACTLY HOW YOU                               
C             WANT THIS TASK TO BE CARRIED OUT.                                 
C                                                                               
C  Y(*) -- THIS ARRAY CONTAINS THE SOLUTION                                     
C          COMPONENTS AT T                                                      
C                                                                               
C  YPRIME(*) -- THIS ARRAY CONTAINS THE DERIVATIVES                             
C               OF THE SOLUTION COMPONENTS AT T                                 
C                                                                               
C                                                                               
C  RTOL,ATOL -- THESE QUANTITIES REPRESENT                                      
C               ABSOLUTE AND RELATIVE ERROR                                     
C               TOLERANCES WHICH YOU PROVIDE TO INDICATE                        
C               HOW ACCURATELY YOU WISH THE SOLUTION                            
C               TO BE COMPUTED. YOU MAY CHOOSE THEM                             
C               TO BE BOTH SCALARS OR ELSE BOTH                                 
C               VECTORS.                                                        
C                                                                               
C  IDID -- THIS SCALAR QUANTITY IS AN INDICATOR REPORTING                       
C          WHAT THE CODE DID. YOU MUST MONITOR THIS                             
C          INTEGER VARIABLE TO DECIDE WHAT ACTION TO                            
C          TAKE NEXT.                                                           
C                                                                               
C  RWORK(*),LRW -- RWORK(*) IS A REAL WORK ARRAY OF                             
C                  LENGTH LRW WHICH PROVIDES THE CODE                           
C                  WITH NEEDED STORAGE SPACE.                                   
C                                                                               
C  IWORK(*),LIW -- IWORK(*) IS AN INTEGER WORK ARRAY                            
C                  OF LENGTH LIW WHICH PROVIDES THE CODE                        
C                  WITH NEEDED STORAGE SPACE.                                   
C                                                                               
C  RPAR,IPAR -- THESE ARE REAL AND INTEGER PARAMETER                            
C               ARRAYS WHICH YOU CAN USE FOR                                    
C               COMMUNICATION BETWEEN YOUR CALLING                              
C               PROGRAM AND THE RES SUBROUTINE                                  
C               (AND THE ODAJAC SUBROUTINE)                                        
C                                                                               
C  JAC -- THIS IS THE NAME OF A SUBROUTINE WHICH YOU                            
C         MAY CHOOSE TO PROVIDE FOR DEFINING                                    
C         A MATRIX OF PARTIAL DERIVATIVES                                       
C         DESCRIBED BELOW.                                                      
C                                                                               
C  QUANTITIES WHICH ARE USED AS INPUT ITEMS ARE                                 
C     NEQ,T,Y(*),YPRIME(*),TOUT,INFO(*),                                        
C     RTOL,ATOL,RWORK(1),RWORK(2),RWORK(3),LRW,IWORK(1),                        
C     IWORK(2),IWORK(3),IWORK(4), IWORK(5) AND LIW.                             
C                                                                               
C  QUANTITIES WHICH MAY BE ALTERED BY THE CODE ARE                              
C     T,Y(*),YPRIME(*),INFO(1),RTOL,ATOL,                                       
C     IDID,RWORK(*) AND IWORK(*)                                                
C                                                                               
C  ----------INPUT-WHAT TO DO ON THE FIRST CALL TO ODASSL---------------        
C                                                                               
C                                                                               
C  THE FIRST CALL OF THE CODE IS DEFINED TO BE THE START OF EACH NEW            
C  PROBLEM. READ THROUGH THE DESCRIPTIONS OF ALL THE FOLLOWING ITEMS,           
C  PROVIDE SUFFICIENT STORAGE SPACE FOR DESIGNATED ARRAYS, SET                  
C  APPROPRIATE VARIABLES FOR THE INITIALIZATION OF THE PROBLEM, AND             
C  GIVE INFORMATION ABOUT HOW YOU WANT THE PROBLEM TO BE SOLVED.                
C                                                                               
C                                                                               
C  RES -- PROVIDE A SUBROUTINE OF THE FORM                                      
C             SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)                   
C         TO DEFINE THE SYSTEM OF DIFFERENTIAL/ALGEBRAIC                        
C         EQUATIONS WHICH IS TO BE SOLVED. FOR THE GIVEN VALUES                 
C         OF T,Y AND YPRIME, THE SUBROUTINE SHOULD                              
C         RETURN THE RESIDUAL OF THE DIFFERENTIAL/ALGEBRAIC                     
C         SYSTEM                                                                
C             DELTA = G(T,Y,YPRIME)                                             
C         (DELTA(*) IS A VECTOR OF LENGTH NEQ WHICH IS                          
C         OUTPUT FOR RES.)                                                      
C                                                                               
C         SUBROUTINE RES MUST NOT ALTER T,Y OR YPRIME.                          
C         YOU MUST DECLARE THE NAME RES IN AN EXTERNAL                          
C         STATEMENT IN YOUR PROGRAM THAT CALLS ODASSL.                          
C         YOU MUST DIMENSION Y,YPRIME AND DELTA IN RES.                         
C                                                                               
C         IRES IS AN INTEGER FLAG WHICH IS ALWAYS EQUAL TO                      
C         ZERO ON INPUT.  SUBROUTINE RES SHOULD ALTER IRES                      
C         ONLY IF IT ENCOUNTERS AN ILLEGAL VALUE OF Y OR                        
C         A STOP CONDITION.  SET IRES = -1 IF AN INPUT VALUE                    
C         IS ILLEGAL, AND ODASSL WILL TRY TO SOLVE THE PROBLEM                  
C         WITHOUT GETTING IRES = -1.  IF IRES = -2, ODASSL                      
C         WILL RETURN CONTROL TO THE CALLING PROGRAM                            
C         WITH IDID = -11.                                                      
C                                                                               
C         RPAR AND IPAR ARE REAL AND INTEGER PARAMETER ARRAYS WHICH             
C         YOU CAN USE FOR COMMUNICATION BETWEEN YOUR CALLING PROGRAM            
C         AND SUBROUTINE RES. THEY ARE NOT ALTERED BY ODASSL. IF YOU            
C         DO NOT NEED RPAR OR IPAR, IGNORE THESE PARAMETERS BY TREAT-           
C         ING THEM AS DUMMY ARGUMENTS. IF YOU DO CHOOSE TO USE THEM,            
C         DIMENSION THEM IN YOUR CALLING PROGRAM AND IN RES AS ARRAYS           
C         OF APPROPRIATE LENGTH.                                                
C                                                                               
C  NEQ -- SET IT TO THE NUMBER OF DIFFERENTIAL EQUATIONS.                       
C         (NEQ .GE. 1)  
C
C  NY --  DIMENSION OF SOLUTION VECTOR Y                                                         
C                                                                               
C         FOR DAES THE NUMBER OF EQUATIONS (NEQ) AND THE DIMENSION              
C         OF Y (NY) IS IDENTICAL. FOR OVERDETERMINED DAES (ODAES)               
C         NY .LT. NEQ HOLDS.                      
C                                                                               
C  T -- SET IT TO THE INITIAL POINT OF THE INTEGRATION.                         
C       T MUST BE DEFINED AS A VARIABLE.                                        
C                                                                               
C  Y(*) -- SET THIS VECTOR TO THE INITIAL VALUES OF THE NY SOLUTION             
C          COMPONENTS AT THE INITIAL POINT. YOU MUST DIMENSION Y OF             
C          LENGTH AT LEAST NY IN YOUR CALLING PROGRAM.                          
C                                                                               
C  YPRIME(*) -- SET THIS VECTOR TO THE INITIAL VALUES OF                        
C               THE NY FIRST DERIVATIVES OF THE SOLUTION                        
C               COMPONENTS AT THE INITIAL POINT. YOU                            
C               MUST DIMENSION YPRIME AT LEAST NY                               
C               IN YOUR CALLING PROGRAM.  IF YOU DO NOT                         
C               KNOW INITIAL VALUES OF SOME OF THE SOLUTION                     
C               COMPONENTS, SEE THE EXPLANATION OF INFO(11).                    
C                                                                               
C  TOUT - SET IT TO THE FIRST POINT AT WHICH A SOLUTION                         
C         IS DESIRED. YOU CAN NOT TAKE TOUT = T.                                
C         INTEGRATION EITHER FORWARD IN T (TOUT .GT. T) OR                      
C         BACKWARD IN T (TOUT .LT. T) IS PERMITTED.                             
C                                                                               
C         THE CODE ADVANCES THE SOLUTION FROM T TO TOUT USING                   
C         STEP SIZES WHICH ARE AUTOMATICALLY SELECTED SO AS TO                  
C         ACHIEVE THE DESIRED ACCURACY. IF YOU WISH, THE CODE WILL              
C         RETURN WITH THE SOLUTION AND ITS DERIVATIVE AT                        
C         INTERMEDIATE STEPS (INTERMEDIATE-OUTPUT MODE) SO THAT                 
C         YOU CAN MONITOR THEM, BUT YOU STILL MUST PROVIDE TOUT IN              
C         ACCORD WITH THE BASIC AIM OF THE CODE.                                
C                                                                               
C         THE FIRST STEP TAKEN BY THE CODE IS A CRITICAL ONE                    
C         BECAUSE IT MUST REFLECT HOW FAST THE SOLUTION CHANGES NEAR            
C         THE INITIAL POINT. THE CODE AUTOMATICALLY SELECTS AN                  
C         INITIAL STEP SIZE WHICH IS PRACTICALLY ALWAYS SUITABLE FOR            
C         THE PROBLEM. BY USING THE FACT THAT THE CODE WILL NOT STEP            
C         PAST TOUT IN THE FIRST STEP, YOU COULD, IF NECESSARY,                 
C         RESTRICT THE LENGTH OF THE INITIAL STEP SIZE.                         
C                                                                               
C         FOR SOME PROBLEMS IT MAY NOT BE PERMISSABLE TO INTEGRATE              
C         PAST A POINT TSTOP BECAUSE A DISCONTINUITY OCCURS THERE               
C         OR THE SOLUTION OR ITS DERIVATIVE IS NOT DEFINED BEYOND               
C         TSTOP. WHEN YOU HAVE DECLARED A TSTOP POINT (SEE INFO(4)              
C         AND RWORK(1)), YOU HAVE TOLD THE CODE NOT TO INTEGRATE                
C         PAST TSTOP. IN THIS CASE ANY TOUT BEYOND TSTOP IS INVALID             
C         INPUT.                                                                
C                                                                               
C  INFO(*) - USE THE INFO ARRAY TO GIVE THE CODE MORE DETAILS ABOUT             
C            HOW YOU WANT YOUR PROBLEM SOLVED. THIS ARRAY SHOULD BE             
C            DIMENSIONED OF LENGTH 15, THOUGH ODASSL USES                       
C            ONLY THE FIRST 14 ENTRIES. YOU MUST RESPOND TO ALL OF              
C            THE FOLLOWING ITEMS WHICH ARE ARRANGED AS QUESTIONS. THE           
C            SIMPLEST USE OF THE CODE IN THE DAE-APPLICATIONS                   
C            CORRESPONDS TO ANSWERING ALL                                       
C            QUESTIONS AS YES ,I.E. SETTING ALL ENTRIES OF INFO TO 0.           
C                                                                               
C       INFO(1) - THIS PARAMETER ENABLES THE CODE TO INITIALIZE                 
C              ITSELF. YOU MUST SET IT TO INDICATE THE START OF EVERY           
C              NEW PROBLEM.                                                     
C                                                                               
C          **** IS THIS THE FIRST CALL FOR THIS PROBLEM ...                     
C                YES - SET INFO(1) = 0                                          
C                 NO - NOT APPLICABLE HERE.                                     
C                      SEE BELOW FOR CONTINUATION CALLS.  ****                  
C                                                                               
C       INFO(2) - HOW MUCH ACCURACY YOU WANT OF YOUR SOLUTION                   
C              IS SPECIFIED BY THE ERROR TOLERANCES RTOL AND ATOL.              
C              THE SIMPLEST USE IS TO TAKE THEM BOTH TO BE SCALARS.             
C              TO OBTAIN MORE FLEXIBILITY, THEY CAN BOTH BE VECTORS.            
C              THE CODE MUST BE TOLD YOUR CHOICE.                               
C                                                                               
C          **** ARE BOTH ERROR TOLERANCES RTOL, ATOL SCALARS ...                
C                YES - SET INFO(2) = 0                                          
C                      AND INPUT SCALARS FOR BOTH RTOL AND ATOL                 
C                 NO - SET INFO(2) = 1                                          
C                      AND INPUT ARRAYS FOR BOTH RTOL AND ATOL ****             
C                                                                               
C       INFO(3) - THE CODE INTEGRATES FROM T IN THE DIRECTION                   
C              OF TOUT BY STEPS. IF YOU WISH, IT WILL RETURN THE                
C              COMPUTED SOLUTION AND DERIVATIVE AT THE NEXT                     
C              INTERMEDIATE STEP (THE INTERMEDIATE-OUTPUT MODE) OR              
C              TOUT, WHICHEVER COMES FIRST. THIS IS A GOOD WAY TO               
C              PROCEED IF YOU WANT TO SEE THE BEHAVIOR OF THE SOLUTION.         
C              IF YOU MUST HAVE SOLUTIONS AT A GREAT MANY SPECIFIC              
C              TOUT POINTS, THIS CODE WILL COMPUTE THEM EFFICIENTLY.            
C                                                                               
C          **** DO YOU WANT THE SOLUTION ONLY AT                                
C                TOUT (AND NOT AT THE NEXT INTERMEDIATE STEP) ...               
C                 YES - SET INFO(3) = 0                                         
C                  NO - SET INFO(3) = 1 ****                                    
C                                                                               
C       INFO(4) - TO HANDLE SOLUTIONS AT A GREAT MANY SPECIFIC                  
C              VALUES TOUT EFFICIENTLY, THIS CODE MAY INTEGRATE PAST            
C              TOUT AND INTERPOLATE TO OBTAIN THE RESULT AT TOUT.               
C              SOMETIMES IT IS NOT POSSIBLE TO INTEGRATE BEYOND SOME            
C              POINT TSTOP BECAUSE THE EQUATION CHANGES THERE OR IT IS          
C              NOT DEFINED PAST TSTOP. THEN YOU MUST TELL THE CODE              
C              NOT TO GO PAST.                                                  
C                                                                               
C           **** CAN THE INTEGRATION BE CARRIED OUT WITHOUT ANY                 
C                RESTRICTIONS ON THE INDEPENDENT VARIABLE T ...                 
C                 YES - SET INFO(4)=0                                           
C                  NO - SET INFO(4)=1                                           
C                       AND DEFINE THE STOPPING POINT TSTOP BY                  
C                       SETTING RWORK(1)=TSTOP ****                             
C                                                                               
C       INFO(5) - TO SOLVE DIFFERENTIAL/ALGEBRAIC PROBLEMS IT IS                
C              NECESSARY TO USE A MATRIX OF PARTIAL DERIVATIVES OF THE          
C              SYSTEM OF DIFFERENTIAL EQUATIONS.  IF YOU DO NOT                 
C              PROVIDE A SUBROUTINE TO EVALUATE IT ANALYTICALLY (SEE            
C              DESCRIPTION OF THE ITEM JAC IN THE CALL LIST), IT WILL           
C              BE APPROXIMATED BY NUMERICAL DIFFERENCING IN THIS CODE.          
C              ALTHOUGH IT IS LESS TROUBLE FOR YOU TO HAVE THE CODE             
C              COMPUTE PARTIAL DERIVATIVES BY NUMERICAL DIFFERENCING,           
C              THE SOLUTION WILL BE MORE RELIABLE IF YOU PROVIDE THE            
C              DERIVATIVES VIA JAC. SOMETIMES NUMERICAL DIFFERENCING            
C              IS CHEAPER THAN EVALUATING DERIVATIVES IN JAC AND                
C              SOMETIMES IT IS NOT - THIS DEPENDS ON YOUR PROBLEM.              
C                                                                               
C           **** DO YOU WANT THE CODE TO EVALUATE THE PARTIAL                   
C                  DERIVATIVES AUTOMATICALLY BY NUMERICAL DIFFERENCES ..        
C                   YES - SET INFO(5)=0                                         
C                    NO - SET INFO(5)=1  (NOT PROVIDED)                         
C                                                                               
C                                                                               
C                                                                               
C       INFO(6) - ODASSL MAKES NOT USE OF THE SPECIAL SPARSITY                  
C              STRUCTURE OF THE JACOBIAN. THIS INPUT VALUE IS                   
C              MEANINGLESS TO ODASSL.                                           
C                                                                               
C                                                                               
C                                                                               
C        INFO(7) -- YOU CAN SPECIFY A MAXIMUM (ABSOLUTE VALUE OF)               
C              STEPSIZE, SO THAT THE CODE                                       
C              WILL AVOID PASSING OVER VERY                                     
C              LARGE REGIONS.                                                   
C                                                                               
C          ****  DO YOU WANT THE CODE TO DECIDE                                 
C                ON ITS OWN MAXIMUM STEPSIZE?                                   
C                YES - SET INFO(7)=0                                            
C                 NO - SET INFO(7)=1                                            
C                      AND DEFINE HMAX BY SETTING                               
C                      RWORK(2)=HMAX ****                                       
C                                                                               
C        INFO(8) -- DIFFERENTIAL/ALGEBRAIC PROBLEMS                             
C              MAY OCCAISIONALLY SUFFER FROM                                    
C              SEVERE SCALING DIFFICULTIES ON THE                               
C              FIRST STEP. IF YOU KNOW A GREAT DEAL                             
C              ABOUT THE SCALING OF YOUR PROBLEM, YOU CAN                       
C              HELP TO ALLEVIATE THIS PROBLEM BY                                
C              SPECIFYING AN INITIAL STEPSIZE HO.                               
C                                                                               
C          ****  DO YOU WANT THE CODE TO DEFINE                                 
C                ITS OWN INITIAL STEPSIZE?                                      
C                YES - SET INFO(8)=0                                            
C                 NO - SET INFO(8)=1                                            
C                      AND DEFINE HO BY SETTING                                 
C                      RWORK(3)=HO ****                                         
C                                                                               
C        INFO(9) -- IF STORAGE IS A SEVERE PROBLEM,                             
C              YOU CAN SAVE SOME LOCATIONS BY                                   
C              RESTRICTING THE MAXIMUM ORDER MAXORD.                            
C              THE DEFAULT VALUE IS 5. FOR EACH                                 
C              ORDER DECREASE BELOW 5, THE CODE                                 
C              REQUIRES NEQ FEWER LOCATIONS, HOWEVER                            
C              IT IS LIKELY TO BE SLOWER. IN ANY                                
C              CASE, YOU MUST HAVE 1 .LE. MAXORD .LE. 5                         
C          ****  DO YOU WANT THE MAXIMUM ORDER TO                               
C                DEFAULT TO 5?                                                  
C                YES - SET INFO(9)=0                                            
C                 NO - SET INFO(9)=1                                            
C                      AND DEFINE MAXORD BY SETTING                             
C                      IWORK(3)=MAXORD ****                                     
C                                                                               
C        INFO(10) --IF YOU KNOW THAT THE SOLUTIONS TO YOUR EQUATIONS WIL        
C               ALWAYS BE NONNEGATIVE, IT MAY HELP TO SET THIS                  
C               PARAMETER.  HOWEVER, IT IS PROBABLY BEST TO                     
C               TRY THE CODE WITHOUT USING THIS OPTION FIRST,                   
C               AND ONLY TO USE THIS OPTION IF THAT DOESN'T                     
C               WORK VERY WELL.                                                 
C           ****  DO YOU WANT THE CODE TO SOLVE THE PROBLEM WITHOUT             
C                 INVOKING ANY SPECIAL NONNEGATIVITY CONSTRAINTS?               
C                  YES - SET INFO(10)=0                                         
C                   NO - SET INFO(10)=1                                         
C                                                                               
C        INFO(11) --ODASSL NORMALLY REQUIRES THE INITIAL T,                     
C               Y, AND YPRIME TO BE CONSISTENT.  THAT IS,                       
C               YOU MUST HAVE G(T,Y,YPRIME) = 0 AT THE INITIAL                  
C               TIME.  IF YOU DO NOT KNOW THE INITIAL                           
C               DERIVATIVE PRECISELY, YOU CAN LET ODASSL TRY                    
C               TO COMPUTE IT.                                                  
C          ****   ARE THE INITIAL T, Y, YPRIME CONSISTENT?                      
C                 YES - SET INFO(11) = 0                                        
C                  NO - (NOT PROVIDED BY ODASSL)                                
C                                                                               
C        INFO(12) -- FOR BETTER RESULTS THE ROWS OF THE ITERATION               
C               MATRIX CORRESPONDING TO THE DIFFERENTIAL PART                   
C               SHOULD BE SCALED.                                               
C          ****   DO YOU WANT SCALING?                                          
C                  NO - SET INFO(12) = 0                                        
C                 YES - SET INFO(12) = 1,                                       
C                       IF INFO(12) = 1 THEN                                    
C                       EQUATION ML UP TO EQUATION MU ARE SCALED.               
C                       YOU MUST SPECIFY ML AND MU IN THE                       
C                       IWORK - ARRAY:                                          
C                       IWORK(1) = ML                                           
C                       IWORK(2) = MU                                           
C                                                                               
C   RTOL, ATOL -- YOU MUST ASSIGN RELATIVE (RTOL) AND ABSOLUTE (ATOL            
C               ERROR TOLERANCES TO TELL THE CODE HOW ACCURATELY YOU WAN        
C               THE SOLUTION TO BE COMPUTED. THEY MUST BE DEFINED AS            
C               VARIABLES BECAUSE THE CODE MAY CHANGE THEM. YOU HAVE TWO        
C               CHOICES --                                                      
C                     BOTH RTOL AND ATOL ARE SCALARS. (INFO(2)=0)               
C                     BOTH RTOL AND ATOL ARE VECTORS. (INFO(2)=1)               
C               IN EITHER CASE ALL COMPONENTS MUST BE NON-NEGATIVE.             
C                                                                               
C               THE TOLERANCES ARE USED BY THE CODE IN A LOCAL ERROR TES        
C               AT EACH STEP WHICH REQUIRES ROUGHLY THAT                        
C                     ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL                    
C               FOR EACH VECTOR COMPONENT.                                      
C               (MORE SPECIFICALLY, A ROOT-MEAN-SQUARE NORM IS USED TO          
C               MEASURE THE SIZE OF VECTORS, AND THE ERROR TEST USES THE        
C               MAGNITUDE OF THE SOLUTION AT THE BEGINNING OF THE STEP.)        
C                                                                               
C               THE TRUE (GLOBAL) ERROR IS THE DIFFERENCE BETWEEN THE TR        
C               SOLUTION OF THE INITIAL VALUE PROBLEM AND THE COMPUTED          
C               APPROXIMATION. PRACTICALLY ALL PRESENT DAY CODES.               
C               INCLUDING THIS ONE, CONTROL THE LOCAL ERROR AT EACH STEP        
C               AND DO NOT EVEN ATTEMPT TO CONTROL THE GLOBAL ERROR             
C               DIRECTLY.                                                       
C               USUALLY, BUT NOT ALWAYS, THE TRUE ACCURACY OF                   
C               THE COMPUTED Y IS COMPARABLE TO THE ERROR TOLERANCES. TH        
C               CODE WILL USUALLY, BUT NOT ALWAYS, DELIVER A MORE ACCURA        
C               SOLUTION IF YOU REDUCE THE TOLERANCES AND INTEGRATE AGAI        
C               BY COMPARING TWO SUCH SOLUTIONS YOU CAN GET A FAIRLY            
C               RELIABLE IDEA OF THE TRUE ERROR IN THE SOLUTION AT THE          
C               BIGGER TOLERANCES.                                              
C                                                                               
C               SETTING ATOL=0. RESULTS IN A PURE RELATIVE ERROR TEST ON        
C               THAT COMPONENT. SETTING RTOL=0. RESULTS IN A PURE ABSOLU        
C               ERROR TEST ON THAT COMPONENT. A MIXED TEST WITH NON-ZERO        
C               RTOL AND ATOL CORRESPONDS ROUGHLY TO A RELATIVE ERROR           
C               TEST WHEN THE SOLUTION COMPONENT IS MUCH BIGGER THAN ATO        
C               AND TO AN ABSOLUTE ERROR TEST WHEN THE SOLUTION COMPONEN        
C               IS SMALLER THAN THE THRESHOLD ATOL.                             
C                                                                               
C               THE CODE WILL NOT ATTEMPT TO COMPUTE A SOLUTION AT AN           
C               ACCURACY UNREASONABLE FOR THE MACHINE BEING USED. IT WIL        
C               ADVISE YOU IF YOU ASK FOR TOO MUCH ACCURACY AND INFORM          
C               YOU AS TO THE MAXIMUM ACCURACY IT BELIEVES POSSIBLE.            
C                                                                               
C  RWORK(*) -- DIMENSION THIS REAL WORK ARRAY OF LENGTH LRW IN YOUR             
C               CALLING PROGRAM.                                                
C                                                                               
C  LRW -- SET IT TO THE DECLARED LENGTH OF THE RWORK ARRAY.                     
C               YOU MUST HAVE                                                   
C                    LRW .GE. 40+(MAXORD+3)*NY + NEQ**2 + 3*NEQ                                                      
C                                                                               
C  IWORK(*) -- DIMENSION THIS INTEGER WORK ARRAY OF LENGTH LIW IN               
C             YOUR CALLING PROGRAM.                                             
C                                                                               
C  LIW -- SET IT TO THE DECLARED LENGTH OF THE IWORK ARRAY.                     
C               YOU MUST HAVE LIW .GE. 22+NEQ                                   
C                                                                               
C  RPAR, IPAR -- THESE ARE PARAMETER ARRAYS, OF REAL AND INTEGER                
C               TYPE, RESPECTIVELY. YOU CAN USE THEM FOR COMMUNICATION          
C               BETWEEN YOUR PROGRAM THAT CALLS ODASSL AND THE                  
C               RES SUBROUTINE (AND THE JAC SUBROUTINE). THEY ARE NOT           
C               ALTERED BY ODASSL. IF YOU DO NOT NEED RPAR OR IPAR, IGNO        
C               THESE PARAMETERS BY TREATING THEM AS DUMMY ARGUMENTS. IF        
C               YOU DO CHOOSE TO USE THEM, DIMENSION THEM IN YOUR CALLIN        
C               PROGRAM AND IN RES (AND IN JAC) AS ARRAYS OF APPROPRIATE        
C               LENGTH.                                                         
C                                                                               
C  JAC -- THIS PARAMETER PARAMETER SHOULD BE IGNORED AND TREATED AS             
C               A DUMMY ARGUMENT.                                               
C                                                                               
C                                                                               
C                                                                               
C  OPTIONALLY REPLACEABLE NORM ROUTINE:                                         
C  ODASSL USES A WEIGHTED NORM DDANRM TO MEASURE THE SIZE                       
C  OF VECTORS SUCH AS THE ESTIMATED ERROR IN EACH STEP.                         
C  A FUNCTION SUBPROGRAM                                                        
C    DOUBLE PRECISION FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR)                       
C    DIMENSION V(NEQ),WT(NEQ)                                                   
C  IS USED TO DEFINE THIS NORM.  HERE, V IS THE VECTOR                          
C  WHOSE NORM IS TO BE COMPUTED, AND WT IS A VECTOR OF                          
C  WEIGHTS.  A DDANRM ROUTINE HAS BEEN INCLUDED WITH ODASSL                     
C  WHICH COMPUTES THE WEIGHTED ROOT-MEAN-SQUARE NORM                            
C  GIVEN BY                                                                     
C    DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)                                    
C  THIS NORM IS SUITABLE FOR MOST PROBLEMS.  IN SOME                            
C  SPECIAL CASES, IT MAY BE MORE CONVENIENT AND/OR                              
C  EFFICIENT TO DEFINE YOUR OWN NORM BY WRITING A FUNCTION                      
C  SUBPROGRAM TO BE CALLED INSTEAD OF DDANRM.  THIS SHOULD                      
C  HOWEVER, BE ATTEMPTED ONLY AFTER CAREFUL THOUGHT AND                         
C  CONSIDERATION.                                                               
C                                                                               
C                                                                               
C -----OUTPUT-AFTER ANY RETURN FROM ODASSL----                                  
C                                                                               
C  THE PRINCIPAL AIM OF THE CODE IS TO RETURN A COMPUTED SOLUTION AT            
C  TOUT, ALTHOUGH IT IS ALSO POSSIBLE TO OBTAIN INTERMEDIATE RESULTS            
C  ALONG THE WAY. TO FIND OUT WHETHER THE CODE ACHIEVED ITS GOAL                
C  OR IF THE INTEGRATION PROCESS WAS INTERRUPTED BEFORE THE TASK WAS            
C  COMPLETED, YOU MUST CHECK THE IDID PARAMETER.                                
C                                                                               
C                                                                               
C   T -- THE SOLUTION WAS SUCCESSFULLY ADVANCED TO THE                          
C               OUTPUT VALUE OF T.                                              
C                                                                               
C   Y(*) -- CONTAINS THE COMPUTED SOLUTION APPROXIMATION AT T.                  
C                                                                               
C   YPRIME(*) -- CONTAINS THE COMPUTED DERIVATIVE                               
C               APPROXIMATION AT T                                              
C                                                                               
C   IDID -- REPORTS WHAT THE CODE DID                                           
C                                                                               
C                     *** TASK COMPLETED ***                                    
C                REPORTED BY POSITIVE VALUES OF IDID                            
C                                                                               
C           IDID = 1 -- A STEP WAS SUCCESSFULLY TAKEN IN THE                    
C                   INTERMEDIATE-OUTPUT MODE. THE CODE HAS NOT                  
C                   YET REACHED TOUT.                                           
C                                                                               
C           IDID = 2 -- THE INTEGRATION TO TOUT WAS SUCCESSFULLY                
C                   COMPLETED (T=TOUT) BY STEPPING EXACTLY TO TOUT.             
C                                                                               
C           IDID = 3 -- THE INTEGRATION TO TOUT WAS SUCCESSFULLY                
C                   COMPLETED (T=TOUT) BY STEPPING PAST TOUT.                   
C                   Y(*) IS OBTAINED BY INTERPOLATION.                          
C                   YPRIME(*) IS OBTAINED BY INTERPOLATION.                     
C                                                                               
C                    *** TASK INTERRUPTED ***                                   
C                REPORTED BY NEGATIVE VALUES OF IDID                            
C                                                                               
C           IDID = -1 -- A LARGE AMOUNT OF WORK HAS BEEN EXPENDED.              
C                   (ABOUT 500 STEPS)                                           
C                                                                               
C           IDID = -2 -- THE ERROR TOLERANCES ARE TOO STRINGENT.                
C                                                                               
C           IDID = -3 -- THE LOCAL ERROR TEST CANNOT BE SATISFIED               
C                   BECAUSE YOU SPECIFIED A ZERO COMPONENT IN ATOL              
C                   AND THE CORRESPONDING COMPUTED SOLUTION                     
C                   COMPONENT IS ZERO. THUS, A PURE RELATIVE ERROR              
C                   TEST IS IMPOSSIBLE FOR THIS COMPONENT.                      
C                                                                               
C           IDID = -6 -- ODASSL HAD REPEATED ERROR TEST                         
C                   FAILURES ON THE LAST ATTEMPTED STEP.                        
C                                                                               
C           IDID = -7 -- THE CORRECTOR COULD NOT CONVERGE.                      
C                                                                               
C           IDID = -8 -- THE MATRIX OF PARTIAL DERIVATIVES                      
C                   IS SINGULAR.                                                
C                                                                               
C           IDID = -9 -- THE CORRECTOR COULD NOT CONVERGE.                      
C                   THERE WERE REPEATED ERROR TEST FAILURES                     
C                   IN THIS STEP.                                               
C                                                                               
C           IDID =-10 -- THE CORRECTOR COULD NOT CONVERGE                       
C                   BECAUSE IRES WAS EQUAL TO MINUS ONE.                        
C                                                                               
C           IDID =-11 -- IRES EQUAL TO -2 WAS ENCOUNTERED                       
C                   AND CONTROL IS BEING RETURNED TO THE                        
C                   CALLING PROGRAM.                                            
C                                                                               
C           IDID =-12 -- ODASSL FAILED TO COMPUTE THE INITIAL                   
C                   YPRIME.                                                     
C                                                                               
C                                                                               
C                                                                               
C           IDID = -13,..,-32 -- NOT APPLICABLE FOR THIS CODE                   
C                                                                               
C                    *** TASK TERMINATED ***                                    
C                REPORTED BY THE VALUE OF IDID=-33                              
C                                                                               
C           IDID = -33 -- THE CODE HAS ENCOUNTERED TROUBLE FROM WHICH           
C                   IT CANNOT RECOVER. A MESSAGE IS PRINTED                     
C                   EXPLAINING THE TROUBLE AND CONTROL IS RETURNED              
C                   TO THE CALLING PROGRAM. FOR EXAMPLE, THIS OCCURS            
C                   WHEN INVALID INPUT IS DETECTED.                             
C                                                                               
C   RTOL, ATOL -- THESE QUANTITIES REMAIN UNCHANGED EXCEPT WHEN                 
C               IDID = -2. IN THIS CASE, THE ERROR TOLERANCES HAVE BEEN         
C               INCREASED BY THE CODE TO VALUES WHICH ARE ESTIMATED TO B        
C               APPROPRIATE FOR CONTINUING THE INTEGRATION. HOWEVER, THE        
C               REPORTED SOLUTION AT T WAS OBTAINED USING THE INPUT VALU        
C               OF RTOL AND ATOL.                                               
C                                                                               
C   RWORK, IWORK -- CONTAIN INFORMATION WHICH IS USUALLY OF NO                  
C               INTEREST TO THE USER BUT NECESSARY FOR SUBSEQUENT CALLS.        
C               HOWEVER, YOU MAY FIND USE FOR                                   
C                                                                               
C               RWORK(3)--WHICH CONTAINS THE STEP SIZE H TO BE                  
C                       ATTEMPTED ON THE NEXT STEP.                             
C                                                                               
C               RWORK(4)--WHICH CONTAINS THE CURRENT VALUE OF THE               
C                       INDEPENDENT VARIABLE, I.E. THE FARTHEST POINT           
C                       INTEGRATION HAS REACHED. THIS WILL BE DIFFERENT         
C                       FROM T ONLY WHEN INTERPOLATION HAS BEEN                 
C                       PERFORMED (IDID=3).                                     
C                                                                               
C               RWORK(7)--WHICH CONTAINS THE STEPSIZE USED                      
C                       ON THE LAST SUCCESSFUL STEP.                            
C                                                                               
C               IWORK(7)--WHICH CONTAINS THE ORDER OF THE METHOD TO             
C                       BE ATTEMPTED ON THE NEXT STEP.                          
C                                                                               
C               IWORK(8)--WHICH CONTAINS THE ORDER OF THE METHOD USED           
C                       ON THE LAST STEP.                                       
C                                                                               
C               IWORK(11)--WHICH CONTAINS THE NUMBER OF STEPS TAKEN SO F        
C                                                                               
C               IWORK(12)--WHICH CONTAINS THE NUMBER OF CALLS TO RES            
C                        SO FAR.                                                
C                                                                               
C               IWORK(13)--WHICH CONTAINS THE NUMBER OF EVALUATIONS OF          
C                        THE MATRIX OF PARTIAL DERIVATIVES NEEDED SO FAR        
C                                                                               
C               IWORK(14)--WHICH CONTAINS THE TOTAL NUMBER                      
C                        OF ERROR TEST FAILURES SO FAR.                         
C                                                                               
C               IWORK(15)--WHICH CONTAINS THE TOTAL NUMBER                      
C                        OF CONVERGENCE TEST FAILURES SO FAR.                   
C                        (INCLUDES SINGULAR ITERATION MATRIX                    
C                        FAILURES.)                                             
C                                                                               
C                                                                               
C                                                                               
C   INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION                             
C            (CALLS AFTER THE FIRST)                **                          
C                                                                               
C     THIS CODE IS ORGANIZED SO THAT SUBSEQUENT CALLS TO CONTINUE THE           
C     INTEGRATION INVOLVE LITTLE (IF ANY) ADDITIONAL EFFORT ON YOUR             
C     PART. YOU MUST MONITOR THE IDID PARAMETER IN ORDER TO DETERMINE           
C     WHAT TO DO NEXT.                                                          
C                                                                               
C     RECALLING THAT THE PRINCIPAL TASK OF THE CODE IS TO INTEGRATE             
C     FROM T TO TOUT (THE INTERVAL MODE), USUALLY ALL YOU WILL NEED             
C     TO DO IS SPECIFY A NEW TOUT UPON REACHING THE CURRENT TOUT.               
C                                                                               
C     DO NOT ALTER ANY QUANTITY NOT SPECIFICALLY PERMITTED BELOW,               
C     IN PARTICULAR DO NOT ALTER NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*)         
C     OR THE DIFFERENTIAL EQUATION IN SUBROUTINE RES. ANY SUCH                  
C     ALTERATION CONSTITUTES A NEW PROBLEM AND MUST BE TREATED AS SUCH,         
C     I.E. YOU MUST START AFRESH.                                               
C                                                                               
C     YOU CANNOT CHANGE FROM VECTOR TO SCALAR ERROR CONTROL OR VICE             
C     VERSA (INFO(2)) BUT YOU CAN CHANGE THE SIZE OF THE ENTRIES OF             
C     RTOL, ATOL. INCREASING A TOLERANCE MAKES THE EQUATION EASIER              
C     TO INTEGRATE. DECREASING A TOLERANCE WILL MAKE THE EQUATION               
C     HARDER TO INTEGRATE AND SHOULD GENERALLY BE AVOIDED.                      
C                                                                               
C     YOU CAN SWITCH FROM THE INTERMEDIATE-OUTPUT MODE TO THE                   
C     INTERVAL MODE (INFO(3)) OR VICE VERSA AT ANY TIME.                        
C                                                                               
C     IF IT HAS BEEN NECESSARY TO PREVENT THE INTEGRATION FROM GOING            
C     PAST A POINT TSTOP (INFO(4), RWORK(1)), KEEP IN MIND THAT THE             
C     CODE WILL NOT INTEGRATE TO ANY TOUT BEYOUND THE CURRENTLY                 
C     SPECIFIED TSTOP. ONCE TSTOP HAS BEEN REACHED YOU MUST CHANGE              
C     THE VALUE OF TSTOP OR SET INFO(4)=0. YOU MAY CHANGE INFO(4)               
C     OR TSTOP AT ANY TIME BUT YOU MUST SUPPLY THE VALUE OF TSTOP IN            
C     RWORK(1) WHENEVER YOU SET INFO(4)=1.                                      
C                                                                               
C     DO NOT CHANGE INFO(5), INFO(6), IWORK(1), OR IWORK(2)                     
C     UNLESS YOU ARE GOING TO RESTART THE CODE.                                 
C                                                                               
C                    *** FOLLOWING A COMPLETED TASK ***                         
C     IF                                                                        
C     IDID = 1, CALL THE CODE AGAIN TO CONTINUE THE INTEGRATION                 
C                  ANOTHER STEP IN THE DIRECTION OF TOUT.                       
C                                                                               
C     IDID = 2 OR 3, DEFINE A NEW TOUT AND CALL THE CODE AGAIN.                 
C                  TOUT MUST BE DIFFERENT FROM T. YOU CANNOT CHANGE             
C                  THE DIRECTION OF INTEGRATION WITHOUT RESTARTING.             
C                                                                               
C                    *** FOLLOWING AN INTERRUPTED TASK ***                      
C                  TO SHOW THE CODE THAT YOU REALIZE THE TASK WAS               
C                  INTERRUPTED AND THAT YOU WANT TO CONTINUE, YOU               
C                  MUST TAKE APPROPRIATE ACTION AND SET INFO(1) = 1             
C     IF                                                                        
C     IDID = -1, THE CODE HAS TAKEN ABOUT 500 STEPS.                            
C                  IF YOU WANT TO CONTINUE, SET INFO(1) = 1 AND                 
C                  CALL THE CODE AGAIN. AN ADDITIONAL 500 STEPS                 
C                  WILL BE ALLOWED.                                             
C                                                                               
C                                                                               
C     IDID = -2, THE ERROR TOLERANCES RTOL, ATOL HAVE BEEN                      
C                  INCREASED TO VALUES THE CODE ESTIMATES APPROPRIATE           
C                  FOR CONTINUING. YOU MAY WANT TO CHANGE THEM                  
C                  YOURSELF. IF YOU ARE SURE YOU WANT TO CONTINUE               
C                  WITH RELAXED ERROR TOLERANCES, SET INFO(1)=1 AND             
C                  CALL THE CODE AGAIN.                                         
C                                                                               
C     IDID = -3, A SOLUTION COMPONENT IS ZERO AND YOU SET THE                   
C                  CORRESPONDING COMPONENT OF ATOL TO ZERO. IF YOU              
C                  ARE SURE YOU WANT TO CONTINUE, YOU MUST FIRST                
C                  ALTER THE ERROR CRITERION TO USE POSITIVE VALUES             
C                  FOR THOSE COMPONENTS OF ATOL CORRESPONDING TO ZERO           
C                  SOLUTION COMPONENTS, THEN SET INFO(1)=1 AND CALL             
C                  THE CODE AGAIN.                                              
C                                                                               
C     IDID = -4,-5  --- CANNOT OCCUR WITH THIS CODE                             
C                                                                               
C     IDID = -6, REPEATED ERROR TEST FAILURES OCCURRED ON THE                   
C                  LAST ATTEMPTED STEP IN ODASSL. A SINGULARITY IN THE          
C                  SOLUTION MAY BE PRESENT. IF YOU ARE ABSOLUTELY               
C                  CERTAIN YOU WANT TO CONTINUE, YOU SHOULD RESTART             
C                  THE INTEGRATION.(PROVIDE INITIAL VALUES OF Y AND             
C                  YPRIME WHICH ARE CONSISTENT)                                 
C                                                                               
C     IDID = -7, REPEATED CONVERGENCE TEST FAILURES OCCURRED                    
C                  ON THE LAST ATTEMPTED STEP IN ODASSL. AN INACCURATE O        
C                  ILLCONDITIONED JACOBIAN MAY BE THE PROBLEM. IF YOU           
C                  ARE ABSOLUTELY CERTAIN YOU WANT TO CONTINUE, YOU             
C                  SHOULD RESTART THE INTEGRATION.                              
C                                                                               
C     IDID = -8, THE MATRIX OF PARTIAL DERIVATIVES IS SINGULAR.                 
C                  SOME OF YOUR EQUATIONS MAY BE REDUNDANT.                     
C                  ODASSL CANNOT SOLVE THE PROBLEM AS STATED.                   
C                  IT IS POSSIBLE THAT THE REDUNDANT EQUATIONS                  
C                  COULD BE REMOVED, AND THEN ODASSL COULD                      
C                  SOLVE THE PROBLEM. IT IS ALSO POSSIBLE                       
C                  THAT A SOLUTION TO YOUR PROBLEM EITHER                       
C                  DOES NOT EXIST OR IS NOT UNIQUE.                             
C                                                                               
C     IDID = -9, ODASSL HAD MULTIPLE CONVERGENCE TEST                           
C                  FAILURES, PRECEEDED BY MULTIPLE ERROR                        
C                  TEST FAILURES, ON THE LAST ATTEMPTED STEP.                   
C                  IT IS POSSIBLE THAT YOUR PROBLEM                             
C                  IS ILL-POSED, AND CANNOT BE SOLVED                           
C                  USING THIS CODE.  OR, THERE MAY BE A                         
C                  DISCONTINUITY OR A SINGULARITY IN THE                        
C                  SOLUTION.  IF YOU ARE ABSOLUTELY CERTAIN                     
C                  YOU WANT TO CONTINUE, YOU SHOULD RESTART                     
C                  THE INTEGRATION.                                             
C                                                                               
C    IDID =-10, ODASSL HAD MULTIPLE CONVERGENCE TEST FAILURES                   
C                  BECAUSE IRES WAS EQUAL TO MINUS ONE.                         
C                  IF YOU ARE ABSOLUTELY CERTAIN YOU WANT                       
C                  TO CONTINUE, YOU SHOULD RESTART THE                          
C                  INTEGRATION.                                                 
C                                                                               
C    IDID =-11, IRES=-2 WAS ENCOUNTERED, AND CONTROL IS BEING                   
C                  RETURNED TO THE CALLING PROGRAM.                             
C                                                                               
C    IDID =-12, ODASSL FAILED TO COMPUTE THE INITIAL YPRIME.                    
C               THIS COULD HAPPEN BECAUSE THE INITIAL                           
C               APPROXIMATION TO YPRIME WAS NOT VERY GOOD, OR                   
C               IF A YPRIME CONSISTENT WITH THE INITIAL Y                       
C               DOES NOT EXIST.  THE PROBLEM COULD ALSO BE CAUSED               
C               BY AN INACCURATE OR SINGULAR ITERATION MATRIX.                  
C                                                                               
C                                                                               
C                                                                               
C     IDID = -13,..,-32 --- CANNOT OCCUR WITH THIS CODE                         
C                                                                               
C                       *** FOLLOWING A TERMINATED TASK ***                     
C     IF IDID= -33, YOU CANNOT CONTINUE THE SOLUTION OF THIS                    
C                  PROBLEM. AN ATTEMPT TO DO SO WILL RESULT IN YOUR             
C                  RUN BEING TERMINATED.                                        
C                                                                               
C  ---------------------------------------------------------------------        
C                                                                               
C***REFERENCES  A DESCRIPTION OF DASSL: A DIFFERENTIAL/ALGEBRAIC                
C                  SYSTEM SOLVER, L. R. PETZOLD, SAND82-8637,                   
C                  SANDIA NATIONAL LABORATORIES, SEPTEMBER 1982.                
C               NUMERICAL SOLUTION OF DIFFERENTIAL-ALGEBRAIC EQUATIONS                 
C                  FOR CONSTRAINED MECHANICAL MOTION                           
C                  C. FUEHRER AND BENEDICT LEIMKUHLER                  
C                  SUBMITTED FOR PUBLICATION TO NUMERISCHE MATHEMATIK 1989            
C                                                                               
C***ROUTINES CALLED  ODASTP,DDAINI,DDANRM,DDAWTS,DDATRP,XERRWV,D1MACH           
C***COMMON BLOCKS    ODA001                                                     
C***END PROLOGUE ODASSL                                                         
C                                                                               
C                                                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      LOGICAL DONE                                                              
      EXTERNAL RES,JAC                                                          
      DIMENSION Y(*),YPRIME(*)                                                  
      DIMENSION INFO(15)                                                        
      DIMENSION RWORK(*),IWORK(*)                                               
      DIMENSION RTOL(*),ATOL(*)                                                 
      DIMENSION RPAR(*),IPAR(*)                                                 
      COMMON/ODA001/NPD,NTEMP,INORM,                                            
     *   LML,LMU,LMXORD,LNY,LM1,LMTYPE,                                         
     *   LNST,LNRE,LNJE,LETF,LCTF,LIPVT 
      SAVE  
      DATA LTSTOP,LHMAX,LH,LTN,                                                 
     *   LCJ,LCJOLD,LHOLD,LS,LROUND,                                            
     *   LALPHA,LBETA,LGAMMA,                                                   
     *   LPSI,LSIGMA,LDELTA                                                     
     *   /1,2,3,4,                                                              
     *   5,6,7,8,9,                                                             
     *   11,17,23,                                                              
     *   29,35,41/  
                                                                     
      IF(INFO(1).NE.0)GO TO 100                                                 
C                                                                               
C-----------------------------------------------------------------------        
C     THIS BLOCK IS EXECUTED FOR THE INITIAL CALL ONLY.                         
C     IT CONTAINS CHECKING OF INPUTS AND INITIALIZATIONS.                       
C-----------------------------------------------------------------------        
C                                                                               
C     FIRST CHECK INFO ARRAY TO MAKE SURE ALL ELEMENTS OF INFO                  
C     ARE EITHER ZERO OR ONE.                                                   
C               
                                                     
      DO 10 I=2,12                                                              
         IF(INFO(I).NE.0.AND.INFO(I).NE.1) GO TO 701                            
10       CONTINUE                                                               
C                                                                               
C                                                                               
C     SET POINTERS INTO IWORK                                                   
      LML=1                                                                     
      LMU=2                                                                     
      LMXORD=3                                                                  
      LNY=4                                                                     
      LM1=5                                                                     
      LMTYPE= 6                                                                 
      LJCALC=17                                                                 
      LPHASE=16                                                                 
      LK=7                                                                      
      LKOLD=8                                                                   
      LNS=9                                                                     
      LNSTL=10                                                                  
      LNST=11                                                                   
      LNRE=12                                                                   
      LNJE=13                                                                   
      LETF=14                                                                   
      LCTF=15                                                                   
      LIPVT=23                                                                  
      LIWM=1                                                                    
C                                                                               
C     CHECK DIMENSIONS                                                  
C                                                                               
C                                                                               
      IF(NEQ.LE.0 .OR. NY .LE. 0 .OR. NY .GT. NEQ) GO TO 702                    
      MXORD=5                                                                   
      IF(INFO(9).EQ.0)GO TO 20                                                  
         MXORD=IWORK(LMXORD)                                                    
         IF(MXORD.LT.1.OR.MXORD.GT.5)GO TO 703                                  
20       IWORK(LMXORD)=MXORD                                                    
C                                                                               
C     COMPUTE MTYPE,LENPD,LENRW.CHECK ML AND MU.   
      IF (INFO(5) .EQ. 0) THEN
         IWORK(LMTYPE) = 2
      ELSE
         IWORK(LMTYPE) = 1
      END IF                             
C     IF(INFO(6).NE.0)GO TO 40  
C        LENCOR: WORKING SPACE REQUIRED BY ODACOR
         LENCOR = NEQ**2 + 2*NEQ                                                         
C                                                                               
         LENRW=40+(IWORK(LMXORD)+3)*NY+NEQ+LENCOR
C                                                                               
C                                                                               
40    INORM = INFO(12)                                                          
      IF(INORM.NE.1) GO TO 60                                                   
      IF(IWORK(LML).LE.0.OR.IWORK(LMU).GT.NEQ) GO TO 717                        
      IF(IWORK(LMU).LT.IWORK(LML)) GO TO 718                                    
C                                                                               
C     CHECK LENGTHS OF RWORK AND IWORK                                          
60    LENIW=22+NEQ                                                              
      IF(LRW.LT.LENRW)GO TO 704                                                 
      IF(LIW.LT.LENIW)GO TO 705                                                 
C                                                                               
C     CHECK TO SEE THAT TOUT IS DIFFERENT FROM T                                
      IF(TOUT .EQ. T)GO TO 719                                                  
C                                                                               
C     CHECK HMAX                                                                
      IF(INFO(7).EQ.0)GO TO 70                                                  
         HMAX=RWORK(LHMAX)                                                      
         IF(HMAX.LE.0.0D0)GO TO 710                                             
70    CONTINUE                                                                  
C                                                                               
C     INITIALIZE COUNTERS                                                       
      IWORK(LNST)=0                                                             
      IWORK(LNRE)=0                                                             
      IWORK(LNJE)=0                                                             
C                                                                               
      IWORK(LNSTL)=0                                                            
      IDID=1                                                                    
      GO TO 200                                                                 
C                                                                               
C ----------------------------------------------------------------------        
C     THIS BLOCK IS FOR CONTINUATION CALLS                                      
C     ONLY. HERE WE CHECK INFO(1),AND IF THE                                    
C     LAST STEP WAS INTERRUPTED WE CHECK WHETHER                                
C     APPROPRIATE ACTION WAS TAKEN.                                             
C-----------------------------------------------------------------------        
C                                                                               
100   CONTINUE                                                                  
      IF(INFO(1).EQ.1)GO TO 110                                                 
      IF(INFO(1).NE.-1)GO TO 701                                                
C     IF WE ARE HERE, THE LAST STEP WAS INTERRUPTED                             
C     BY AN ERROR CONDITION FROM ODASTP,AND                                     
C     APPROPRIATE ACTION WAS NOT TAKEN. THIS                                    
C     IS A FATAL ERROR.                                                         
      CALL XERRWV(                                                              
     *49HDASSL--  THE LAST STEP TERMINATED WITH A NEGATIVE,                     
     *49,201,0,0,0,0,0,0.0D0,0.0D0)                                             
      CALL XERRWV(                                                              
     *47HDASSL--  VALUE (=I1) OF IDID AND NO APPROPRIATE,                       
     *47,202,0,1,IDID,0,0,0.0D0,0.0D0)                                          
      CALL XERRWV(                                                              
     *41HDASSL--  ACTION WAS TAKEN. RUN TERMINATED,                             
     *41,203,1,0,0,0,0,0.0D0,0.0D0)                                             
      RETURN                                                                    
110   CONTINUE                                                                  
      IWORK(LNSTL)=IWORK(LNST)                                                  
C                                                                               
C-----------------------------------------------------------------------        
C     THIS BLOCK IS EXECUTED ON ALL CALLS.                                      
C     THE ERROR TOLERANCE PARAMETERS ARE                                        
C     CHECKED, AND THE WORK ARRAY POINTERS                                      
C     ARE SET.                                                                  
C-----------------------------------------------------------------------        
C                                                                               
200   CONTINUE                                                                  
C     CHECK RTOL,ATOL                                                           
      NZFLG=0                                                                   
      RTOLI=RTOL(1)                                                             
      ATOLI=ATOL(1)                                                             
      DO 210 I=1,NY                                                             
         IF(INFO(2).EQ.1)RTOLI=RTOL(I)                                          
         IF(INFO(2).EQ.1)ATOLI=ATOL(I)                                          
         IF(RTOLI.GT.0.0D0.OR.ATOLI.GT.0.0D0)NZFLG=1                            
         IF(RTOLI.LT.0.0D0)GO TO 706                                            
         IF(ATOLI.LT.0.0D0)GO TO 707                                            
210      CONTINUE                                                               
      IF(NZFLG.EQ.0)GO TO 708                                                   
C                                                                               
C     SET UP RWORK STORAGE.IWORK STORAGE IS FIXED                               
C     IN DATA STATEMENT.                                                        
      LE=LDELTA+NEQ                                                             
      LWT=LE+NY                                                                 
      LPHI=LWT+NY                                                               
      LPD=LPHI+(IWORK(LMXORD)+1)*NY                                             
      LWM=LPD                                                                   
      NPD=1                                                                     
      NTEMP=NPD+LENCOR                                                           
      IF(INFO(1).EQ.1)GO TO 400                                                 
C                                                                               
C-----------------------------------------------------------------------        
C     THIS BLOCK IS EXECUTED ON THE INITIAL CALL                                
C     ONLY. SET THE INITIAL STEP SIZE, AND                                      
C     THE ERROR WEIGHT VECTOR, AND PHI.                                         
C     COMPUTE INITIAL YPRIME, IF NECESSARY.                                     
C-----------------------------------------------------------------------        
C                                                                               
300   CONTINUE                                                                  
      TN=T                                                                      
      IDID=1                                                                    
C                                                                               
C     SET ERROR WEIGHT VECTOR WT                                                
      CALL DDAWTS(NY,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)                  
      DO 305 I = 1,NY                                                           
         IF(RWORK(LWT+I-1).LE.0.0D0) GO TO 713                                  
305      CONTINUE                                                               
C                                                                               
C     COMPUTE UNIT ROUNDOFF AND HMIN                                            
      UROUND = D1MACH(4) 
      RWORK(LROUND) = UROUND                                                    
      HMIN = 4.0D0*UROUND*DMAX1(DABS(T),DABS(TOUT))                             
C                                                                               
C     CHECK INITIAL INTERVAL TO SEE THAT IT IS LONG ENOUGH                      
      TDIST = DABS(TOUT - T)                                                    
      IF(TDIST .LT. HMIN) GO TO 714                                             
C                                                                               
C     CHECK HO, IF THIS WAS INPUT                                               
      IF (INFO(8) .EQ. 0) GO TO 310                                             
         HO = RWORK(LH)                                                         
         IF ((TOUT - T)*HO .LT. 0.0D0) GO TO 711                                
         IF (HO .EQ. 0.0D0) GO TO 712                                           
         GO TO 320                                                              
310    CONTINUE                                                                 
C                                                                               
C     COMPUTE INITIAL STEPSIZE, TO BE USED BY EITHER                            
C     ODASTP OR DDAINI, DEPENDING ON INFO(11)                                   
      HO = 0.001D0*TDIST                                                        
      YPNORM = DDANRM(NY,YPRIME,RWORK(LWT),RPAR,IPAR)                           
      IF (YPNORM .GT. 0.5D0/HO) HO = 0.5D0/YPNORM                               
      HO = DSIGN(HO,TOUT-T)                                                     
C     ADJUST HO IF NECESSARY TO MEET HMAX BOUND                                 
320   IF (INFO(7) .EQ. 0) GO TO 330                                             
         RH = DABS(HO)/HMAX                                                     
         IF (RH .GT. 1.0D0) HO = HO/RH                                          
C     COMPUTE TSTOP, IF APPLICABLE                                              
330   IF (INFO(4) .EQ. 0) GO TO 340                                             
         TSTOP = RWORK(LTSTOP)                                                  
         IF ((TSTOP - T)*HO .LT. 0.0D0) GO TO 715                               
         IF ((T + HO - TSTOP)*HO .GT. 0.0D0) HO = TSTOP - T                     
         IF ((TSTOP - TOUT)*HO .LT. 0.0D0) GO TO 709                            
C                                                                               
C     COMPUTE INITIAL DERIVATIVE, IF APPLICABLE                                 
340   IF (INFO(11) .EQ. 0) GO TO 350                                            
C                 DAINI NOT AVAILABLE                                           
C     CALL DDAINI(T,Y,YPRIME,NEQ,                                               
C    *  RES,JAC,HO,RWORK(LWT),IDID,RPAR,IPAR,                                   
C    *  RWORK(LPHI),RWORK(LDELTA),RWORK(LE),                                    
C    *  RWORK(LWM),IWORK(LIWM),HMIN,RWORK(LROUND),INFO(10))                     
      IF (IDID .LT. 0) GO TO 390                                                
C                                                                               
C     LOAD H WITH HO.  STORE H IN RWORK(LH)                                     
350   H = HO                                                                    
      RWORK(LH) = H                                                             
C                                                                               
C     LOAD Y AND H*YPRIME INTO PHI(*,1) AND PHI(*,2)                            
360   ITEMP = LPHI + NY                                                         
      DO 370 I = 1,NY                                                           
         RWORK(LPHI + I - 1) = Y(I)                                             
370      RWORK(ITEMP + I - 1) = H*YPRIME(I)                                     
C                                                                               
390   GO TO 500                                                                 
C                                                                               
C-------------------------------------------------------                        
C     THIS BLOCK IS FOR CONTINUATION CALLS ONLY. ITS                            
C     PURPOSE IS TO CHECK STOP CONDITIONS BEFORE                                
C     TAKING A STEP.                                                            
C     ADJUST H IF NECESSARY TO MEET HMAX BOUND                                  
C-------------------------------------------------------                        
C                                                                               
400   CONTINUE                                                                  
      DONE = .FALSE.                                                            
      TN=RWORK(LTN)                                                             
      H=RWORK(LH)                                                               
      IF(INFO(7) .EQ. 0) GO TO 410                                              
         RH = DABS(H)/HMAX                                                      
         IF(RH .GT. 1.0D0) H = H/RH                                             
410   CONTINUE                                                                  
      IF(T .EQ. TOUT) GO TO 719                                                 
      IF((T - TOUT)*H .GT. 0.0D0) GO TO 711                                     
      IF(INFO(4) .EQ. 1) GO TO 430                                              
      IF(INFO(3) .EQ. 1) GO TO 420                                              
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 490                                         
      CALL DDATRP(TN,TOUT,Y,YPRIME,NY,IWORK(LKOLD),                             
     *  RWORK(LPHI),RWORK(LPSI))                                                
      T=TOUT                                                                    
      IDID = 3                                                                  
      DONE = .TRUE.                                                             
      GO TO 490                                                                 
420   IF((TN-T)*H .LE. 0.0D0) GO TO 490                                         
      IF((TN - TOUT)*H .GT. 0.0D0) GO TO 425                                    
      CALL DDATRP(TN,TN,Y,YPRIME,NY,IWORK(LKOLD),                               
     *  RWORK(LPHI),RWORK(LPSI))                                                
      T = TN                                                                    
      IDID = 1                                                                  
      DONE = .TRUE.                                                             
      GO TO 490                                                                 
425   CONTINUE                                                                  
      CALL DDATRP(TN,TOUT,Y,YPRIME,NY,IWORK(LKOLD),                             
     *  RWORK(LPHI),RWORK(LPSI))                                                
      T = TOUT                                                                  
      IDID = 3                                                                  
      DONE = .TRUE.                                                             
      GO TO 490                                                                 
430   IF(INFO(3) .EQ. 1) GO TO 440                                              
      TSTOP=RWORK(LTSTOP)                                                       
      IF((TN-TSTOP)*H.GT.0.0D0) GO TO 715                                       
      IF((TSTOP-TOUT)*H.LT.0.0D0)GO TO 709                                      
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 450                                         
      CALL DDATRP(TN,TOUT,Y,YPRIME,NY,IWORK(LKOLD),                             
     *   RWORK(LPHI),RWORK(LPSI))                                               
      T=TOUT                                                                    
      IDID = 3                                                                  
      DONE = .TRUE.                                                             
      GO TO 490                                                                 
440   TSTOP = RWORK(LTSTOP)                                                     
      IF((TN-TSTOP)*H .GT. 0.0D0) GO TO 715                                     
      IF((TSTOP-TOUT)*H .LT. 0.0D0) GO TO 709                                   
      IF((TN-T)*H .LE. 0.0D0) GO TO 450                                         
      IF((TN - TOUT)*H .GT. 0.0D0) GO TO 445                                    
      CALL DDATRP(TN,TN,Y,YPRIME,NY,IWORK(LKOLD),                               
     *  RWORK(LPHI),RWORK(LPSI))                                                
      T = TN                                                                    
      IDID = 1                                                                  
      DONE = .TRUE.                                                             
      GO TO 490                                                                 
445   CONTINUE                                                                  
      CALL DDATRP(TN,TOUT,Y,YPRIME,NY,IWORK(LKOLD),                             
     *  RWORK(LPHI),RWORK(LPSI))                                                
      T = TOUT                                                                  
      IDID = 3                                                                  
      DONE = .TRUE.                                                             
      GO TO 490                                                                 
450   CONTINUE                                                                  
C     CHECK WHETHER WE ARE WITH IN ROUNDOFF OF TSTOP                            
      IF(DABS(TN-TSTOP).GT.100.0D0*UROUND*                                      
     *   (DABS(TN)+DABS(H)))GO TO 460                                           
      IDID=2                                                                    
      T=TSTOP                                                                   
      DONE = .TRUE.                                                             
      GO TO 490                                                                 
460   TNEXT=TN+H*(1.0D0+4.0D0*UROUND)                                           
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 490                                     
      H=(TSTOP-TN)*(1.0D0-4.0D0*UROUND)                                         
      RWORK(LH)=H                                                               
C                                                                               
490   IF (DONE) GO TO 590                                                       
C                                                                               
C-------------------------------------------------------                        
C     THE NEXT BLOCK CONTAINS THE CALL TO THE                                   
C     ONE-STEP INTEGRATOR ODASTP.                                               
C     THIS IS A LOOPING POINT FOR THE INTEGRATION                               
C     STEPS.                                                                    
C     CHECK FOR TOO MANY STEPS.                                                 
C     UPDATE WT.                                                                
C     CHECK FOR TOO MUCH ACCURACY REQUESTED.                                    
C     COMPUTE MINIMUM STEPSIZE.                                                 
C-------------------------------------------------------                        
C                                                                               
500   CONTINUE                                                                  
C     CHECK FOR FAILURE TO COMPUTE INITIAL YPRIME                               
      IF (IDID .EQ. -12) GO TO 527                                              
C                                                                               
C     CHECK FOR TOO MANY STEPS                                                  
      IF((IWORK(LNST)-IWORK(LNSTL)).LT.5000)                                    
     *   GO TO 510                                                              
           IDID=-1                                                              
           GO TO 527                                                            
C                                                                               
C     UPDATE WT                                                                 
510   CALL DDAWTS(NY,INFO(2),RTOL,ATOL,RWORK(LPHI),                             
     *  RWORK(LWT),RPAR,IPAR)                                                   
      DO 520 I=1,NY                                                             
         IF(RWORK(I+LWT-1).GT.0.0D0)GO TO 520                                   
           IDID=-3                                                              
           GO TO 527                                                            
520   CONTINUE                                                                  
C                                                                               
C     TEST FOR TOO MUCH ACCURACY REQUESTED.                                     
      R=DDANRM(NY,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*                            
     *   100.0D0*UROUND                                                         
      IF(R.LE.1.0D0)GO TO 525                                                   
C     MULTIPLY RTOL AND ATOL BY R AND RETURN                                    
      IF(INFO(2).EQ.1)GO TO 523                                                 
           RTOL(1)=R*RTOL(1)                                                    
           ATOL(1)=R*ATOL(1)                                                    
           IDID=-2                                                              
           GO TO 527                                                            
523   DO 524 I=1,NY                                                             
           RTOL(I)=R*RTOL(I)                                                    
524        ATOL(I)=R*ATOL(I)                                                    
      IDID=-2                                                                   
      GO TO 527                                                                 
525   CONTINUE                                                                  
C                                                                               
C     COMPUTE MINIMUM STEPSIZE                                                  
      HMIN=4.0D0*UROUND*DMAX1(DABS(TN),DABS(TOUT))                              
C                                                                               
      CALL ODASTP(TN,Y,YPRIME,NEQ,NY,                                      
     *   RES,JAC,H,RWORK(LWT),INFO(1),IDID,RPAR,IPAR,                           
     *   RWORK(LPHI),RWORK(LDELTA),RWORK(LE),                                   
     *   RWORK(LWM),IWORK(LIWM),                                                
     *   RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),                              
     *   RWORK(LPSI),RWORK(LSIGMA),                                             
     *   RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),                                 
     *   RWORK(LS),HMIN,RWORK(LROUND),                                          
     *   IWORK(LPHASE),IWORK(LJCALC),IWORK(LK),                                 
     *   IWORK(LKOLD),IWORK(LNS))    
                                              
527   IF(IDID.LT.0)GO TO 600                                                    
C                                                                               
C------------------------------------------------------                         
C     THIS BLOCK HANDLES THE CASE OF A SUCCESSFUL                               
C     RETURN FROM ODASTP (IDID=1) TEST FOR                                      
C     STOP CONDITIONS.                                                          
C--------------------------------------------------------                       
C                                                                               
      IF(INFO(4).NE.0)GO TO 540                                                 
           IF(INFO(3).NE.0)GO TO 530                                            
             IF((TN-TOUT)*H.LT.0.0D0)GO TO 500                                  
             CALL DDATRP(TN,TOUT,Y,YPRIME,NY,                                   
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))                            
             IDID=3                                                             
             T=TOUT                                                             
             GO TO 580                                                          
530          IF((TN-TOUT)*H.GE.0.0D0)GO TO 535                                  
             T=TN                                                               
             IDID=1                                                             
             GO TO 580                                                          
535          CALL DDATRP(TN,TOUT,Y,YPRIME,NY,                                   
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))                            
             IDID=3                                                             
             T=TOUT                                                             
             GO TO 580                                                          
540   IF(INFO(3).NE.0)GO TO 550                                                 
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 542                                         
         CALL DDATRP(TN,TOUT,Y,YPRIME,NY,                                       
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))                                
         T=TOUT                                                                 
         IDID=3                                                                 
         GO TO 580                                                              
542   IF(DABS(TN-TSTOP).LE.100.0D0*UROUND*                                      
     *   (DABS(TN)+DABS(H)))GO TO 545                                           
      TNEXT=TN+H*(1.0D0+4.0D0*UROUND)                                           
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 500                                     
      H=(TSTOP-TN)*(1.0D0-4.0D0*UROUND)                                         
      GO TO 500                                                                 
545   IDID=2                                                                    
      T=TSTOP                                                                   
      GO TO 580                                                                 
550   IF((TN-TOUT)*H.GE.0.0D0)GO TO 555                                         
      IF(DABS(TN-TSTOP).LE.100.0D0*UROUND*(DABS(TN)+DABS(H)))GO TO 552          
      T=TN                                                                      
      IDID=1                                                                    
      GO TO 580                                                                 
552   IDID=2                                                                    
      T=TSTOP                                                                   
      GO TO 580                                                                 
555   CALL DDATRP(TN,TOUT,Y,YPRIME,NY,                                          
     *   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))                                  
      T=TOUT                                                                    
      IDID=3                                                                    
580   CONTINUE                                                                  
C                                                                               
C--------------------------------------------------------                       
C     ALL SUCCESSFUL RETURNS FROM ODASSL ARE MADE FROM                          
C     THIS BLOCK.                                                               
C--------------------------------------------------------                       
C                                                                               
590   CONTINUE                                                                  
      RWORK(LTN)=TN                                                             
      RWORK(LH)=H
C      open(20,file='test.txt',status='unknown')
C      write(20,*) t                                                                 
      RETURN                                                                    
C                                                                               
C-----------------------------------------------------------------------        
C     THIS BLOCK HANDLES ALL UNSUCCESSFUL                                       
C     RETURNS OTHER THAN FOR ILLEGAL INPUT.                                     
C-----------------------------------------------------------------------        
C                                                                               
600   CONTINUE                                                                  
      ITEMP=-IDID                                                               
      GO TO (610,620,630,690,690,640,650,660,670,675,                           
     *  680,685), ITEMP                                                         
C                                                                               
C     THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE                              
C     REACHING TOUT                                                             
610   CALL XERRWV(                                                              
     *38HDASSL--  AT CURRENT T (=R1)  500 STEPS,                                
     *38,610,0,0,0,0,1,TN,0.0D0)                                                
      CALL XERRWV(48HDASSL--  TAKEN ON THIS CALL BEFORE REACHING TOUT,          
     *48,611,0,0,0,0,0,0.0D0,0.0D0)                                             
      GO TO 690                                                                 
C                                                                               
C     TOO MUCH ACCURACY FOR MACHINE PRECISION                                   
620   CALL XERRWV(                                                              
     *47HDASSL--  AT T (=R1) TOO MUCH ACCURACY REQUESTED,                       
     *47,620,0,0,0,0,1,TN,0.0D0)                                                
      CALL XERRWV(                                                              
     *48HDASSL--  FOR PRECISION OF MACHINE. RTOL AND ATOL,                      
     *48,621,0,0,0,0,0,0.0D0,0.0D0)                                             
      CALL XERRWV(                                                              
     *45HDASSL--  WERE INCREASED TO APPROPRIATE VALUES,                         
     *45,622,0,0,0,0,0,0.0D0,0.0D0)                                             
C                                                                               
      GO TO 690                                                                 
C     WT(I) .LE. 0.0D0 FOR SOME I (NOT AT START OF PROBLEM)                     
630   CALL XERRWV(                                                              
     *38HDASSL--  AT T (=R1) SOME ELEMENT OF WT,                                
     *38,630,0,0,0,0,1,TN,0.0D0)                                                
      CALL XERRWV(28HDASSL--  HAS BECOME .LE. 0.0,                              
     *28,631,0,0,0,0,0,0.0D0,0.0D0)                                             
      GO TO 690                                                                 
C                                                                               
C     ERROR TEST FAILED REPEATEDLY OR WITH H=HMIN                               
640   CALL XERRWV(                                                              
     *44HDASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE,                          
     *44,640,0,0,0,0,2,TN,H)                                                    
      CALL XERRWV(                                                              
     *57HDASSL--  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN,             
     *57,641,0,0,0,0,0,0.0D0,0.0D0)                                             
      GO TO 690                                                                 
C                                                                               
C     CORRECTOR CONVERGENCE FAILED REPEATEDLY OR WITH H=HMIN                    
650   CALL XERRWV(                                                              
     *44HDASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE,                          
     *44,650,0,0,0,0,2,TN,H)                                                    
      CALL XERRWV(                                                              
     *48HDASSL--  CORRECTOR FAILED TO CONVERGE REPEATEDLY,                      
     *48,651,0,0,0,0,0,0.0D0,0.0D0)                                             
      CALL XERRWV(                                                              
     *28HDASSL--  OR WITH ABS(H)=HMIN,                                          
     *28,652,0,0,0,0,0,0.0D0,0.0D0)                                             
      GO TO 690                                                                 
C                                                                               
C     THE ITERATION MATRIX IS SINGULAR                                          
660   CALL XERRWV(                                                              
     *44HDASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE,                          
     *44,660,0,0,0,0,2,TN,H)                                                    
      CALL XERRWV(                                                              
     *37HDASSL--  ITERATION MATRIX IS SINGULAR,                                 
     *37,661,0,0,0,0,0,0.0D0,0.0D0)                                             
      GO TO 690                                                                 
C                                                                               
C     CORRECTOR FAILURE PRECEEDED BY ERROR TEST FAILURES.                       
670   CALL XERRWV(                                                              
     *44HDASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE,                          
     *44,670,0,0,0,0,2,TN,H)                                                    
      CALL XERRWV(                                                              
     *49HDASSL--  CORRECTOR COULD NOT CONVERGE.  ALSO, THE,                     
     *49,671,0,0,0,0,0,0.0D0,0.0D0)                                             
      CALL XERRWV(                                                              
     *38HDASSL--  ERROR TEST FAILED REPEATEDLY.,                                
     *38,672,0,0,0,0,0,0.0D0,0.0D0)                                             
      GO TO 690                                                                 
C                                                                               
C     CORRECTOR FAILURE BECAUSE IRES = -1                                       
675   CALL XERRWV(                                                              
     *44HDASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE,                          
     *44,675,0,0,0,0,2,TN,H)                                                    
      CALL XERRWV(                                                              
     *45HDASSL--  CORRECTOR COULD NOT CONVERGE BECAUSE,                         
     *455,676,0,0,0,0,0,0.0D0,0.0D0)                                            
      CALL XERRWV(                                                              
     *36HDASSL--  IRES WAS EQUAL TO MINUS ONE,                                  
     *36,677,0,0,0,0,0,0.0D0,0.0D0)                                             
      GO TO 690                                                                 
C                                                                               
C     FAILURE BECAUSE IRES = -2                                                 
680   CALL XERRWV(                                                              
     *40HDASSL--  AT T (=R1) AND STEPSIZE H (=R2),                              
     *40,680,0,0,0,0,2,TN,H)                                                    
      CALL XERRWV(                                                              
     *36HDASSL--  IRES WAS EQUAL TO MINUS TWO,                                  
     *36,681,0,0,0,0,0,0.0D0,0.0D0)                                             
      GO TO 690                                                                 
C                                                                               
C     FAILED TO COMPUTE INITIAL YPRIME                                          
685   CALL XERRWV(                                                              
     *44HDASSL--  AT T (=R1) AND STEPSIZE H (=R2) THE,                          
     *44,685,0,0,0,0,2,TN,HO)                                                   
      CALL XERRWV(                                                              
     *45HDASSL--  INITIAL YPRIME COULD NOT BE COMPUTED,                         
     *45,686,0,0,0,0,0,0.0D0,0.0D0)                                             
      GO TO 690                                                                 
690   CONTINUE                                                                  
      INFO(1)=-1                                                                
      T=TN                                                                      
      RWORK(LTN)=TN                                                             
      RWORK(LH)=H                                                               
      RETURN                                                                    
C-----------------------------------------------------------------------        
C     THIS BLOCK HANDLES ALL ERROR RETURNS DUE                                  
C     TO ILLEGAL INPUT, AS DETECTED BEFORE CALLING                              
C     ODASTP. FIRST THE ERROR MESSAGE ROUTINE IS                                
C     CALLED. IF THIS HAPPENS TWICE IN                                          
C     SUCCESSION, EXECUTION IS TERMINATED                                       
C                                                                               
C-----------------------------------------------------------------------        
701   CALL XERRWV(                                                              
     *55HDASSL--  SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE,               
     *55,1,0,0,0,0,0,0.0D0,0.0D0)                                               
      GO TO 750                                                                 
702   CALL XERRWV(25HDASSL--  NEQ (=I1) .LE. 0,                                 
     *25,2,0,1,NEQ,0,0,0.0D0,0.0D0)                                             
      CALL XERRWV(25HOR   --  NY  (=I1) .LE. 0,                                 
     *25,2,0,1,NY ,0,0,0.0D0,0.0D0)                                             
      CALL XERRWV(33HOR   --  NY  (=I1) .GT. NEQ (=I2),                         
     *25,2,0,1,NY , NEQ,0,0.0D0,0.0D0)                                          
      GO TO 750                                                                 
703   CALL XERRWV(34HDASSL--  MAXORD (=I1) NOT IN RANGE,                        
     *34,3,0,1,MXORD,0,0,0.0D0,0.0D0)                                           
      GO TO 750                                                                 
704   CALL XERRWV(                                                              
     *60HDASSL--  RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2),          
     *60,4,0,2,LENRW,LRW,0,0.0D0,0.0D0)                                         
      GO TO 750                                                                 
705   CALL XERRWV(                                                              
     *60HDASSL--  IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2),          
     *60,5,0,2,LENIW,LIW,0,0.0D0,0.0D0)                                         
      GO TO 750                                                                 
706   CALL XERRWV(                                                              
     *39HDASSL--  SOME ELEMENT OF RTOL IS .LT. 0,                               
     *39,6,0,0,0,0,0,0.0D0,0.0D0)                                               
      GO TO 750                                                                 
707   CALL XERRWV(                                                              
     *39HDASSL--  SOME ELEMENT OF ATOL IS .LT. 0,                               
     *39,7,0,0,0,0,0,0.0D0,0.0D0)                                               
      GO TO 750                                                                 
708   CALL XERRWV(                                                              
     *47HDASSL--  ALL ELEMENTS OF RTOL AND ATOL ARE ZERO,                       
     *47,8,0,0,0,0,0,0.0D0,0.0D0)                                               
      GO TO 750                                                                 
709   CALL XERRWV(                                                              
     *54HDASSL--  INFO(4) = 1 AND TSTOP (=R1) BEHIND TOUT (=R2),                
     *54,9,0,0,0,0,2,TSTOP,TOUT)                                                
      GO TO 750                                                                 
710   CALL XERRWV(28HDASSL--  HMAX (=R1) .LT. 0.0,                              
     *28,10,0,0,0,0,1,HMAX,0.0D0)                                               
      GO TO 750                                                                 
711   CALL XERRWV(34HDASSL--  TOUT (=R1) BEHIND T (=R2),                        
     *34,11,0,0,0,0,2,TOUT,T)                                                   
      GO TO 750                                                                 
712   CALL XERRWV(29HDASSL--  INFO(8)=1 AND H0=0.0,                             
     *29,12,0,0,0,0,0,0.0D0,0.0D0)                                              
      GO TO 750                                                                 
713   CALL XERRWV(39HDASSL--  SOME ELEMENT OF WT IS .LE. 0.0,                   
     *39,13,0,0,0,0,0,0.0D0,0.0D0)                                              
      GO TO 750                                                                 
714   CALL XERRWV(                                                              
     *61HDASSL--  TOUT (=R1) TOO CLOSE TO T (=R2) TO START INTEGRATION,         
     *61,14,0,0,0,0,2,TOUT,T)                                                   
      GO TO 750                                                                 
715   CALL XERRWV(                                                              
     *49HDASSL--  INFO(4)=1 AND TSTOP (=R1) BEHIND T (=R2),                     
     *49,15,0,0,0,0,2,TSTOP,T)                                                  
      GO TO 750                                                                 
C                                                                               
717   CALL XERRWV(                                                              
     *52HDASSL--  ML (=I1) .LE. 0 .OR. MU (=I2) .GT. NEQ     ,                  
     *52,17,0,2,IWORK(LML),IWORK(LMU),0,0.0D0,0.0D0)                            
      GO TO 750                                                                 
C                                                                               
718   CALL XERRWV(                                                              
     *52HDASSL--  MU (=I1) .LT. ML (=I2)                     ,                  
     *52,18,0,2,IWORK(LML),IWORK(LMU),0,0.0D0,0.0D0)                            
      GO TO 750                                                                 
719   CALL XERRWV(                                                              
     *39HDASSL--  TOUT (=R1) IS EQUAL TO T (=R2),                               
     *39,19,0,0,0,0,2,TOUT,T)                                                   
      GO TO 750                                                                 
750   IF(INFO(1).EQ.-1) GO TO 760                                               
      INFO(1)=-1                                                                
      IDID=-33                                                                  
      RETURN                                                                    
760   CALL XERRWV(                                                              
     *46HDASSL--  REPEATED OCCURRENCES OF ILLEGAL INPUT,                        
     *46,801,0,0,0,0,0,0.0D0,0.0D0)                                             
770   CALL XERRWV(                                                              
     *47HDASSL--  RUN TERMINATED. APPARENT INFINITE LOOP,                       
     *47,802,1,0,0,0,0,0.0D0,0.0D0)                                             
      RETURN                                                                    
C-----------END OF SUBROUTINE ODASSL------------------------------------        
      END                                                                       
