/*
 * This file has been modified by Modelon AB, changes are Copyright
 * (c) 2014 Modelon AB.
 * 
 * Original copyright notice:
 * 
 * Copyright (c) 2002, The Regents of the University of California. 
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by S.D. Cohen, A.C. Hindmarsh, R. Serban, 
 *            D. Shumaker, and A.G. Taylor.
 * UCRL-CODE-155951    (CVODE)
 * UCRL-CODE-155950    (CVODES)
 * UCRL-CODE-155952    (IDA)
 * UCRL-CODE-237203    (IDAS)
 * UCRL-CODE-155953    (KINSOL)
 * All rights reserved. 

 * This file is part of SUNDIALS.

 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the disclaimer below.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the disclaimer (as noted below)
 * in the documentation and/or other materials provided with the
 * distribution.
 * 
 * 3. Neither the name of the UC/LLNL nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * REGENTS OF THE UNIVERSITY OF CALIFORNIA, THE U.S. DEPARTMENT OF ENERGY
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * Additional BSD Notice
 * ---------------------
 * 1. This notice is required to be provided under our contract with
 * the U.S. Department of Energy (DOE). This work was produced at the
 * University of California, Lawrence Livermore National Laboratory
 * under Contract No. W-7405-ENG-48 with the DOE.
 * 
 * 2. Neither the United States Government nor the University of
 * California nor any of their employees, makes any warranty, express
 * or implied, or assumes any liability or responsibility for the
 * accuracy, completeness, or usefulness of any information, apparatus,
 * product, or process disclosed, or represents that its use would not
 * infringe privately-owned rights.
 * 
 * 3. Also, reference herein to any specific commercial products,
 * process, or services by trade name, trademark, manufacturer or
 * otherwise does not necessarily constitute or imply its endorsement,
 * recommendation, or favoring by the United States Government or the
 * University of California. The views and opinions of authors expressed
 * herein do not necessarily state or reflect those of the United States
 * Government or the University of California, and shall not be used for
 * advertising or product endorsement purposes.
 */


#ifndef _KINJMOD_H
#define _KINJMOD_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <sundials/sundials_direct.h>
#include <sundials/sundials_nvector.h>

/*
 * =================================================================
 *              K I N P I N V    C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * KINPINV return values 
 * -----------------------------------------------------------------
 */

#define KINPINV_SUCCESS           0
#define KINPINV_MEM_NULL         -1
#define KINPINV_LMEM_NULL        -2
#define KINPINV_ILL_INPUT        -3
#define KINPINV_MEM_FAIL         -4

/* Additional last_flag values */

#define KINPINV_JACFUNC_UNRECVR  -5
#define KINPINV_JACFUNC_RECVR    -6

/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type: KINPinvJacFn
 * -----------------------------------------------------------------
 *
 * A dense Jacobian approximation function Jac must be of type 
 * KINDlsDenseJacFn. Its parameters are:
 *
 * N        - problem size.
 *
 * u        - current iterate (unscaled) [input]
 *
 * fu       - vector (type N_Vector) containing result of nonlinear
 *            system function evaluated at current iterate:
 *            fu = F(u) [input]
 *
 * J        - dense matrix (of type DlsMat) that will be loaded
 *            by a KINDlsDenseJacFn with an approximation to the
 *            Jacobian matrix J = (dF_i/dy_j).
 *
 * user_data   - pointer to user data - the same as the user_data
 *            parameter passed to KINSetFdata.
 *
 * tmp1, tmp2 - available scratch vectors (volatile storage)
 *
 * A KINPinvJacFn should return 0 if successful, a positive 
 * value if a recoverable error occurred, and a negative value if 
 * an unrecoverable error occurred.
 *
 * -----------------------------------------------------------------
 *
 * NOTE: The following are two efficient ways to load a dense Jac:         
 * (1) (with macros - no explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = DENSE_COL(Jac,j);                                 
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * (2) (without macros - explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = (Jac->data)[j];                                   
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * A third way, using the DENSE_ELEM(A,i,j) macro, is much less   
 * efficient in general.  It is only appropriate for use in small 
 * problems in which efficiency of access is NOT a major concern. 
 *                                                                
 * -----------------------------------------------------------------
 */
  
  
typedef int (*KINPinvJacFn)(int N,
			    N_Vector u, N_Vector fu, 
			    DlsMat J, void *user_data,
			    N_Vector tmp1, N_Vector tmp2);
  

/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the KINPinv linear solver
 * -----------------------------------------------------------------
 *
 * KINDlsSetDenseJacFn specifies the dense Jacobian approximation
 * routine to be used for a direct dense linear solver.
 *
 * By default, a difference quotient approximation, supplied with
 * the solver is used.
 *
 * The return value is one of:
 *    KINPINV_SUCCESS   if successful
 *    KINPINV_MEM_NULL  if the KINSOL memory was NULL
 *    KINPINV_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINPinvSetJacFn(void *kinmem, KINPinvJacFn jac);

/*
  Set regularization parameter
  */
SUNDIALS_EXPORT int KINPinvSetRegParam(void *kinmem, realtype reg_p);

/*
 * -----------------------------------------------------------------
 * Optional outputs from a KINDLS linear solver
 * -----------------------------------------------------------------
 *
 * KINPinvGetWorkSpace    returns the real and integer workspace used
 *                       by the KINDLS linear solver.
 * KINPinvGetNumJacEvals  returns the number of calls made to the
 *                       Jacobian evaluation routine.
 * KINPinvGetNumFuncEvals returns the number of calls to the user's F
 *                       routine due to finite difference Jacobian
 *                       evaluation.
 * KINPinvGetLastFlag     returns the last error flag set by any of
 *                       the KINDLS interface functions.
 * KINPinvGetReturnFlagName returns the name of the constant
 *                       associated with a KINDLS return flag
 *
 * The return value of KINPinvGet* is one of:
 *    KINPINV_SUCCESS   if successful
 *    KINPINV_MEM_NULL  if the KINSOL memory was NULL
 *    KINPINV_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */


SUNDIALS_EXPORT int KINPinvGetWorkSpace(void *kinmem, long int *lenrwB, long int *leniwB);
SUNDIALS_EXPORT int KINPinvGetNumJacEvals(void *kinmem, long int *njevalsB);
SUNDIALS_EXPORT int KINPinvGetNumFuncEvals(void *kinmem, long int *nfevalsB);
SUNDIALS_EXPORT int KINPinvGetLastFlag(void *kinmem, int *flag);
SUNDIALS_EXPORT char *KINPinvGetReturnFlagName(int flag);

#ifdef __cplusplus
}
#endif

#endif
