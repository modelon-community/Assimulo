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

#ifndef _KINJMOD_IMPL_H
#define _KINJMOD_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "kinsol_jmod_wSLU.h"
#include "slu_ddefs.h"


/*
 * -----------------------------------------------------------------
 * Types: KINPinvMemRec, KINPinvMem                             
 * -----------------------------------------------------------------
 * The type KINPinvMem is pointer to a KINPinvMemRec.
 * This structure contains KINPinv/KINSLUG solver-specific data. 
 * -----------------------------------------------------------------
 */

typedef struct KINPinvMemRec {

  int d_type;              /* SUNDIALS_DENSE or SUNDIALS_BAND              */

  int d_n;                 /* problem dimension                            */

  int d_ml;                /* lower bandwidth of Jacobian                  */
  int d_mu;                /* upper bandwidth of Jacobian                  */ 
  int d_smu;               /* upper bandwith of M = MIN(N-1,d_mu+d_ml)     */

  booleantype d_jacDQ;     /* TRUE if using internal DQ Jacobian approx.   */
  KINPinvJacFn d_djac;     /* dense Jacobian routine to be called          */
  KINSLUGJacFn d_spjac;    /* sparse Jacobian routine to be called         */

  void *d_J_data;          /* J_data is passed to djac or bjac             */
    
  DlsMat d_J;              /* problem Jacobian                             */

  SuperMatrix *sp_J;       /* sparse problem jacobian                      */
  SuperMatrix *sp_L;       /* L in the sparse LU factorization             */
  SuperMatrix *sp_U;       /* U in the sparse LU factorization             */
  SuperMatrix *sp_B;       /* sparse right hand side                       */
  int *sp_perm_r;          /* row permutations from partial pivoting       */
  int *sp_perm_c;          /* column permutation vector                    */

  SuperMatrix *sp_JTJ;     /* Matrix needed for regularisation             */

  superlu_options_t *sp_options; /* options struct for SuperLU             */
  SuperLUStat_t *sp_stat;  /* statistcis struct for SuperLU                */

  int *d_pivots;           /* pivot array for PM = LU                      */
  realtype *d_beta;
  realtype d_reg_param;    /* Regularization parameter                     */
    
  long int d_nje;          /* no. of calls to jac                          */
    
  long int d_nfeDQ;        /* no. of calls to F due to DQ Jacobian approx. */
    
  int d_last_flag;         /* last error return flag                       */

  DlsMat d_JTJ;
  booleantype d_regularized; /* Boolean set to true if problem is regularized*/
  booleantype d_redojac;
    
} *KINPinvMem;


/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

int kinPinvDQJac(int N,
		 N_Vector u, N_Vector fu,
		 DlsMat Jac, void *data,
		 N_Vector tmp1, N_Vector tmp2);

/*
 * -----------------------------------------------------------------
 * Error Messages
 * -----------------------------------------------------------------
 */

#define MSGD_KINMEM_NULL "KINSOL memory is NULL."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_MEM_FAIL    "A memory request failed."
#define MSGD_LMEM_NULL   "Linear solver memory is NULL."
#define MSGD_BAD_SIZES   "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#ifdef __cplusplus
}
#endif

#endif
