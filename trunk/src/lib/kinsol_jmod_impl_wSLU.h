/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2008/04/18 19:42:43 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Common implementation header file for the KINDLS linear solvers.
 * -----------------------------------------------------------------
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
 * Types: KINDlsMemRec, KINDlsMem                             
 * -----------------------------------------------------------------
 * The type KINDlsMem is pointer to a KINDlsMemRec.
 * This structure contains KINDLS solver-specific data. 
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
