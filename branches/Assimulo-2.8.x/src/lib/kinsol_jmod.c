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

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "kinsol/kinsol_impl.h"
#include "kinsol_jmod_impl.h"
#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT RCONST(1000.0)

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define lrw1           (kin_mem->kin_lrw1)
#define liw1           (kin_mem->kin_liw1)
#define uround         (kin_mem->kin_uround)
#define func           (kin_mem->kin_func)
#define user_data      (kin_mem->kin_user_data)
#define printfl        (kin_mem->kin_printfl)
#define linit          (kin_mem->kin_linit)
#define lsetup         (kin_mem->kin_lsetup)
#define lsolve         (kin_mem->kin_lsolve)
#define lfree          (kin_mem->kin_lfree)
#define lmem           (kin_mem->kin_lmem)
#define inexact_ls     (kin_mem->kin_inexact_ls)
#define uu             (kin_mem->kin_uu)
#define fval           (kin_mem->kin_fval)
#define uscale         (kin_mem->kin_uscale)
#define fscale         (kin_mem->kin_fscale)
#define sqrt_relfunc   (kin_mem->kin_sqrt_relfunc)
#define sJpnorm        (kin_mem->kin_sJpnorm)
#define sfdotJp        (kin_mem->kin_sfdotJp)
#define errfp          (kin_mem->kin_errfp)
#define infofp         (kin_mem->kin_infofp)
#define setupNonNull   (kin_mem->kin_setupNonNull)
#define vtemp1         (kin_mem->kin_vtemp1)
#define vec_tmpl       (kin_mem->kin_vtemp1)
#define vtemp2         (kin_mem->kin_vtemp2)

#define mtype          (kinpinv_mem->d_type)
#define n              (kinpinv_mem->d_n)
#define ml             (kinpinv_mem->d_ml)
#define mu             (kinpinv_mem->d_mu)
#define smu            (kinpinv_mem->d_smu)
#define jacDQ          (kinpinv_mem->d_jacDQ)
#define djac           (kinpinv_mem->d_djac)
#define bjac           (kinpinv_mem->d_bjac)
#define J              (kinpinv_mem->d_J)
#define pivots         (kinpinv_mem->d_pivots)
#define nje            (kinpinv_mem->d_nje)
#define nfeDQ          (kinpinv_mem->d_nfeDQ)
#define J_data         (kinpinv_mem->d_J_data)
#define last_flag      (kinpinv_mem->d_last_flag)
#define RTR            (kinpinv_mem->d_RTR)
#define reg_param      (kinpinv_mem->d_reg_param)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * KINPinvSetJacFn
 * -----------------------------------------------------------------
 */

int KINPinvSetJacFn(void *kinmem, KINPinvJacFn jac)
{
  KINMem kin_mem;
  KINPinvMem kinpinv_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINPINV_MEM_NULL, "KINPINV", "KINPinvSetJacFn", MSGD_KINMEM_NULL);
    return(KINPINV_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINPINV_LMEM_NULL, "KINPINV", "KINPinvSetJacFn", MSGD_LMEM_NULL);
    return(KINPINV_LMEM_NULL);
  }
  kinpinv_mem = (KINPinvMem) lmem;

  if (jac != NULL) {
    jacDQ = FALSE;
    djac = jac;
  } else {
    jacDQ = TRUE;
  }

  return(KINPINV_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINPinvSetRegParam
 * -----------------------------------------------------------------
 */

int KINPinvSetRegParam(void *kinmem, realtype reg_p)
{
  KINMem kin_mem;
  KINPinvMem kinpinv_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINPINV_MEM_NULL, "KINPINV", "KINPinvSetJacFn", MSGD_KINMEM_NULL);
    return(KINPINV_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINPINV_LMEM_NULL, "KINPINV", "KINPinvSetJacFn", MSGD_LMEM_NULL);
    return(KINPINV_LMEM_NULL);
  }
  kinpinv_mem = (KINPinvMem) lmem;

  reg_param = reg_p;

  return(KINPINV_SUCCESS);
}


/*
 * -----------------------------------------------------------------
 * KINPinvGetWorkSpace
 * -----------------------------------------------------------------
 */

int KINPinvGetWorkSpace(void *kinmem, long int *lenrwLS, long int *leniwLS)
{
  KINMem kin_mem;
  KINPinvMem kinpinv_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINPINV_MEM_NULL, "KINPINV", "KINPinvGetWorkSpace", MSGD_KINMEM_NULL);
    return(KINPINV_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINPINV_LMEM_NULL, "KINPINV", "KINPinvGetWorkSpace", MSGD_LMEM_NULL);
    return(KINPINV_LMEM_NULL);
  }
  kinpinv_mem = (KINPinvMem) lmem;

  if (mtype == SUNDIALS_DENSE) {
    *lenrwLS = n*n;
    *leniwLS = n;
  } else if (mtype == SUNDIALS_BAND) {
    *lenrwLS = n*(smu + mu + 2*ml + 2);
    *leniwLS = n;
  }

  return(KINPINV_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINPinvGetNumJacEvals
 * -----------------------------------------------------------------
 */

int KINPinvGetNumJacEvals(void *kinmem, long int *njevals)
{
  KINMem kin_mem;
  KINPinvMem kinpinv_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINPINV_MEM_NULL, "KINPINV", "KINPinvGetNumJacEvals", MSGD_KINMEM_NULL);
    return(KINPINV_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINPINV_LMEM_NULL, "KINPINV", "KINPinvGetNumJacEvals", MSGD_LMEM_NULL);
    return(KINPINV_LMEM_NULL);
  }
  kinpinv_mem = (KINPinvMem) lmem;

  *njevals = nje;

  return(KINPINV_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINPinvGetNumFuncEvals
 * -----------------------------------------------------------------
 */

int KINPinvGetNumFuncEvals(void *kinmem, long int *nfevalsLS)
{
  KINMem kin_mem;
  KINPinvMem kinpinv_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINPINV_MEM_NULL, "KINPINV", "KINPinvGetNumFuncEvals", MSGD_KINMEM_NULL);
    return(KINPINV_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINPINV_LMEM_NULL, "KINPINV", "KINPinvGetNumGuncEvals", MSGD_LMEM_NULL);
    return(KINPINV_LMEM_NULL);
  }
  kinpinv_mem = (KINPinvMem) lmem;

  *nfevalsLS = nfeDQ;

  return(KINPINV_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINPinvGetLastFlag
 * -----------------------------------------------------------------
 */

int KINPinvGetLastFlag(void *kinmem, int *flag)
{
  KINMem kin_mem;
  KINPinvMem kinpinv_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINPINV_MEM_NULL, "KINPINV", "KINPinvGetLastFlag", MSGD_KINMEM_NULL);
    return(KINPINV_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  if (lmem == NULL) {
    KINProcessError(kin_mem, KINPINV_LMEM_NULL, "KINPINV", "KINPinvGetLastFlag", MSGD_LMEM_NULL);
    return(KINPINV_LMEM_NULL);
  }
  kinpinv_mem = (KINPinvMem) lmem;

  *flag = last_flag;

  return(KINPINV_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * KINPinvGetReturnFlagName
 * -----------------------------------------------------------------
 */

char *KINPinvGetReturnFlagName(int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case KINPINV_SUCCESS:
    sprintf(name, "KINPINV_SUCCESS");
    break;
  case KINPINV_MEM_NULL:
    sprintf(name, "KINPINV_MEM_NULL");
    break;
  case KINPINV_LMEM_NULL:
    sprintf(name, "KINPINV_LMEM_NULL");
    break;
  case KINPINV_ILL_INPUT:
    sprintf(name, "KINPINV_ILL_INPUT");
    break;
  case KINPINV_MEM_FAIL:
    sprintf(name, "KINPINV_MEM_FAIL");
    break;
  default:
    sprintf(name, "NONE");
  }

  return(name);
}


/* 
 * =================================================================
 * DQ JACOBIAN APPROXIMATIONS
 * =================================================================
 */




/*
 * -----------------------------------------------------------------
 * kinPinvDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of F(u). It assumes that a dense matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 *
 * The increment used in the finite-difference approximation
 *   J_ij = ( F_i(u+sigma_j * e_j) - F_i(u)  ) / sigma_j
 * is
 *  sigma_j = max{|u_j|, |1/uscale_j|} * sqrt(uround)
 *
 * Note: uscale_j = 1/typ(u_j)
 *
 * NOTE: Any type of failure of the system function her leads to an
 *       unrecoverable failure of the Jacobian function and thus
 *       of the linear solver setup function, stopping KINSOL.
 * -----------------------------------------------------------------
 */

int kinPinvDQJac(int N,
                 N_Vector u, N_Vector fu,
                 DlsMat Jac, void *data,
                 N_Vector tmp1, N_Vector tmp2)
{
  realtype inc, inc_inv, ujsaved, ujscale, sign;
  realtype *tmp2_data, *u_data, *uscale_data;
  N_Vector ftemp, jthCol;
  long int j;
  int retval = 0;

  KINMem kin_mem;
  KINPinvMem  kinpinv_mem;

  /* data points to kin_mem */
  kin_mem = (KINMem) data;
  kinpinv_mem = (KINPinvMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for u and uscale */
  u_data   = N_VGetArrayPointer(u);
  uscale_data = N_VGetArrayPointer(uscale);

  /* This is the only for loop for 0..N-1 in KINSOL */

  for (j = 0; j < N; j++) {

    /* Generate the jth col of Jac(u) */

    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);

    ujsaved = u_data[j];
    ujscale = ONE/uscale_data[j];
    sign = (ujsaved >= ZERO) ? ONE : -ONE;
#if SUNDIALS_26
    inc = sqrt_relfunc*SUNMAX(SUNRabs(ujsaved), ujscale)*sign;
#else 
    inc = sqrt_relfunc*MAX(ABS(ujsaved), ujscale)*sign;
#endif
    u_data[j] += inc;

    retval = func(u, ftemp, user_data);
    nfeDQ++;
    if(retval > 0) {
        /* try to recover by stepping in the opposite direction */
        inc = -inc;
        u_data[j] = ujsaved + inc;
    
        retval = func(u, ftemp, user_data);
        nfeDQ++;
    }
    if (retval != 0) break; 

    u_data[j] = ujsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fu, jthCol);

  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}
