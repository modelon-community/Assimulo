/*
 * Copyright (C) 2010 Modelon AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3 of the License.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.

 * You should have received a copy of the GNU Lesser General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*
 * This is a linear solver using Tikhonov regularization and the sparse
 * solver SuperLU
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "slu_ddefs.h"
#include "supermatrix.h"
#include "reg_routines.h"
#include "kinsol_jmod_impl_wSLU.h"
#include "kinslug.h"
#include "kinsol/kinsol_impl.h"

#include <nvector/nvector_serial.h>
#include <sundials/sundials_direct.h>
#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

#define DEBUG

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */


/* help functions */
void DenSLU2Vector(N_Vector res, SuperMatrix *x);
void Vector2DenSLU(N_Vector x, SuperMatrix *res);


/* KINPinv linit, lsetup, lsolve, and lfree routines */

static int kinSLUGInit(KINMem kin_mem);
static int kinSLUGSetup(KINMem kin_mem);
static int kinSLUGSolve(KINMem kin_mem, N_Vector x, N_Vector b,
                         realtype *res_norm);
static void kinSLUGFree(KINMem kin_mem);


/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define lrw1           (kin_mem->kin_lrw1)
#define liw1           (kin_mem->kin_liw1)
#define func           (kin_mem->kin_func)
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
#define spjac          (kinpinv_mem->d_spjac)
#define J              (kinpinv_mem->d_J)
#define spJ            (kinpinv_mem->sp_J)
#define L              (kinpinv_mem->sp_L)
#define U              (kinpinv_mem->sp_U)
#define B              (kinpinv_mem->sp_B)
#define perm_r         (kinpinv_mem->sp_perm_r)
#define perm_c         (kinpinv_mem->sp_perm_c)
#define spJTJ          (kinpinv_mem->sp_JTJ)
#define options        (kinpinv_mem->sp_options)
#define stat           (kinpinv_mem->sp_stat)
#define pivots         (kinpinv_mem->d_pivots)
#define nje            (kinpinv_mem->d_nje)
#define nfeDQ          (kinpinv_mem->d_nfeDQ)
#define J_data         (kinpinv_mem->d_J_data)
#define last_flag      (kinpinv_mem->d_last_flag)
#define JTJ            (kinpinv_mem->d_JTJ)
#define regularized    (kinpinv_mem->d_regularized)
#define redojac        (kinpinv_mem->d_redojac)
#define beta           (kinpinv_mem->d_beta)
#define reg_param      (kinpinv_mem->d_reg_param)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */
             
/*
 * -----------------------------------------------------------------
 * KINSLUG
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the sparse regularization solver module. 
 * KINSLUG sets the kin_linit, kin_lsetup, kin_lsolve, kin_lfree fields 
 * in *kinmem to be kinSLUGInit, kinSLUGSetup, kinSLUGSolve, and 
 * kinSLUGFree, respectively.  
 * It allocates memory for a structure of type KINPinvMemRec and sets 
 * the kin_lmem field in *kinmem to the address of this structure.  
 * It sets setupNonNull in *kinmem to TRUE, and the djac field to the 
 * default kinPinvDenseDQJac.
 * Finally, it allocates memory for J and pivots.
 *
 * NOTE: The sparse regularized linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINSLUG will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int KINSLUG(void *kinmem, int N)
{
  KINMem kin_mem;
  KINPinvMem kinpinv_mem;
  /* Check if kinmem is different from NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINPINV_MEM_NULL, "KINPINV", "KINPinv", MSGD_KINMEM_NULL);
    return(KINPINV_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is present */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    KINProcessError(kin_mem, KINPINV_ILL_INPUT, "KINPINV", "KINPinv", MSGD_BAD_NVECTOR);
    return(KINPINV_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */
  linit  = kinSLUGInit;
  lsetup = kinSLUGSetup;
  lsolve = kinSLUGSolve;
  lfree  = kinSLUGFree;

  /* Get memory for KINDlsMemRec */
  kinpinv_mem = NULL;
  kinpinv_mem = (KINPinvMem) malloc(sizeof(struct KINPinvMemRec));
  if (kinpinv_mem == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    return(KINPINV_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;  

  /* Set default Jacobian routine and Jacobian data */
  jacDQ  = TRUE;
  spjac   = NULL;
  J_data = NULL;
  last_flag = KINPINV_SUCCESS;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for J,JTJ and pivots */
  
  spJ = NULL;
  spJ = calloc(1,sizeof(SuperMatrix)) ;
  /*spJ = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );*/
  if (spJ == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  L = NULL;
  L = calloc(1,sizeof(SuperMatrix)) ;
  if (L == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    free(spJ);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  U = NULL;
  U = calloc(1,sizeof(SuperMatrix)) ;
  if (U == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    free(spJ);
    free(L);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  B = NULL;
  B = calloc(1,sizeof(SuperMatrix)) ;
  if (B == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    free(spJ);
    free(L);
    free(U);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }
  
  perm_r = NULL;
  perm_r = intMalloc(N);
  if (perm_r == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    free(spJ);
    free(L);
    free(U);
    free(B);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  perm_c = NULL;
  perm_c = calloc(N,sizeof(int));
  if (perm_c == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    free(spJ);
    free(L);
    free(U);
    free(B);
    free(perm_r);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  options = NULL;
  options = (superlu_options_t *) calloc(1,sizeof(superlu_options_t));
  if (options == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    free(spJ);
    free(L);
    free(U);
    free(B);
    free(perm_r);
    free(perm_c);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  stat = NULL;
  stat = (SuperLUStat_t *) calloc(1,sizeof(SuperLUStat_t));
  if (stat == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    free(spJ);
    free(L);
    free(U);
    free(B);
    free(perm_r);
    free(perm_c);
    free(options);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  spJTJ = NULL;
  spJTJ= (SuperMatrix *) calloc(1,sizeof(SuperMatrix));
  if (spJTJ == NULL) {
    KINProcessError(kin_mem, KINPINV_MEM_FAIL, "KINSLUG", "KINSLUG", MSGD_MEM_FAIL);
    free(spJ);
    free(L);
    free(U);
    free(B);
    free(perm_r);
    free(perm_c);
    free(options);
    free(stat);
    free(kinpinv_mem); kinpinv_mem = NULL;
    return(KINPINV_MEM_FAIL);
  }

  if (printfl>0) printf("Sparse data structures created \n");

  /* Set the default input options. */
  set_default_options(options);
  options->ColPerm = NATURAL;
  if (printfl>0)  printf("Sparse options set \n");

  /* Initialize the statistics variables. */
  StatInit(stat);
  if (printfl>0) printf("Sparse stat initialized \n");

  /* This is a direct linear solver */
  inexact_ls = FALSE;

  /* Attach linear solver memory to integrator memory */
  lmem = kinpinv_mem;

  /* Set reg_param 'not set' */
  reg_param = 0;
  

  return(KINPINV_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * kinSLUGInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the Sparse
 * linear solver.
 * -----------------------------------------------------------------
 */

static int kinSLUGInit(KINMem kin_mem)
{
  KINPinvMem kinpinv_mem;

  kinpinv_mem = (KINPinvMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  
  /*
   * Set where to find jacobian data
   * if jacDQ=True it will be calculated by finite differenses
   * otherwise a user-supplied jacobian will be used
   */


  if (jacDQ) {
    djac = kinPinvDQJac;
    J_data = kin_mem;
  } else {
    J_data = kin_mem->kin_user_data;
  }

  /* Set regularization parameter */
  if (reg_param == 0){
    reg_param = 1;
  } 
  last_flag = KINPINV_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinSLUGSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the linear solver and
 * prepares J transpose J  + h^2 I if necessary for regularization.
 * -----------------------------------------------------------------
 */

static int kinSLUGSetup(KINMem kin_mem)
{
  KINPinvMem kinpinv_mem;
  SuperMatrix spJC  ; /* Matrix postmultiplied by P_c */
  int retval,relax,panel_size,*etree,permc_spec;
  double rp = 0;

  kinpinv_mem = (KINPinvMem) lmem;

  regularized = FALSE;

  /* Calculate value of jacobian */
  nje++;
  
  retval = spjac(n, uu, fval, spJ, J_data, vtemp1, vtemp2);
  if (retval != 0) {
    last_flag = -1;
    return(-1);
  }

  /* Set flag that the matrix needs to be factorized */

  options->Fact = DOFACT;

  
  /* Get column permutation vector perm_c */
  permc_spec = options->ColPerm;
  get_perm_c(permc_spec,spJ,perm_c);

  /* Permute Jacobian and calculate and calculate etree */
  etree = calloc(n,sizeof(int));

  sp_preorder(options, spJ, perm_c, etree, &spJC);

  /* Get sizes for supernodes and panels from SUperLU */
  panel_size = sp_ienv(1);
  relax = sp_ienv(2);
  

  /* Try to do a LU factorization of spJ */
  dgstrf(options, &spJC, relax, panel_size, etree,
	 NULL, 0,perm_c, perm_r, L, U, stat, &retval);
  if (retval != 0){
    /* Factorization failed, determine reason */
    if (retval < 0) {
      printf("Argument %i has illegal value\n",(-retval));
    } else if (retval > n) {
      printf("%i nb of bytes allocated when memory allocation failed \n",(retval+n));
    } else {
      /* Jacobian is singular, perform regularization */
      if (printfl > 0) printf("Jacobian is singular at %i, regularizing \n ",retval);

      
      /* Recalculate jacobian */
      
      Destroy_CompCol_Matrix(spJ);
      nje++;
      retval = spjac(n, uu, fval, spJ, J_data, vtemp1, vtemp2);
      if (retval != 0) {
	last_flag = -1;
	free(etree);
	return(-1);
      }

      /* Calculate regularization parameter */
      /* Put rhs in B2 */
      B = calloc(1,sizeof(SuperMatrix));
      Vector2DenSLU(fval,B);
      
      rp = getRegParam(spJ,B);
      
      printf("rp: %f \n",rp*rp);
      if (rp > 1) {
	reg_param = 1;
      } else {
	reg_param = rp;
      }
      if (printfl > 0) printf("Regparam: %f, \n ",reg_param);
      /* Calculate regularized matrix */      
      spJTJ = regSparseMatrix(spJ,sqrt(reg_param));
      
      /* Get column permutation vector perm_c */
      permc_spec = options->ColPerm;
      get_perm_c(permc_spec,spJTJ,perm_c);
      
      /* Do a preorder on the column etree */
      sp_preorder(options, spJTJ, perm_c, etree, &spJC);

      /* Get sizes for supernodes and panels from SUperLU */
      panel_size = sp_ienv(1);
      relax = sp_ienv(2);
  

      /* Try to do a LU factorization of spJ */
      dgstrf(options, &spJC, relax, panel_size, etree,
	     NULL, 0,perm_c, perm_r, L, U, stat, &retval);

      regularized = TRUE;
      if (retval == 0) {
	if (printfl > 0) printf("Jacobian regularized \n");
	/* set flag that the matrix is factorized */
	options->Fact = FACTORED;
      } else {
	if (printfl >0) {
	  if (retval > n) {
	    printf("Error in Jacobian regularization, %d nb of bytes allocated when allocation failed. n: %d",retval,n);
	  } else {
	    printf("Error in Jacobian regularization, retval: %d, n: %d,reg_param: %f",retval,n,reg_param);
	  }
	} 
      }
    }
  
    
  } else {
    /* set flag that the matrix is factorized */
    options->Fact = FACTORED;
    fflush(stdout);
  }

  redojac = FALSE;
  /* Free allocated memory and return flags */
  SUPERLU_FREE(etree);
  Destroy_CompCol_Permuted(&spJC);

  last_flag = retval;
  if (retval > 0) return(-1);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * Vector2SDenSLU
 * -----------------------------------------------------------------
 * Routine used to get the information from a NVector into the structure
 * needed by SuperLU.
 * -----------------------------------------------------------------
 */
void Vector2DenSLU(N_Vector x, SuperMatrix *res) {
  realtype *data;
  int length = NV_LENGTH_S(x);
  Stype_t S = SLU_DN ;
  Dtype_t D = SLU_D ;
  Mtype_t M = SLU_GE;

  /* Get data from the N_Vector */
  data = N_VGetArrayPointer(x);
  
  /* Create dense Supermatrix */
  dCreate_Dense_Matrix(res,length,1,data,length,S,D,M);

}

/*
 * -----------------------------------------------------------------
 * DenSLU2Vector
 * -----------------------------------------------------------------
 * Routine used to get the information from a SuperMatrix into a
 * N_Vector.
 * -----------------------------------------------------------------
 */
void DenSLU2Vector(N_Vector res, SuperMatrix *x) {
  realtype *data;
  DNformat *storage;
  
  /* Get data pointer from SuperMatrix */
  storage = x->Store;
  data = storage->nzval;
  
  /* Set data pointer in N_Vector */
  N_VSetArrayPointer(data,res);

} 

/*
 * -----------------------------------------------------------------
 * kinSLUGSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the sparse linear solver
 * by calling the sparse backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int kinSLUGSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINPinvMem kinpinv_mem;
  SuperMatrix *B2;
  int info;
  kinpinv_mem = (KINPinvMem) lmem;

  if (redojac) {
    return 1;
  }
  /* Copy the right-hand side into x*/
  N_VScale(ONE, b, x);

  
  if (regularized) {

    /* Put rhs in B2 */
    B2 = calloc(1,sizeof(SuperMatrix));
    Vector2DenSLU(x,B2);
    
    /* Get regularized rhs */
    B = getRegRHS(spJ,B2);

    regularized = FALSE;
    redojac = TRUE;
    
    /* Take out the garbage */
    free(B2);
  } else {
    /* put rhs in B */
    Vector2DenSLU(x,B);

  }

  /* Back-solve, get solution in B and propagate it to x*/
  dgstrs (NOTRANS, L, U, perm_c, perm_r, B, stat, &info);

  if (info < 0) {
    printf("Solving failed with flag %i \n",info);
    last_flag = info;
      
  } else {
    DenSLU2Vector(x,B);
  }
  
  
  

  /* Compute the terms Jpnorm and sfdotJp for use in the global strategy
     routines and in KINForcingTerm. Both of these terms are subsequently
     corrected if the step is reduced by constraints or the line search.

     sJpnorm is the norm of the scaled product (scaled by fscale) of
     the current Jacobian matrix J and the step vector p.

     sfdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale. */

  sJpnorm = N_VWL2Norm(b,fscale);
  N_VProd(b, fscale, b);
  N_VProd(b, fscale, b);
  sfdotJp = N_VDotProd(fval, b);
  
  
  last_flag = KINPINV_SUCCESS;
  
  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinSLUGFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void kinSLUGFree(KINMem kin_mem)
{
  KINPinvMem  kinpinv_mem;

  kinpinv_mem = (KINPinvMem) lmem;
  /* De-allocate storage */

  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(spJ);
  Destroy_Dense_Matrix(B);
  Destroy_SuperNode_Matrix(L);
  Destroy_CompCol_Matrix(U);
  Destroy_CompCol_Matrix(spJTJ);
  StatFree(stat);
  free(stat);
  free(options);
  free(spJ);
  free(B);
  free(L);
  free(U);
  free(spJTJ);

  free(kinpinv_mem); kinpinv_mem = NULL;
}

