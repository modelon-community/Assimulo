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
 * Johan Ylikiiskilä - johan.ylikiiskila@gmail.com
 *
 */

#include "reg_routines.h"
#include <math.h>

SuperMatrix* regSparseMatrix(SuperMatrix *jac, double h)
{
  SuperMatrix *res = NULL;
  int i=0,j=0,k=0,eli=0,elj=0,rowi=0,rowj=0,nnz=0,jnnz=0,jn=0,is_zero=0,max_nz=0;
  NCformat *Jstore = NULL;
  double  *nzt = NULL, *jnzv=NULL,*nzval=NULL;
  double reg_param = h*h;
  int *rit = NULL, *colptr=NULL, *jri=NULL, *jcp=NULL,*rowind=NULL;
  Stype_t S = SLU_NC;
  Dtype_t D = SLU_D;
  Mtype_t M = SLU_GE;

  double *tmp_nzv = NULL;
  int *tmp_ri = NULL;


  /* Extract data */
  Jstore = jac->Store;
  jn = jac->nrow;
  jnnz = Jstore->nnz;
  jnzv = Jstore->nzval;
  jri = Jstore->rowind;
  jcp = Jstore->colptr;

  /* Initialize data and allocate memory for result */
  max_nz = 2*jnnz;
  res = calloc(1,sizeof(SuperMatrix));
  nnz = 0;
  nzt = calloc(max_nz,sizeof(double));
  rit = calloc(max_nz,sizeof(int));
  colptr = calloc(jn+1,sizeof(int));
  
  /* Calculate nnz,nzval,rowind and colptr */
  colptr[0] = 0;
  
  /* Iterate over the columns of both J^T and J */
  for (i = 0;i < jn;i++) {
    for (j = 0;j < jn;j++) {
      /* Get first index of first elt in both columns */
      rowi = jcp[i];
      rowj = jcp[j];

      if (i==j) {
	/* Diagonal elt, perform regularization */
	is_zero = 0;
	rit[nnz] = j;
	nzt[nnz++] += reg_param;
      } else {
	is_zero = 1;
      }

      /* Chack for matching and perform multiplication */
      while ((rowi < jcp[i+1])&&(rowj < jcp[j+1])) {
	/* Get row numbers of the two current elts */
	eli = jri[rowi];
	elj = jri[rowj];
	if (eli == elj) {
	  /* Match! Non-zero elt to be inserted! */
	  /* If it is the first time, add to the non_zero structure */
	  if (is_zero) {
	    is_zero = 0;
	    rit[nnz++] = j;
	  }
	  nzt[nnz-1] += jnzv[rowi]*jnzv[rowj];
	  rowi++;
	  rowj++;
	} else if (eli < elj) {
	  rowi++;
	} else {
	  rowj++;
	}
      }
    }
    /* New column, add info to colptr */
    colptr[i+1] = nnz;

    /* 
     * Check if we need to expand the allocated memory 
     * and expand if necessary
     */
    
    if (nnz > (max_nz-jn)) {
      tmp_nzv = calloc(2*max_nz,sizeof(double));
      tmp_ri = calloc(2*max_nz,sizeof(int));
      
      for (k = 0; k < nnz; k++) {
	tmp_nzv[k] = nzt[k];
	tmp_ri[k] = rit[k];
      }

      free(nzt);
      free(rit);
      
      nzt = tmp_nzv;
      rit = tmp_ri;

      tmp_nzv = NULL;
      tmp_ri = NULL;
      
      max_nz = 2*max_nz;
    }
  }
  
  /* Reallocate to the needed size of storage */
  nzval = realloc(nzt,nnz*sizeof(double));
  rowind = realloc(rit,nnz*sizeof(int));
  
  /* Now create the matrix */
  dCreate_CompCol_Matrix(res, jn, jn, nnz, nzval, rowind, colptr, S, D, M);
  return res ;
}

SuperMatrix* getRegRHS( SuperMatrix *jac, SuperMatrix *B) {
  SuperMatrix *res;
  NCformat *Jstore;
  DNformat *Bstore;
  double *data, *jac_data, *rhs_data;
  int *jac_ri, *jac_cp;
  int n,i,j,k ;
  
  /* Extract data from Jacobian */
  Jstore = jac->Store;
  jac_data = Jstore->nzval;
  jac_ri = Jstore->rowind;
  jac_cp = Jstore->colptr;
  
  /* Extract data from right hand side */
  n = B->nrow;
  Bstore = B->Store;
  rhs_data = Bstore->nzval;

  /* Allocate data for result */
  data = calloc(n,sizeof(double));

  /* Calculate data for result */
  for (i = 0; i< n; i++) {
    for (j = jac_cp[i]; j < jac_cp[i+1]; j++) {
      k = jac_ri[j] ;
      data[i] += jac_data[j]*rhs_data[k];
    }
  }
  
  /* Create the Matrix in SuperMatrix format */
  res = calloc(1,sizeof(SuperMatrix));
  dCreate_Dense_Matrix(res,n,1,data,n,SLU_DN,SLU_D,SLU_GE);

  return res;
}

double getRegParam(SuperMatrix *jac, SuperMatrix *B) {
  double res = 0;
  NCformat *Jstore;
  DNformat *Bstore;
  double data, *jac_data, *rhs_data;
  int *jac_ri, *jac_cp;
  int n,i,j,k ;
  
  /* Extract data from Jacobian */
  Jstore = jac->Store;
  jac_data = Jstore->nzval;
  jac_ri = Jstore->rowind;
  jac_cp = Jstore->colptr;
  
  /* Extract data from right hand side */
  n = B->nrow;
  Bstore = B->Store;
  rhs_data = Bstore->nzval;

  /* Calculate data for result */
  for (i = 0; i< n; i++) {
    data = 0;
    for (j = jac_cp[i]; j < jac_cp[i+1]; j++) {
      k = jac_ri[j] ;
      data += jac_data[j]*rhs_data[k];
    }
    res += data*data;
  }

  return sqrt(res);
}
