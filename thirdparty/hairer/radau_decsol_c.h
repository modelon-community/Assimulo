#ifndef RADAU_DECSOL_C_H
#define RADAU_DECSOL_C_H

#include <stdint.h>
#include "radau5_superlu_double.h"
#include "radau5_superlu_complex.h"

#define TRUE_ (1)
#define FALSE_ (0)
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

typedef int integer;
typedef double doublereal;
typedef int logical;

// FP_CB = FunctionPointer_CallBack
typedef int (*FP_CB_f)(integer, doublereal*, doublereal*, doublereal*,
                       doublereal*, integer*, void*);
typedef int (*FP_CB_jac)(integer, doublereal*, doublereal*, doublereal*,
                         integer*, doublereal*, integer*, void*);
typedef int (*FP_CB_mas)(integer, doublereal*, integer*, doublereal*,
                         integer*, void*);
typedef int (*FP_CB_solout)(integer*, doublereal*, doublereal*, doublereal*,
                            doublereal*, doublereal*, integer*, integer*,
                            doublereal*, integer*, integer*, void*);

typedef int (*FP_CB_jac_sparse)(int, double*, double*, int*, double*, int*, int*, doublereal*, integer*, void*);

int radau5_c(integer n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *y,
            doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *atol,
            integer *itol, FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, integer *ijac,
            integer *mljac, integer *mujac, integer *imas, integer *mlmas,
            integer *mumas, FP_CB_solout solout, void* solout_PY, integer *iout, doublereal *work,
            integer *lwork, integer *iwork, integer *liwork, doublereal *rpar, integer *ipar, integer *idid,
            double* jac_data, int* jac_indices, int* jac_indptr,
	        SuperLU_aux_d* slu_aux_d, SuperLU_aux_z* slu_aux_z);

doublereal contr5_c(integer *i__, doublereal *x, doublereal *cont, integer * lrc);

int radau_sparse_aux_init(double**, int**, int**, int, int);
int radau_sparse_aux_finalize(double**, int**, int**);
                   
#endif
