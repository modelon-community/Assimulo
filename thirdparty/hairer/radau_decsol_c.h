#ifndef RADAU_DECSOL_C_H
#define RADAU_DECSOL_C_H

#include <stdint.h>
#include "radau5_superlu_double.h"
#include "radau5_superlu_complex.h"

#define TRUE_ (1)
#define FALSE_ (0)
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
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

typedef struct {double r, i;} doublecomplex;
typedef int (*FP_CB_jac_sparse)(int, double*, double*, int*, double*, int*, int*, doublereal*, integer*, void*);
// typedef int (*FP_CB_assemble_sys_d)(int, double, int *, double *, int *, int *, double *, int *, int *, int, double*);
// typedef int (*FP_CB_assemble_sys_z)(int, double, double, int *, double *, int *, int *, doublecomplex *, int *, int *, int, double*);

// int superlu_setup_z(SuperLU_aux_z *, double, double, double *, int *, int *, int, double*, FP_CB_assemble_sys_z);
int superlu_setup_z(SuperLU_aux_z *, double, double, double *, int *, int *, int, double*, int);

// int radau5_c(integer n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *y,
//             doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *atol,
//             integer *itol, FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, integer *ijac,
//             integer *mljac, integer *mujac, FP_CB_mas mas, void* mas_PY, integer *imas, integer *mlmas,
//             integer *mumas, FP_CB_solout solout, void* solout_PY, integer *iout, doublereal *work,
//             integer *lwork, integer *iwork, integer *liwork, doublereal *rpar, integer *ipar, integer *idid,
//             FP_CB_assemble_sys_d sys_d, FP_CB_assemble_sys_z sys_z, integer);
int radau5_c(integer n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *y,
            doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *atol,
            integer *itol, FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, integer *ijac,
            integer *mljac, integer *mujac, FP_CB_mas mas, void* mas_PY, integer *imas, integer *mlmas,
            integer *mumas, FP_CB_solout solout, void* solout_PY, integer *iout, doublereal *work,
            integer *lwork, integer *iwork, integer *liwork, doublereal *rpar, integer *ipar, integer *idid,
            integer num_threads);

doublereal contr5_c(integer *i__, doublereal *x, doublereal *cont, integer * lrc);
                   
#endif
