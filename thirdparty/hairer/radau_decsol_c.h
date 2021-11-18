#ifndef RADAU_DECSOL_C_H
#define RADAU_DECSOL_C_H

#include <stdint.h>

#define TRUE_ (1)
#define FALSE_ (0)

typedef int64_t integer;
typedef double doublereal;
typedef int64_t logical;

// FP_CB = FunctionPointer_CallBack
typedef int (*FP_CB_f)(integer*, doublereal*, doublereal*, doublereal*,
                       doublereal*, integer*, void*);
typedef int (*FP_CB_jac)(integer*, doublereal*, doublereal*, doublereal*,
                         integer*, doublereal*, integer*, void*);
typedef int (*FP_CB_mas)(integer*, doublereal*, integer*, doublereal*,
                         integer*, void*);
typedef int (*FP_CB_solout)(integer*, doublereal*, doublereal*, doublereal*,
                            doublereal*, doublereal*, integer*, integer*,
                            doublereal*, integer*, integer*, void*);

int radau5_c(integer *n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *
            y, doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *
            atol, integer *itol, FP_CB_jac jac, void* jac_PY, integer *ijac, integer *mljac, integer 
            *mujac, FP_CB_mas mas, void* mas_PY, integer *imas, integer *mlmas, integer *mumas, FP_CB_solout 
            solout, void* solout_PY, integer *iout, doublereal *work, integer *lwork, integer *
            iwork, integer *liwork, doublereal *rpar, integer *ipar, integer *
            idid);

doublereal contr5_c(integer *i__, doublereal *x, doublereal *cont, integer * lrc);
                   
#endif