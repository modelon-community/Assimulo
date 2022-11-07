#ifndef _RADAU5_C_H
#define _RADAU5_C_H

#include "radau5_impl.h"

/* FP_CB = FunctionPointer_CallBack */
typedef int (*FP_CB_f)(int, double, double*, double*, void*);
typedef int (*FP_CB_jac)(int, double, double*, double*, void*);
typedef int (*FP_CB_solout)(int, double, double, double*, double*, int, int*, void*);
typedef int (*FP_CB_jac_sparse)(int, double, double*, int*, double*, int*, int*, void*);

int radau5_solve(void *radau_mem, FP_CB_f fcn, void *fcn_EXT,
				 double *x, double *y, double *xend, double *h__,
				 double *rtol, double *atol,
				 FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void *jac_EXT, int *ijac,
				 FP_CB_solout solout, void *solout_EXT, int *iout, int *idid);

int radau_get_cont_output_single(void *radau_mem, int i, double x, double *out);
int radau_get_cont_output(void *radau_mem, double x, double *out);

#endif /*_RADAU5_C_H*/
