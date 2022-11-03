#ifndef _RADAU5_C_H
#define _RADAU5_C_H

#include "radau5_impl.h"

/* FP_CB = FunctionPointer_CallBack */
typedef int (*FP_CB_f)(int, double, double*, double*, void*);
typedef int (*FP_CB_jac)(int, double, double*, double*, void*);
typedef int (*FP_CB_solout)(int, double, double, double*, double*, int*, int*, void*);
typedef int (*FP_CB_jac_sparse)(int, double, double*, int*, double*, int*, int*, void*);

int radau5_c(void *radau_mem, FP_CB_f fcn, void *fcn_EXT,
			 double *x, double *y, double *xend, double *h__,
			 double *rtol, double *atol,
			 FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void *jac_EXT, int *ijac,
			 FP_CB_solout solout, void *solout_EXT, int *iout,
			 double *work, int *idid);

int radcor_(radau_mem_t *rmem, int n, FP_CB_f fcn, void *fcn_EXT,
			double *x, double *y, double *xend, double *hmax, double *h__,
			double *rtol, double *atol,
			FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void *jac_EXT, int *ijac,
			FP_CB_solout solout, void *solout_EXT, int *iout, int *idid,
			double *uround, double *thet,
			double *fnewt, double *quot1, double *quot2,
			double *facl, double *facr,
			double *z1, double *z2, double *z3,
			double *y0, double *scal,
			double *f1, double *f2, double *f3,
			double *fjac, double *e1, double *e2r, double *e2i,
			double *cont, double *werr);

int radau_get_cont_output_single(void *radau_mem, int i, double x, double *out);
int radau_get_cont_output(void *radau_mem, double x, double *out);

int dec_(int n, double *a, int *ip, int *ier);
int sol_(int n, double *a, double *b, int *ip);

int decc_(int n, double *ar, double *ai, int *ip, int *ier);
int solc_(int n, double *ar, double *ai, double *br, double *bi, int *ip);

int decomr_(radau_linsol_mem_t *mem, int n, double *fjac, double fac1, double *e1, int *ier);
int decomc_(radau_linsol_mem_t *mem, int n, double *fjac, double alphn, double betan,
	double *e2r, double *e2i, int *ier);

int slvrad_(radau_mem_t *rmem, int n, double fac1, double alphn, double betan, 
	double *e1, double *e2r, double *e2i, 
	double *z1, double *z2, double *z3,
	double *f1, double *f2, double *f3);

int estrad_(radau_mem_t *rmem, int n, double h,
	double dd1, double dd2, double dd3,
	FP_CB_f fcn, void *fcn_EXT,
	double *y0, double *y,
	double x,double *e1,
	double *z1, double *z2, double *z3,
	double *cont, double *werr,
	double *f1, double *f2,
	double *err, int first, int reject);

#endif /*_RADAU5_C_H*/
