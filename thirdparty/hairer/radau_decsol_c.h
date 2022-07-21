#ifndef RADAU_DECSOL_C_H
#define RADAU_DECSOL_C_H

#include <stdint.h>
#include "radau5_superlu_double.h"
#include "radau5_superlu_complex.h"

#define TRUE_ (1)
#define FALSE_ (0)
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (double)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

// FP_CB = FunctionPointer_CallBack
typedef int (*FP_CB_f)(int, double*, double*, double*, void*);
typedef int (*FP_CB_jac)(int, double*, double*, double*, void*);
typedef int (*FP_CB_solout)(int*, double*, double*, double*,
                            double*, double*, int*, int*,
                            int*, void*);
typedef int (*FP_CB_jac_sparse)(int, double*, double*, int*, double*, int*, int*, void*);

int radau5_c(int n, FP_CB_f fcn, void* fcn_PY,
			 double *x, double *y, double *xend, double *h__,
			 double *rtol, double *atol, int *itol,
			 FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, int *ijac,
			 FP_CB_solout solout, void* solout_PY, int *iout,
			 double *work, int *lwork, int *iwork, int *liwork, int *idid,
			 double* jac_data, int* jac_indices, int* jac_indptr,
			 SuperLU_aux_d* slu_aux_d, SuperLU_aux_z* slu_aux_z);

int radcor_(int n, FP_CB_f fcn, void* fcn_PY,
			double *x, double *y, double *xend, double *hmax, double *h__,
			double *rtol, double *atol, int *itol,
			FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, int *ijac, 
			FP_CB_solout solout, void* solout_PY, int *iout, int *idid,
			int *nmax, double *uround, double *safe, double *thet,
			double *fnewt, double *quot1, double *quot2, int *nit,
			int *ijob, int *startn, int *pred,
			double *facl, double *facr,
			double *z1, double *z2, double *z3,
			double *y0, double *scal,
			double *f1, double *f2, double *f3,
			double *fjac, double *e1, double *e2r, double *e2i,
			int *ip1, int *ip2, double *cont,
			int *nfcn, int *njac, int *nstep, int *naccpt,
			int *nrejct, int *ndec, int *nsol,
			double *werr, int nnz,
			double* jac_data, int* jac_indices, int* jac_indptr,
			SuperLU_aux_d* slu_aux_d, SuperLU_aux_z* slu_aux_z);

double contr5_c(int *i__, double *x, double *cont, int * lrc);

int dec_(int n, double *a, int *ip, int *ier);
int sol_(int n, double *a, double *b, int *ip);

int decc_(int n, double *ar, double *ai, int *ip, int *ier);
int solc_(int n, double *ar, double *ai, double *br, double *bi, int *ip);

int decomr_(int n, double *fjac,double *fac1, double *e1,
	int *ip1, int *ier, int *ijob, SuperLU_aux_d* slu_aux,
	double* jac_data, int* jac_indices, int* jac_indptr, int fresh_jacobian, int jac_nnz);
int decomc_(int n, double *fjac, double *alphn, double *betan,
	double *e2r, double *e2i, int *ip2,
	int *ier, int *ijob, SuperLU_aux_z* slu_aux,
	double* jac_data, int* jac_indices, int* jac_indptr, int fresh_jacobian, int jac_nnz);

int slvrad_(int n, double *fac1, double *alphn, double *betan, 
	double *e1, double *e2r, double *e2i, 
	double *z1, double *z2, double *z3,
	double *f1, double *f2, double *f3,
	int *ip1, int *ip2, int *ier, int *ijob,
	SuperLU_aux_d* slu_aux_d, SuperLU_aux_z* slu_aux_z);

int estrad_(int n, double *h__,
	double *dd1, double *dd2, double *dd3,
	FP_CB_f fcn, void* fcn_PY, int *nfcn,
	double *y0, double *y, int *ijob,
	double *x,double *e1,
	double *z1, double *z2, double *z3,
	double *cont, double *werr,
	double *f1, double *f2, int *ip1,
	double *scal, double *err, int *first, int *reject, 
	SuperLU_aux_d* slu_aux_d, int *ier);

int radau_sparse_aux_init(double**, int**, int**, int, int);
int radau_sparse_aux_finalize(double**, int**, int**);
                   
#endif
