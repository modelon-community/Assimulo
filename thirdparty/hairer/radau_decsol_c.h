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
typedef int (*FP_CB_f)(integer, doublereal*, doublereal*, doublereal*, void*);
typedef int (*FP_CB_jac)(integer, doublereal*, doublereal*, doublereal*, void*);
typedef int (*FP_CB_solout)(integer*, doublereal*, doublereal*, doublereal*,
                            doublereal*, doublereal*, integer*, integer*,
                            integer*, void*);
typedef int (*FP_CB_jac_sparse)(int, double*, double*, int*, double*, int*, int*, void*);

int radau5_c(integer n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *y,
            doublereal *xend, doublereal *h__, doublereal *rtol, doublereal *atol,
            integer *itol, FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, integer *ijac,
            integer *mljac, integer *mujac,
            FP_CB_solout solout, void* solout_PY, integer *iout, doublereal *work,
            integer *lwork, integer *iwork, integer *liwork, integer *idid,
            double* jac_data, int* jac_indices, int* jac_indptr,
	        SuperLU_aux_d* slu_aux_d, SuperLU_aux_z* slu_aux_z);

int radcor_(integer n, FP_CB_f fcn, void* fcn_PY, doublereal *x, doublereal *
	y, doublereal *xend, doublereal *hmax, doublereal *h__, doublereal *
	rtol, doublereal *atol, integer *itol, FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, integer *ijac, 
	FP_CB_solout solout, void* solout_PY, integer *iout, integer *idid, integer *nmax, 
	doublereal *uround, doublereal *safe, doublereal *thet, doublereal *
	fnewt, doublereal *quot1, doublereal *quot2, integer *nit, integer *
	ijob, logical *startn, integer *nind1, integer *nind2, integer *nind3,
	logical *pred, doublereal *facl, doublereal *facr, integer *m1, 
	integer *ldjac,
	integer *lde1, doublereal *z1, doublereal *z2, 
	doublereal *z3, doublereal *y0, doublereal *scal, doublereal *f1, 
	doublereal *f2, doublereal *f3, doublereal *fjac, doublereal *e1, 
	doublereal *e2r, doublereal *e2i, integer *ip1, 
	integer *ip2, integer *iphes, doublereal *cont, integer *nfcn, 
	integer *njac, integer *nstep, integer *naccpt, integer *nrejct, 
	integer *ndec, integer *nsol,
	doublereal *werr, integer nnz,
	double* jac_data, int* jac_indices, int* jac_indptr,
	SuperLU_aux_d* slu_aux_d, SuperLU_aux_z* slu_aux_z);

doublereal contr5_c(integer *i__, doublereal *x, doublereal *cont, integer * lrc);

int dec_(integer n, integer *ndim, doublereal *a, integer *ip, integer *ier);
int sol_(integer n, integer *ndim, doublereal *a, doublereal *b, integer *ip);

int decc_(integer n, integer *ndim, doublereal *ar, doublereal *ai, integer *ip, integer *ier);
int solc_(integer n, integer *ndim, doublereal *ar, doublereal *ai, doublereal *br, doublereal *bi, integer *ip);

int decomr_(integer n, doublereal *fjac, integer *ldjac, 
	doublereal *fac1, doublereal *e1,
	integer *lde1, integer *ip1, integer *ier, integer *ijob,
	integer *iphes, SuperLU_aux_d* slu_aux,
	double* jac_data, int* jac_indices, int* jac_indptr, int fresh_jacobian, int jac_nnz);
int decomc_(integer n, doublereal *fjac, integer *ldjac, 
	doublereal *alphn, doublereal *betan,
	doublereal *e2r, doublereal *e2i, integer *lde1, integer *ip2,
	integer *ier, integer *ijob, SuperLU_aux_z* slu_aux,
	double* jac_data, int* jac_indices, int* jac_indptr, int fresh_jacobian, int jac_nnz);

int slvrad_(integer n, doublereal *fac1, doublereal *alphn, doublereal *betan, 
	doublereal *e1, doublereal *e2r, doublereal *e2i, integer *lde1, 
	doublereal *z1, doublereal *z2, doublereal *z3, doublereal *f1, 
	doublereal *f2, doublereal *f3, integer *ip1, 
	integer *ip2, integer *iphes, integer *ier, integer *ijob,
	SuperLU_aux_d* slu_aux_d, SuperLU_aux_z* slu_aux_z);

int estrad_(integer n,
	doublereal *h__, doublereal *dd1, 
	doublereal *dd2, doublereal *dd3, FP_CB_f fcn, void* fcn_PY, integer *nfcn, doublereal 
	*y0, doublereal *y, integer *ijob, doublereal *x,
	doublereal *e1, integer *lde1, doublereal *
	z1, doublereal *z2, doublereal *z3, doublereal *cont, doublereal *
	werr, doublereal *f1, doublereal *f2, integer *ip1,
	doublereal *scal, doublereal *err, logical *first, logical *reject, 
	SuperLU_aux_d* slu_aux_d, integer *ier);

int radau_sparse_aux_init(double**, int**, int**, int, int);
int radau_sparse_aux_finalize(double**, int**, int**);
                   
#endif
