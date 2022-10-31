#ifndef _RADAU5_C_H
#define _RADAU5_C_H

#ifdef __RADAU5_WITH_SUPERLU
#include "superlu_double.h"
#include "superlu_complex.h"
#include "superlu_util.h"
#endif /*__RADAU5_WITH_SUPERLU*/

#define TRUE_ (1)
#define FALSE_ (0)
#define radau5_abs(x) ((x) >= 0 ? (x) : -(x))
#define radau_min(a,b) ((a) <= (b) ? (a) : (b))
#define radau_max(a,b) ((a) >= (b) ? (a) : (b))
#define copysign(a,b) (((a < 0 && b > 0) || (a > 0 && b < 0)) ? (-a) : (a))

/* Radau return flags */
#define RADAU_OK 0 

/* setup */
#define RADAU_SETUP_INVALID_INPUT        -1
#define RADAU_SETUP_MALLOC_FAILURE       -2
#define RADAU_SETUP_SUPERLU_NOT_ENABLED  -100

#define RADAU_SUCCESS                                1
#define RADAU_SUCCESS_SOLOUT_INTERRUPT               2 
#define RADAU_ERROR_INCONSISTENT_INPUT              -1
#define RADAU_ERROR_NMAX_TOO_SMALL                  -2
#define RADAU_ERROR_STEPSIZE_TOO_SMALL              -3
#define RADAU_ERROR_JAC_SINGULAR                    -4
#define RADAU_ERROR_REP_STEP_REJECT                 -5
#define RADAU_ERROR_NNZ_TOO_SMALL                   -6
#define RADAU_ERROR_WRONG_SPARSE_JAC_FORMAT         -7
#define RADAU_ERROR_UNEXPECTED_SUPERLU_FAILURE      -8
#define RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE       -9
#define RADAU_ERROR_UNRECOVERABLE_CALLBACK_ERROR    -10

#define RADAU_CALLBACK_OK                        0
#define RADAU_CALLBACK_ERROR_RECOVERABLE        -1
#define RADAU_CALLBACK_ERROR_NONRECOVERABLE     -2
#define RADAU_CALLBACK_ERROR_INVALID_JAC_FORMAT -3
/* this one always has to have the smaller number among all RADAU_CALLBACK_ERROR_x */
#define RADAU_CALLBACK_ERROR_INVALID_NNZ        -10

#define RADAU_SUPERLU_INVALID_INPUT_N             -1
#define RADAU_SUPERLU_INVALID_INPUT_NNZ           -2
#define RADAU_SUPERLU_INVALID_INPUT_NNZ_TOO_LARGE -3
#define RADAU_SUPERLU_INVALID_INPUT_NPROC         -4

/* Struct for linear solver related memory info */
struct radau_linsol_mem_t{
	int n, fresh_jacobian;
	double *jac; /* both dense and sparse */

	/* DENSE */
	double *e1, *e2r, *e2i; /* dense LU */
	int *ip1, *ip2; /* dense LU pivots */

	/* sparse LU with SUPERLU */
	int nnz, nproc, nnz_actual;

    int *jac_indices, *jac_indptr;

	/* superlu auxiliary structs */
	void *slu_aux_d, *slu_aux_z; 
};
typedef struct radau_linsol_mem_t radau_linsol_mem_t;

struct radau_mem_t{
	double *work; /* base work parameters; TODO */
	double *werr; /* local error estimate*/
	double *z1, *z2, *z3; /* transformed state vector */
	double *y0;
	double *scal;
	double *f1, *f2, *f3; /* newton rhs */
	double *con; /* interpolation*/

	radau_linsol_mem_t *lin_sol;

	int *iwork; /* base iwork parameters; TODO */
};
typedef struct radau_mem_t radau_mem_t;

int setup_radau_mem(int n, int sparseLU, int nprocs, int nnz, void **mem_out);
int setup_radau_linsol_mem(int n, int sparseLU, int nprocs, int nnz, radau_linsol_mem_t **mem_out);

/* FP_CB = FunctionPointer_CallBack */
typedef int (*FP_CB_f)(int, double*, double*, double*, void*);
typedef int (*FP_CB_jac)(int, double*, double*, double*, void*);
typedef int (*FP_CB_solout)(int*, double*, double*, double*,
                            double*, double*, int*, int*,
                            int*, void*);
typedef int (*FP_CB_jac_sparse)(int, double*, double*, int*, double*, int*, int*, void*);

int radau5_c(void* radau_mem, int n, FP_CB_f fcn, void* fcn_PY,
			 double *x, double *y, double *xend, double *h__,
			 double *rtol, double *atol, int *itol,
			 FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, int *ijac, int sparse_LU,
			 FP_CB_solout solout, void* solout_PY, int *iout,
			 double *work, int *lwork, int *iwork, int *liwork, int *idid);

int radcor_(radau_mem_t *rmem, int n, FP_CB_f fcn, void* fcn_PY,
			double *x, double *y, double *xend, double *hmax, double *h__,
			double *rtol, double *atol, int *itol,
			FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, int *ijac,
			FP_CB_solout solout, void* solout_PY, int *iout, int *idid,
			int *nmax, double *uround, double *safe, double *thet,
			double *fnewt, double *quot1, double *quot2, int *nit,
			int sparse_LU, int *startn, int *pred,
			double *facl, double *facr,
			double *z1, double *z2, double *z3,
			double *y0, double *scal,
			double *f1, double *f2, double *f3,
			double *fjac, double *e1, double *e2r, double *e2i,
			int *ip1, int *ip2, double *cont,
			int *nfcn, int *njac, int *nstep, int *naccpt,
			int *nrejct, int *ndec, int *nsol,
			double *werr);

double contr5_c(int *i__, double *x, double *cont, int * lrc);

int dec_(int n, double *a, int *ip, int *ier);
int sol_(int n, double *a, double *b, int *ip);

int decc_(int n, double *ar, double *ai, int *ip, int *ier);
int solc_(int n, double *ar, double *ai, double *br, double *bi, int *ip);

int decomr_(radau_linsol_mem_t *mem, int n, double *fjac,double *fac1, double *e1,
	int *ip1, int *ier, int sparse_LU);
int decomc_(radau_linsol_mem_t *mem, int n, double *fjac, double *alphn, double *betan,
	double *e2r, double *e2i, int *ip2,
	int *ier, int sparse_LU);

int slvrad_(radau_mem_t *rmem, int n, double *fac1, double *alphn, double *betan, 
	double *e1, double *e2r, double *e2i, 
	double *z1, double *z2, double *z3,
	double *f1, double *f2, double *f3,
	int *ip1, int *ip2, int sparse_LU);

int estrad_(radau_mem_t *rmem, int n, double *h__,
	double *dd1, double *dd2, double *dd3,
	FP_CB_f fcn, void* fcn_PY, int *nfcn,
	double *y0, double *y, int sparse_LU,
	double *x,double *e1,
	double *z1, double *z2, double *z3,
	double *cont, double *werr,
	double *f1, double *f2, int *ip1,
	double *scal, double *err, int *first, int *reject,  int *ier);

void free_radau_mem(void **radau_mem);
void free_radau_linsol_mem(radau_linsol_mem_t **mem);

#endif /*_RADAU5_C_H*/
