#ifndef RADAU_DECSOL_C_H
#define RADAU_DECSOL_C_H

#include "radau5_superlu_double.h"
#include "radau5_superlu_complex.h"
#include "superlu_util.h"

#define TRUE_ (1)
#define FALSE_ (0)
#define radau5_abs(x) ((x) >= 0 ? (x) : -(x))
#define radau_min(a,b) ((a) <= (b) ? (a) : (b))
#define radau_max(a,b) ((a) >= (b) ? (a) : (b))
#define copysign(a,b) (((a < 0 && b > 0) || (a > 0 && b < 0)) ? (-a) : (a))

/* Radau return flags */
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

#define RADAU_SUPERLU_INVALID_INPUT_N             -1
#define RADAU_SUPERLU_INVALID_INPUT_NNZ           -2
#define RADAU_SUPERLU_INVALID_INPUT_NNZ_TOO_LARGE -3
#define RADAU_SUPERLU_INVALID_INPUT_NPROC         -4

#define RADAU_CALLBACK_OK                        0
#define RADAU_CALLBACK_ERROR_RECOVERABLE        -1
#define RADAU_CALLBACK_ERROR_NONRECOVERABLE     -2
#define RADAU_CALLBACK_ERROR_INVALID_JAC_FORMAT -3
/* this one always has to have the smaller number among all RADAU_CALLBACK_ERROR_x */
#define RADAU_CALLBACK_ERROR_INVALID_NNZ        -10

struct Radau_SuperLU_aux{
    int n, nnz, nprocs; 

    int fresh_jacobian;
    int nnz_actual;

    /* Jacobian data */
    double *jac_data;
    int *jac_indices;
    int *jac_indptr;

    SuperLU_aux_d *slu_aux_d;
    SuperLU_aux_z *slu_aux_z;
};
typedef struct Radau_SuperLU_aux Radau_SuperLU_aux;

/* FP_CB = FunctionPointer_CallBack */
typedef int (*FP_CB_f)(int, double*, double*, double*, void*);
typedef int (*FP_CB_jac)(int, double*, double*, double*, void*);
typedef int (*FP_CB_solout)(int*, double*, double*, double*,
                            double*, double*, int*, int*,
                            int*, void*);
typedef int (*FP_CB_jac_sparse)(int, double*, double*, int*, double*, int*, int*, void*);

int radau5_c(int n, FP_CB_f fcn, void* fcn_PY,
			 double *x, double *y, double *xend, double *h__,
			 double *rtol, double *atol, int *itol,
			 FP_CB_jac jac, FP_CB_jac_sparse jac_sparse, void* jac_PY, int *ijac, int sparse_LU,
			 FP_CB_solout solout, void* solout_PY, int *iout,
			 double *work, int *lwork, int *iwork, int *liwork, int *idid,
			 Radau_SuperLU_aux *radau_slu_aux);

int radcor_(int n, FP_CB_f fcn, void* fcn_PY,
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
			double *werr, Radau_SuperLU_aux *radau_slu_aux);

double contr5_c(int *i__, double *x, double *cont, int * lrc);

int dec_(int n, double *a, int *ip, int *ier);
int sol_(int n, double *a, double *b, int *ip);

int decc_(int n, double *ar, double *ai, int *ip, int *ier);
int solc_(int n, double *ar, double *ai, double *br, double *bi, int *ip);

int decomr_(int n, double *fjac,double *fac1, double *e1,
	int *ip1, int *ier, int sparse_LU, Radau_SuperLU_aux *radau_slu_aux);
int decomc_(int n, double *fjac, double *alphn, double *betan,
	double *e2r, double *e2i, int *ip2,
	int *ier, int sparse_LU, Radau_SuperLU_aux *radau_slu_aux);

int slvrad_(int n, double *fac1, double *alphn, double *betan, 
	double *e1, double *e2r, double *e2i, 
	double *z1, double *z2, double *z3,
	double *f1, double *f2, double *f3,
	int *ip1, int *ip2, int sparse_LU,
	Radau_SuperLU_aux *radau_slu_aux);

int estrad_(int n, double *h__,
	double *dd1, double *dd2, double *dd3,
	FP_CB_f fcn, void* fcn_PY, int *nfcn,
	double *y0, double *y, int sparse_LU,
	double *x,double *e1,
	double *z1, double *z2, double *z3,
	double *cont, double *werr,
	double *f1, double *f2, int *ip1,
	double *scal, double *err, int *first, int *reject, 
	Radau_SuperLU_aux *radau_slu_aux, int *ier);

Radau_SuperLU_aux* radau_superlu_aux_setup(int, int, int, int*);
int radau_superlu_aux_finalize(Radau_SuperLU_aux*);

#endif
