#ifndef _RADAU5_IMPL_H
#define _RADAU5_IMPL_H

/* Radau return flags */
#define RADAU_OK 0 

/* setup */
#define RADAU_SETUP_INVALID_INPUT        -1
#define RADAU_SETUP_MALLOC_FAILURE       -2
#define RADAU_SETUP_SUPERLU_NOT_ENABLED  -100

/* Setting parameters */
#define RADAU_PARA_RADAU_MEM_NULL       -1
#define RADAU_PARA_INCONSISTENT_INPUT	-2

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

/* CONTAINS BASIC DATA STRUCTURES */

/* Struct for linear solver related memory info */
struct radau_linsol_mem_t{
	int n; /* problem size */
	int LU_with_fresh_jac; /* flag to superLU, if jacobian is fresh */
	int sparseLU; /* flag if using sparse solver */
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

/* structure for statistics logged when running Radau5 */
struct radau_stats_t{
	int nfcn; /* number of rhs function evals */
	int njac; /* jacobian evaluations */
	int nsteps; /* steps */
	int naccpt; /* accepted steps */
	int nreject; /* rejected steps */

	/* one LU decomp/sol includes both the real&complex LU*/
	int ludecomps; /* LU decompositions */
	int lusolves; /* LU solves */
};
typedef struct radau_stats_t radau_stats_t;

/* structure for parameters, both internal and input, for Radau5 */
struct radau_parameters_t{
	int use_jac_external; /* switch for internal/external jacobian */
	int callback_enabled; /* switch for calling callback function after each successful step */
	int atol_as_vec; /* switch if atol is a vector or scalar */

	int nmax; /* maximum number of timesteps */
	int nmax_newton; /* maximum number of newton iterations */

	/* switch for newton starting values strategy */
	/* newton_start_zero = 0 : EXTRAPOLATED COLLOCATION SOLUTION */
	/* newton_start_zero != 0 : ZERO STARTING VALUES, RECOMMENDED IF NEWTON'S METHOD HAS DIFFICULTIES WITH CONVERGENCE*/
	int newton_start_zero; 
	/* switch for predictive step size strategy */
	/* pred_step_control = 0; : MOD. PREDICTIVE CONTROLLER (GUSTAFSSON) */
	/* pred_step_control != 0; : CLASSICAL STEP SIZE CONTROL */
	/* = 0 is considered safer, while != 0 may often yield slightly faster runs for simple problems*/
	int pred_step_control;

	double step_size_safety; /* safety factor for stepsize control */

	/* INTERNAL PARAMETERS */

	int jac_is_fresh; /* flag if computed jacobian is fresh, used in fallbacks */
	double h_old; /* previous step-size, used in predictive controller */
	double dynold; /* newton related parameter */
	double erracc; /* bound in predictive controller */
	double thqold; /* newton related parameter */

	int new_jac_req; /* flag if new jacobian is required */
};
typedef struct radau_parameters_t radau_parameters_t;

/* structure of (computed) mathematical constants inside Radau5 */
struct radau_math_const_t{
	double expm;
	double expm_inv;

	double c1, c2, c1mc2;
	double dd1, dd2, dd3;
	double u1, alph, beta;

	double c1m1, c2m1;
};
typedef struct radau_math_const_t radau_math_const_t;
struct radau_mem_t{
	int n; /* problem size */

	double *work; /* base work parameters; TODO */
	double *werr; /* local error estimate*/
	double *z1, *z2, *z3; /* transformed state vector */
	double *y0;
	double *scal;
	double *f1, *f2, *f3; /* newton rhs */
	double *cont; /* interpolation*/

	double *rtol, *atol; /* internal rtol and atol vectors */

	/* derived mathematical constants */
	radau_math_const_t *mconst;

	/* linear solver related data */
	radau_linsol_mem_t *lin_sol;

	/* parameters relevant for solver control */
	radau_parameters_t *para;

	/* solver statistics */
	radau_stats_t *stats;

	/* dense output */
	double xsol;  /* latest solution point */
	double hsol;  /* latest successful stepsize */
};
typedef struct radau_mem_t radau_mem_t;

#endif /*_RADAU5_IMPL_H*/
