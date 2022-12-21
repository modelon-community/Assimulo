#ifndef _RADAU5_IMPL_H
#define _RADAU5_IMPL_H

/* Radau return flags */
#define RADAU_OK 0

#define RADAU_SUCCESS_SOLOUT_INTERRUPT 1

#define RADAU_ERROR_MEM_NULL                    -1
#define RADAU_ERROR_INCONSISTENT_INPUT          -2
#define RADAU_ERROR_SUPERLU_NOT_ENABLED         -3
#define RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE   -4
#define RADAU_ERROR_NMAX_TOO_SMALL              -5
#define RADAU_ERROR_STEPSIZE_TOO_SMALL          -6
#define RADAU_ERROR_JAC_SINGULAR                -7
#define RADAU_ERROR_REP_STEP_REJECT             -8
#define RADAU_ERROR_NNZ_TOO_SMALL               -9
#define RADAU_ERROR_CALLBACK_UNRECOVERABLE      -10
#define RADAU_ERROR_CALLBACK_JAC_FORMAT         -11
#define RADAU_ERROR_CALLBACK_INVALID_NNZ        -12
#define RADAU_ERROR_UNEXPECTED_SUPERLU_FAILURE  -13
#define RADAU_ERROR_DENSE_CALLBACK              -14

/* FP_CB = FunctionPointer_CallBack */
/* See comments on radau5_solve in radau5.c for guide on the inputs */
typedef int (*FP_CB_f)(int, double, double*, double*, void*);
typedef int (*FP_CB_jac)(int, double, double*, double*, void*);
typedef int (*FP_CB_solout)(int, double, double*, double*, double*, int, void*);
typedef int (*FP_CB_jac_sparse)(int, double, double*, int*, double*, int*, int*, void*);

/* forward declarations of data structures */
struct radau_mem_t;
struct radau_linsol_mem_t;
struct radau_stats_t;
struct radau_inputs_t;
struct radau_math_const_t;

/* shorthands */
typedef struct radau_linsol_mem_t radau_linsol_mem_t;
typedef struct radau_stats_t radau_stats_t;
typedef struct radau_inputs_t radau_inputs_t;
typedef struct radau_math_const_t radau_math_const_t;
typedef struct radau_mem_t radau_mem_t;

/* Main Radau data structure */
struct radau_mem_t{
	int n; /* problem size */

	double *werr; /* local error estimate*/
	double *z1, *z2, *z3; /* transformed state vector */
	double *y0; /* auxiliary array for rhs eval at full time-points */
	double *scal; /* norm scaling factors */
	double *f1, *f2, *f3; /* newton rhs */
	double *cont; /* interpolation*/
	double *rtol, *atol; /* internal rtol and atol vectors */

	char err_log[256]; /* logging error message in case of failure */

	/* function calls */
	FP_CB_solout solout; /* callback function, called after each successful timestep */
	void *solout_ext; /* extra input to callback function */

	/* derived mathematical constants */
	radau_math_const_t *mconst;

	/* linear solver related data */
	radau_linsol_mem_t *lin_sol;

	/* input parameters relevant for solver control */
	radau_inputs_t *input;

	/* solver statistics */
	radau_stats_t *stats;

	/* dense output */
	int _dense_output_valid; /* switch if solver is in callback, since only then dense output calls are allowed */
	double xsol;  /* latest solution point */
	double hsol;  /* latest successful stepsize */

    /* internal parameters */
    /* these can affect the solver of several radau_solve() calls */
	int jac_is_fresh; /* flag if computed jacobian is fresh, used in fallbacks */
	double h_old; /* previous step-size, used in predictive controller */
	double erracc; /* bound in predictive controller */
    double faccon; /* factor tracked over several timestep to track jacobian recomputes in newton */

	int new_jac_req; /* flag if new jacobian is required */

	/* need to be stored in case LU factorization from previous solve calls is used */
	double fac1; /* diagonal factor for last real LU factorization */
	double alphn;  /* real diagonal factor for last complex LU factorization */
	double betan; /* complex diagonal factor for last complex LU factorization */
};

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
    int *jac_indices, *jac_indptr; /* sparse structure info */

	/* superlu auxiliary structs */
	void *slu_aux_d, *slu_aux_z; 
};


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


/* structure for input parameters, both internal and input, for Radau5 */
struct radau_inputs_t{
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
	int hmax_set; /* flag if hmax has been set manually */
	int fnewt_set; /* flag if fnewt has been set manually */

	double step_size_safety; /* safety factor for stepsize control */
	double uround; /* machine epsilon; 1.0D0+UROUND>1.0D0 */
	double theta_jac_recomp; /* theta; factor for jacobian recomputation */
	double fnewt; /* fnewt; newton cancellation factor: kappa*tol */
	double quot1, quot2; /* quot1; if quot1 < HNEW/HOLD < quot2, stepsize is not changed */
	double hmax; /* maximal step-size */
	double fac_lower, fac_upper; /* maximal limit for step-size (in|de)crease */
};

/* structure of (computed) mathematical constants inside Radau5 */
struct radau_math_const_t{
	double expm;
	double expm_inv;

	double c1, c2, c1mc2;
	double dd1, dd2, dd3;
	double u1, alph, beta;

	double c1m1, c2m1;
};

/* Error messages used in multiple places */
#define MSG_MALLOC_FAIL "Unexpected malloc failure."
#define MSG_MEM_NULL    "Unexpected NULL pointer for Radau memory structure."

#endif /*_RADAU5_IMPL_H*/
