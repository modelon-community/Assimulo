#include "radau5_io.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef __RADAU5_WITH_SUPERLU
#include "superlu_double.h"
#include "superlu_complex.h"
#endif /*__RADAU5_WITH_SUPERLU*/

#define TRUE_ (1)
#define FALSE_ (0)

/* forward declarations of private functions */
static int  _radau_setup_linsol_mem(radau_mem_t *rmem, int n, int sparseLU, int nprocs, int nnz);
static void _radau_setup_math_consts(radau_math_const_t* mconst);
static void _radau_reset_stats(radau_stats_t *mem);
static int  _radau_set_default_inputs(radau_inputs_t **mem_out);

static void _free_radau_linsol_mem(radau_linsol_mem_t **lin_sol_mem);
static void _free_radau_stats_mem(radau_stats_t **mem);
static void free_radau_inputs_mem(radau_inputs_t **mem);

/* setup radau memory structure with inputs that are required to be fixed. */
int radau_setup_mem(int n, int sparseLU, int nprocs, int nnz, void **mem_out){
	radau_mem_t *rmem;
	int ret = RADAU_OK;
	rmem = (radau_mem_t*)malloc(sizeof(radau_mem_t));
	if (!rmem){return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;}

	/* input sanity check */
	if (n < 1){
		sprintf(rmem->err_log, "Problem size must be positive integer, received n = %i", n);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}

	rmem->n = n; /* problem size */

	rmem->werr = (double*)malloc(n*sizeof(double));
	rmem->z1   = (double*)malloc(n*sizeof(double));
	rmem->z2   = (double*)malloc(n*sizeof(double));
	rmem->z3   = (double*)malloc(n*sizeof(double));
	rmem->y0   = (double*)malloc(n*sizeof(double));
	rmem->scal = (double*)malloc(n*sizeof(double));
	rmem->f1   = (double*)malloc(n*sizeof(double));
	rmem->f2   = (double*)malloc(n*sizeof(double));
	rmem->f3   = (double*)malloc(n*sizeof(double));
	rmem->cont = (double*)malloc(4*n*sizeof(double));
	rmem->rtol = (double*)malloc(n*sizeof(double));
	rmem->atol = (double*)malloc(n*sizeof(double));

	if(!rmem->werr || !rmem->z1 || !rmem->z3 || !rmem->y0 || !rmem->scal || !rmem->f1 || !rmem->f2 || !rmem->f3 || !rmem->cont){
		sprintf(rmem->err_log, MSG_MALLOC_FAIL);
		return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;
	}

	/* computed maths constants */
	rmem->mconst = (radau_math_const_t*)malloc(sizeof(radau_math_const_t));
	if(!rmem->mconst){
		sprintf(rmem->err_log, MSG_MALLOC_FAIL);
		return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;
	}
	_radau_setup_math_consts(rmem->mconst);

	/* Setup linear solver */
	ret = _radau_setup_linsol_mem(rmem, n, sparseLU, nprocs, nnz);
	if (ret < 0){ return ret;}

	/* Setup stats */
	rmem->stats = (radau_stats_t*)malloc(sizeof(radau_stats_t));
	if(!rmem->stats){
		sprintf(rmem->err_log, MSG_MALLOC_FAIL);
		return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;
	}
	/* reset - here: initialize - stats*/
	_radau_reset_stats(rmem->stats);

	/* Setup input by their default values */
	ret = _radau_set_default_inputs(&rmem->input);
	if (ret < 0){ return ret;}

    ret = radau_reinit((void*)rmem);
    if (ret < 0){ return ret;}

	*mem_out = (void*)rmem;
	return RADAU_OK;
} /* radau_setup_mem */

/* re-initializes internal radau_mem, affects internal parameters & stats, but not inputs */
int radau_reinit(void *radau_mem){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	/* default error message */
	sprintf(rmem->err_log, "No issues logged.");

	_radau_reset_stats(rmem->stats);

	/* dense output */
	rmem->_dense_output_valid = FALSE_;
	rmem->xsol = 0;  /* latest solution point */
	rmem->hsol = 0;  /* latest successful stepsize */

    /* internal parameters */
	rmem->jac_is_fresh = FALSE_;
	rmem->h_old = 0.; /* initial will be set in code */
	rmem->erracc = .01;
    rmem->faccon = 1;

	rmem->new_jac_req = TRUE_; /* by default: require new Jacobian  */

	rmem->fac1 = 0;
	rmem->alphn = 0;
	rmem->betan = 0;

	return RADAU_OK;
} /* radau_reinit */

/* returns all solver statistics, e.g., as number of function evaluations */
int radau_get_stats(void *radau_mem, int *nfcn, int *njac, int *nsteps, int *naccpt, int *nreject, int * ludecomps, int *lusolves){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	*nfcn = rmem->stats->nfcn;
	*njac = rmem->stats->njac;
	*nsteps = rmem->stats->nsteps;
	*naccpt = rmem->stats->naccpt;
	*nreject = rmem->stats->nreject;
	*ludecomps = rmem->stats->ludecomps;
	*lusolves = rmem->stats->lusolves;

	return RADAU_OK;
} /* radau_get_stats */

/* Get a detailed error message */
char *radau_get_err_msg(void *radau_mem){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){
		return "Unexpected NULL pointer for Radau memory structure when retrieving error message.";
	}else{
		return rmem->err_log;
	}
} /* radau_get_err_msg */

/* Set nmax parameter; maximal number of steps */
int radau_set_nmax(void *radau_mem, int val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	if (val <= 0){
		sprintf(rmem->err_log, "Input for nmax must be nonnegative, received nmax = %i.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->nmax = val;
	}
	return RADAU_OK;
} /* radau_set_nmax */

/* Set nmax_newton parameter; maximal number of newton steps */
int radau_set_nmax_newton(void *radau_mem, int val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	if (val <= 0){
		sprintf(rmem->err_log, "Input for nmax_newton must be nonnegative, received nmax = %i.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->nmax_newton = val;
	}
	return RADAU_OK;
} /* radau_set_nmax_newton */

/* Set newton_start_zero parameter; newton starting strategy */
int radau_set_newton_start_zero(void *radau_mem, int val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	rmem->input->newton_start_zero = (val != 0) ? TRUE_ : FALSE_;
	return RADAU_OK;
} /* radau_set_newton_start_zero */

/* Set pred_step_control parameter; step size control strategy */
int radau_set_pred_step_control(void *radau_mem, int val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	rmem->input->pred_step_control = (val != 0) ? TRUE_ : FALSE_;
	return RADAU_OK;
} /* radau_set_pred_step_control */

/* Set safety factor in timestep control */
int radau_set_step_size_safety(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	if (val <= .001 || val >= 1.){
		sprintf(rmem->err_log, "Input for step_size_safety must be between 0.001 and 1. received = %g.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->step_size_safety = val;
	}
	return RADAU_OK;
} /* radau_set_step_size_safety */

/* Set machine epsilon */
int radau_set_uround(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}
    rmem->input->_checked = FALSE_;

	if (val <= 1e-19 || val >= 1.){
		sprintf(rmem->err_log, "Input for uround must be between 1e-19 and 1. received = %e.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->uround = val;
	}
	return RADAU_OK;
} /* radau_set_uround */

/* Set theta; factor for jacobian recomputation */
int radau_set_theta_jac_recomp(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	if (val >= 1.){
		sprintf(rmem->err_log, "Input for theta_jac must be smaller 1. received = %e.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->theta_jac_recomp = val;
	}
	return RADAU_OK;
} /* radau_set_theta_jac_recomp */

/* Set fnewt; newton cancellation factor: kappa*tol */
int radau_set_fnewt(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}
    rmem->input->_checked = FALSE_;

	if (val <= 0.){
		sprintf(rmem->err_log, "Input for fnewt must be nonnegative, received = %g.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->fnewt = val;
		rmem->input->fnewt_set = TRUE_;
	}

	return RADAU_OK;
} /* radau_set_fnewt */

/* Set quot1; if quot1 < HNEW/HOLD < quot2, stepsize is not changed */
int radau_set_quot1(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	if (val > 1.){
		sprintf(rmem->err_log, "Input for quot1 must be larger 1, received = %g.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->quot1 = val;
	}
	return RADAU_OK;
} /* radau_set_quot1 */

/* Set quot2; if quot1 < HNEW/HOLD < quot2, stepsize is not changed */
int radau_set_quot2(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	if (val < 1.){
		sprintf(rmem->err_log, "Input for quot2 must be smaller 1, received = %g.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->quot2 = val;
	}
	return RADAU_OK;
} /* radau_set_quot2 */

/* Set maximal step-size */
int radau_set_hmax(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	if (val <= 0){
		sprintf(rmem->err_log, "Input for hmax must be positive, received = %g.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->hmax = val;
		rmem->input->hmax_set = TRUE_;
	}
	return RADAU_OK;
} /* radau_set_hmax */

/* Set maximal factor for step-size decrease */
int radau_set_fac_lower(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	if (val >= 1.){
		sprintf(rmem->err_log, "Input for fac_lower must be smaller 1, received = %g.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->fac_lower = 1./val;
	}
	return RADAU_OK;
} /* radau_set_fac_lower */

/* Set maximal factor for step-size increase */
int radau_set_fac_upper(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_ERROR_MEM_NULL;}

	if (val <= 1.){
		sprintf(rmem->err_log, "Input for fac_upper must be larger 1, received = %g.", val);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}else{
		rmem->input->fac_upper = 1./val;
	}
	return RADAU_OK;
} /* radau_set_fac_upper */

/* free all memory and delete structure */
void radau_free_mem(void **radau_mem){
	radau_mem_t *rmem = (radau_mem_t*) *radau_mem;
	if(!rmem){ return;}
	free(rmem->werr);
	free(rmem->z1);
	free(rmem->z2);
	free(rmem->z3);
	free(rmem->y0);
	free(rmem->scal);
	free(rmem->f1);
	free(rmem->f2);
	free(rmem->f3);
	free(rmem->cont);

	free(rmem->rtol);
	free(rmem->atol);

	free(rmem->mconst);
	_free_radau_linsol_mem(&rmem->lin_sol);
	_free_radau_stats_mem(&rmem->stats);
	free_radau_inputs_mem(&rmem->input);

	free(rmem);
} /* radau_free_mem */

/* setting up various computed mathematical constants used */
static void _radau_setup_math_consts(radau_math_const_t* mconst){
	double sq6 = sqrt(6.);
	double p1 = pow(81., .33333333333333331);
	double p2 = pow(9, .33333333333333331);
	double cno;

	mconst->c1 = (4. - sq6) / 10.;
    mconst->c2 = (sq6 + 4.) / 10.;

    mconst->c1mc2 = mconst->c1 - mconst->c2;
	mconst->c1m1 = mconst->c1 - 1;
	mconst->c2m1 = mconst->c2 - 1;

    mconst->dd1 = -(sq6 * 7. + 13.) / 3.;
    mconst->dd2 = (sq6 * 7. - 13.) / 3.;
    mconst->dd3 = -.33333333333333331;

    mconst->u1 = 1./((p1 + 6. - p2) / 30.);
    mconst->alph = (12. - p1 + p2) / 60.;
    mconst->beta = (p1 + p2) * sqrt(3.) / 60.;
    cno = mconst->alph * mconst->alph + mconst->beta * mconst->beta;
    mconst->alph /= cno;
    mconst->beta /= cno;
} /* _radau_setup_math_consts */

/* Setup linear solver related memory */
static int _radau_setup_linsol_mem(radau_mem_t *rmem, int n, int sparseLU, int nprocs, int nnz){
	int n_sq;
	rmem->lin_sol = (radau_linsol_mem_t*)malloc(sizeof(radau_linsol_mem_t));
	if (!rmem->lin_sol){
		sprintf(rmem->err_log, MSG_MALLOC_FAIL);
		return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;
	}
	
	/* input sanity check */
	if (n < 1){
		sprintf(rmem->err_log, "Problem size must be positive integer, received n = %i", n);
		return RADAU_ERROR_INCONSISTENT_INPUT;
	}
	n_sq = n*n;

	rmem->lin_sol->n = n;
	rmem->lin_sol->sparseLU = sparseLU;

	/* initialize all pointers with 0, since we do not use all */
	rmem->lin_sol->jac = 0;
	rmem->lin_sol->e1 = 0;
	rmem->lin_sol->e2r = 0;
	rmem->lin_sol->e2i = 0;
	rmem->lin_sol->ip1 = 0;
	rmem->lin_sol->ip2 = 0;
	rmem->lin_sol->jac_indices = 0;
	rmem->lin_sol->jac_indptr = 0;

	rmem->lin_sol->slu_aux_d = 0;
	rmem->lin_sol->slu_aux_z = 0;
	
	if(sparseLU){
		#ifdef __RADAU5_WITH_SUPERLU
			if (nnz < 0){
				sprintf(rmem->err_log, "Input nnz must be nonnegative, received nnz = %i", nnz);
				return RADAU_ERROR_INCONSISTENT_INPUT;
			}
			if (nprocs < 1){
				sprintf(rmem->err_log, "Input nprocs must be positive, received nprocs = %i", nprocs);
				return RADAU_ERROR_INCONSISTENT_INPUT;
			}
			rmem->lin_sol->nnz = nnz > n*n ? n*n : nnz; /* automatically cap at n*n */
			rmem->lin_sol->nnz_actual = rmem->lin_sol->nnz;
			rmem->lin_sol->nproc = nprocs;
			rmem->lin_sol->LU_with_fresh_jac = FALSE_;

			/* allocate memory for sparse Jacobian structure */
			rmem->lin_sol->jac = (double*)malloc((nnz + n)*sizeof(double));
			rmem->lin_sol->jac_indices = (int*)malloc((nnz + n)*sizeof(int));
			rmem->lin_sol->jac_indptr  = (int*)malloc((n + 1)*sizeof(int));

			/* Create auxiliary superLU structures */
			rmem->lin_sol->slu_aux_d = (void*)superlu_init_d(nprocs, n, nnz);
			rmem->lin_sol->slu_aux_z = (void*)superlu_init_z(nprocs, n, nnz);
			if(!rmem->lin_sol->jac || !rmem->lin_sol->jac_indices || !rmem->lin_sol->jac_indptr || !rmem->lin_sol->slu_aux_d || !rmem->lin_sol->slu_aux_z){
				sprintf(rmem->err_log, MSG_MALLOC_FAIL);
				return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;
			}
		#else
			sprintf(rmem->err_log, "sparseLU set to true, but code has not been compiled with '__RADAU5_WITH_SUPERLU' flag.");
			return RADAU_ERROR_SUPERLU_NOT_ENABLED;
		#endif /*__RADAU5_WITH_SUPERLU*/

	}else{ /* DENSE */
		rmem->lin_sol->jac = (double*)malloc(n_sq*sizeof(double));
		rmem->lin_sol->e1  = (double*)malloc(n_sq*sizeof(double));
		rmem->lin_sol->e2r = (double*)malloc(n_sq*sizeof(double));
		rmem->lin_sol->e2i = (double*)malloc(n_sq*sizeof(double));

		rmem->lin_sol->ip1 = (int*)malloc(n*sizeof(int));
		rmem->lin_sol->ip2 = (int*)malloc(n*sizeof(int));

		if(!rmem->lin_sol->jac || !rmem->lin_sol->e1 || !rmem->lin_sol->e2r || !rmem->lin_sol->e2i || !rmem->lin_sol->ip1 || !rmem->lin_sol->ip2){
			sprintf(rmem->err_log, MSG_MALLOC_FAIL);
			return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;
		}
	}
	return RADAU_OK;
} /* _radau_setup_linsol_mem */

/* Reset all runtime related statistics */
static void _radau_reset_stats(radau_stats_t *rmem){
	rmem->nfcn = 0;
	rmem->njac = 0;
	rmem->nsteps = 0;
	rmem->naccpt = 0;
	rmem->nreject = 0;
	rmem->ludecomps = 0;
	rmem->lusolves = 0;
} /* _radau_reset_stats */

/* Initialize inputs structure and assign default values */
static int _radau_set_default_inputs(radau_inputs_t **input_out){
	radau_inputs_t *mem = (radau_inputs_t*)malloc(sizeof(radau_inputs_t));
	if(!mem){ return RADAU_ERROR_UNEXPECTED_MALLOC_FAILURE;}

    mem->_checked = FALSE_;

	/* external parameters */
	mem->nmax = 100000;
	mem->nmax_newton = 7;

	mem->newton_start_zero= FALSE_;
	mem->pred_step_control = TRUE_;

	mem->hmax_set = FALSE_;
	mem->fnewt_set = FALSE_;

	mem->step_size_safety = .9;

	mem->uround = 1.e-16;
	mem->theta_jac_recomp = .001;
	mem->fnewt = .001; /* radau_max(uround * 10 / tolst, radau_min(.03, pow(tolst, .5)));*/
	mem->quot1 = 1.;
	mem->quot2 = 1.2;
	mem->hmax = 0;
	mem->fac_lower = 5.;
	mem->fac_upper = 0.125;

	*input_out = mem;
	return RADAU_OK;
} /* _radau_set_default_inputs */

/* Free all linear solver related memory */
static void _free_radau_linsol_mem(radau_linsol_mem_t **lin_sol_mem){
	radau_linsol_mem_t *mem = (radau_linsol_mem_t*) *lin_sol_mem;
	if(!mem){ return;}

	free(mem->jac);
	free(mem->e1);
	free(mem->e2r);
	free(mem->e2i);
	free(mem->ip1);
	free(mem->ip2);

	free(mem->jac_indices);
	free(mem->jac_indptr);

	#ifdef __RADAU5_WITH_SUPERLU
		superlu_finalize_d((SuperLU_aux_d*)mem->slu_aux_d);
		superlu_finalize_z((SuperLU_aux_z*)mem->slu_aux_z);
	#endif /* __RADAU5_WITH_SUPERLU */

	free(mem);
} /* _free_radau_linsol_mem */

/* Free stats related memorys */
static void _free_radau_stats_mem(radau_stats_t **stats_mem){
	radau_stats_t *mem = (radau_stats_t*) *stats_mem;
	if(!mem){ return;}

	free(mem);
} /* _free_radau_stats_mem */

/* Free inputs related memory */
static void free_radau_inputs_mem(radau_inputs_t **para_mem){
	radau_inputs_t *mem = (radau_inputs_t*) *para_mem;
	if(!mem){ return;}

	free(mem);
} /* free_radau_inputs_mem */
