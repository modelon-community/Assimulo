#include "radau5_io.h"
#include <stdlib.h>
#include <math.h>

#ifdef __RADAU5_WITH_SUPERLU
#include "superlu_double.h"
#include "superlu_complex.h"
#include "superlu_util.h"
#endif /*__RADAU5_WITH_SUPERLU*/

#define TRUE_ (1)
#define FALSE_ (0)

/* setting up various computed mathematical constants used */
static void setup_math_consts(radau_math_const_t* mconst){
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
}

/* set default values for internal parameters */
static void set_radau_para_default_internal(radau_parameters_t* mem){
	/* internal parameters */
	mem->jac_is_fresh = FALSE_;
	mem->h_old = 0.; /* will be set in code */
	mem->dynold = 0.;
	mem->erracc = .01;
	mem->thqold = 0;

	mem->new_jac_req = TRUE_; /* by default: require new Jacobian  */
}

int setup_radau_mem(int n, int sparseLU, int nprocs, int nnz, void **mem_out){
	radau_mem_t *rmem;
	int ret = RADAU_OK;
	
	/* input sanity checks */
	if (n < 1){return RADAU_SETUP_INVALID_INPUT;}
	rmem = (radau_mem_t*)malloc(sizeof(radau_mem_t));
	if (!rmem){return RADAU_SETUP_MALLOC_FAILURE;}

	rmem->n = n; /* problem size */

	rmem->work = (double*)malloc(20*sizeof(double)); /* TODO */
	rmem->werr = (double*)malloc(n*sizeof(double));
	rmem->z1 = (double*)malloc(n*sizeof(double));
	rmem->z2 = (double*)malloc(n*sizeof(double));
	rmem->z3 = (double*)malloc(n*sizeof(double));
	rmem->y0 = (double*)malloc(n*sizeof(double));
	rmem->scal = (double*)malloc(n*sizeof(double));
	rmem->f1 = (double*)malloc(n*sizeof(double));
	rmem->f2 = (double*)malloc(n*sizeof(double));
	rmem->f3 = (double*)malloc(n*sizeof(double));
	rmem->cont = (double*)malloc(4*n*sizeof(double));

	rmem->rtol = (double*)malloc(n*sizeof(double));
	rmem->atol = (double*)malloc(n*sizeof(double));

	if(!rmem->work || !rmem->werr || !rmem->z1 || !rmem->z3 || !rmem->y0 || !rmem->scal || !rmem->f1 || !rmem->f2 || !rmem->f3 || !rmem->cont){
		return RADAU_SETUP_MALLOC_FAILURE;
	}

	/* maths constants */
	rmem->mconst = (radau_math_const_t*)malloc(sizeof(radau_math_const_t));
	if(!rmem->mconst){ return RADAU_SETUP_MALLOC_FAILURE;}
	setup_math_consts(rmem->mconst);

	/* Setup linear solver */
	ret = setup_radau_linsol_mem(n, sparseLU, nprocs, nnz, &rmem->lin_sol);
	if (ret < 0){ return ret;}

	/* Setup stats */
	rmem->stats = (radau_stats_t*)malloc(sizeof(radau_stats_t));
	if(!rmem->stats){ return RADAU_SETUP_MALLOC_FAILURE;}
	/* reset - here: initialize - stats*/
	ret = reset_radau_stats((void*)rmem);
	if (ret < 0){ return ret;}

	/* Setup parameters, by default values */
	ret = setup_radau_para_default(&rmem->para);
	if (ret < 0){ return ret;}

	/* dense output */
	rmem->xsol = 0;  /* latest solution point */
	rmem->hsol = 0;  /* latest successful stepsize */

	*mem_out = (void*)rmem;
	return RADAU_OK;
}

int reset_radau_internal_mem(void *radau_mem){
	int ret;
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}

	ret = reset_radau_stats(radau_mem);
	if (ret < 0){ return RADAU_PARA_RADAU_MEM_NULL;}

	set_radau_para_default_internal(rmem->para);

	return RADAU_OK;
}


int setup_radau_linsol_mem(int n, int sparseLU, int nprocs, int nnz, radau_linsol_mem_t **mem_out){
	radau_linsol_mem_t *mem = (radau_linsol_mem_t*)malloc(sizeof(radau_linsol_mem_t));
	int n_sq = n*n;

	/* input sanity check */
	if (n < 1){return RADAU_SETUP_INVALID_INPUT;}

	mem->n = n;
	mem->sparseLU = sparseLU;

	/* initialize all pointers with 0, since we do not use all */
	mem->jac = 0;
	mem->e1 = 0;
	mem->e2r = 0;
	mem->e2i = 0;
	mem->ip1 = 0;
	mem->ip2 = 0;
	mem->jac_indices = 0;
	mem->jac_indptr = 0;

	mem->slu_aux_d = 0;
	mem->slu_aux_z = 0;
	
	if(sparseLU){
		#ifdef __RADAU5_WITH_SUPERLU
			if ((nnz < 0) || (nnz > n_sq) || (nprocs < 1)) {return RADAU_SETUP_INVALID_INPUT;}
			mem->nnz = nnz;
			mem->nnz_actual = nnz;
			mem->nproc = nprocs;
			mem->LU_with_fresh_jac = FALSE_;

			/* allocate memory for sparse Jacobian structure */
			mem->jac = (double*)malloc((nnz + n)*sizeof(double));
			mem->jac_indices = (int*)malloc((nnz + n)*sizeof(int));
			mem->jac_indptr = (int*)malloc((n + 1)*sizeof(int));

			/* Create auxiliary superLU structures */
			mem->slu_aux_d = (void*)superlu_init_d(nprocs, n, nnz);
			mem->slu_aux_z = (void*)superlu_init_z(nprocs, n, nnz);
			if(!mem->jac || !mem->jac_indices || !mem->jac_indptr || !mem->slu_aux_d || !mem->slu_aux_z){
				return RADAU_SETUP_MALLOC_FAILURE;
			}
		#else
			return RADAU_SETUP_SUPERLU_NOT_ENABLED;
		#endif /*__RADAU5_WITH_SUPERLU*/

	}else{ /* DENSE */
		mem->jac = (double*)malloc(n_sq*sizeof(double));
		mem->e1 = (double*)malloc(n_sq*sizeof(double));
		mem->e2r = (double*)malloc(n_sq*sizeof(double));
		mem->e2i = (double*)malloc(n_sq*sizeof(double));

		mem->ip1 = (int*)malloc(n*sizeof(int));
		mem->ip2 = (int*)malloc(n*sizeof(int));

		if(!mem->jac || !mem->e1 || !mem->e2r || !mem->e2i || !mem->ip1 || !mem->ip2){
			return RADAU_SETUP_MALLOC_FAILURE;
		}
	}
	*mem_out = mem;
	return RADAU_OK;
}


int reset_radau_stats(void *radau_mem){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}

	rmem->stats->nfcn = 0;
	rmem->stats->njac = 0;
	rmem->stats->nsteps = 0;
	rmem->stats->naccpt = 0;
	rmem->stats->nreject = 0;
	rmem->stats->ludecomps = 0;
	rmem->stats->lusolves = 0;

	return RADAU_OK;
}

/* returns all runtime stats */
int radau_get_stats(void *radau_mem, int *nfcn, int *njac, int *nsteps, int *naccpt, int *nreject, int * ludecomps, int *lusolves){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}

	*nfcn = rmem->stats->nfcn;
	*njac = rmem->stats->njac;
	*nsteps = rmem->stats->nsteps;
	*naccpt = rmem->stats->naccpt;
	*nreject = rmem->stats->nreject;
	*ludecomps = rmem->stats->ludecomps;
	*lusolves = rmem->stats->lusolves;

	return RADAU_OK;
}

int setup_radau_para_default(radau_parameters_t **mem_out){
	radau_parameters_t *mem = (radau_parameters_t*)malloc(sizeof(radau_parameters_t));
	if(!mem){ return RADAU_SETUP_MALLOC_FAILURE;}

	/* external parameters */
	mem->nmax = 100000;
	mem->nmax_newton = 7;

	mem->newton_start_zero= FALSE_;
	mem->pred_step_control = TRUE_;

	mem->step_size_safety = .9;

	set_radau_para_default_internal(mem);

	*mem_out = mem;
	return RADAU_OK;
}

/* Set nmax parameter; maximal number of steps */
int radau_set_para_nmax(void *radau_mem, int val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}

	if (val <= 0){
		return RADAU_PARA_INCONSISTENT_INPUT;
	}else{
		rmem->para->nmax = val;
	}
	return RADAU_OK;
}

/* Set nmax_newton parameter; maximal number of newton steps */
int radau_set_para_nmax_newton(void *radau_mem, int val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}

	if (val <= 0){
		return RADAU_PARA_INCONSISTENT_INPUT;
	}else{
		rmem->para->nmax_newton = val;
	}
	return RADAU_OK;
}

/* Set newton_start_zero parameter; newton starting strategy */
int radau_set_para_newton_start_zero(void *radau_mem, int val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}

	rmem->para->newton_start_zero = (val != 0) ? TRUE_ : FALSE_;
	return RADAU_OK;
}

/* Set pred_step_control parameter; step size control strategy */
int radau_set_para_pred_step_control(void *radau_mem, int val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}

	rmem->para->pred_step_control = (val != 0) ? TRUE_ : FALSE_;
	return RADAU_OK;
}

int radau_set_para_step_size_safety(void *radau_mem, double val){
	radau_mem_t *rmem = (radau_mem_t*)radau_mem;
	if (!rmem){ return RADAU_PARA_RADAU_MEM_NULL;}

	if (val <= .001 || val >= 1.){
		return RADAU_PARA_INCONSISTENT_INPUT;
	}else{
		rmem->para->step_size_safety = val;
	}
	return RADAU_OK;
}

void free_radau_mem(void **radau_mem){
	radau_mem_t *rmem = (radau_mem_t*) *radau_mem;
	if(!rmem){ return;}
	free(rmem->work);
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
	free_radau_linsol_mem(&rmem->lin_sol);
	free_radau_stats_mem(&rmem->stats);
	free_radau_parameters_mem(&rmem->para);

	free(rmem);
}


void free_radau_linsol_mem(radau_linsol_mem_t **lin_sol_mem){
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
}

void free_radau_stats_mem(radau_stats_t **stats_mem){
	radau_stats_t *mem = (radau_stats_t*) *stats_mem;
	if(!mem){ return;}

	free(mem);
}

void free_radau_parameters_mem(radau_parameters_t **para_mem){
	radau_parameters_t *mem = (radau_parameters_t*) *para_mem;
	if(!mem){ return;}

	free(mem);
}
