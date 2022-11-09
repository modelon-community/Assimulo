#ifndef _RADAU5_IO_H
#define _RADAU5_IO_H

#include "radau5_impl.h"

/* setup radau memory structure with inputs that are required to be fixed. */
int radau_setup_mem(int n, int sparseLU, int nprocs, int nnz, void **mem_out);
/* re-initializes internal radau_mem, affects internal parameters & stats, but not inputs */ 
int radau_reinit(void *radau_mem); 
/* returns all solver statistics, e.g., as number of function evaluations */
int radau_get_stats(void *radau_mem, int *nfcn, int *njac, int *nsteps, int *naccpt, int *nreject, int * ludecomps, int *lusolves);
/* Get a detailed error message */
char *radau_get_err_msg(void *radau_mem);

/* INPUT PARAMETER SETTING */
/* Some parameters have immediate errors checks */
/* Others depend on each other and are only checked when solver is run */
int radau_set_nmax              (void *radau_mem, int val); /* maximal number of steps */
int radau_set_nmax_newton       (void *radau_mem, int val); /* max number of newton steps */
int radau_set_newton_startn     (void *radau_mem, int val); /* newton starting strategy switch */
int radau_set_pred_step_control (void *radau_mem, int val); /* predictive step-size control switch */

int radau_set_step_size_safety  (void *radau_mem, double val); /* safety factor in step-size control */
int radau_set_uround            (void *radau_mem, double val); /* machine epsilon */
int radau_set_theta_jac_recomp  (void *radau_mem, double val); /* theta; factor for jacobian recomputation */
int radau_set_fnewt             (void *radau_mem, double val); /* fnewt; newton cancellation factor: kappa*tol */
int radau_set_quot1             (void *radau_mem, double val); /* quot1; if quot1 < HNEW/HOLD < quot2, stepsize is not changed */ 
int radau_set_quot2             (void *radau_mem, double val); /* quot2; if quot1 < HNEW/HOLD < quot2, stepsize is not changed */
int radau_set_hmax              (void *radau_mem, double val); /* maximal step-size */
int radau_set_fac_lower         (void *radau_mem, double val); /* maximal factor for step-size decrease */
int radau_set_fac_upper         (void *radau_mem, double val); /* maximal factor for step-size increase */

/* TODO: Currently not functioning */
int radau_set_solout(void *radau_mem, FP_CB_solout solout, void *solout_ext); /* set callback function, this automatically enables callbacks */

/* free all memory and delete structure */
void radau_free_mem(void **radau_mem);

#endif /*_RADAU5_IO_H*/
