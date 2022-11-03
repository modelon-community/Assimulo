#ifndef _RADAU5_IO_H
#define _RADAU5_IO_H

#include "radau5_impl.h"

int setup_radau_mem(int n, int sparseLU, int nprocs, int nnz, void **mem_out);
int reset_radau_internal_mem(void *radau_mem); /* reset internal radau_mem, affects internal parameters & stats */ 
int setup_radau_linsol_mem(int n, int sparseLU, int nprocs, int nnz, radau_linsol_mem_t **mem_out);

int reset_radau_stats(void *radau_mem);
int radau_get_stats(void *radau_mem, int *nfcn, int *njac, int *nsteps, int *naccpt, int *nreject, int * ludecomps, int *lusolves);

int setup_radau_para_default(radau_parameters_t **mem_out);
int radau_set_para_nmax(void *radau_mem, int val);
int radau_set_para_nmax_newton(void *radau_mem, int val);
int radau_set_para_newton_startn(void *radau_mem, int val);
int radau_set_para_pred_step_control(void *radau_mem, int val);

int radau_set_para_step_size_safety(void *radau_mem, double val);


void free_radau_mem(void **radau_mem);
void free_radau_linsol_mem(radau_linsol_mem_t **mem);
void free_radau_stats_mem(radau_stats_t **mem);
void free_radau_parameters_mem(radau_parameters_t **mem);

#endif /*_RADAU5_IO_H*/
