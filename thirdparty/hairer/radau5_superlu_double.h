#ifndef _TEST_PROXY_1_H_
#define _TEST_PROXY_1_H_

struct SuperLU_aux_d;
typedef struct SuperLU_aux_d SuperLU_aux_d;

typedef int (*CB_assemble_sys_d)(int, double, int *, double *, int *, int *, double *, int *, int *, int, double*);

SuperLU_aux_d* superlu_init_d(int, int, int);
int superlu_setup_d(SuperLU_aux_d *, double, double *, int *, int *, int, double*, CB_assemble_sys_d);
int superlu_factorize_d(SuperLU_aux_d *, int);
int superlu_solve_d(SuperLU_aux_d *, double *);
int superlu_finalize_d(SuperLU_aux_d *);

#endif
