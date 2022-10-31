#ifndef _SUPERLU_DOUBLE_H_
#define _SUPERLU_DOUBLE_H_

struct SuperLU_aux_d;
typedef struct SuperLU_aux_d SuperLU_aux_d;

#define min(a,b) ((a) <= (b) ? (a) : (b))

SuperLU_aux_d* superlu_init_d(int, int, int);
int superlu_setup_d(SuperLU_aux_d*, double, double *, int *, int *, int, int);
int superlu_factorize_d(SuperLU_aux_d*);
int superlu_solve_d(SuperLU_aux_d*, double *);
int superlu_finalize_d(SuperLU_aux_d*);

#endif
