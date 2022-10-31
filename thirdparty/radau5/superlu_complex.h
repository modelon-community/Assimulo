#ifndef _SUPERLU_COMPLEX_H_
#define _SUPERLU_COMPLEX_H_

struct SuperLU_aux_z;
typedef struct SuperLU_aux_z SuperLU_aux_z;

#define min(a,b) ((a) <= (b) ? (a) : (b))

SuperLU_aux_z* superlu_init_z(int, int, int);
int superlu_setup_z(SuperLU_aux_z*, double, double, double *, int *, int *, int, int);
int superlu_factorize_z(SuperLU_aux_z*);
int superlu_solve_z(SuperLU_aux_z*, double *, double*);
int superlu_finalize_z(SuperLU_aux_z*);

#endif
