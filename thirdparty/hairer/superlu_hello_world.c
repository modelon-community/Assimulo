// Based on superlu.c example from userguide

#include "slu_mt_ddefs.h"
#include <stdio.h>

main(int argc, char *argv[])
{
    // TODO: particular printf option for int_t?
    // Misc
    int i; 

    SuperMatrix A, L, U, B;
    // User input matrix and rhs vector
    double *a, *rhs, *sol;
    // Some constants for setting up the matrix
    double s, u, p, e, r, l;
    // sizes for a and rhs
    int_t m = 5, n = 5; 
    // nnz = Number Non Zero
    int_t nnz = 12; 
    // required for storing matrix in compressed column format
    int_t *asub, *xa;

    /* ---------- ALLOCATE SPACE FOR MATRIX ---------- */
    if ( !(a = doubleMalloc(nnz)) ) SUPERLU_ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) SUPERLU_ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) SUPERLU_ABORT("Malloc fails for xa[].");

    /* ---------- INITIALIZE MATRIX ---------- */
    // set the various constants
    s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
    // initialize matrix in compressed column format
    a[0] = s; a[1] = l; a[2] = l; a[3] = u; a[4] = l; a[5] = l;
    a[6] = u; a[7] = p; a[8] = u; a[9] = e; a[10]= u; a[11]= r;
    asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
    asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
    asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
    xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;

    /* ---------- Create matrix A in the format expected by SuperLU.  ----------*/
    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

    /* ---------- ALLOCATE SPACE FOR RHS ---------- */
    int nrhs = 1; // size of rhs, here: simple vector
    if ( !(rhs = doubleMalloc(m * nrhs)) ) SUPERLU_ABORT("Malloc fails for rhs[].");

    /* ---------- INITIALIZE RHS ---------- */
    for (i = 0; i < m; ++i){
        rhs[i] = 1.0;
    } 
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    int_t *perm_r, *perm_c; // row & column permutation vectors
    if (!(perm_r = intMalloc(A.nrow))) SUPERLU_ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(A.ncol))) SUPERLU_ABORT("Malloc fails for perm_c[].");
    int_t permc_spec = 3; // approximate minimum degree for unsymmetric matrices
    get_perm_c(permc_spec, &A, perm_c);

    Gstat_t Gstat;
    superlumt_options_t superlumt_options;
    SuperMatrix AC;

    int_t info;
    int_t nprocs = 1;
    // Other option: EQUIBRILATE: Scale row/colums to unit norm; good if matrix is poorly scaled?
    fact_t fact = DOFACT; // if factorized matrix is being supplied, if not: how to factorize
    trans_t trans = NOTRANS; // whether to solve transposed system or nott
    yes_no_t refact = NO; // NO for first time, YES for re-factorization
    int_t panel_size = sp_ienv(1); // Tuning parameter, system specific
    int_t relax = sp_ienv(2); // Tuning parameter, system specific
    double diag_pivot_thresh = 1.0; // Default
    yes_no_t usepr = NO; // Whether the pivoting will use perm_r specified by the user, NO = it becomes output of pdgstrf function
    double drop_tol = 0.0; // Default, not implemented for pdgstrf, TODO: Can skip?
    int_t lwork = 0; // flag; work-space allocated internally
    void *work = NULL; // internal work space; not referenced due to lwork = 0

    printf("Allocating stats ... \n");
    StatAlloc((int_t) n, nprocs, panel_size, relax, &Gstat);
    printf("Allocating stats DONE \n");

    printf("Initializing stats ... \n");
    StatInit((int_t) n, nprocs, &Gstat); // initializing stats
    printf("Initializing stats DONE \n");

    printf("Initializing Factorization ... \n");
    pdgstrf_init(nprocs, fact, trans, refact, panel_size, relax,
	            diag_pivot_thresh, usepr, drop_tol, perm_c, perm_r,
	            work, lwork, &A, &AC, &superlumt_options, &Gstat); // initialize options
    printf("Initializing Factorization DONE \n");

    printf("Factorizing ... \n");
    pdgstrf(&superlumt_options, &AC, perm_r, &L, &U, &Gstat, &info); // Factorization
    printf("Factorizing DONE \n");
    
    printf("Solving ... \n");
    dgstrs(trans, &L, &U, perm_r, perm_c, &B, &Gstat, &info); // Solve
    printf("Solving DONE \n");

    // solution = (-0.031250  0.065476  0.013393  0.062500  0.032738)
    printf("Solution: (");
    for(i = 0; i < m; ++i){
        printf(" %f ", rhs[i]);
    }
    printf(" )");
    printf("   info : %i \n", info);

    int j = 0;
    int mmax = 10;
    int fac = 1;
    for(j = 0; j < mmax; j++){
        fac *= 2;
        for (i = 0; i < m; ++i){
            rhs[i] = fac;
        }
        dgstrs(trans, &L, &U, perm_r, perm_c, &B, &Gstat, &info); // Solve
        printf("Solution: (");
        for(i = 0; i < m; ++i){
            printf(" %f ", rhs[i]);
        }
        printf(" )");
        printf("   info : %i \n", info);
    }

    // /* De-allocate storage */
    SUPERLU_FREE (rhs);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    Destroy_CompCol_Matrix(&AC);
    StatFree(&Gstat);

    SUPERLU_FREE (a);
    SUPERLU_FREE (xa);
    SUPERLU_FREE (asub);
}
