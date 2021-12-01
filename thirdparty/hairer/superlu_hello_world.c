// Based on superlu.c example from userguide

#include "slu_mt_ddefs.h"
#include <stdio.h>

struct SuperLU_aux{
    int_t nprocs;

    int_t n, nnz;

    int_t panel_size, relax; // System specific tuning parameters
    fact_t fact; // if factorized matrix is being supplied, if not: how to factorize
    trans_t trans; // whether to solve transposed system or not
    yes_no_t refact; // NO for first time fac., YES for re-factorization
    yes_no_t usepr; // Whether the pivoting will use perm_r specified by the user, NO = it becomes output of pdgstrf function

    double diag_pivot_thresh, drop_tol;
    int_t lwork;
    void *work;

    int_t *perm_r, *perm_c; // row & column permutation vectors
    Gstat_t *Gstat;
    superlumt_options_t *slu_options;

    SuperMatrix *A, *B, *AC, *L, *U;
};

// TODO: where should nnz go; init or setup?
int superlu_init(int_t nprocs, struct SuperLU_aux *slu_aux,
                int_t n, int_t nnz);

int superlu_setup(struct SuperLU_aux *slu_aux,
                double *A_data, int_t *A_asub, int_t *A_xa);

int_t superlu_factorize(struct SuperLU_aux *slu_aux, int first);

int_t superlu_solve(struct SuperLU_aux *slu_aux, double *rhs);

int superlu_finalize(struct SuperLU_aux *slu_aux, double *rhs);

int main(int argc, char *argv[])
{
    // TODO: particular printf option for int_t?
    // Misc
    int i; 
    int_t info;

    // User input matrix and rhs vector
    double *a, *rhs;
    // Some constants for setting up the matrix
    double s, u, p, e, r, l;
    int_t n = 5; // size of (square) matrix and rhs vector
    int_t nnz = 12; // nnz = Number Non Zero
    int_t *asub, *xa; // indices & indptr

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

    /* ---------- ALLOCATE SPACE & INIT RHS ---------- */
    if ( !(rhs = doubleMalloc(n)) ) SUPERLU_ABORT("Malloc fails for rhs[].");
    for (i = 0; i < n; ++i){
        rhs[i] = 1.0;
    } 

    struct SuperLU_aux slu_aux;
    int_t nprocs = 1;
    // initialize various parameters and data structures
    superlu_init(nprocs, &slu_aux, n, nnz);

    // set up matrices and rhs
    superlu_setup(&slu_aux, a, asub, xa);

    superlu_factorize(&slu_aux, 1);
    
    info = superlu_solve(&slu_aux, rhs);
    if (info != 0){
        return (int) info;
    }

    printf("Expected solution: (-0.031250  0.065476  0.013393  0.062500  0.032738) \n");
    printf("Solution: (");
    for(i = 0; i < n; ++i){
        printf(" %f ", rhs[i]);
    }
    printf(" )\n");

    int j = 0;
    int mmax = 10;
    int fac = 1;
    for(j = 0; j < mmax; j++){
        fac *= 2;
        for (i = 0; i < n; ++i){
            rhs[i] = fac;
        }
        info = superlu_solve(&slu_aux, rhs);
        if (info != 0){
            return (int) info;
        }
        printf("Solution: (");
        for(i = 0; i < n; ++i){
            printf(" %f ", rhs[i]);
        }
        printf(" ) \n");
    }

    superlu_finalize(&slu_aux, rhs);

    return 0;
}

int superlu_init(int_t nprocs, struct SuperLU_aux *slu_aux,
                 int_t n, int_t nnz){
    // TODO: Are any failures in here recoverable?
    slu_aux->nprocs = nprocs;
    slu_aux->n = n;
    slu_aux->nnz = nnz;

    slu_aux->panel_size = sp_ienv(1);
    slu_aux->relax = sp_ienv(2);

    if (!(slu_aux->A = (SuperMatrix *)malloc(sizeof(SuperMatrix)))) SUPERLU_ABORT("Malloc failed for A.");
    if (!(slu_aux->B = (SuperMatrix *)malloc(sizeof(SuperMatrix)))) SUPERLU_ABORT("Malloc failed for B.");
    if (!(slu_aux->AC = (SuperMatrix *)malloc(sizeof(SuperMatrix)))) SUPERLU_ABORT("Malloc failed for AC.");
    if (!(slu_aux->L = (SuperMatrix *)malloc(sizeof(SuperMatrix)))) SUPERLU_ABORT("Malloc failed for L.");
    if (!(slu_aux->U = (SuperMatrix *)malloc(sizeof(SuperMatrix)))) SUPERLU_ABORT("Malloc failed for U.");

    slu_aux->Gstat = (Gstat_t *)malloc(sizeof(Gstat_t));
    slu_aux->slu_options = (superlumt_options_t *)malloc(sizeof(superlumt_options_t));

    slu_aux->fact = DOFACT; // if factorized matrix is being supplied, if not: how to factorize
    // other option: EQUIBRILATE: Scale row/colums to unit norm; good if matrix is poorly scaled?
    slu_aux->trans = NOTRANS; // whether to solve transposed system or not
    slu_aux->refact = NO; // NO for first time, YES for re-factorization
    slu_aux->diag_pivot_thresh = 1.0; // Default
    slu_aux->usepr = NO; // Whether the pivoting will use perm_r specified by the user, NO = it becomes output of pdgstrf function
    slu_aux->drop_tol = 0.0; // Default, not implemented for pdgstrf, TODO: Can skip?
    slu_aux->lwork = 0; // flag; work-space allocated internally
    slu_aux->work = NULL; // internal work space; not referenced due to lwork = 0

    StatAlloc(slu_aux->n, slu_aux->nprocs, slu_aux->panel_size, slu_aux->relax, slu_aux->Gstat);
    StatInit(slu_aux->n, slu_aux->nprocs, slu_aux->Gstat);

    if (!(slu_aux->perm_r = intMalloc(slu_aux->n))) SUPERLU_ABORT("Malloc failed for perm_r[].");
    if (!(slu_aux->perm_c = intMalloc(slu_aux->n))) SUPERLU_ABORT("Malloc failed for perm_c[].");

    dCreate_Dense_Matrix(slu_aux->B, slu_aux->n, 1, NULL, slu_aux->n, SLU_DN, SLU_D, SLU_GE);

    return 0;
}

int superlu_setup(struct SuperLU_aux *slu_aux,
                  double *A_data, int_t *A_asub, int_t *A_xa){
        // TODO: Need to clear storage in case of repeated setup?
        // if (slu_aux->A->Store){
        //     SUPERLU_FREE(slu_aux->A->Store);
        // }
        dCreate_CompCol_Matrix(slu_aux->A, slu_aux->n, slu_aux->n, slu_aux->nnz,
                               A_data, A_asub, A_xa, SLU_NC, SLU_D, SLU_GE);
        get_perm_c(3, slu_aux->A, slu_aux->perm_c); // 3 = approximate minimum degree for unsymmetric matrices
        return 0;
}

int_t superlu_factorize(struct SuperLU_aux *slu_aux, int first){
    int_t info;
    // TODO: enable functionality for re-factorization
    if (first){ // first time factorization
        pdgstrf_init(slu_aux->nprocs, slu_aux->fact, slu_aux->trans, slu_aux->refact, slu_aux->panel_size, slu_aux->relax,
                    slu_aux->diag_pivot_thresh, slu_aux->usepr, slu_aux->drop_tol, slu_aux->perm_c, slu_aux->perm_r,
                    slu_aux->work, slu_aux->lwork, slu_aux->A, slu_aux->AC, slu_aux->slu_options, slu_aux->Gstat); // initialize options

        pdgstrf(slu_aux->slu_options, slu_aux->AC, slu_aux->perm_r, slu_aux->L, slu_aux->U, slu_aux->Gstat, &info); // Factorization
    } else { // re-facatorization
        ;
    }
    return info;
}

int_t superlu_solve(struct SuperLU_aux *slu_aux, double *rhs){
    int_t info;
    DNformat *Bstore;

    Bstore = (DNformat *) slu_aux->B->Store;
    Bstore->nzval = rhs;

    dgstrs(slu_aux->trans, slu_aux->L, slu_aux->U, slu_aux->perm_r, slu_aux->perm_c, slu_aux->B, slu_aux->Gstat, &info); // Solve
    return info;
}

int superlu_finalize(struct SuperLU_aux *slu_aux, double *rhs){
    SUPERLU_FREE (rhs);
    SUPERLU_FREE (slu_aux->perm_r);
    SUPERLU_FREE (slu_aux->perm_c);
    Destroy_SuperNode_Matrix(slu_aux->L);
    Destroy_CompCol_Matrix(slu_aux->U);
    Destroy_CompCol_Matrix(slu_aux->A);
    Destroy_SuperMatrix_Store(slu_aux->B);
    StatFree(slu_aux->Gstat);

    SUPERLU_FREE(slu_aux->slu_options->etree);
    SUPERLU_FREE(slu_aux->slu_options->colcnt_h);
    SUPERLU_FREE(slu_aux->slu_options->part_super_h);
    Destroy_CompCol_Permuted(slu_aux->AC);
    return 0;
}