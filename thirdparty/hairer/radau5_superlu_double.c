#include "slu_mt_ddefs.h"
#include "radau5_superlu_double.h"

struct SuperLU_aux_d{
    int nprocs, n, nnz_jac;
    int setup_done, fact_done; // flags for which memory to free in the end

    double *data_sys;
    int *indices_sys, *indptr_sys;

    int panel_size, relax; // System specific tuning parameters
    fact_t fact; // if factorized matrix is being supplied, if not: how to factorize
    trans_t trans; // whether to solve transposed system or not
    yes_no_t refact; // NO for first time fac., YES for re-factorization
    yes_no_t usepr; // Whether the pivoting will use perm_r specified by the user, NO = it becomes output of pdgstrf function

    double diag_pivot_thresh, drop_tol;
    int lwork;
    void *work;

    int *perm_r, *perm_c; // row & column permutation vectors
    Gstat_t *Gstat;
    superlumt_options_t *slu_options;

    SuperMatrix *A, *B, *AC, *L, *U;
};

// Initialization of required data structures for SuperLU
SuperLU_aux_d* superlu_init_d(int nprocs, int n, int nnz){
    SuperLU_aux_d *slu_aux = (SuperLU_aux_d *)malloc(sizeof(SuperLU_aux_d));
    if (!slu_aux) {SUPERLU_ABORT("Malloc failed for double slu_aux.");}
    slu_aux->setup_done = 0;
    slu_aux->fact_done = 0;

    slu_aux->nprocs = nprocs;
    slu_aux->nnz_jac = nnz;
    slu_aux->n = n;

    slu_aux->panel_size = sp_ienv(1);
    slu_aux->relax = sp_ienv(2);

    slu_aux->A =  (SuperMatrix *)malloc(sizeof(SuperMatrix));
    slu_aux->B =  (SuperMatrix *)malloc(sizeof(SuperMatrix));
    slu_aux->AC = (SuperMatrix *)malloc(sizeof(SuperMatrix));
    slu_aux->L =  (SuperMatrix *)malloc(sizeof(SuperMatrix));
    slu_aux->U =  (SuperMatrix *)malloc(sizeof(SuperMatrix));

    if (!slu_aux->A)  {SUPERLU_ABORT("Malloc failed for double A.");}
    if (!slu_aux->B)  {SUPERLU_ABORT("Malloc failed for double B.");}
    if (!slu_aux->AC) {SUPERLU_ABORT("Malloc failed for double AC.");}
    if (!slu_aux->L)  {SUPERLU_ABORT("Malloc failed for double L.");}
    if (!slu_aux->U)  {SUPERLU_ABORT("Malloc failed for double U.");}
    slu_aux->A->Store = NULL;

    slu_aux->Gstat = (Gstat_t *)malloc(sizeof(Gstat_t));
    slu_aux->slu_options = (superlumt_options_t *)malloc(sizeof(superlumt_options_t));

    if (!slu_aux->Gstat) {SUPERLU_ABORT("Malloc failed for double Gstat.");}
    if (!slu_aux->slu_options) {SUPERLU_ABORT("Malloc failed for double slu_options.");}

    slu_aux->fact = DOFACT; // if factorized matrix is being supplied, if not: how to factorize
    // other option: EQUIBRILATE: Scale row/colums to unit norm; good if matrix is poorly scaled?
    slu_aux->trans = NOTRANS; // whether to solve transposed system or not
    slu_aux->refact = NO; // NO for first time, YES for re-factorization
    slu_aux->diag_pivot_thresh = 1.0; // Default
    slu_aux->usepr = NO; // Whether the pivoting will use perm_r specified by the user, NO = it becomes output of pdgstrf function
    slu_aux->drop_tol = 0.0; // Default, not implemented for pdgstrf
    slu_aux->lwork = 0; // flag; work-space allocated internally
    slu_aux->work = NULL; // internal work space; not referenced due to lwork = 0

    StatAlloc(slu_aux->n, slu_aux->nprocs, slu_aux->panel_size, slu_aux->relax, slu_aux->Gstat);
    StatInit(slu_aux->n, slu_aux->nprocs, slu_aux->Gstat);

    slu_aux->perm_r = intMalloc(slu_aux->n);
    slu_aux->perm_c = intMalloc(slu_aux->n);

    if (!slu_aux->perm_r) {SUPERLU_ABORT("Malloc failed for double perm_r[].");}
    if (!slu_aux->perm_c) {SUPERLU_ABORT("Malloc failed for double perm_c[].");}

    dCreate_Dense_Matrix(slu_aux->B, slu_aux->n, 1, NULL, slu_aux->n, SLU_DN, SLU_D, SLU_GE);

    // allocate memory for storing matrix of linear system
    // min(nnz_jac + n, n*n) is upper bound on storage requirement of linear system
    slu_aux->data_sys = doubleMalloc(min(slu_aux->nnz_jac + slu_aux->n, n*n));
    if (!slu_aux->data_sys)    {SUPERLU_ABORT("Malloc fails for double data_sys[].");}

    return slu_aux;
}

// Setting up the matrix to be factorized
int superlu_setup_d(SuperLU_aux_d *slu_aux, double scale,
                    double *data_J, int *indices_J, int *indptr_J,
                    int fresh_jacobian, int jac_nnz){
    NCformat *AStore = slu_aux->A->Store;
    SUPERLU_FREE(AStore);

    // number of non-zero elements maz have changed during recent jacobian evaluation
    slu_aux -> nnz_jac = jac_nnz;

    int current_idx = 0;
    int i, j;

    // build system matrix scale * I - JAC
    // Copy jacobian data to slu_aux struct
    slu_aux->indices_sys = indices_J;
    slu_aux->indptr_sys = indptr_J;

    for (i = 0; i < slu_aux-> n; i++){
        for (j = indptr_J[i]; j < indptr_J[i+1]; j++){
            slu_aux->data_sys[current_idx] = -data_J[current_idx];
            if (i == indices_J[current_idx]){
                slu_aux->data_sys[current_idx] = slu_aux->data_sys[current_idx] + scale;
            }
            current_idx++;
        }
    }

    dCreate_CompCol_Matrix(slu_aux->A, slu_aux->n, slu_aux->n, slu_aux->nnz_jac,
                           slu_aux->data_sys, slu_aux->indices_sys, slu_aux->indptr_sys,
                           SLU_NC, SLU_D, SLU_GE);

    if (fresh_jacobian){
        get_perm_c(3, slu_aux->A, slu_aux->perm_c); // 3 = approximate minimum degree for unsymmetrical matrices
        slu_aux->refact = NO; // new jacobian, do new factorization
    }else{
        slu_aux->refact = YES; // same jacobian structure, re-factorization 
    }
    slu_aux->setup_done = 1;
    return 0;
}

// Factorize matrix
int superlu_factorize_d(SuperLU_aux_d *slu_aux){
    int info;
    // clean up memory in case of re-factorization
    if (slu_aux->fact_done){
        NCPformat *ACstore = slu_aux->AC->Store;
        SUPERLU_FREE(ACstore->colend);
        SUPERLU_FREE(ACstore->colbeg);
        SUPERLU_FREE(ACstore);
        if (slu_aux->refact == NO){
            SCPformat *LStore = slu_aux->L->Store;
            SUPERLU_FREE(LStore->col_to_sup);
            SUPERLU_FREE(LStore->sup_to_colbeg);
            SUPERLU_FREE(LStore->sup_to_colend);
            Destroy_SuperNode_Matrix(slu_aux->L);

            NCPformat *UStore = slu_aux->U->Store;
            SUPERLU_FREE(UStore->colend);
            Destroy_CompCol_Matrix(slu_aux->U);

            SUPERLU_FREE(slu_aux->slu_options->etree);
            SUPERLU_FREE(slu_aux->slu_options->colcnt_h);
            SUPERLU_FREE(slu_aux->slu_options->part_super_h);
        }
    }
    // initialize options
    pdgstrf_init(slu_aux->nprocs, slu_aux->fact, slu_aux->trans, slu_aux->refact, slu_aux->panel_size, slu_aux->relax,
                 slu_aux->diag_pivot_thresh, slu_aux->usepr, slu_aux->drop_tol, slu_aux->perm_c, slu_aux->perm_r,
                 slu_aux->work, slu_aux->lwork, slu_aux->A, slu_aux->AC, slu_aux->slu_options, slu_aux->Gstat);
    // Factorization
    pdgstrf(slu_aux->slu_options, slu_aux->AC, slu_aux->perm_r, slu_aux->L, slu_aux->U, slu_aux->Gstat, &info);
    slu_aux->refact = YES;
    slu_aux->fact_done = 1;
    return info;
}

// Solve linear system based on previous factorization
int superlu_solve_d(SuperLU_aux_d *slu_aux, double *rhs){
    int info;
    DNformat *Bstore;

    Bstore = (DNformat *) slu_aux->B->Store;
    Bstore->nzval = rhs;
    // Solve
    dgstrs(slu_aux->trans, slu_aux->L, slu_aux->U, slu_aux->perm_r, slu_aux->perm_c, slu_aux->B, slu_aux->Gstat, &info);
    return info;
}

// de-allocate memory
int superlu_finalize_d(SuperLU_aux_d *slu_aux){
    SUPERLU_FREE(slu_aux->perm_r);
    SUPERLU_FREE(slu_aux->perm_c);

    Destroy_SuperMatrix_Store(slu_aux->B);
    StatFree(slu_aux->Gstat);

    if (slu_aux->setup_done){
        free(slu_aux->A->Store);
    }
    SUPERLU_FREE(slu_aux->data_sys);

    if (slu_aux->fact_done){
        SCPformat *LStore = slu_aux->L->Store;
        SUPERLU_FREE(LStore->col_to_sup);
        SUPERLU_FREE(LStore->sup_to_colbeg);
        SUPERLU_FREE(LStore->sup_to_colend);
        Destroy_SuperNode_Matrix(slu_aux->L);

        NCPformat *UStore = slu_aux->U->Store;
        SUPERLU_FREE(UStore->colend);
        Destroy_CompCol_Matrix(slu_aux->U);
        
        Destroy_CompCol_Permuted(slu_aux->AC);
        SUPERLU_FREE(slu_aux->slu_options->etree);
        SUPERLU_FREE(slu_aux->slu_options->colcnt_h);
        SUPERLU_FREE(slu_aux->slu_options->part_super_h);
    }

    free(slu_aux->A);
    free(slu_aux->B);
    free(slu_aux->L);
    free(slu_aux->U);
    free(slu_aux->AC);
    free(slu_aux->Gstat);
    free(slu_aux->slu_options);
    free(slu_aux->work);
    free(slu_aux);
    return 0;
}
