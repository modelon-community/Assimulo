#include "superlu_util.h"

int sparse_csc_add_diagonal(int n, int* nnz, double* jac_data, int* jac_indices, int* jac_indptr){
    // Takes a sparse (CSC) matrix and adds empty diagonal entries, if not already present.
    // Works under the assumption that corresponding arrays have sufficient memory allocated.
    int i, j, diag_found, diag_added, indptr_incr, intptr_base, idx_data;

    // Quick version: Scan through jacobian and see if any diagonal entries are missing
    for(i = 0; i < n; i++){
        diag_found = 0;
        intptr_base = jac_indptr[i];
        for(j = 0; j < jac_indptr[i+1] - intptr_base; j++){
            if (jac_indices[intptr_base + j] ==  i){ // diagonal entry found
                diag_found = 1;
                break;
            }
        }

        if (!diag_found){ // diagonal not found for a particular column
            break;
        }
    }

    if (diag_found){ // above loop ended with all diagonals found, diagonal already part of matrix
        return 0;
    }
    // else: diagonal needs to be added
    
    // Copy over original jacobian data; inputs also serve as outputs
    double* jac_data_original = (double*) malloc(*nnz * sizeof(double));
    int* jac_indices_original = (int*) malloc(*nnz * sizeof(int));
    int* jac_indptr_original = (int*) malloc((n + 1) * sizeof(int));

    if ((!jac_data_original) || (!jac_indices_original) || (!jac_indptr)) {
        return MALLOC_FAILURE;
    }
    
    for(i = 0; i < *nnz; i++){
        jac_data_original[i] = jac_data[i];
    }
    for(i = 0; i < *nnz; i++){
        jac_indices_original[i] = jac_indices[i];
    }
    for(i = 0; i < n + 1; i++){
        jac_indptr_original[i] = jac_indptr[i];
    }

    // iterate through matrix; copy over data and possibly insert missing diagonals
    indptr_incr = 0; // increment of indptr from previous diagonal insertions
    idx_data = 0; // idx to track output positions for jac_data, jac_indices
    for(i = 0; i < n; i++){
        diag_found = 0;
        diag_added = 0;
        intptr_base = jac_indptr_original[i];
        for(j = 0; j < jac_indptr_original[i+1] - intptr_base; j++){
            if ((!diag_found) && (jac_indices_original[intptr_base + j] == i)){
                diag_found = 1;
            }
            if ((!diag_found) && (!diag_added) && (jac_indices_original[intptr_base + j] > i)){ // diagonal position exceeded; insert diagonal
                jac_data[idx_data] = 0.;
                jac_indices[idx_data] = i;
                idx_data++;
                indptr_incr++;
                diag_added = 1;
            }
            jac_data[idx_data] = jac_data_original[intptr_base + j];
            jac_indices[idx_data] = jac_indices_original[intptr_base + j];
            idx_data++;
        }

        if ((!diag_found) && (!diag_added)){ // column done, but diagonal neither found nor added
            jac_data[idx_data] = 0.;
            jac_indices[idx_data] = i;
            idx_data++;
            indptr_incr++;
            diag_added = 1;
        }
        jac_indptr[i+1] = jac_indptr_original[i+1] + indptr_incr;
    }
    *nnz += indptr_incr; // increment nnz by number of added diagonal elements
    free(jac_data_original);
    free(jac_indices_original);
    free(jac_indptr_original);
    return 0;
}
