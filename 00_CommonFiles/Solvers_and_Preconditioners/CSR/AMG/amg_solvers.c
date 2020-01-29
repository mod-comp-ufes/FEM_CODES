#include "amg_matrix.h"

void SOR_relax(matrix *A, double *f, double *u, double omega, int iter) {
    int i, j, k, n = A->m;
    double sum;
    for (k = 0; k < iter; k++) {
        for (i = 0; i < n; i++) {
            sum = 0.0;
            for (j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
                if (A->col_ind[j] != i) {
                    sum += A->val[j] * u[A->col_ind[j]];
                }
            }
            u[i] += omega * ((f[i] - sum) / A->diag[i] - u[i]);
        }
    }
};

void SOR_relax_rev(matrix *A, double *f, double *u, double omega, int iter) {
    int i, j, k, n = A->m;
    double sum;
    for (k = 0; k < iter; k++) {
        for (i = n - 1; i >= 0; i--) {
            sum = 0.0;
            for (j = A->row_ptr[i]; j < A->row_ptr[i + 1]; j++) {
                if (A->col_ind[j] != i) {
                    sum += A->val[j] * u[A->col_ind[j]];
                }
            }
            u[i] += omega * ((f[i] - sum) / A->diag[i] - u[i]);
        }
    }
};
