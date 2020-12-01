#include <stdlib.h>
#include <math.h>
#include "amg_matrix.h"
# include "../../amg_precond.h"
# include "../../../Allocation_Operations/allocations.h"
# include "../../../BLAS_Operations/ourBLAS.h"

void precond_sol(double *q, double *p, int precond, AMG_precond_data *data) { // solve Mp = q

    precondAMG *amg_s = data->AMG_data;
    precondDPA *amg_d = data->DPA_data;
    switch (precond) { // preconditioning
        case 1: // AMG
            data->f[0] = q;
            data->u[0] = p;
            AMG_Vcycle_precond(data->NCL, data->A, data->f, data->u, data->r, amg_s->I_cf, amg_s->I_cf_t, data->omega, data->nr);
            break;
        case 2: // DPA AMG
            data->f[0] = q;
            data->u[0] = p;
            dpa_Vcycle_precond(data->NCL, data->A, data->f, data->u, data->r, data->omega, data->nr, amg_d->match);
            break;
        case 3: // XDPA
            data->f[0] = q;
            data->u[0] = p;
            xdpa_precond(data->NCL, data->A, data->f, data->u, data->r, data->omega, data->nr, amg_d->match);
            break;
            
    }
};

void matvec_product(double *p, matrix *A, double *v, int precond, AMG_precond_data *data) {
    double *q;
    if (precond) {
        q = (double *) mycalloc("q of matvec_product",A->m , sizeof (double));
        mat_vec(q, A, v);
        precond_sol(q, p, precond, data);
        myfree(q);
    } else {
        mat_vec(p, A, v);
    }
};

void residual(double *r, double *f, matrix *A, double *u, int precond, AMG_precond_data *data) {
    int i, n = A->m;
    matvec_product(r, A, u, precond, data);
    for (i = 0; i < n; i++) r[i] = f[i] - r[i];
};

double norm_inf(double *v, int n) {
    int i;
    double max = fabs(v[0]);
    for (i = 1; i < n; i++) {
        if (fabs(v[i]) > max) max = fabs(v[i]);
    }
    return max;
};

double norm_euclid(double *v, int n) {
    int i;
    double sum = 0.0;
    for (i = 0; i < n; i++) sum += v[i] * v[i];
    return sqrt(sum);
};

int GMRES(matrix *A, double *f, double *u, int k, double tol, int lmax) {
    int ind, i, j, l, n = A->m;
    double **U, **H, *e, *y, *c, *s, tol2, delta, aux, r;
    U = (double **) malloc((k + 1) * sizeof (double *));
    for (ind = 0; ind < k + 1; ind++) U[ind] = (double *) calloc(n, sizeof (double));
    H = (double **) malloc((k + 1) * sizeof (double *));
    for (ind = 0; ind < k + 1; ind++) H[ind] = (double *) malloc(k * sizeof (double));
    e = (double *) malloc((k + 1) * sizeof (double));
    y = (double *) malloc(k * sizeof (double));
    c = (double *) malloc(k * sizeof (double));
    s = (double *) malloc(k * sizeof (double));

    tol2 = tol * norm_euclid(f, n);
    l = 0;
    do {
		i=0;
        residual(U[i], f, A, u,0,0);
        e[i] = norm_euclid(U[i], n);
        for (aux = 1.0 / e[i], ind = 0; ind < n; ind++) U[i][ind] *= aux;
        delta = e[i];

        for (; i < k && delta > tol2; i++, l++) {

            mat_vec(U[i + 1], A, U[i]);
            // Gram-Schmidt orthogonalization
            for (j = 0; j <= i; j++) {
                H[j][i] = ddot(n, U[i + 1], U[j]);
                daxpy(n, -H[j][i], U[j], U[i + 1]);
            }
            H[i + 1][i] = norm_euclid(U[i + 1], n);
            for (aux = 1.0 / H[i + 1][i], ind = 0; ind < n; ind++) U[i + 1][ind] *= aux;
            // QR algorithm
            for (j = 0; j <= i - 1; j++) {
                aux = H[j][i];
                H[j][i] = c[j] * H[j][i] + s[j] * H[j + 1][i];
                H[j + 1][i] = -s[j] * aux + c[j] * H[j + 1][i];
            }
            r = sqrt(H[i][i] * H[i][i] + H[i + 1][i] * H[i + 1][i]);
            c[i] = H[i][i] / r;
            s[i] = H[i + 1][i] / r;
            H[i][i] = r; // H[i+1][i] = 0.0;
            e[i + 1] = -s[i] * e[i];
            e[i] = c[i] * e[i];
            delta = fabs(e[i + 1]);
        }
        i--;
        for (j = i; j >= 0; j--) {
            for (aux = 0.0, ind = j + 1; ind <= i; ind++) aux += H[j][ind] * y[ind];
            y[j] = (e[j] - aux) / H[j][j];
        }
        for (ind = 0; ind < n; ind++) {
            for (j = 0; j <= i; j++) u[ind] += y[j] * U[j][ind];
        }
    } while (l < lmax && delta > tol2);


    free(s);
    free(c);
    free(y);
    free(e);
    for (ind = 0; ind < k + 1; ind++) free(H[ind]);
    free(H);
    for (ind = 0; ind < k + 1; ind++) free(U[ind]);
    free(U);
    return (delta <= tol2) ? 1 : 0;
};