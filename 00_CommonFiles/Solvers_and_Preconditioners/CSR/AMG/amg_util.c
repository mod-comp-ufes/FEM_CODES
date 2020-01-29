#include <stdlib.h>
#include <math.h>
#include "amg_matrix.h"
# include "../../amg_precond.h"
# include "../../../Allocation_Operations/allocations.h"

void precond_sol(double *q, double *p, int precond, AMG_precond_data *data) { // solve Mp = q

    precondAMG *amg_s = data->AMG_data;
    precondDPA *amg_d = data->DPA_data;
    int i, j;
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
