# ifndef util_h
# define util_h

# include "amg_matrix.h"
# include "../../amg_precond.h"

void precond_sol (double *q, double *p, int precond, AMG_precond_data *data);
void matvec_product (double *p, matrix *A, double *v, int precond, AMG_precond_data *data);
void residual (double *r, double *f, matrix *A, double *u, int precond, AMG_precond_data *data);
double norm_inf (double *v, int n);
double norm_euclid(double *v, int n);

# endif
