# ifndef amg_h
# define amg_h
# include "amg_matrix.h"

typedef struct {
    int aggress, refinePass, truncq;
    double trunc_fact;
    matrix **I_cf, **I_cf_t;
} precondAMG;

void AMG_setup (matrix *A_o, double *f_o, double *u_o, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, matrix ***A, double ***f, double ***u, double **r, matrix ***I_cf, matrix ***I_cf_t);
void AMG_destroy (int NCL, matrix **A, double **f, double **u, double *r, matrix **I_cf, matrix **I_cf_t);
void AMG_Vcycle_precond (int NCL, matrix **A, double **f, double **u, double *r, matrix **I_cf, matrix **I_cf_t, double omega, int nr);
int AMG (matrix *A, double *f, double *u, int NCL, double str_thr, int aggress, int refinePass, int truncq, double trunc_fact, double omega, int nr, double tol, int lmax);

# endif
