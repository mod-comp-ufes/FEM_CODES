# ifndef amg_dpa_h
# define amg_dpa_h

#include "amg_matrix.h"

typedef struct {
	int **match;
} precondDPA;

void dpa_setup(matrix *A_o, double *f_o, double *u_o, int NCL, matrix ***A, double ***f, double ***u, double **r, int ***match, double beta);
int    dpa_AMG (matrix *A_o, double *f_o, double *u_o, int NCL, double omega, int nr, double tol, int lmax, double beta);
void   dpa_Vcycle_precond(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match);
void   xdpa_precond(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match);
void   dpa_destroy(int NCL, matrix **A, double **f, double **u, double *r, int **match);

# endif
