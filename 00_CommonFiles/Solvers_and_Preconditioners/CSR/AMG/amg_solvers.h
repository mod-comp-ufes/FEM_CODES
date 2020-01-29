# ifndef amg_sor_relax_h
# define amg_sor_relax_h
# include "amg_matrix.h"

void SOR_relax (matrix *A, double *f, double *u, double omega, int iter);
void SOR_relax_rev (matrix *A, double *f, double *u, double omega, int iter);

# endif
