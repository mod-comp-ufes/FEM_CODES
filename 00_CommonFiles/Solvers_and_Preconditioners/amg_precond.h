# ifndef amg_precond_h
# define amg_precond_h

# include "CSR/AMG/amg.h"
# include "CSR/AMG/amg_dpa.h"

typedef struct {
    int     precond;
    int     NCL;
    matrix  **A;
    double  **f;
    double  **u;
    double  *r;
    double  str_thr;
    double  omega;
    int     nr;
    precondAMG  *AMG_data;
    precondDPA  *DPA_data;
} AMG_precond_data;

/* amg_precond.c prototypes are defined on preconditioners.h */

void AMG_precond_data_destroy(AMG_precond_data *data);

# endif
