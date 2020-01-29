# ifndef matrix_h
# define matrix_h

typedef struct { // compressed row storage
	int m, n, nnz, *col_ind, *row_ptr;
	double *val, *diag;
} matrix;

matrix * create_m (int m, int n, int nnz);
matrix * copy_m(matrix * A);
void destroy_m (matrix *A);
void print_m (matrix *A);
double get_Aij (matrix *A, int i, int j);
matrix * csr_transpose (matrix *A);
matrix * sum(matrix *A, matrix *A_t);
matrix * sum_abs(matrix *A, matrix *A_t);
matrix * sum_upper(matrix *A, matrix *A_t);
double maxCoeff_abs (matrix *A);
matrix * matmat_sparseproduct (matrix *A, matrix *B);
void mat_vec (double *p, matrix *A, double *v);
void mat_vec_plus (double *p, matrix *A, double *v);

# endif
