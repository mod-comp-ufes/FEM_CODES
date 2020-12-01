#include "amg_util.h"
#include "../../solvers.h"

int AMG_GMRES(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
			FemFunctionsType *FemFunctions,	double *B, double *X) {
    matrix *A 		= MatrixData->amg_precond_data->A[0];
    double *f 		= B;
    double *u 		= X;
    int    k 		= Parameters->KrylovBasisVectorsQuantity;
    double tol 		= Parameters->SolverTolerance;
    int    lmax		= Parameters->LinearMaxIter;
    int    precond 	= MatrixData->amg_precond_data->precond;

    #ifdef SSTranspEquation2D
	Parameters->gmres++;
    #endif

    int ind, i, j, l, n = A->m;
    double **U, **H, *e, *y, *c, *s, *f2, tol2, delta, aux, r;
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
        i = 0;
        residual(U[i], f, A, u, precond, MatrixData->amg_precond_data);
        e[i] = norm_euclid(U[i], n);
        for (aux = 1.0 / e[i], ind = 0; ind < n; ind++) U[i][ind] *= aux;
        delta = e[i];

        for (i=0; i < k && delta > tol2; i++, l++) {

            matvec_product(U[i + 1], A, U[i], precond, MatrixData->amg_precond_data);
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
    AMG_precond_data_destroy(MatrixData->amg_precond_data);
    return (delta <= tol2) ? 1 : 0;
};
