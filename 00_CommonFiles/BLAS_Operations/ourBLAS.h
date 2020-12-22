#ifndef _ourBLAS_h_
#define _ourBLAS_h_

#include <stdio.h>
#include <math.h>
#include <string.h>


// y = ax + y
int daxpy(int n, double a, double *x, double *y);
// y = x
int dcopy(int n, double *x, double *y);
// r = x dot y, scalar product
double ddot(int n, double *x, double *y);
// x = ax
int dscal(int n, double a, double *x);
// x = 0.0
int dzero(int n, double *x);
// v = 0
int izero(int n, int *v);
// Ax = b, A upper triangular
int dtrsvUP (int n, double **A, double *b, double *x);

// calculate minor of matrix OR build new matrix : k-had = minor
void minor_matrix(double **b, double **a, int i, int n);
// calculate determinte of matrix
double determinant(double **a, int n);
// calculate inverse of matrix
void inverse_matrix(double **Mat, double **Inv, int n);

// max(v)
double dmax(int n, double *v);

// matriz multiplication
void matrix_multiplication(double **Mat1, int lin1, int col1, double **Mat2, int lin2, int col2, double **Sol);

#endif
