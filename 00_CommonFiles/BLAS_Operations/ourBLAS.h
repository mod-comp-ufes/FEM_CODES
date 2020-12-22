#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <conio.h>

int daxpy(int, double, double *, double *);
int dcopy(int, double *, double *);
double ddot(int, double *, double *);
int dscal(int, double, double *);
int dzero(int, double *);
int izero(int, int *);
int dtrsvUP (int, double **, double *, double *);

void minor_matrix(double **, double **, int , int );
double determinant(double **, int );
void inverse_matrix(double **, double **, int );

double dmax(int, double *);

void matrix_multiplication(double **Mat1, int lin1, int col1, double **Mat2, int lin2, int col2, double **Sol);
