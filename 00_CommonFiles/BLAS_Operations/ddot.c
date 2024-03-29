#include "ourBLAS.h"

// r = x dot y -> produto escalar
double ddot(int n, double *x, double *y)
{
	long int i, m;
	double stemp;

	stemp = 0.0;
	m = n-4;

 	for (i = 0; i < m; i += 5)
 		stemp += x[i] * y[i] + x[i+1] * y[i+1] + x[i+2] * y[i+2] + x[i+3] * y[i+3] + x[i+4] * y[i+4];

	for ( ; i < n; i++)        /* clean-up loop */
        	stemp += x[i] * y[i];

	return stemp;
} 
