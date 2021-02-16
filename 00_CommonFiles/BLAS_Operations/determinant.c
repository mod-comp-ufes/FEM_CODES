#include "ourBLAS.h"
#include "../Allocation_Operations/allocations.h"

// calculate determinte of matrix
double determinant(double **a, int n) {
	int i;
	double **b, sum = 0.0;
	
	b = (double**) mycalloc("b line of 'determinant'", n, sizeof(double*));
	for (i = 0; i < n; i++)
		b[i] = (double*) mycalloc("b colum of 'determinant'", n, sizeof(double));
	
	/*printf("\n Matriz dentro de determinant \n");
	int j;
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("a[%d][%d] = %lf\n", i, j, a[i][j]);
		}
	}
	getchar();*/
	
	if (n == 1){
		return a[0][0];
	} else if(n == 2){
		return (a[0][0]*a[1][1]-a[0][1]*a[1][0]);
	} else {
		for(i = 0; i < n; i++){
			minor_matrix(b,a,i,n);
			sum = (double)(sum + a[0][i]*pow(-1,i)*determinant(b,(n-1))); // sum = determinte matrix
		}
	}
	
	free(b);
	
	return sum;
	
}// end function
