#include "ourBLAS.h"

// calculate minor of matrix OR build new matrix : k-had = minor
void minor_matrix(double **b, double **a, int i, int n) {
	int j, l, h = 0, k = 0;
	
	/*printf("\nDentro de minor_matrix\n");
	for(l = 0; l < n; l++){
		for(j = 0; j < n; j++){
			printf("am[%d][%d] = %lf\n", l, j, a[l][j]);
		}
	}*/
	
	for(l = 1; l < n; l++){
		for(j = 0; j < n; j++){
			if(j == i)
				continue;
			b[h][k] = a[l][j];
			k++;
			if(k == (n-1)){
				h++;
				k=0;
			}
		}
	}
	
}// end function
