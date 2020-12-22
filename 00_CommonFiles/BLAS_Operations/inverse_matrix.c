#include "ourBLAS.h"
#include "../Allocation_Operations/allocations.h"

//	calculate inverse of matrix
void inverse_matrix(double **Mat, double **Inv, int n){
	int l, h, m, k, i, j;
	double **M, **I, **B, **C, **Caux, det;
	
	M = (double**) mycalloc("M line of 'inverse_matrix'", n, sizeof(double*));
	I = (double**) mycalloc("I line of 'inverse_matrix'", n, sizeof(double*));
	B = (double**) mycalloc("B line of 'inverse_matrix'", n, sizeof(double*));
	C = (double**) mycalloc("C line of 'inverse_matrix'", n, sizeof(double*));
	Caux = (double**) mycalloc("Caux line of 'inverse_matrix'", n, sizeof(double*));
	for (i = 0; i < n; i++){
		M[i] = (double*) mycalloc("M colum of 'inverse_matrix'", n, sizeof(double));
		I[i] = (double*) mycalloc("I colum of 'inverse_matrix'", n, sizeof(double));
		B[i] = (double*) mycalloc("B colum of 'inverse_matrix'", n, sizeof(double));
		C[i] = (double*) mycalloc("C colum of 'inverse_matrix'", n, sizeof(double));
		Caux[i] = (double*) mycalloc("Caux colum of 'inverse_matrix'", n, sizeof(double));
	}
	
	for(i=0; i<n; i++){
		for(j=0;j<n;j++){
			M[i][j] = Mat[i][j];
		}
	}
	
	/*printf("\n Dentro de inverse_matrix \n");
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("M[%d][%d] = %lf\n", i, j, M[i][j]);
		}
	}
	getchar();*/
	
	// calculating the determinant
	det = determinant(M,n);
	
	/*printf("det M = %.15lf\n", det);
	getchar();*/
	
	if(det < 1e-40){
		printf("\n Inverse of Entered Matrix is not possible \n");
	} else if(n == 1){
		I[0][0] = 1.0/M[0][0];
	} else {
		// calculating the cofactor
		for (h = 0; h < n; h++){
			for (l = 0; l < n; l++){
				m = 0;
				k = 0;
				for (i = 0; i < n; i++){
					for (j = 0; j < n; j++){
						if (i != h && j != l){
							B[m][k] = M[i][j];
							if (k < (n-2)){
								k++;
							}else{
								k=0;
								m++;
							}
						}
					}
				}
				C[h][l] = pow(-1,(h+l))*determinant(B,(n-1));	// C = cofactor Matrix
			}
		}
		// calculating the inverse
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++){
				Caux[i][j] = C[j][i];
			}
		}
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++){
				I[i][j] = Caux[i][j]/det;	// I = inverse matrix
			}
		}
	}
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			Inv[i][j] = I[i][j];
		}
	}
	
	/*printf("\n Inversa Dentro de inverse_matrix \n");
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			printf("I[%d][%d] = %lf\n", i, j, I[i][j]);
		}
	}
	getchar();*/
	
	free(M);
	free(I);
	free(B);
	free(C);
	free(Caux);
	
}// end function
