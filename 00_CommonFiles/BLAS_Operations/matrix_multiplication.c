#include "ourBLAS.h"

void matrix_multiplication(double **Mat1, int lin1, int col1, double **Mat2, int lin2, int col2, double **Sol) {
	//Definição de variaveis
	int i, j, x;
	double aux = 0.0;

	if(col1 == lin2) {
		
		//printf("To no if\n");
		
		//Processamento e saida em tela  =  PRODUTO DAS MATRIZES
		for(i = 0; i < lin1; i++) {
			
			//printf("To no for 1i: %d\n",i);
			
			for(j = 0; j < col2; j++) {
			
				//printf("to no for j: %d\n",j);
			
				Sol[i][j] = 0;
				for(x = 0; x < lin2; x++) {
					
					//printf("to no for x: %d\n",x);
					
					aux +=  Mat1[i][x] * Mat2[x][j];
				}

				Sol[i][j] = aux;
				aux = 0;
			}
		}
	}else {
		printf("\n\n Nao ha com multiplicar as matrizes dadas ");
	}

}
