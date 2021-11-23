#include "SSTransportEquation3D.h"

int AssemblyGlobalMatrix(MatrixDataType *MatrixData, FemStructsType *FemStructs, int e, double (*Me)[4], int neq){
	
	int size = NNOEL*NDOF;
	int i, j, J, eq[size], ind[size]; // eq' armazena o número da equação e i' a posição onde estava essa equação
	int **lm;
	double **N;
	
	N = MatrixData->G;
	lm = FemStructs->lm;
		
	J = 0;
	for(i = 0; i < size; i++){
		if(lm[e][i] != neq){
			eq[J] = lm[e][i];
			ind[J] = i;
			J++;
		}
	}

/*	printf("\n Elem = %d\n", e);
	printf("J = %d\n", J);
	for(i = 0; i < J; i++){
		printf("eq = %d \t ind = %d \n",eq[i],ind[i]);
	}
	getchar();*/
	
	for(i = 0; i < J; i++){
		for(j = 0; j < J; j++){
			N[eq[i]][eq[j]] = N[eq[i]][eq[j]] + Me[ind[i]][ind[j]];
		}
	}

	return 0;
	
	
}
