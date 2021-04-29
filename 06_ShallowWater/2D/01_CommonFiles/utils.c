#include "ShalowWater.h"


void printU(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, eq1, eq2, eq3, nnodes = Parameters->nnodes;
	double X, Y, *U = FemStructs->u;

    for(I=0; I<nnodes; I++){
		eq1 = FemStructs->Node[I].id[0];
		eq2 = FemStructs->Node[I].id[1];
		eq3 = FemStructs->Node[I].id[2];
		X = FemStructs->Node[I].x;
		Y = FemStructs->Node[I].y;
		printf("(%.2f\t%.2f)\t", FemStructs->Node[I].x, FemStructs->Node[I].y);

		if(eq1>=0)
			printf("%.3lf\t", U[eq1]);
		else
			printf("%.3lf\t", FemFunctions->hpresc(X, Y));

		if(eq2>=0)
			printf("%.3lf\t", U[eq2]);
		else
			printf("%.3lf\t", FemFunctions->qxpresc(X, Y));

	   	if(eq3>=0)
	   		printf("%.3lf\t", U[eq3]);
	   	else
			printf("%.3lf\t", FemFunctions->qypresc(X, Y));

		printf("%d\t%d\t%d\n", eq1, eq2, eq3);
	}
	getchar();
}

void checknull(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, eq1, eq2, eq3, nnodes = Parameters->nnodes;
	double X, Y, *U = FemStructs->u;
	
    for(I=0; I<nnodes; I++){
		eq1 = FemStructs->Node[I].id[0];
		eq2 = FemStructs->Node[I].id[1];
		eq3 = FemStructs->Node[I].id[2];
		X = FemStructs->Node[I].x;
		Y = FemStructs->Node[I].y;

		if(U[eq1] > 10.0) {
			printf("(%.2f\t%.2f)\n", FemStructs->Node[I].x, FemStructs->Node[I].y);
			getchar();
		}
		if(isnan(U[eq1]) || isinf(U[eq1]) || 
		   isnan(U[eq2]) || isinf(U[eq2]) || 
		   isnan(U[eq3]) || isinf(U[eq3])) {
			printf("(%.2f\t%.2f)\n", FemStructs->Node[I].x, FemStructs->Node[I].y);

			if(eq1>=0)
				printf("%.3lf\t", U[eq1]);

			if(eq2>=0)
				printf("%.3lf\t", U[eq2]);

	   		if(eq3>=0)
	   			printf("%.3lf\t", U[eq3]);
			getchar();
			}
	}
}

int AssemblyGlobalMatrix(MatrixDataType *MatrixData, FemStructsType *FemStructs, int e, double (*Me)[3], int neq)
{	
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