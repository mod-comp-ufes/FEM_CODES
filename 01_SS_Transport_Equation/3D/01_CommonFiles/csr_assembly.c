#include "SSTransportEquation3D.h"

void csr_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*ke)[4])
{
	double *K;
	int *CSR_by_Element;
	int nnzero = Parameters->nnzero;
	 
	K = MatrixData->AA;	
	CSR_by_Element = MatrixData->Scheme_by_Element[E];
	
	K[CSR_by_Element[0]] += ke[0][0];
	K[CSR_by_Element[1]] += ke[0][1];
	K[CSR_by_Element[2]] += ke[0][2];
	K[CSR_by_Element[3]] += ke[0][3];
	K[CSR_by_Element[4]] += ke[1][0];
	K[CSR_by_Element[5]] += ke[1][1];
	K[CSR_by_Element[6]] += ke[1][2];
	K[CSR_by_Element[7]] += ke[1][3];
	K[CSR_by_Element[8]] += ke[2][0];
	K[CSR_by_Element[9]] += ke[2][1];
	K[CSR_by_Element[10]] += ke[2][2];
	K[CSR_by_Element[11]] += ke[2][3];
	K[CSR_by_Element[12]] += ke[3][0];
	K[CSR_by_Element[13]] += ke[3][1];
	K[CSR_by_Element[14]] += ke[3][2];
	K[CSR_by_Element[15]] += ke[3][3];
	
	K[nnzero] = 0;

/*	printf("neq: %d \t nnzero: %d\n\n", Parameters->neq, nnzero);

	if (E < 10){
		printf("--- CSR_by_Element \t %d ---\n", E);
		for(int i = 0; i< 16; i++)
			printf("%d\t",CSR_by_Element[i]);
		printf("\n");

		printf("\n--- Matriz LM %d ---\n", E);
		printf("%d \t %d \t %d \t %d \n", FemStructs->lm[E][0], FemStructs->lm[E][1], FemStructs->lm[E][2], FemStructs->lm[E][3]);

		getchar();
	} */
}




