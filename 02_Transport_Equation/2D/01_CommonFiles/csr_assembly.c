#include "TranspEquation.h"

void csr_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*ke)[3])
{
	int nnzero = Parameters->nnzero;
	double *K;
	int *CSR_by_Element;
	 
		
	K = MatrixData->AA;	
	CSR_by_Element = MatrixData->Scheme_by_Element[E];
	
	K[CSR_by_Element[0]] += ke[0][0];
	K[CSR_by_Element[1]] += ke[0][1];
	K[CSR_by_Element[2]] += ke[0][2];
	K[CSR_by_Element[3]] += ke[1][0];
	K[CSR_by_Element[4]] += ke[1][1];
	K[CSR_by_Element[5]] += ke[1][2];
	K[CSR_by_Element[6]] += ke[2][0];
	K[CSR_by_Element[7]] += ke[2][1];
	K[CSR_by_Element[8]] += ke[2][2];
 
	K[nnzero] = 0;

}




