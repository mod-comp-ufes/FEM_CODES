#include "EulerEquations.h"

void csr_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*Me)[12])
{
	int I, J, K;
	double *M;
	int *CSR_by_Element;
	int nnzero = Parameters->nnzero;
			
	M = MatrixData->AA;
	CSR_by_Element = MatrixData->Scheme_by_Element[E];

	for (I=0, K=0; I<12; I++)
		for (J=0; J<12; J++)
			M[CSR_by_Element[K++]] += Me[I][J];

	M[nnzero] = 0;

}





