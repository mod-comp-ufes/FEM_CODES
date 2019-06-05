#include "SSNavierStokesEquations.h"

int setzeros(ParametersType *Parameters, MatrixDataType *MatrixData)
{
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3)==0){
		int size = NDOF*NNOEL;
		int size2 = size*size;

		memset(MatrixData->Aaux,0,Parameters->nel*size2*sizeof(double));
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){

		memset(MatrixData->AA,0,(Parameters->nnzero)*sizeof(double));
		
	}
	return 0;
}

