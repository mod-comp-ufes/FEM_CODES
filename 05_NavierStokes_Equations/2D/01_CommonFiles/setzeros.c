#include "NavierStokesEquations.h"

int setzeros(ParametersType *Parameters, MatrixDataType *MatrixData)
{
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3)==0){
		int size = NDOF*NNOEL;
		int size2 = size*size;

		memset(MatrixData->Aaux,0,Parameters->nel*size2*sizeof(double));
	}
	else if (strncmp(Parameters->MatrixVectorProductScheme,"EDE",3)==0){

		memset(MatrixData->Aaux,0,((Parameters->nedge+1)*NDOF*NDOF*(NNOEL+1))*sizeof(double));
	}	
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){

		memset(MatrixData->AA,0,(Parameters->nnzero)*sizeof(double));
		
	}
	return 0;
}

