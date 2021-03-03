#include "ShalowWater.h"


int setzeros(ParametersType *Parameters, MatrixDataType *MatrixData)
{
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") == 0)
		memset(MatrixData->AA,0,(Parameters->nnzero)*sizeof(double));

	return 0;
}
