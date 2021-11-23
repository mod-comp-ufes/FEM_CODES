#include "SSTransportEquation3D.h"

void ebe_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*Me)[4])
{
	int i, j, J = 0;
	int mDim = NNOEL*NDOF;
	double **M;
	
	M = MatrixData->A;

	for (i = 0; i < mDim; i++){
		for (j = 0; j < mDim; j++, J++){
			M[E][J] = Me[i][j];
		}
	}

}


