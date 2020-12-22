#include "EulerEquations.h"

void ebe_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*Me)[12])
{
	int i, j, J = 0;
	double **M;

	M = MatrixData->A;

	for (i = 0; i < 12; i++){
		for (j = 0; j < 12; j++, J++){
			M[E][J] = Me[i][j];
		}
	}

}


