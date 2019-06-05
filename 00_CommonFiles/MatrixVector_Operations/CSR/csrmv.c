#include "../matvec.h"

int csrmv(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *P, double *Q)
{
	double *AA;
	int i, j, *IA, *JA;
	int neq;

	neq = Parameters->neq;

	AA = MatrixData->AA;
	IA = MatrixData->IA;
	JA = MatrixData->JA;

	for (i = 0; i < neq; i++)
		Q[i] = 0;

	for (i = 0; i < neq; i++)
		for (j = IA[i]; j < IA[i+1]; j++)
			Q[i] += AA[j]*P[JA[j]];

	return 0;

}






