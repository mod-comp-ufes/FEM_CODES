#include "../preconditioners.h"

int SGS_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1, lm2; //auxiliar
	double z0, z1, z2; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **A = MatrixData->A; //matrix data

	z[neq] = 0;

	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];

		z0 = z[lm0];
		z1 = z[lm1];
		z2 = z[lm2];

		z1 += - A[I][3]*z0;
		z2 += - A[I][6]*z0 - A[I][7]*z1;

		z[lm1] = z1;
		z[lm2] = z2;

		z[neq] = 0.0;

	}

	return 0;
}

int SGS_precondR_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1, lm2; //auxiliar
	double z0, z1, z2; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **A = MatrixData->A; //matrix data

	z[neq] = 0;

	for (I = nel-1; I >= 0; I--){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];

		z0 = z[lm0];
		z1 = z[lm1];
		z2 = z[lm2];

		z1 += - A[I][5]*z2;
		z0 += - A[I][1]*z1 - A[I][2]*z2;

		z[lm0] = z0;
		z[lm1] = z1;

		z[neq] = 0.0;

	}

	return 0;
}


