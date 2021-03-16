#include "../preconditioners.h"
/* p = A z
 * A = (D + w U)^{-1} (D + w L)^{-1}
 * (D + w L) * p = w
 * (D + w U) * w = z
 * z out
*/
int SSOR_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1, lm2; //auxiliar
	double z0, z1, z2; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **A = MatrixData->A; //matrix data

	double omega;
	omega = atof(&(Parameters->Preconditioner[3])); // In SORw, w is the weight coefficient

	z[neq] = 0;

	// OPTION 1: (D + wU) (D +wL)

	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];

		z0 = z[lm0];
		z1 = z[lm1];
		z2 = z[lm2];

		z1 = (z1 - omega*A[I][3]*z0);
		z2 = (z2 - omega*A[I][6]*z0 - omega*A[I][7]*z1);

		z[lm1] = z1;
		z[lm2] = z2;

		z[neq] = 0.0;

	}

	return 0;
}

int SSOR_precondR_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1, lm2; //auxiliar
	double z0, z1, z2; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **A = MatrixData->A; //matrix data

	double omega;
	omega = atof(&(Parameters->Preconditioner[3])); // In SORw, w is the weight coefficient

	z[neq] = 0;

	// OPTION 1: (D + wU) (D +wL)

	for (I = nel-1; I >= 0; I--){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];

		z0 = z[lm0];
		z1 = z[lm1];
		z2 = z[lm2];

		z1 = z1 - omega*A[I][5]*z2;
		z0 = z0 - omega*A[I][1]*z1 - omega*A[I][2]*z2;

		z[lm0] = z0;
		z[lm1] = z1;

		z[neq] = 0.0;

	}

	return 0;
}


