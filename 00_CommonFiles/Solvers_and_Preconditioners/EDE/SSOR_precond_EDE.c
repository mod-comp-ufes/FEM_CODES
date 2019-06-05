#include "../preconditioners.h"
/* p = A z
 * A = (D + w U)^{-1} (D + w L)^{-1}
 * (D + w L) * p = w
 * (D + w U) * w = z
 * z out
*/
int SSOR_precond_EDE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1; //auxiliar
	double z0, z1; //auxiliar
	int nedge = Parameters->nedge; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm2; //array LM
	double **A = MatrixData->A; //matrix data

	double omega;
	omega = atof(&(Parameters->Preconditioner[3])); // In SORw, w is the weight coefficient

	z[neq] = 0;

	// OPTION 1: (D + wU) (D +wL)

	for (I = 0; I < nedge; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];

		z0 = z[lm0];
		z1 = z[lm1];

		z1 += - omega*A[I][2]*z0;

		z[lm1] = z1;

		z[neq] = 0.0;

	}

	return 0;
}

int SSOR_precondR_EDE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1; //auxiliar
	double z0, z1; //auxiliar
	int nedge = Parameters->nedge; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm2; //array LM
	double **A = MatrixData->A; //matrix data

	double omega;
	omega = atof(&(Parameters->Preconditioner[3])); // In SORw, w is the weight coefficient

	z[neq] = 0;

	// OPTION 1: (D + wU) (D +wL)

	for (I = nedge-1; I >= 0; I--){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];

		z0 = z[lm0];
		z1 = z[lm1];

		z0 += - omega*A[I][1]*z1;

		z[lm0] = z0;

		z[neq] = 0.0;

	}

	return 0;
}


