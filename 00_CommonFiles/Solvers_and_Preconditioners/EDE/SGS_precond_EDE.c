#include "../preconditioners.h"
/* z = A p
 * A = (D + U)^{-1} (D + L)^{-1}
 * (L + D) * p = w
 * (D + U) * z = p
 * z out
*/
int SGS_precond_EDE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1; //auxiliar
	double z0, z1; //auxiliar
	int nedge = Parameters->nedge; //number of edges
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm2; //array LM
	double **A = MatrixData->A;

	z[neq] = 0;

	for (I = 0; I < nedge; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];

		z0 = z[lm0];
		z1 = z[lm1];

		z1 += - A[I][2]*z0;

		z[lm1] = z1;

		z[neq] = 0.0;
	}

	return 0;
}

int SGS_precondR_EDE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1; //auxiliar
	double z0, z1; //auxiliar
	int nedge = Parameters->nedge; //number of edges
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm2; //array LM
	double **A = MatrixData->A;

	z[neq] = 0;

	for (I = nedge-1; I >=0; I--){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];

		z0 = z[lm0];
		z1 = z[lm1];

		z0 += - A[I][1]*z1;

		z[lm0] = z0;

		z[neq] = 0.0;
	}

	
	return 0;

}

