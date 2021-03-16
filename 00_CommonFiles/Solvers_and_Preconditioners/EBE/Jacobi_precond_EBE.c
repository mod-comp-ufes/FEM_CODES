#include "../preconditioners.h"
/* Define preconditioner action on matrix-vector product */
/* p = A z
 * A = D^{-1}
 * D * p = z
 * z out
*/
int Jacobi_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1, lm2; //auxiliar
	double z0, z1, z2; //auxiliar
	double p0, p1, p2; //auxiliar
	int nel = Parameters->nel;  //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **inv = MatrixData->invDe; //  1.0 / De ou 1.0 / (1.0 + De)

	for (I = 0; I < neq; I++)
		z[I] = 0;

	z[neq] = 0;

	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];

		p0 = p[lm0];
		p1 = p[lm1];
		p2 = p[lm2];

		/* D * p = z */
		z0 = p0 * inv[I][0];
		z1 = p1 * inv[I][1];
		z2 = p2 * inv[I][2];

		z[lm0] += z0;
		z[lm1] += z1;
		z[lm2] += z2;

		z[neq] = 0.0;

	}
	
	return 0;
}
