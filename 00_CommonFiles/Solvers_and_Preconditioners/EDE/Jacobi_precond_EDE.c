#include "../preconditioners.h"
/* Define preconditioner action on matrix-vector product */
/* z = A p
 * A = D^{-1}
 * D * z = p
 * z out
*/
int Jacobi_precond_EDE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1; //auxiliar
	double z0, z1; //auxiliar
	double p0, p1; //auxiliar
	int nedge = Parameters->nedge;  //number of edges
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm2; //array LM
	double **inv = MatrixData->invDe; // 1.0 / De

	for (I=0; I<neq; I++)
		z[I] = 0;
	
	z[neq] = 0;

	for (I = 0; I < nedge; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];

		p0 = p[lm0];
		p1 = p[lm1];

		/* D * z = p */
		z0 = p0 * inv[I][0];
		z1 = p1 * inv[I][1];

		z[lm0] += z0;
		z[lm1] += z1;

		z[neq] = 0.0;
	}

	return 0;
}
