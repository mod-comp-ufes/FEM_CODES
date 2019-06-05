#include "../matvec.h"
#include "../../Solvers_and_Preconditioners/preconditioners.h"

int edemv(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *P, double *Q)
{
	int I, i1, i2;
	double p1, p2;
	double **A;
	int neq, nedge;
	int **lm = FemStructs->lm2;

	neq = Parameters->neq;
	nedge = Parameters->nedge;

	A = MatrixData->A;

	for (I = 0; I < neq; I++)
		Q[I] = 0;

	P[neq] = 0;
	Q[neq] = 0;

	for(I = 0; I < nedge; I++){

		i1 = lm[I][0];
		i2 = lm[I][1];

		p1 = P[i1];
		p2 = P[i2];
				
		Q[i1] += A[I][0]*p1 + A[I][1]*p2;
		Q[i2] += A[I][2]*p1 + A[I][3]*p2;
	
		Q[neq] = 0.0;

	}

	return 0;
}







