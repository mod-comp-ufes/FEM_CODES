#include "../matvec.h"
#include "../../Solvers_and_Preconditioners/preconditioners.h"

int ebemv(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *P, double *Q)
{

	int I, i1, i2, i3;
	double p1, p2, p3;
	double **A;
	int neq, nel;
	int **lm = FemStructs->lm;

	neq = Parameters->neq;
	nel = Parameters->nel;

	A = MatrixData->A;

	for (I = 0; I < neq; I++)
		Q[I] = 0;

	P[neq] = 0;
	Q[neq] = 0;

	for(I = 0; I < nel; I++){

		i1 = lm[I][0];
		i2 = lm[I][1];
		i3 = lm[I][2];

		p1 = P[i1];
		p2 = P[i2];
		p3 = P[i3];	
				
		Q[i1] += A[I][0]*p1 + A[I][1]*p2 + A[I][2]*p3;
		Q[i2] += A[I][3]*p1 + A[I][4]*p2 + A[I][5]*p3;
		Q[i3] += A[I][6]*p1 + A[I][7]*p2 + A[I][8]*p3;
	
		Q[neq] = 0.0;

	}

	return 0;
}







