#include "../matvec.h"

int ebemv3Dndof1(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *P, double *Q)
{

	int I, i1, i2, i3, i4;
	double p1, p2, p3, p4;
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
		i4 = lm[I][3];

		p1 = P[i1];
		p2 = P[i2];
		p3 = P[i3];	
		p4 = P[i4];
		
		Q[i1] += A[I][0]*p1 + A[I][1]*p2 + A[I][2]*p3 + A[I][3]*p4;
		Q[i2] += A[I][4]*p1 + A[I][5]*p2 + A[I][6]*p3 + A[I][7]*p4;
		Q[i3] += A[I][8]*p1 + A[I][9]*p2 + A[I][10]*p3 + A[I][11]*p4;
		Q[i4] += A[I][12]*p1 + A[I][13]*p2 + A[I][14]*p3 + A[I][15]*p4;
	
		Q[neq] = 0.0;

	}

	return 0;
}
