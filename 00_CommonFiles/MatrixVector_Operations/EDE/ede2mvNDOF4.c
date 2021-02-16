#include "../matvec.h"
#include "../../Solvers_and_Preconditioners/preconditioners.h"

int ede2mvNDOF4(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *P, double *Q)
{
	int I, i1, i2, i3, i4, i5, i6, i7, i8;
	double p1,p2,p3,p4,p5,p6,p7,p8;
	double **A;
	int neq, nedge;
	int **lm = FemStructs->lm2;

	neq = Parameters->neq;
	nedge = Parameters->nedge;

	A = MatrixData->A;

	for (I=0; I<neq; I++)
		Q[I]=0;

	P[neq] = 0;
	Q[neq] = 0;

	for (I=0; I< nedge; I++){
		i1 = lm[I][0];
		i2 = lm[I][1];
		i3 = lm[I][2];
		i4 = lm[I][3];
		i5 = lm[I][4];
		i6 = lm[I][5];
		i7 = lm[I][6];
		i8 = lm[I][7];
		
		p1 = P[i1];
		p2 = P[i2];
		p3 = P[i3];
		p4 = P[i4];
		p5 = P[i5];
		p6 = P[i6];
		p7 = P[i7];
		p8 = P[i8];

		Q[i1] += A[I][4 ]*p5 + A[I][5 ]*p6 + A[I][6 ]*p7 + A[I][7 ]*p8;
		Q[i2] += A[I][12]*p5 + A[I][13]*p6 + A[I][14]*p7 + A[I][15]*p8;
		Q[i3] += A[I][20]*p5 + A[I][21]*p6 + A[I][22]*p7 + A[I][23]*p8;
		Q[i4] += A[I][28]*p5 + A[I][29]*p6 + A[I][30]*p7 + A[I][31]*p8;

		Q[i5] += A[I][32]*p1 + A[I][33]*p2 + A[I][34]*p3 + A[I][35]*p4;
		Q[i6] += A[I][40]*p1 + A[I][41]*p2 + A[I][42]*p3 + A[I][43]*p4;
		Q[i7] += A[I][48]*p1 + A[I][49]*p2 + A[I][50]*p3 + A[I][51]*p4;
		Q[i8] += A[I][56]*p1 + A[I][57]*p2 + A[I][58]*p3 + A[I][59]*p4;

		Q[neq] = 0;
	}

/*	if (strcasecmp(Parameters->Preconditioner,"NOT")==0)
		NOBlockDiag2_precond(Parameters, MatrixData, FemStructs, P, Q);
	else*/
		BlockDiag2_precond(Parameters, MatrixData, FemStructs, P, Q);

	return 0;
}
