#include "../preconditioners.h"
/* p = A z
 * A = D^{-1}
 * D * p = z
 * z out
*/
int JacobiDOF4_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{
	int I;
	int lm0, lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11; //auxiliar
	double z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11; //auxiliar
	double p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11; //auxiliar
	int nel = Parameters->nel; //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM


	for (I = 0; I < neq; I++){
		z[I] = 0;
	}
	z[neq] = 0;

	//OPTION 1: 12X12
	double **inv = MatrixData->invDe; // 1.0 / De ou 1.0 / (1.0 + De)
	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
		lm1 = lm[I][1];
		lm2 = lm[I][2];
		lm3 = lm[I][3];
		lm4 = lm[I][4];
		lm5 = lm[I][5];
		lm6 = lm[I][6];
		lm7 = lm[I][7];
		lm8 = lm[I][8];
		lm9 = lm[I][9];
		lm10 = lm[I][10];
		lm11 = lm[I][11];

		p0 = p[lm0];
		p1 = p[lm1];
		p2 = p[lm2];
		p3 = p[lm3];
		p4 = p[lm4];
		p5 = p[lm5];
		p6 = p[lm6];
		p7 = p[lm7];
		p8 = p[lm8];
		p9 = p[lm9];
		p10 = p[lm10];
		p11 = p[lm11];

		// D * p = z
		z0 = p0 * inv[I][0];
		z1 = p1 * inv[I][1];
		z2 = p2 * inv[I][2];
		z3 = p3 * inv[I][3];
		z4 = p4 * inv[I][4];
		z5 = p5 * inv[I][5];
		z6 = p6 * inv[I][6];
		z7 = p7 * inv[I][7];
		z8 = p8 * inv[I][8];
		z9 = p9 * inv[I][9];
		z10 = p10 * inv[I][10];
		z11 = p11 * inv[I][11];

		z[lm0] += z0;
		z[lm1] += z1;
		z[lm2] += z2;
		z[lm3] += z3;
		z[lm4] += z4;
		z[lm5] += z5;
		z[lm6] += z6;
		z[lm7] += z7;
		z[lm8] += z8;
		z[lm9] += z9;
		z[lm10] += z10;
		z[lm11] += z11;
	}

	z[neq] = 0.0;

  return 0;
}
