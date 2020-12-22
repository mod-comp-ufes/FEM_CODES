	#include "preconditioners.h"
/* Define preconditioner action on matrix-vector product */
/* bpre = M^{-1}*b
 * M*bpre = b
 * M = LU
 * L*z = b
 * U*brep = z
 * bpre out
*/
int LU_precond_EBE_NNOEL4 (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *b, double *bpre){
	int I;
	int lm0, lm1, lm2, lm3; //auxiliar
	double b0, b1, b2, b3, bp0, bp1, bp2, bp3, z0, z1, z2, z3; //auxiliar
	int nel = Parameters->nel;  //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **LUe = MatrixData->LUe;
	
/*	int i;
	printf("\nVetor b dentro do precondicionador\n");
	for(i = 0; i <= neq; i++){
		printf("%lf\n",b[i]);
	}
	getchar(); */

	b[neq] = 0;

	for (I = 0; I < nel; I++){

		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];
		lm3 = lm[I][3];

		// dados do vetor que correspondem ao elemento
		b0 = b[lm0];
		b1 = b[lm1];
		b2 = b[lm2];
		b3 = b[lm3];

		// L * z = b (como o sistema é triangular, resolve por substituição)
		z0 = b0;
		z1 = b1 - LUe[I][4]*z0;
		z2 = b2 - LUe[I][8]*z0 - LUe[I][9]*z1;
		z3 = b3 - LUe[I][12]*z0 - LUe[I][13]*z1 - LUe[I][14]*z2;		

		// U * bpred = z (como o sistema é triangula, resolve por substituição)
		bp3 = z3/LUe[I][15];
		bp2 = (z2 - LUe[I][11]*bp3)/LUe[I][10];
		bp1 = (z1 - LUe[I][6]*bp2 - LUe[I][7]*bp3)/LUe[I][5];
		bp0 = (z0 - LUe[I][1]*bp1 - LUe[I][2]*bp2 - LUe[I][3]*bp3)/LUe[I][0];

		bpre[lm0] += bp0;
		bpre[lm1] += bp1;
		bpre[lm2] += bp2;
		bpre[lm3] += bp3;

		bpre[neq] = 0.0;
		
	}
	
/*	printf("\nVetor bpre dentro do precondicionador\n");
	for(i = 0; i <= neq; i++){
		printf("%lf\n",bpre[i]);
	}
	getchar();*/
	
	return 0;
}
