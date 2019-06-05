	#include "../preconditioners.h"
/* Define preconditioner action on matrix-vector product */
/* p = A z
 * A = (LU)^{-1}
 * L * p = w
 * U * w = z
 * z out
*/
int LU_precond_EBE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z){
	int I;
	int lm0, lm1, lm2; //auxiliar
	double z0, z1, z2; //auxiliar
	int nel = Parameters->nel;  //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm; //array LM
	double **LUe = MatrixData->LUe;

	z[neq] = 0;

	for (I = 0; I < nel; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];
		lm2 = lm[I][2];

		z0 = z[lm0];
		z1 = z[lm1];
		z2 = z[lm2];

	/*        printf("LM%d->",I);
                int j;
                for (j=0;j<3;j++)
                       printf("%d ",lm[I][j]);
                printf("\n ");
                printf("z0=%lf z1=%lf z2=%lf\n",z0,z1,z2);
*/


		// L * p = w
//		z1 = (z1 - LUe[I][3] * z0)*LUe[I][4]; // LUe[I][4] is store as 1/LUe[I][4]
//		z2 = (z2 - LUe[I][6] * z0 - LUe[I][7] * z1)*LUe[I][8]; // LUe[I][8] is stored as 1/LUe[8]		
		z1 += - LUe[I][3] * z0; // LUe[I][4] is store as 1/LUe[I][4]
		z2 += - LUe[I][6] * z0 - LUe[I][7] * z1; // LUe[I][8] is stored as 1/LUe[8]		

		// U * w = z
		z2 *= LUe[I][8];
		z1 = (z1 - LUe[I][5] * z2) * LUe[I][4];
		z0 += - LUe[I][1] * z1 - LUe[I][2] * z2;

		z[lm0] = z0;
		z[lm1] = z1;
		z[lm2] = z2;

		z[neq] = 0.0;

/*		printf("e=%d\n",I);
		int j;
		for (j=0;j<9;j++)
			printf("LUe[%d]=%lf Ae[%d]=%lf\n",j, LUe[I][j], j, MatrixData->A[I][j]);*/
	}
	
	return 0;
}


