#include "../preconditioners.h"
/* Define preconditioner action on matrix-vector product */
/* p = A z
 * A = (LU)^{-1}
 * L * p = w
 * U * w = z
 * z out
*/
int LU_precond_EDE (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z){
	int I;
	int lm0, lm1; //auxiliar
	double z0, z1; //auxiliar
	int nedge = Parameters->nedge;  //number of elements
	int neq = Parameters->neq; //number of equations
	int **lm = FemStructs->lm2; //array LM
	double **LUe = MatrixData->LUe;

	z[neq] = 0;

	for (I = 0; I < nedge; I++){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];

		z0 = z[lm0];
		z1 = z[lm1];

		// L * p = w
//		z1 = (z1 - LUe[I][2] * z0)*LUe[I][3]; // LUe[I][3] is store as 1/LUe[I][4]
		z1 += - LUe[I][2] * z0; // LUe[I][3] is store as 1/LUe[I][4]

		// U * w = z
//		z1 *= LUe[I][3];
//		z0 += - LUe[I][1] * z1;
	
//		z[lm0] = z0;
		z[lm1] = z1;

		z[neq] = 0.0;
	}

	for (I = nedge-1; I >=0 ; I--){
		lm0 = lm[I][0];
	 	lm1 = lm[I][1];

		z0 = z[lm0];
		z1 = z[lm1];

		// L * p = w
//		z1 = (z1 - LUe[I][2] * z0)*LUe[I][3]; // LUe[I][3] is store as 1/LUe[I][4]
//		z1 += - LUe[I][2] * z0; // LUe[I][3] is store as 1/LUe[I][4]

		// U * w = z
		z1 *= LUe[I][3];
		z0 += - LUe[I][1] * z1;
	
		z[lm0] = z0;
		z[lm1] = z1;

		z[neq] = 0.0;
	}	




	return 0;
}


