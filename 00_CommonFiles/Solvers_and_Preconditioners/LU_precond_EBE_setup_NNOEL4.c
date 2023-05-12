#include "preconditioners.h"
#include "../Allocation_Operations/allocations.h"

int LU_precond_EBE_setup_NNOEL4 (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int nel = Parameters->nel;
	double **LUe, **A;

	if (tag == 1){
		MatrixData->LUeaux = mycalloc("Ueaux of 'LU_precond_EBE_setup'", 16*nel,sizeof(double));
		MatrixData->LUe = mycalloc("Ue of 'LU_precond_EBE_setup'",nel,sizeof(double*));
		for (I = 0; I < nel; I++){
			MatrixData->LUe[I] = &(MatrixData->LUeaux[16*I]);
		}
	}
	
	memset(MatrixData->LUeaux,0,16*nel*sizeof(double));// zerando a matriz LUe
	A = MatrixData->A;
	LUe = MatrixData->LUe;

	for (I = 0; I < nel; I++){
	/*	A[I][0] = 2;
		A[I][1] = 1;
		A[I][2] = 1;
		A[I][3] = 0;
		A[I][4] = 4;
		A[I][5] = 3;
		A[I][6] = 3;
		A[I][7] = 1;
		A[I][8] = 8;
		A[I][9] = 7;
		A[I][10] = 9;
		A[I][11] = 5;
		A[I][12] = 6;
		A[I][13] = 7;
		A[I][14] = 9;
		A[I][15] = 8; */
		
		// line 1
		LUe[I][0] = A[I][0]; // u11
		LUe[I][1] = A[I][1]; // u12
		LUe[I][2] = A[I][2]; // u13
		LUe[I][3] = A[I][3]; // u14
		// line 2
		LUe[I][4] = A[I][4]/A[I][0]; // multiplicador l21
		LUe[I][5] = A[I][5] - LUe[I][4]*A[I][1]; // u22 
		LUe[I][6] = A[I][6] - LUe[I][4]*A[I][2]; // u23
		LUe[I][7] = A[I][7] - LUe[I][4]*A[I][3]; // u24
		// line 3
		LUe[I][8]  = A[I][8]/A[I][0]; // multiplicador l31
		LUe[I][9]  = (A[I][9] - LUe[I][8]*A[I][1])/LUe[I][5]; // multiplicador l32 
		LUe[I][10] = (A[I][10] - LUe[I][8]*A[I][2]) - LUe[I][9]*LUe[I][6]; // u33 
		LUe[I][11] = (A[I][11] - LUe[I][8]*A[I][3]) - LUe[I][9]*LUe[I][7]; // u34
		// line 4
		LUe[I][12] = A[I][12]/A[I][0]; // multiplicador l41
		LUe[I][13] = (A[I][13] - LUe[I][12]*A[I][1])/LUe[I][5]; // multiplicador l42 
		LUe[I][14] = ((A[I][14] - LUe[I][12]*A[I][2]) - LUe[I][13]*LUe[I][6])/LUe[I][10]; //multiplicador 43
		LUe[I][15] = ((A[I][15] - LUe[I][12]*A[I][3]) - LUe[I][13]*LUe[I][7]) - LUe[I][14]*LUe[I][11]; // u44
		
		/*printf("Matriz LU elemento %d\n", I);
		printf("%lf\t%lf\t%lf\t%lf\n", LUe[I][0], LUe[I][1], LUe[I][2], LUe[I][3]);
		printf("%lf\t%lf\t%lf\t%lf\n", LUe[I][4], LUe[I][5], LUe[I][6], LUe[I][7]);
		printf("%lf\t%lf\t%lf\t%lf\n", LUe[I][8], LUe[I][9], LUe[I][10], LUe[I][11]);
		printf("%lf\t%lf\t%lf\t%lf\n", LUe[I][12], LUe[I][13], LUe[I][14], LUe[I][15]);
		getchar();*/
	}

	return 0;

}

