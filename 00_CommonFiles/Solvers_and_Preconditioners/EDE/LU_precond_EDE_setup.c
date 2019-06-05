#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int LU_precond_EDE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int nedge = Parameters->nedge;
	double **LUe, **A;

	if (tag == 1){
		MatrixData->LUeaux = mycalloc("Ueaux of 'LU_precond_EBE_setup'", 4*nedge,sizeof(double));
		MatrixData->LUe = mycalloc("Ue of 'LU_precond_EBE_setup'",nedge,sizeof(double*));
		for (I=0; I<nedge; I++){
			MatrixData->LUe[I] = &(MatrixData->LUeaux[4*I]);
		}
	}
	A = MatrixData->A;
	LUe = MatrixData->LUe;

	for (I = 0; I < nedge; I++){
		LUe[I][0] = 1.0;
		LUe[I][1] = A[I][1];
		LUe[I][2] = A[I][2];
	//	LUe[I][3]= 1./sqrt(1.0 - A[I][1]*A[I][2]); //store 1/Lue[4] to improve LU linear solver
		LUe[I][3] = 1./(1.0 - A[I][1]*A[I][2]);
	}

	LU_precond_EBE (Parameters, MatrixData, FemStructs, F, F);

	return 0;

}

