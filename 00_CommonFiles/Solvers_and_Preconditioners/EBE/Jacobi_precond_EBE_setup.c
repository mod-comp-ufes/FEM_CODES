#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int Jacobi_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int neq = Parameters->neq;
	int nel = Parameters->nel;

	if (tag==1){
		MatrixData->invDeaux = mycalloc("invDeaux of 'Jacobi_precond_EBE_setup'",NNOEL*nel,sizeof(double));
		MatrixData->invDe = mycalloc("invDe of 'Jacobi_precond_EBE_setup'",nel,sizeof(double*));
		for (I=0; I<nel; I++)
			MatrixData->invDe[I] = &(MatrixData->invDeaux[NNOEL*I]);
	}

	for (I=0; I<nel; I++){
		MatrixData->invDe[I][0] = 1.0 / (MatrixData->A[I][0]);
		MatrixData->invDe[I][1] = 1.0 / (MatrixData->A[I][4]);
		MatrixData->invDe[I][2] = 1.0 / (MatrixData->A[I][8]);

		/*
		MatrixData->invDe[I][0] = 1.0 / (1.0 + MatrixData->A[I][0]);
		MatrixData->invDe[I][1] = 1.0 / (1.0 + MatrixData->A[I][4]);
		MatrixData->invDe[I][2] = 1.0 / (1.0 + MatrixData->A[I][8]);
		*/
	}

	/* F preconditioning */
	double *faux = calloc((neq + 1), sizeof(double));
	for (I = 0; I < neq; I++){
		faux[I] = F[I];
	}
	faux[neq] = 0.0;
	F[neq] = 0.0;

	Jacobi_precond_EBE (Parameters, MatrixData, FemStructs, faux, F);

	free(faux);

	return 0;

}
