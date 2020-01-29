#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int SGS_precond_EDE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I;
	int neq = Parameters->neq;
	int nedge = Parameters->nedge;

	if (tag==1){
		MatrixData->invDeaux = mycalloc("invDeaux of 'SGS_precond_EDE_setup'",2*nedge,sizeof(double));
		MatrixData->invDe = mycalloc("invDe of 'SGS_precond_EDE_setup'",nedge,sizeof(double*));
		for (I=0; I<nedge; I++)
			MatrixData->invDe[I] = &(MatrixData->invDeaux[2*I]);
	}

	for (I=0; I<nedge; I++){
		MatrixData->invDe[I][0] = 1.0/(1.0 + MatrixData->A[I][0]);
		MatrixData->invDe[I][1] = 1.0/(1.0 + MatrixData->A[I][3]);
	}
	/* F preconditioning */
	double *faux = mycalloc("faux of SGS_precond_EDE_setup",(neq + 1), sizeof(double));
	for (I = 0; I < neq; I++){
		faux[I] = F[I];
	}
	faux[neq] = 0.0;
	F[neq] = 0.0;

	SGS_precond_EDE (Parameters, MatrixData, FemStructs, faux, F);

	myfree(faux);

	return 0;
}
