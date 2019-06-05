#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int SSOR_precond_EBE_setup2 (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I, E, K;
	int **lm = FemStructs->lm;
	int neq = Parameters->neq;
	int nel = Parameters->nel;
	int size = NNOEL;
	double *Diag, *invDiag, **A;

	if (tag==1){
		MatrixData->invDeaux = mycalloc("invDeaux of 'SOR_precond_EBE_setup'",NNOEL*nel,sizeof(double));
		MatrixData->invDe = mycalloc("invDe of 'SOR_precond_EBE_setup'",nel,sizeof(double*));
		for (I=0; I<nel; I++)
			MatrixData->invDe[I] = &(MatrixData->invDeaux[NNOEL*I]);
	}

	A = MatrixData->A;
	Diag = MatrixData->Diag;
	invDiag = MatrixData->invDiag;

	for (I=0; I<neq; I++)
		Diag[I] = 0;

	for (E=0; E<nel; E++){
		for (I=0, K=0; I<size; I++){
			Diag[lm[E][I]] += A[E][K];
			K += size+1;
		}
		Diag[neq] = 0;
	}

	for(I=0; I<neq; I++)
		invDiag[I] = 1.0/Diag[I];

	for(I=0; I<neq; I++)
		F[I] *= invDiag[I];

	invDiag[neq] = 0.0;
	for (I=0; I<nel; I++){
		MatrixData->invDe[I][0] = 1.0/(1.0 + invDiag[lm[I][0]]*A[I][0]);
		MatrixData->invDe[I][1] = 1.0/(1.0 + invDiag[lm[I][1]]*A[I][4]);
		MatrixData->invDe[I][2] = 1.0/(1.0 + invDiag[lm[I][2]]*A[I][8]);
	}

	/* F preconditioning */
	double *faux = calloc((neq + 1), sizeof(double));
	for (I = 0; I < neq; I++){
		faux[I] = F[I];
	}
	faux[neq] = 0.0;
	F[neq] = 0.0;

	SSOR_precond_EBE (Parameters, MatrixData, FemStructs, faux, F);

	free(faux);

	return 0;
}
