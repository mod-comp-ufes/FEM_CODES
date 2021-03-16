#include "preconditioners.h"

int Diag_precond_EDE_setup(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int E, I, K, nedge, neq, size=2*NDOF; 
	double *Diag, *invDiag, **A;
	int **lm = FemStructs->lm2;

	A = MatrixData->A;
	Diag = MatrixData->Diag;
	invDiag = MatrixData->invDiag;
	neq = Parameters-> neq;
	nedge = Parameters-> nedge;
	
	for (I=0; I<neq; I++)
		Diag[I] = 0;

	for (E=0; E<nedge; E++){
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


	return 0;
}


