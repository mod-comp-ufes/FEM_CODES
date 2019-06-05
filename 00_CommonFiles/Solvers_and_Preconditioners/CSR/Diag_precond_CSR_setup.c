#include "../preconditioners.h"

int Diag_precond_CSR_setup(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I, J, neq, *IA, *JA;
	double *Diag, *invDiag, *AA;

	AA = MatrixData->AA;
	JA = MatrixData->JA;
	IA = MatrixData->IA;
	Diag = MatrixData->Diag;
	invDiag = MatrixData->invDiag;
	neq = Parameters->neq;

	for (I=0; I<neq; I++){
		for (J=IA[I]; J<IA[I+1]; J++){
			if (JA[J] == I)
				Diag[I] = AA[J];	
		}
	}

	for(I=0; I<neq; I++)
		invDiag[I] = 1.0/Diag[I];	
	
	for(I=0; I<neq; I++)
		F[I] *= invDiag[I];
		
	return 0;
}



