#include "../preconditioners.h"

int SGSScaling_precond_EBE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I, neq = Parameters->neq;
	double *Diag = MatrixData->Diag;
	
	Diag_precond_EBE_setup(Parameters, MatrixData, FemStructs, 1, F);

	for (I=0; I<neq; I++)
		F[I]*=sqrt(Diag[I]);
	
	SGSScaling_precond_EBE (Parameters, MatrixData, FemStructs, F, F);

	return 0;
}


