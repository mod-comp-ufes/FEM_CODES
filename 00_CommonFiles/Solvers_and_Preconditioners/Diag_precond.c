#include "preconditioners.h"

int Diag_precond(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z)
{		
	int I, neq = Parameters->neq;
	double *invDiag = MatrixData->invDiag;


	for (I=0; I<neq; I++)
		z[I] = p[I]*invDiag[I];


	return 0;
}
	
