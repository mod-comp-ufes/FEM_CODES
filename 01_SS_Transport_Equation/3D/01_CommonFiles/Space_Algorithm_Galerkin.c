#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

void Space_Algorithm_Galerkin(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	
	double *u, *F;
	
	u = FemStructs->u;
	F = FemStructs->F;

	FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
	
	FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, 1, F);
	
	FemOtherFunctions->Solver(Parameters, MatrixData, FemStructs, FemFunctions, F, u);
	printf("IterGMRES: %d\n", Parameters->ContGMRES);
	
	if (strcasecmp(Parameters->ExactSolution,"YES")==0){
		Calculating_Errors(Parameters, FemStructs, FemFunctions);
	}
	
}	
	
