#include "ShalowWater.h"
#include "../../../00_CommonFiles/MatrixVector_Operations/matvec.h"


int NO_scaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs)
{
	return 0;
}

int NO_unscaling(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *x)
{
	return 0;
} 

int setMatrixVectorProductType(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	FemFunctions->assembly = csr_assembly;
	FemFunctions->ProductMatrixVector = csrmv;
	//FemFunctions->mv = csrmv;
	FemFunctions->unscaling = NO_unscaling;
	FemFunctions->scaling = NO_scaling;

	return 0;
}
