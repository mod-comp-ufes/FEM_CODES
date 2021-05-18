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
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") == 0) {
		FemFunctions->assembly = csr_assembly;
		FemFunctions->ProductMatrixVector = csrmv;
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE") == 0) {
		FemFunctions->assembly = ebe_assembly;
		FemFunctions->ProductMatrixVector = ebemvNDOF3;
	}
	FemFunctions->unscaling = NO_unscaling;
	FemFunctions->scaling = NO_scaling;

	return 0;
}
