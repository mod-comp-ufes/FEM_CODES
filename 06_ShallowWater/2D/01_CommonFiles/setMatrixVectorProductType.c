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
	if(strcasecmp(Parameters->reordering,"NOT") != 0 && strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") !=0)
	{
		printf("Reordering is applied only with CSR scheme!\n");
		exit(1);
	}
	else if(strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0)
	{
		FemFunctions->assembly = csr_assembly;
		FemFunctions->mv = csrmv;
	}
	else
	{
		printf("\nMatrix vector product scheme is not defined correctly!\n\n");
		exit(1);
	}

	FemFunctions->unscaling = NO_unscaling;
	FemFunctions->scaling = NO_scaling;

	return 0;
}
