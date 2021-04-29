#include "ShalowWater.h"
#include "../../../00_CommonFiles/Time_Integration/time_integration.h"


int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int (*TimeIntegration)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);
	
	setProblem(Parameters, FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters,FemOtherFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions, &TimeIntegration);
	setStopCriteria(Parameters, FemFunctions);
	//TimeIntegration(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
	PredictorMulticorrector(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
	
	return 0;
}
