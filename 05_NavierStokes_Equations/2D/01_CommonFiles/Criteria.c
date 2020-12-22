#include "NavierStokesEquations.h"

void setStopCriteria(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if(strcasecmp(Parameters->StopMulticorrection,"ITERATION") == 0)
		FemFunctions->StopCriteria = StopByIterations;
	else if(strcasecmp(Parameters->StopMulticorrection,"NORM") == 0)
		FemFunctions->StopCriteria = StopByNorm;
	else{
		printf("Stop criteria not defined!\n");
		exit(1);
	}
}


int StopByIterations(ParametersType *Parameters, double norm_a, double norm_Da, int i)
{
	if (i>=Parameters->NumberCorrection)
		return 1;
	else
		return 0;
}



int StopByNorm(ParametersType *Parameters, double norm_a, double norm_Da, int i)
{
	if (norm_Da < (Parameters->CorrectionTolerance)*norm_a || i>=Parameters->NumberCorrection)
		return 1;
	else
		return 0;
}


