#include "ShalowWater.h"


int setStabilizationForm(ParametersType *Parameters,FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions,
						int (**TimeIntegration)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *))
{
	if (strcasecmp(Parameters->StabilizationForm,"SUPG")==0)
	{
		FemOtherFunctions->Build = Build_M_D_F_SUPG;
		if (strcasecmp(Parameters->TimeIntegration,"Predictor")==0) {
			*TimeIntegration = PredictorMulticorrector;
			setStopCriteria(Parameters, FemFunctions);
		}
		else
		{
			printf("Time integration method is not defined correctly!\n");
			exit(1);
		}

		if (strcasecmp(Parameters->ShockCapture,"CAU")==0)
			FemFunctions->ShockCapture = deltaCAU;
		else{
			printf("Shock capture is not defined correctly!\n");
			exit(1);
		}
	}
	else
	{
		printf("Stabilizantion form is not defined!\n");
		exit(1);
	}

	return 0;
}


void setStopCriteria(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if(strcasecmp(Parameters->StopMulticorrection,"ITERATION") == 0)
		FemFunctions->StopCriteria = StopByIterations;
	else if(strcasecmp(Parameters->StopMulticorrection,"NORM") == 0)
		FemFunctions->StopCriteria = StopByNorm;
	else{
		printf("Stop criteria is not defined!\n");
		exit(1);
	}

	if(strcasecmp(Parameters->StopAtSteadyState,"YES") == 0)
		FemFunctions->StopTimeIntegration = StopBySteadyState;
	else if(strcasecmp(Parameters->StopAtSteadyState,"NOT") == 0)
		FemFunctions->StopTimeIntegration = StopByTime;
	else{
		printf("Stop time integration is not defined!\n");
		exit(1);
	}

}