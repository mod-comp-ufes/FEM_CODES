#include "ShalowWater.h"


int setStabilizationForm(ParametersType *Parameters,FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions,
						int (**TimeIntegration)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *))
{
	
	if (strcasecmp(Parameters->StabilizationForm,"SUPG")==0)
	{
		FemOtherFunctions->Build = Build_M_D_F_SUPG;
		if (strcasecmp(Parameters->TimeIntegration,"Predictor")==0){		
			*TimeIntegration = PredictorMulticorrector;
		}
		else{
			printf("Time integration method is not defined correctly!\n");
			exit(1);
		}

		if (strcasecmp(Parameters->ShockCapture,"delta91-MOD")==0){
			FemFunctions->ShockCapture = delta91_MOD;
		}
		else{
			printf("Shock capture is not defined correctly!\n");
			exit(1);
		}
	}
	else {
		printf("Stabilizantion form is not defined!\n");
		exit(1);
	}

	return 0;
}
