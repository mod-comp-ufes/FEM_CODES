#include "ShalowWater.h"
#include "../../../00_CommonFiles/Time_Integration/time_integration.h"


int setStabilizationForm(ParametersType *Parameters,FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions,
						int (**Predictor)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *))
{
	
	if (strcasecmp(Parameters->StabilizationForm,"SUPG")==0)
	{
		FemOtherFunctions->Build = Build_M_D_F_SUPG;
		if (strcasecmp(Parameters->TimeIntegration,"Predictor1")==0){		
			*Predictor = Predictor_Old;
		}
		else if (strcasecmp(Parameters->TimeIntegration,"Predictor1_BDF")==0){		
			*Predictor = Predictor_Old_BDF;
		}
		else if (strcasecmp(Parameters->TimeIntegration,"Predictor1_TRBDF2")==0){		
			*Predictor = Predictor_Old_TRBDF2;
		}
		else{
			printf("Time integration method is not defined correctly!\n");
			exit(1);
		}

		if (strcasecmp(Parameters->ShockCapture,"CAU")==0){
			FemFunctions->ShockCapture = Delta_CAU;
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
