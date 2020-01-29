#include "SSNavierStokesEquations.h"

int setStabilizationForm(ParametersType *Parameters,FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	
	if (strcasecmp(Parameters->StabilizationForm,"SUPG")==0){
			FemOtherFunctions->Build = Build_K_F_SUPG;
	}		
	else if (strcasecmp(Parameters->StabilizationForm,"MS")==0){
			FemOtherFunctions->Build = Build_K_F_MS;
	}
	else if (strcasecmp(Parameters->StabilizationForm,"MS_Time")==0){
			FemOtherFunctions->Build = Build_K_F_MS_Time;
	}		
	else {
		printf("Stabilizantion form not defined!\n");
		exit(1);
	}

	return 0;
}

