#include "SSTransportEquation3D.h"

int setStabilizationForm(ParametersType *Parameters, FemOtherFunctionsType *FemOtherFunctions){
	
	if (strcasecmp(Parameters->StabilizationForm,"DD")==0){
	
		FemOtherFunctions->Build = Build_DD;

	}else if (strcasecmp(Parameters->StabilizationForm,"CAU")==0){
	
		FemOtherFunctions->Build = Build_CAU;

	}else if (strcasecmp(Parameters->StabilizationForm,"SUPG")==0){
	
		FemOtherFunctions->Build = Build_SUPG;
		
	}else if (strcasecmp(Parameters->StabilizationForm,"VMS")==0){
	
		FemOtherFunctions->Build = Build_VMS;

	}else if (strcasecmp(Parameters->StabilizationForm,"NSGS")==0){
	
		FemOtherFunctions->Build = Build_NSGS;

	}else if (strcasecmp(Parameters->StabilizationForm,"Galerkin")==0){
	
		FemOtherFunctions->Build = Build;

	}else {
		printf("Stabilizantion form is not defined correctly!\n");
		exit(1);
	}

	return 0;
}

