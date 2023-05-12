#include "SSTranspEquation.h"

int setStabilizationForm(ParametersType *Parameters, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	if (strcasecmp(Parameters->StabilizationForm,"NOT")==0) {
		FemOtherFunctions->Build = Build_Galerkin;
	}
	else if (strcasecmp(Parameters->StabilizationForm,"SUPG")==0){
	
		FemOtherFunctions->Build = Build_K_F_SUPG;

		if (strcasecmp(Parameters->ShockCapture,"CAU")==0){
			FemFunctions->ShockCapture = CAU_ShockCapture;
		}
		else if (strcasecmp(Parameters->ShockCapture,"YZBeta")==0){
			FemFunctions->ShockCapture = YZBeta_ShockCapture;
		}	
		else{
			printf("Shock capture is not defined correctly!\n");
			exit(1);
		}
		
	}
	else if (strcasecmp(Parameters->StabilizationForm,"DD")==0){
	
		FemOtherFunctions->Build = Build_K_F_DD;
		
		if (strcasecmp(Parameters->ShockCapture,"CAUDD")==0){
			FemFunctions->ShockCapture = CAU_DD_ShockCapture;
		}
		else if (strcasecmp(Parameters->ShockCapture,"CAU")==0){
			FemFunctions->ShockCapture = CAU_ShockCapture;
		}
		else if (strcasecmp(Parameters->ShockCapture,"YZBeta")==0){
			FemFunctions->ShockCapture = YZBeta_ShockCapture;
		}	
		else{
			printf("Shock capture is not defined correctly!\n");
			exit(1);
		}	
	}
	else {
		printf("Stabilizantion form is not defined correctly!\n");
		exit(1);
	}

	//Set h shock parameter
	if (strcasecmp(Parameters->h_Shock,"NOT")==0){ 
		
	}else if (strcasecmp(Parameters->h_Shock,"sqrtArea")==0){ 
		FemFunctions->h_shock = h_shock_sqrtArea;
	}else if (strcasecmp(Parameters->h_Shock,"2sqrtArea")==0){ 
		FemFunctions->h_shock = h_shock_2sqrtArea;
	}
	else if (strcasecmp(Parameters->h_Shock,"Option1")==0){  
		FemFunctions->h_shock = h_shock_Option1;
	}
	else if (strcasecmp(Parameters->h_Shock,"Option2")==0){ 
		FemFunctions->h_shock = h_shock_Option2;
	}
 	else{
		printf("h_shock is not defined correctly!\n");
		exit(1);
	}
	return 0;
}


