#include "SSTranspEquation.h"

int set_h_Shock(ParametersType *Parameters, double (**h_shock_calculation)(double, double, double, double, double, double, double, double, double, double, double, double))
{
	if (strcasecmp(Parameters->h_Shock,"2sqrtArea")==0){ 
		*h_shock_calculation = h_shock_2sqrtArea;
	}
	else if (strcasecmp(Parameters->h_Shock,"Option1")==0){  
		*h_shock_calculation = h_shock_Option1;
	}
	else if (strcasecmp(Parameters->h_Shock,"Option2")==0){ 
		*h_shock_calculation = h_shock_Option2;
	}
 	else{
		printf("h_shock not defined!\n");
		exit(1);
	}

	return 0;
}


