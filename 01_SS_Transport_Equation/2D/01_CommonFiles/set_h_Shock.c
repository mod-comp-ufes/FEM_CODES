#include "SSTranspEquation.h"

double h_shock_not(ParametersType *Parameters, ElementType *Element, int e, double Be_x, double Be_y, double u1, double u2, double u3, double y23, double y31, double y12, double x32, double x13, double x21, double Area);


int set_h_Shock(ParametersType *Parameters, double (**h_shock_calculation)(double, double, double, double, double, double, double, double, double, double, double, double))
{
	if (strcasecmp(Parameters->h_Shock,"NOT")==0){ 
		*h_shock_calculation = h_shock_not;
	}
	else if (strcasecmp(Parameters->h_Shock,"2sqrtArea")==0){ 
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


double h_shock_not(ParametersType *Parameters, ElementType *Element, int e, double Be_x, double Be_y, double u1, double u2, double u3, double y23, double y31, double y12, double x32, double x13, double x21, double Area)
{
	return 0.0;
}


