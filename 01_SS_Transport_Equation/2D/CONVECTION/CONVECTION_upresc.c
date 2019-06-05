#include "../01_CommonFiles/SSTranspEquation.h"
#include "convection.h" 

double CONVECTION_upresc(double X, double Y)
{
	double u;

	//Pudim
	if ((fabs(X)<=1e-14 && Y<0.25) || fabs(Y)<1e-14)
		u = 0.0;
	else
		u = 1.0;

	return u;

}






