#include "../01_CommonFiles/SSTranspEquation.h"
#include "rampa2.h" 

double RAMPA2_upresc(double X, double Y)
{
	double u;

	if (fabs(Y)<1e-14 && fabs(X) >=0.5)
		u = 1.0;
	else
		u = 0.0;

	return u;

}






