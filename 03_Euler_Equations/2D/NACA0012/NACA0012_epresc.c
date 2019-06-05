#include "naca0012.h"

double NACA0012_epresc(double x, double y)
{
	double e;
	
	e = 0.9464; // rho * e = 7.64286, para c_v = 7.14286
	
	return e;
}

 
