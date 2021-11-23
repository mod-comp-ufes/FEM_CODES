#include "hemker.h" 
#include <math.h>

double HEMKER_upresc(double x, double y, double z, double A, double Lx)
{
	double u = 0.0;
	
	if (y == 0.0){
		u = 0.0;
	}else if ( pow(y - 0.5,2) + pow(z - 0.5,2) - 0.04 < 1e-6){
		u = 1.0;
	}

	return u;

}
