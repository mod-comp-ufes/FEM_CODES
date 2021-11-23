#include "high.h" 
#include <math.h>

double HIGH_upresc(double x, double y, double z, double Re, double Lx)
{
	double u=0.0;
	
	if((x == 0.0) || (y == 0.0) || (z == 0.0)){
		u = 0.0;
	}else if((x == 1.0) || (y == 1.0) || (z == 1.0)){
		u = 1.0;
	}
	
	return u;

}
