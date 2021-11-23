#include "apli6.h" 
#include <math.h>

double APLI6_upresc(double x, double y, double z, double A, double Lx)
{
	double u;
	
	u = exp(-A*x) + exp(-A*y) + exp(-A*z);
	
	return u;

}
