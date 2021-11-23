#include "difrea.h" 
#include <math.h>

double DIFREA_DuDy(double x, double y, double z, double Cst)
{
	double u;
	
	u = PI*sin(PI*x)*cos(PI*y)*sin(PI*z);
	
	return u;

}
