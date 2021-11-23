#include "difrea.h" 
#include <math.h>

double DIFREA_DuDx(double x, double y, double z, double Cst)
{
	double u;
	
	u = PI*cos(PI*x)*sin(PI*y)*sin(PI*z);
	
	return u;

}
