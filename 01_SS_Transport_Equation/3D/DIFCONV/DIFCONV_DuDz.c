#include "difconv.h" 
#include <math.h>

double DIFCONV_DuDz(double x, double y, double z, double Cst)
{
	double u;
	
	u = PI*sin(PI*x)*sin(PI*y)*cos(PI*z);
	
	return u;

}
