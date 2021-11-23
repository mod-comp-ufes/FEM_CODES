#include "difconv.h" 
#include <math.h>

double DIFCONV_f(double x, double y, double z, double A)
{
	double f;

	f = A*3.0*PI*PI*sin(PI*x)*sin(PI*y)*sin(PI*z) + PI*cos(PI*x)*sin(PI*y)*sin(PI*z) + PI*sin(PI*x)*cos(PI*y)*sin(PI*z) + PI*sin(PI*x)*sin(PI*y)*cos(PI*z);
	
	return f;

}






