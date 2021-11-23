#include "difrea.h" 
#include <math.h>

double DIFREA_f(double x, double y, double z, double A)
{
	double f;

	f = (A*3.0*PI*PI + 1)*sin(PI*x)*sin(PI*y)*sin(PI*z);
	
	return f;

}






