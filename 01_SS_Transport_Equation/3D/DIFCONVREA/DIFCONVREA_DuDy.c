
#include "difconvrea.h" 

double DIFCONVREA_DuDy(double x, double y, double z, double Cst)
{
	double u;
	
	u = 100*(1-x)*x*(1-2*y)*(1-z)*z;
	
	return u;

}
