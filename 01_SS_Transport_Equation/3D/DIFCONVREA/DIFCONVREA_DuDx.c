
#include "difconvrea.h" 

double DIFCONVREA_DuDx(double x, double y, double z, double Cst)
{
	double u;
	
	u = 100*(1-2*x)*(1-y)*y*(1-z)*z;
	
	return u;

}
