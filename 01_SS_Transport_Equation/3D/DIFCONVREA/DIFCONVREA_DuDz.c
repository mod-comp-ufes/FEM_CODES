
#include "difconvrea.h" 

double DIFCONVREA_DuDz(double x, double y, double z, double Cst)
{
	double u;
	
	u = 100*(1-x)*x*(1-y)*y*(1-2*z);
	
	return u;

}
