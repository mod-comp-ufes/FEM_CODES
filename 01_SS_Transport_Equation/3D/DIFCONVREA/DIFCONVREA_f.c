
#include "difconvrea.h" 

double DIFCONVREA_f(double x, double y, double z, double k)
{
	double f;
	
	// f = -k(d^2u/dx^2 + d^2u/dy^2 + d^2u/dz^2) + 1du/dx + 1du/du + 1du/dz + 1u
	f = - k*( (-200*(1-y)*y*(1-z)*z) + (-200*(1-x)*x*(1-z)*z) + (-200*(1-x)*x*(1-y)*y) ) + 100*(1-2*x)*(1-y)*y*(1-z)*z + 100*(1-x)*x*(1-2*y)*(1-z)*z + 100*(1-x)*x*(1-y)*y*(1-2*z) + 100*x*y*z*(1-x)*(1-y)*(1-z);
	
	return f;

}






