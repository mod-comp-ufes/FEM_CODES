#include "exata.h" 
#include <math.h>

double EXATA_f2ext(double t, double x, double y)
{
	double u1txy, u2txy, f;

	u1txy = pow(1+t,2)*pow(x,2)*pow(1-x,2)*(2*y-6*pow(y,2)+4*pow(y,3));
	u2txy = pow(1+t,2)*pow(y,2)*pow(1-y,2)*(-2*x+6*pow(x,2)-4*pow(x,3));		
		
	f = 2*(t+1)*pow(y,2)*pow(1-y,2)*(-2*x+6*pow(x,2)-4*pow(x,3)) - u1txy*pow(t+1,2)*pow(y,2)*(-2+12*x-12*pow(x,2))*pow(1-y,2) + u2txy*pow(t+1,2)*(2*y-6*pow(y,2)+4*pow(y,3))*(-2*x+6*pow(x,2)-4*pow(x,3)) + 2*pow(t+1,2)*(y+0.02*((1-6*(y-pow(y,2)))*(-x+3*pow(x,2)-2*pow(x,3)) + pow(1-y,2)*pow(y,2)*(3-6*x)));
	
	return f;

}






