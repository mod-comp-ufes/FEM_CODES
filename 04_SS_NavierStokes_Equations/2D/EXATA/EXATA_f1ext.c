#include "exata.h" 
#include <math.h>

double EXATA_f1ext(double x, double y)
{
	double f;

	//f = 0.1;
	// mu = 0.01
	f = 2*x + pow(x,4)*(0.12 - 0.24*y) + pow(x,3)*(-0.24 + 0.48*y) + y*(-0.04 + 0.12*y - 0.08*pow(y,2)) + 8*pow((-1 + x),3)*pow(x,3)*(-1 + 2*x)*pow(y,2)*pow((1 - 3*y + 2*pow(y,2)),2) + x*y*(0.24 - 0.72*y + 0.48*pow(y,2)) - 4*pow((-1 + x),2)*pow(x,3)*(1 - 3*x + 2*pow(x,2))*pow((-1 + y),2)*pow(y,2)*(1 - 6*y + 6*pow(y,2)) + pow(x,2)*(0.12 - 0.48*y + 0.72*pow(y,2) - 0.48*pow(y,3));
	//
	// mu = 0.1
	//f = 2*x + pow(x,4)*(1.2 - 2.4*y) + pow(x,3)*(-2.4 + 4.8*y) + y*(-0.4 + 1.2*y - 0.8*pow(y,2)) + 8*pow((-1 + x),3)*pow(x,3)*(-1 + 2*x)*pow(y,2)*pow((1 - 3*y + 2*pow(y,2)),2) + x*y*(2.4 - 7.2*y + 4.8*pow(y,2)) - 4*pow((-1 + x),2)*pow(x,3)*(1 - 3*x + 2*pow(x,2))*pow((-1 + y),2)*pow(y,2)*(1 - 6*y + 6*pow(y,2)) + pow(x,2)*(1.2 - 4.8*y + 7.2*pow(y,2) - 4.8*pow(y,3));

	return f;

}






