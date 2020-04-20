#include "cylinder.h" 
#include <math.h>

double CYLINDER_v1presc(double X, double Y)
{
	double v;
	//volker Jhon
	if (fabs(X-0.0)<=1e-15 || fabs(X-2.2)<=1e-15)
		v = 5.948839976*(1.2*Y*(0.41-Y));
	
	//Erturk
	//if (fabs(X+5.0) <= 1e-15 || fabs(Y-5.0) <= 1e-15)
	//	v = 1.0;
	else
		v = 0.0;

	return v;
// fabs( ((Node[I].x-4.5)*(Node[I].x-4.5) + (Node[I].y-4.5)*(Node[I].y-4.5)) - 0.25)<=1e-6)
}






