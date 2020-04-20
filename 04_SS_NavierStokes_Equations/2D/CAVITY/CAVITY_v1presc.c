#include "cavity.h" 
#include <math.h>

double CAVITY_v1presc(double X, double Y)
{
	double v, pi=3.14159265359, eps1, eps2;
	eps1 = 0.1;
	eps2 = 1.0 - eps1;
		
	if (fabs(Y-1.0)<=1e-15){
		if(X<eps1)
			v = 1.0-0.25*(1-cos((eps1-X)*pi/eps1))*(1-cos((eps1-X)*pi/eps1));
		else if(X>eps2)
			v = 1.0-0.25*(1-cos((X-eps2)*pi/eps1))*(1-cos((X-eps2)*pi/eps1));
		else 
		v = 1.0;
	}else
		v = 0.0;

	return v;

}






