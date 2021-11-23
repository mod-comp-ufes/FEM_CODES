#include "parede.h" 
#include <math.h>

double PAREDE_upresc(double x, double y, double z, double Re, double Lx)
{
	double u = 0.0;
	
	if (fabs(y) < 1e-8){ // y = 0

		u = 0.0;
		if (x - 1.0 >= 0.0){ // y = 0 e x > 1
			u = 1.0;
		}

	}
	//else if (fabs(x - 2.0) < 1e-8 || fabs(x) < 1e-8 || fabs(y - 2.0) < 1e-8){ // x = 2, x = 0, y = 2

	//	u = 0.0;

	//}
	
	return u;

}
