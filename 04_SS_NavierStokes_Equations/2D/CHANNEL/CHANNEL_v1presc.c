#include "channel.h" 
#include <math.h>

double CHANNEL_v1presc(double X, double Y)
{
	double v;
	//Old
	/*if (fabs(X)<=1e-15)
		v = 1.0;
	else
		v = 0.0;
	*/
	//Erturk
	if (fabs(X-0.0) <= 1e-15)
		v = -10.66666667*Y*Y + 24.0*Y - 12.0;
	else
	
		v = 0.0;

	return v;

}






