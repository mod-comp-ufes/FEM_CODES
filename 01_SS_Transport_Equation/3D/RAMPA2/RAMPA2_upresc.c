#include "Rampa2.h" 
#include <math.h>

double RAMPA2_upresc(double x, double y, double z, double Re, double Lx)
{
	double u;
	
	u = 0.0;
	
	if(y == 0.0){
		if(x > 0.5){
			u = 1.0;
		}else{
			u = 0.0;
		}
	}
	
	return u;

}
