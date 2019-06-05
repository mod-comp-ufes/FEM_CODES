#include "cylinder.h"

double CYLINDER_v1presc(double x, double y){
	double v1;

	if ( (fabs(x)-0.5) <= 1e-15 && fabs(y) <= 1e-15 ){
		v1 = 0.0;
	}else{
		v1 = 1.0;
	}
	
	return v1;

}
 
