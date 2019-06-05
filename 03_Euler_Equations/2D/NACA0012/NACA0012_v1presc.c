#include "naca0012.h"

double NACA0012_v1presc(double x, double y){
	double v1;
	
	if ( ( fabs(x) <= 1e-15 && fabs(y) <= 1e-15 ) || ( fabs(x - 1.0) <= 1e-15 && fabs(y) <= 1e-15 ) ){
		v1 = 0.0;
	}else{
		v1 = 1.0;
	}
	
	return v1;

}
 
