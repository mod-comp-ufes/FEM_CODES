#include "ourBLAS.h"

double dmax(int n, double *v) {
	int i;
	double max = v[0];
	
	for(i = 1; i < n; i++){
		if(max < v[i]){
			max = v[i];
		}
	}

	return max;
}
