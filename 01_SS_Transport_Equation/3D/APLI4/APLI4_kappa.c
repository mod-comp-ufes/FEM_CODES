#include "apli4.h" 
#include <math.h>

void APLI4_kappa(double *kx, double *ky, double *kz, double x, double y, double z, double A)
{
	*kx = 1.0/A;
	*ky = 1.0/A;
	*kz = 1.0/A;
}
