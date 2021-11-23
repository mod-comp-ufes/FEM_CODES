#include "difconv.h" 
#include <math.h>

void DIFCONV_kappa(double *kx, double *ky, double *kz, double x, double y, double z, double A)
{
	*kx = A;
	*ky = A;
	*kz = A;
}
