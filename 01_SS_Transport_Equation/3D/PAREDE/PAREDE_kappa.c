#include "parede.h" 
#include <math.h>

void PAREDE_kappa(double *kx, double *ky, double *kz, double x, double y, double z, double Cst)
{
	*kx = Cst;
	*ky = Cst;
	*kz = Cst;
}
