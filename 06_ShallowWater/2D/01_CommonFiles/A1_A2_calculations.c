#include "ShalowWater.h"


void A1_A2_calculations(double Ub[3], double A1[3][3], double A2[3][3], double g)
{
	double h = Ub[0];
	double u = Ub[1]/h;
	double v = Ub[2]/h;

	// *** A1 coefficients
	A1[0][0] = 0.0;
	A1[0][1] = 1.0;
	A1[0][2] = 0.0;
	A1[1][0] = g*h - u*u;
	A1[1][1] = 2.0*u;
	A1[1][2] = 0.0;
	A1[2][0] = -u*v;
	A1[2][1] = v;
	A1[2][2] = u;

	// *** A2 coefficients
	A2[0][0] = 0.0;
	A2[0][1] = 0.0;
	A2[0][2] = 1.0;
	A2[1][0] = -u*v;
	A2[1][1] = v;
	A2[1][2] = u;
	A2[2][0] = g*h - v*v;
	A2[2][1] = 0.0;
	A2[2][2] = 2*v;

}
