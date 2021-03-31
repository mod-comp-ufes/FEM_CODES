#include <math.h>
#include "ShalowWater.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"


double delta91_MOD(double Ub[3], double gradUx[3], double gradUy[3], double A1[3][3], double A2[3][3], double Sb[3], 
                   double y23, double y31, double y12, double x32, double x13, double x21, double twoArea, 
				   double tau, double g)
{
	double h = Ub[0], u = Ub[1]/Ub[0], v = Ub[2]/Ub[0];
	double delta = 0.0, delta_91, delta_tau;
	double A1gradUx[3], A2gradUy[3], gradksiUx[3], gradksiUy[3], R[3], Rn, gradUn, gradksiUn;
	
	for (int i = 0; i < 3; i++)
	{
		gradksiUx[i] = (y31*gradUx[i] + x13*gradUy[i]) / twoArea;
		gradksiUy[i] = (y23*gradUx[i] + x32*gradUy[i]) / twoArea;

		A1gradUx[i] = A1[i][0]*gradUx[0] + A1[i][1]*gradUx[1] + A1[i][2]*gradUx[2];
		A2gradUy[i] = A2[i][0]*gradUy[0] + A2[i][1]*gradUy[1] + A2[i][2]*gradUy[2];

		// Calculo de R = AxgradUx + AygradUy - S em A0(-1)
		R[i] = A1gradUx[i] + A2gradUy[i] - Sb[i];
	}

	double a11 = g*h + u*u + v*v, 
	       a12 = -u, 
		   a13 = -v;

	Rn = sqrt(R[0]*(R[0]*a11 + R[1]*a12 + R[2]*a13) +
	     	  R[1]*(R[0]*a12 + R[1]) +
		 	  R[2]*(R[0]*a13 + R[2]));

	gradUn = (gradUx[0] + gradUy[0])*((gradUx[0] + gradUy[0])*a11 + (gradUx[1] + gradUy[1])*a12 + (gradUx[2] + gradUy[2])*a13) +
	     	 (gradUx[1] + gradUy[1])*((gradUx[0] + gradUy[0])*a12 + (gradUx[1] + gradUy[1])) +
		     (gradUx[2] + gradUy[2])*((gradUx[0] + gradUy[0])*a13 + (gradUx[2] + gradUy[2]));

	gradksiUn = sqrt((y31*gradUx[0] + x13*gradUy[0])*((y31*gradUx[0] + x13*gradUy[0])*a11 + (y31*gradUx[1] + x13*gradUy[1])*a12 + (y31*gradUx[2] + x13*gradUy[2])*a13) +
	     	         (y31*gradUx[1] + x13*gradUy[1])*((y31*gradUx[0] + x13*gradUy[0])*a12 + (y31*gradUx[1] + x13*gradUy[1])) +
		             (y31*gradUx[2] + x13*gradUy[2])*((y31*gradUx[0] + x13*gradUy[0])*a13 + (y31*gradUx[2] + x13*gradUy[2])) +
				     (y12*gradUx[0] + x21*gradUy[0])*((y12*gradUx[0] + x21*gradUy[0])*a11 + (y12*gradUx[1] + x21*gradUy[1])*a12 + (y12*gradUx[2] + x21*gradUy[2])*a13) +
	     	         (y12*gradUx[1] + x21*gradUy[1])*((y12*gradUx[0] + x21*gradUy[0])*a12 + (y12*gradUx[1] + x21*gradUy[1])) +
		             (y12*gradUx[2] + x21*gradUy[2])*((y12*gradUx[0] + x21*gradUy[0])*a13 + (y12*gradUx[2] + x21*gradUy[2])));

	delta_91 = Rn/gradksiUn;
	delta_tau = (Rn*tau)/gradUn;

	delta = max(0, delta_91 - delta_tau);

	return delta;
}
