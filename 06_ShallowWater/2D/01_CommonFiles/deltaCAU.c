#include "ShalowWater.h"


double normA0invtilde(double *R, double h, double u, double v)
{
	double r;
	r = 1/h*(R[0]*(R[0]*(GRAVITY*h + u*u + v*v) - R[1]*u - R[2]*v) +
	         R[1]*(-R[0]*u + R[1]) +
		     R[2]*(-R[0]*v + R[2]));

	return r;
}

double deltaCAU(double *Ub, double *Sb, double *gradEta, double *gradUx, double *gradUy,
                double (*A1)[3], double (*A2)[3],
                double y23, double y31, double y12, double x32, double x13, double x21, double twoArea, 
                double tau, ParametersType *Parameters)
{
	int i;
	double h = Ub[0], u = Ub[1]/Ub[0], v = Ub[2]/Ub[0];
	double delta = 0.0, delta_91, delta_tau;
	double gradksiUx[3], gradksiUy[3], R[3], AgradU[3], normR, normgradU, normgradksiU;

	for(i=0; i<3; i++)
	{
		// R = A1gradUx + A2gradUy - S
		AgradU[i] = A1[i][0]*gradUx[0] + A1[i][1]*gradUx[1] + A1[i][2]*gradUx[2] + 
			        A2[i][0]*gradUy[0] + A2[i][1]*gradUy[1] + A2[i][2]*gradUy[2];
		R[i] = AgradU[i] - Sb[i];

		gradksiUx[i] = (y31*gradUx[i] + x13*gradUy[i])/twoArea;
		gradksiUy[i] = (y12*gradUx[i] + x21*gradUy[i])/twoArea;
	}

	normgradksiU = sqrt(normA0invtilde(gradksiUx, h, u, v) + normA0invtilde(gradksiUy, h, u, v));
	normgradU = normA0invtilde(gradUx, h, u, v) + normA0invtilde(gradUy, h, u, v);
	normR = sqrt(normA0invtilde(R, h, u, v));

	delta_91 = (fabs(normgradksiU) >= TOL) ? normR/normgradksiU : 0.0;
	delta_tau = (fabs(normgradU) >= TOL) ? tau*normR/normgradU : 0.0;

	if(normgradU <= TOL)
		delta = 0.0;
	else
		delta = max(0.0, delta_91 - delta_tau);

	return delta;
}
