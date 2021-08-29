#include "ShalowWater.h"


double deltaYZB(double *Ub, double *Sb, double *gradEta, double *gradUx, double *gradUy,
                double (*A1)[3], double (*A2)[3],
                double y23, double y31, double y12, double x32, double x13, double x21, double twoArea, 
                double tau, ParametersType *Parameters)
{
	int i;
	double AgradU[3], hshock, deltashock1, deltashock2, delta;
    double invhRef = 1.0/Parameters->href, 
           invqxRef = 1.0/Parameters->qxref,
           invqyRef = 1.0/Parameters->qyref;

    hshock = fabs(gradEta[0]*y23 + gradEta[1]*x32) +
		     fabs(gradEta[0]*y31 + gradEta[1]*x13) +
			 fabs(gradEta[0]*y12 + gradEta[1]*x21);
    hshock = (fabs(hshock) >= TOL) ? 2*twoArea/hshock : 0.0;
	
    for(i=0; i<3; i++)
	{
		// A1gradUx + A2gradUy
		AgradU[i] = A1[i][0]*gradUx[0] + A1[i][1]*gradUx[1] + A1[i][2]*gradUx[2] + 
			        A2[i][0]*gradUy[0] + A2[i][1]*gradUy[1] + A2[i][2]*gradUy[2];
	}

    deltashock1 = sqrt(invhRef *AgradU[0]*invhRef *AgradU[0] + 
                       invqxRef*AgradU[1]*invqxRef*AgradU[1] + 
                       invqyRef*AgradU[2]*invqyRef*AgradU[2])*
                  1.0/sqrt((invhRef *gradUx[0]*invhRef *gradUx[0] + invhRef *gradUy[0]*invhRef *gradUy[0] +
                            invqxRef*gradUx[1]*invqxRef*gradUx[1] + invqxRef*gradUy[1]*invqxRef*gradUy[1] +
                            invqyRef*gradUx[2]*invqyRef*gradUx[2] + invqyRef*gradUy[2]*invqyRef*gradUy[2]))*hshock/2.0;
    deltashock1 = (isnan(deltashock1) || isinf(deltashock1)) ? 0.0 : deltashock1;

    deltashock2 = sqrt(invhRef*AgradU[0]*invhRef*AgradU[0] + 
                       invqxRef*AgradU[1]*invqxRef*AgradU[1] + 
                       invqyRef*AgradU[2]*invqyRef*AgradU[2])*
                  1.0/sqrt(invhRef *Ub[0]*invhRef *Ub[0] + 
                           invqxRef*Ub[1]*invqxRef*Ub[1] + 
                           invqyRef*Ub[2]*invqyRef*Ub[2])*
                      hshock*hshock/4.0;
    deltashock2 = (isnan(deltashock2) || isinf(deltashock2)) ? 0.0 : deltashock2;

    delta = (deltashock1 + deltashock2)/2.0;

    return delta;
}
