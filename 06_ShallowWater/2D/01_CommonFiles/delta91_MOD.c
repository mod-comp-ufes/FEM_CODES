#include <math.h>
#include "ShalowWater.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"


double delta91_MOD(double* Ue, double *gradUx, double *gradUy, double *Re, double g)
{
	double delta = 0.0, delta_91, delta_tau;
	double normR;
	double h[3], u[3], v[3];

	h[0] = Ue[0];
	h[1] = Ue[3];
	h[2] = Ue[6];

	u[0] = Ue[1]/Ue[0];
	u[1] = Ue[4]/Ue[3];
	u[2] = Ue[7]/Ue[6];

	v[0] = Ue[2]/Ue[0];
	v[1] = Ue[5]/Ue[3];
	v[2] = Ue[8]/Ue[6];

	normR = sqrt(1.0/h[0]*(Re[0]*(Re[0]*(g*h[0] + u[0]*u[0] + v[0]*v[0]) - Re[1]*Re[1]*u[0] - Re[2]*Re[2]*v[0]) +
			               Re[1]*(-Re[0]*u[0] + Re[1]*Re[1]) +
			               Re[2]*(-Re[0]*v[0] + Re[2]*Re[2])) +

			     1.0/h[1]*(Re[3]*(Re[3]*(g*h[1] + u[1]*u[1] + v[1]*v[1]) - Re[4]*Re[4]*u[1] - Re[5]*Re[5]*v[1]) +
			               Re[4]*(-Re[3]*u[1] + Re[4]*Re[4]) +
			               Re[5]*(-Re[3]*v[1] + Re[5]*Re[5])) +

			     1.0/h[2]*(Re[6]*(Re[6]*(g*h[2] + u[2]*u[2] + v[2]*v[2]) - Re[7]*Re[7]*u[2] - Re[8]*Re[8]*v[2]) +
			               Re[7]*(-Re[6]*u[2] + Re[7]*Re[7]) +
			               Re[8]*(-Re[6]*v[2] + Re[8]*Re[8])));


	delta_tau = 1.0/h[0]*(gradUx[0]*(gradUx[0]*(g*h[0] + u[0]*u[0] + v[0]*v[0]) - gradUx[1]*gradUx[1]*u[0] - gradUx[2]*gradUx[2]*v[0]) +
			              gradUx[1]*(-gradUx[0]*u[0] + gradUx[1]*gradUx[1]) +
			              gradUx[2]*(-gradUx[0]*v[0] + gradUx[2]*gradUx[2])) +

			    1.0/h[1]*(gradUx[3]*(gradUx[3]*(g*h[1] + u[1]*u[1] + v[1]*v[1]) - gradUx[4]*gradUx[4]*u[1] - gradUx[5]*gradUx[5]*v[1]) +
			              gradUx[4]*(-gradUx[3]*u[1] + gradUx[4]*gradUx[4]) +
			              gradUx[5]*(-gradUx[3]*v[1] + gradUx[5]*gradUx[5])) +

			    1.0/h[2]*(gradUx[6]*(gradUx[6]*(g*h[2] + u[2]*u[2] + v[2]*v[2]) - gradUx[7]*gradUx[7]*u[2] - gradUx[8]*gradUx[8]*v[2]) +
			              gradUx[7]*(-gradUx[6]*u[2] + gradUx[7]*gradUx[7]) +
			              gradUx[8]*(-gradUx[6]*v[2] + gradUx[8]*gradUx[8])) +

				1.0/h[0]*(gradUy[0]*(gradUy[0]*(g*h[0] + u[0]*u[0] + v[0]*v[0]) - gradUy[1]*gradUy[1]*u[0] - gradUy[2]*gradUy[2]*v[0]) +
			              gradUy[1]*(-gradUy[0]*u[0] + gradUy[1]*gradUy[1]) +
			              gradUy[2]*(-gradUy[0]*v[0] + gradUy[2]*gradUy[2])) +

			    1.0/h[1]*(gradUy[3]*(gradUy[3]*(g*h[1] + u[1]*u[1] + v[1]*v[1]) - gradUy[4]*gradUy[4]*u[1] - gradUy[5]*gradUy[5]*v[1]) +
			              gradUy[4]*(-gradUy[3]*u[1] + gradUy[4]*gradUy[4]) +
			              gradUy[5]*(-gradUy[3]*v[1] + gradUy[5]*gradUy[5])) +

			    1.0/h[2]*(gradUy[6]*(gradUy[6]*(g*h[2] + u[2]*u[2] + v[2]*v[2]) - gradUy[7]*gradUy[7]*u[2] - gradUy[8]*gradUy[8]*v[2]) +
			              gradUy[7]*(-gradUy[6]*u[2] + gradUy[7]*gradUy[7]) +
			              gradUy[8]*(-gradUy[6]*v[2] + gradUy[8]*gradUy[8]));
	delta_tau = normR/delta_tau;

	delta = max(0, delta_91 - delta_tau);

	return delta;
}
