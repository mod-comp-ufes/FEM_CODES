#include "EulerEquations.h"

void rotation(int tag, double theta, double M[12][12], double F[12])
{
	int I, J, K;
	double cs = cos(theta);
	double sn = sin(theta);
	double b,c,e,f,g,h,i,j,k,l,n,o;

	for (K=0; K<3; K++){
		I=tag*NDOF;
		J=K*NDOF;
		b = M[I][J+1];
		c = M[I][J+2];
		e = M[I+1][J];
		f = M[I+1][J+1];
		g = M[I+1][J+2];
		h = M[I+1][J+3];
		i = M[I+2][J];
		j = M[I+2][J+1];
		k = M[I+2][J+2];
		l = M[I+2][J+3];
		n = M[I+3][J+1];
		o = M[I+3][J+2];
		M[I][J+1] = b*cs - c*sn;
		M[I][J+2] = b*sn + c*cs;
		M[I+1][J] = e*cs - i*sn;
		M[I+1][J+1] = cs*(f*cs - j*sn) - sn*(g*cs - k*sn);
		M[I+1][J+2] = cs*(g*cs - k*sn) + sn*(f*cs - j*sn);
		M[I+1][J+3] = h*cs - l*sn;
		M[I+2][J] = e*sn + i*cs;
		M[I+2][J+1] = cs*(f*sn + j*cs) - sn*(g*sn + k*cs);
		M[I+2][J+2] = cs*(g*sn + k*cs) + sn*(f*sn + j*cs);
		M[I+2][J+3] = h*sn + l*cs;
		M[I+3][J+1] = n*cs - o*sn;
		M[I+3][J+2] = n*sn + o*cs;
	}
	
	b = F[I+1];
	c = F[I+2];
	F[I+1] = b*cs - c*sn;
	F[I+2] = b*sn + c*cs;

/* Pt*Ke*P =
[        a,                     b*x - c*y,                     b*y + c*x,         d],
[e*x - i*y, x*(f*x - j*y) - y*(g*x - k*y), x*(g*x - k*y) + y*(f*x - j*y), h*x - l*y],
[e*y + i*x, x*(f*y + j*x) - y*(g*y + k*x), x*(g*y + k*x) + y*(f*y + j*x), h*y + l*x],
[        m,                     n*x - o*y,                     n*y + o*x,         p]])


Pt*F =
[        a],
[b*x - c*y],
[b*y + c*x],
[        d]])*/


	return;	

}




