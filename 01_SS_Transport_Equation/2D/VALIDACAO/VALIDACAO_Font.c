#include "validacao.h"
#include <math.h>

double VALIDACAO_Font(double X, double Y, double k, double gamma, double Be_x, double Be_y)
{
	double f, u, u_x, u_y, u_x2, u_y2;

	u = 10*X*Y*(1-X)*(1-Y)*exp(pow(X, 4.5));
  	u_x = 10*Y*(1-Y)*(exp(pow(X, 4.5)) - 2*exp(pow(X, 4.5))*X + 4.5*exp(pow(X, 4.5))*pow(X, 4.5) - 4.5*exp(pow(X, 4.5))*pow(X, 5.5));
  	u_x2 = 10*Y*(1-Y)*(-20.25*exp(pow(X, 4.5))*pow(X, 9) + 20.25*exp(pow(X, 4.5))*pow(X, 8) - 33.75*exp(pow(X, 4.5))*pow(X, 4.5) + 24.75*exp(pow(X, 4.5))*pow(X, 3.5) - 2*exp(pow(X, 4.5)));
  	u_y = 10*exp(pow(X, 4.5))*X*(1-X)*(1-2*Y);
  	u_y2 = -20*exp(pow(X, 4.5))*X*(1-X);

  	f = -k*(u_x2 + u_y2) + Be_x*u_x + Be_y*u_y + gamma*u;

	return f;

}
