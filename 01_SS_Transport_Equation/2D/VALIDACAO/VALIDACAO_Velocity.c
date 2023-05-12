#include "validacao.h" 

void VALIDACAO_Velocity(double X, double Y, double Beta[2])
{
	double Be_x, Be_y;

	Be_x = 1;
	Be_y = 20*Y;
	Beta[0] = Be_x;
	Beta[1] = Be_y;	

}
