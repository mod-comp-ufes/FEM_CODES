#include "difconv.h" 
#include <math.h>

// Retorna o valor exato em um dado ponto
double DIFCONV_ExactSolution(double x, double y, double z, double A){
	
	double u;
		
	u = sin(PI*x)*sin(PI*y)*sin(PI*z);
		
	return u;
}

// Preenche o vetor u com as soluções em todos os pontos
void DIFCONV_ExactSolutionAllPoints(NodeType *Node, double nnodes, double *u, double A){
	
	int i;
	double x,y,z;
		
	for(i = 0; i < nnodes; i++){
		x = Node[i].x;
		y = Node[i].y;
		z = Node[i].z;
		
		u[i] = sin(PI*x)*sin(PI*y)*sin(PI*z);
		
	}
}
