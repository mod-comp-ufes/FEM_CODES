
#include "difconvrea.h" 

// Retorna o valor exato em um dado ponto
double DIFCONVREA_ExactSolution(double x, double y, double z, double Cst){
	
	double u;
		
	u = 100*x*y*z*(1-x)*(1-y)*(1-z);
		
	return u;
}

// Preenche o vetor u com as soluções em todos os pontos
void DIFCONVREA_ExactSolutionAllPoints(NodeType *Node, double nnodes, double *u, double A){
	
	int i;
	double x,y,z;
		
	for(i = 0; i < nnodes; i++){
		x = Node[i].x;
		y = Node[i].y;
		z = Node[i].z;
		
		u[i] = 100*x*y*z*(1-x)*(1-y)*(1-z);
		
	}
}
