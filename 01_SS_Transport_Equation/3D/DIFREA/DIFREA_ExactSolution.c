#include "difrea.h" 
#include <math.h>

double DIFREA_ExactSolution(double x, double y, double z, double Cst){
	
	double u;
		
	u = sin(PI*x)*sin(PI*y)*sin(PI*z);
		
	return u;
}


void DIFREA_ExactSolutionAllPoints(NodeType *Node, double nnodes, double *u, double A){
	
	int i;
	double x,y,z;
		
	for(i = 0; i < nnodes; i++){
		x = Node[i].x;
		y = Node[i].y;
		z = Node[i].z;
		
		u[i] = sin(PI*x)*sin(PI*y)*sin(PI*z);
		
	}
}
