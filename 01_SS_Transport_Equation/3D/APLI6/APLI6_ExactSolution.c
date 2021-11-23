#include "apli6.h" 
#include <math.h>

double APLI6_ExactSolution(double x, double y, double z, double A){
	
	double u;
		
	u = exp(-A*x) + exp(-A*y) + exp(-A*z);
	
	return u;	
	
}

void APLI6_ExactSolutionAllPoints(NodeType *Node, double nnodes, double *u, double A){
	
	int i;
	double x,y,z;
		
	for(i = 0; i < nnodes; i++){
		x = Node[i].x;
		y = Node[i].y;
		z = Node[i].z;
		
		u[i] = exp(-A*x) + exp(-A*y) + exp(-A*z);
		
	}
}

double APLI6_DuDx(double x, double y, double z, double A){
	
	return (-A*exp(-A*x));

}

double APLI6_DuDy(double x, double y, double z, double A){
	
	return (-A*exp(-A*y));

}

double APLI6_DuDz(double x, double y, double z, double A){
	
	return (-A*exp(-A*z));

}
