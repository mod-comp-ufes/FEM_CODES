# include "naca0012.h"

int NACA0012_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u){

	int I, nnodes;
	double theta;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++)
	{
		if (Node[I].id[0] >= 0){ 
			u[Node[I].id[0]] = 1.0; 
		}
		if (Node[I].id[1] >= 0){ 
			if (Node[I].v1Type<0){
				theta = NACA0012_theta(Node[I].x, Node[I].y);
				u[Node[I].id[1]] = cos(theta); //ajustando a condição inicial para a malha do naca			
			}
			else	
				u[Node[I].id[1]] = 1.0;
		}
		if (Node[I].id[2] >= 0){ 
			u[Node[I].id[2]] = 0.0; 
		}
		if (Node[I].id[3] >= 0){ 
			u[Node[I].id[3]] = 0.9464;  // 0.9464 para Mach=2 // e = 7.64286; para definir mach = 0.5
		}
	}
	return 0;
}


