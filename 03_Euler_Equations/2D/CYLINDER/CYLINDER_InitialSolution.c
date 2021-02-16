# include "cylinder.h"

int CYLINDER_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u){

	int I, nnodes;
	double theta;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++)
	{
		if (Node[I].id[0] >= 0){ 
			u[Node[I].id[0]] = 1.0; 
		}
		if (Node[I].id[1] >= 0){ 
			if (Node[I].v1Type==-1){
				theta = CYLINDER_theta(Node[I].x,Node[I].y);			
				u[Node[I].id[1]] = cos(theta); //ajustando a condição inicial para a malha do cilindro
			}
			else	
				u[Node[I].id[1]] = 1.0; 
		}
		if (Node[I].id[2] >= 0){ 
			u[Node[I].id[2]] = 0.0; 
		}
		if (Node[I].id[3] >= 0){ 
			u[Node[I].id[3]] = 0.9464; 
		}
	}

	return 0;
}


