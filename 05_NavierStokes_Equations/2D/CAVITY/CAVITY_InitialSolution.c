# include "cavity.h"

double CAVITY_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u){
	int I, nnodes;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++)
	{
		if (Node[I].id[0] >= 0){ u[Node[I].id[0]] = 0.0; };
		if (Node[I].id[1] >= 0){ u[Node[I].id[1]] = 0.0; }; 
		if (Node[I].id[2] >= 0){ u[Node[I].id[2]] = 0.0; };
	}

	return 0;
}


