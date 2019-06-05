#include "refletido.h"

int REFLETIDO_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u){
	int I, nnodes;

	nnodes = Parameters->nnodes;

	for(I = 0; I < nnodes; I++)
	{
		if (Node[I].id[0] >= 0){ u[Node[I].id[0]] = 1.0; }
		if (Node[I].id[1] >= 0){ u[Node[I].id[1]] = 2.9; }
		if (Node[I].id[2] >= 0){ u[Node[I].id[2]] = 0.0; }
		if (Node[I].id[3] >= 0){ u[Node[I].id[3]] = 5.990715; }
	}

	return 0;
}


