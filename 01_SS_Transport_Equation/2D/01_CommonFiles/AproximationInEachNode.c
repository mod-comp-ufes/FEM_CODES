#include "SSTranspEquation.h"

/*----- vectors containing the value of the solution (prescribed or calculated) at each node -----*/
int AproximationInEachNode(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *u_out){

	int nnodes = Parameters->nnodes;
	int i, eq;
	double X, Y;
	double *u = FemStructs->u;
	NodeType *Node = FemStructs->Node;

	for (i = 0; i < nnodes; i++){
		eq = Node[i].id;
		
		X = Node[i].x;
		Y = Node[i].y;
		
		if (eq >= 0)
			u_out[i] = u[eq];
		else
			u_out[i] = FemFunctions->upresc(X, Y);
	}

	return 0;
}
