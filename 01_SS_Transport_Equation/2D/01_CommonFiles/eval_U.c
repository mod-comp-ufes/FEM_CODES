#include "SSTranspEquation.h"

void eval_U(ParametersType *Parameters,FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *U)
{
	int I;
	double x, y;
	int nnodes = Parameters->nnodes;
	double *u = FemStructs->u;
	NodeType *Node = FemStructs->Node;


	for (I=0; I<nnodes; I++){

		x = Node[I].x;
		y = Node[I].y;
	
		if (Node[I].id >= 0)
			U[I] = u[Node[I].id];
		else
			U[I] = FemFunctions->upresc(x,y);
		
	}

}
