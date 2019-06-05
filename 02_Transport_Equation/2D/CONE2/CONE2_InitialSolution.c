#include "cone2.h"

int CONE2_InitialSolution(ParametersType *Parameters, NodeType* Node, double *u)
{
	int I, nnodes;
	double X, Y;
	double R;

	nnodes = Parameters->nnodes;

	for(I=0;I<nnodes;I++){
		if (Node[I].id>=0){
			X = Node[I].x;
			Y = Node[I].y;
			R = sqrt((X+2.5)*(X+2.5)+ Y*Y);
			if (R<=1.0)
				u[Node[I].id] = 7.5*(1.0 + cos(R*PI));
			else
				u[Node[I].id] = 0.0;
				
		}
	}

	return 0;
}


