# include "exata.h"

double EXATA_InitialSolution(ParametersType *Parameters, NodeType *Node, double *u){
	int I, nnodes;
	double x, y;

	nnodes = Parameters->nnodes;
	
	for(I = 0; I < nnodes; I++)
	{
		x = Node[I].x;
		y = Node[I].y;
				
		if (Node[I].id[0] >= 0){ u[Node[I].id[0]] = pow(x,2)*pow(1-x,2)*(2*y-6*pow(y,2)+4*pow(y,3)); };
		if (Node[I].id[1] >= 0){ u[Node[I].id[1]] = pow(y,2)*pow(1-y,2)*(-2*x+6*pow(x,2)-4*pow(x,3)); };
		if (Node[I].id[2] >= 0){ u[Node[I].id[2]] = pow(x,2) - pow(y,2); };
		//printf("\n Pasei aqui SOLUCAO INICIAL x = %f, y= %f, u= %f, v=%f, p=%f . \n", x, y, u[Node[I].id[0]], u[Node[I].id[1]], u[Node[I].id[2]]);
			
	}
	return 0;
}


