#include "NavierStokesEquations.h"

int Fill_ID(int *neq_out, int *neqpress_out, NodeType *Node, int nnodes)
{
	int I, neq = 0, neqpress = 0;
	
	for(I = 0; I < nnodes; I++){
		if (Node[I].v1Type == 1){
			Node[I].id[0] = neq++;
		}
		else
			Node[I].id[0] = -1;

		if (Node[I].v2Type == 1){
			Node[I].id[1] = neq++;
		}
		else
			Node[I].id[1] = -1;
		
		if (Node[I].v3Type == 1){
			Node[I].id[2] = neq++;
		}
		else
			Node[I].id[2] = -1;
	
		if (Node[I].pType == 1){
			Node[I].id[3] = neq++;
			neqpress++;
		}	
		else
			Node[I].id[3] = -1;

	}

	*neq_out = neq;
	*neqpress_out = neqpress;

	return 0;
}

