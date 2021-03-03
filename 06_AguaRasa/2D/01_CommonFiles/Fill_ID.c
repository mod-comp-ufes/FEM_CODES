#include "ShalowWater.h"


int Fill_ID(int *neq_out, NodeType *Node, int nnodes)
{
	int I, neq = 0;
	
	for(I = 0; I < nnodes; I++)
	{
		if (Node[I].hType != 0)
			Node[I].id[0] = neq++;
		else
			Node[I].id[0] = -1;

		if (Node[I].qxType != 0)
			Node[I].id[1] = neq++;
		else
			Node[I].id[1] = -1;

		if (Node[I].qyType != 0)
			Node[I].id[2] = neq++;
		else
			Node[I].id[2] = -1;
	}

	*neq_out = neq;

	return 0;
}
