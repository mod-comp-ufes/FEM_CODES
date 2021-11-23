#include "SSTransportEquation3D.h"

// Node[I].uType armazena 1 se para esse nó u é incognita.
int Fill_ID(int *neq_out, NodeType *Node, int nnodes)
{
	int I, neq = 0;

	for(I = 0; I < nnodes; I++){
		if (Node[I].uType == 1){
			Node[I].id = neq;
			neq = neq + 1;
		}
		else
			Node[I].id = -1;

	}

	*neq_out = neq;

	return 0;
}

