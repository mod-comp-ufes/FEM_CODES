#include "EulerEquations.h"

int Fill_ID(int *neq_out, int *neqrho_out, NodeType *Node, int nnodes)
{
	int I, neq = 0, neqrho = 0;
	
	for(I = 0; I < nnodes; I++){
		if (Node[I].rhoType != 0){
			Node[I].id[0] = neq++;
			neqrho++;	
		}
		else
			Node[I].id[0] = -1;

		if (Node[I].v1Type != 0){
			Node[I].id[1] = neq++;
		}
		else
			Node[I].id[1] = -1;

		if (Node[I].v2Type != 0){
			Node[I].id[2] = neq++;
		}
		else
			Node[I].id[2] = -1;
	
		if (Node[I].eType != 0){
			Node[I].id[3] = neq++;
		}	
		else
			Node[I].id[3] = -1;

	}

	*neq_out = neq;
	*neqrho_out = neqrho;

	return 0;
}

