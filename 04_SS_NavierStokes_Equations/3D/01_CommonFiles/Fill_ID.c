#include "SS_NavierStokesEquations3D.h"

// Node[I].v1Type armazena 1 se para esse nó v na direção x é incognita. O mesmo se aplica às demais grandezas.
// As equações são distribuidas primeiro para velocidade na direção x, depois y e z e depois para a pressão.
int Fill_ID(int *neq_out, int *neqv1_out, int *neqv2_out, int *neqv3_out, int *neqpress_out, NodeType *Node, int nnodes)
{
	int I, neq = 0, neqpress = 0, neqv1 = 0, neqv2 = 0, neqv3 = 0;

	for(I = 0; I < nnodes; I++){
		if (Node[I].v1Type == 1){
			Node[I].id[0] = neq++;
			neqv1++;
		}
		else
			Node[I].id[0] = -1;
		
		if (Node[I].v2Type == 1){
			Node[I].id[1] = neq++;
			neqv2++;
		}
		else
			Node[I].id[1] = -1;	
			
		if (Node[I].v3Type == 1){
			Node[I].id[2] = neq++;
			neqv3++;
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


/*	for(I = 0; I < nnodes; I++){
		if (Node[I].v1Type == 1){
			Node[I].id[0] = neq++;
			neqv1++;
		}
		else
			Node[I].id[0] = -1;
	}

	for(I = 0; I < nnodes; I++){
		if (Node[I].v2Type == 1){
			Node[I].id[1] = neq++;
			neqv2++;
		}
		else
			Node[I].id[1] = -1;
	}

	for(I = 0; I < nnodes; I++){
		if (Node[I].v3Type == 1){
			Node[I].id[2] = neq++;
			neqv3++;
		}
		else
			Node[I].id[2] = -1;
	}

	for(I = 0; I < nnodes; I++){
		if (Node[I].pType == 1){
			Node[I].id[3] = neq++;
			neqpress++;
		}	
		else
			Node[I].id[3] = -1;
	} */

	*neq_out = neq;
	*neqv1_out = neqv1;
	*neqv2_out = neqv2;
	*neqv3_out = neqv3;
	*neqpress_out = neqpress;

	return 0;
}

