#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/Reordering/reordering.h"


void csr_List_insertA(NodeListType **CSR_List, int I, int J, int *nnzero_out)
{
	int nnzero;
	NodeListType *current, *previous, *new;

	nnzero = *nnzero_out;
	nnzero++;
	new = mycalloc("new of 'csr_List_insertA'", 1, sizeof(NodeListType));
	new->J = J;
	new->next = NULL;
	//only to initialization
	previous = new;
	previous->next = NULL;
	
	//printf("\n I = %d \t J = %d\n", I, J);

	if (CSR_List[I]==NULL){
		//printf("\n Entrei no CSR_List[I]==NULL: I = %d\n", I);
		CSR_List[I] = new;
	}else if (CSR_List[I]->J > J){
		//printf("\n Entrei no CSR_List[I]->J > J: I = %d e J = %d\n", I, J);
		new->next = CSR_List[I];
		CSR_List[I] = new;
	}
	else {
		//printf("\n Entrei no else: I = %d\n", I);
		current = CSR_List[I];		
		while (current != NULL){
			//printf("\n Entrei no While %d\n", ++i);
			if (J <= current->J)			
				break;
			previous = current;	
			current = current->next;
		}
		if (current == NULL){
			//printf("\n Entrei no current == NULL\n");
			previous->next = new;
		}else if (current->J == J){
			//printf("\n Entrei no current->J == J: J = %d\n", J);
			nnzero--;
			free(new);
		}
		else {
			//printf("\n Entrei no else do else\n");
			new->next = current;
			previous->next = new;
		}
	}
	*nnzero_out = nnzero;
	
/*	printf("\nCSR_List[%d]->J = %d\n", I, CSR_List[I]->J);
	printf("CSR_List[%d]->next = %p\n", I, CSR_List[I]->next);
	for (NodeListType *p = CSR_List[I]; p != NULL; p = p->next)
		printf("p->J = %d \t p->next = %p\n", p->J, p->next);
	getchar(); */
}

int csr_search(int I, int J, NodeListType *current)
{
	int position;

	position = -1;

	while (current != NULL){
		position++;
		if (current->J == J)
			break;	 
		current = current->next;

	}

	return position;
}

void csr_Initialization(ParametersType *Parameters, NodeType *Node, int **JA_out, int **IA_out, int **perm_out, int  **invperm_out,
			int ***lm_out, int **lmaux_out, int ***CSR_by_Element_out)
{
	int neq, nel, nnodes, nnzero, count, pos1, pos2, pos3, pos4,  K, I, I1, I2, I3, I4;
	int **CSR_by_Element, **lm, *lmaux, *IA, *JA;
	NodeListType *current, *temp, **CSR_List;

	nel  = Parameters->nel;
	neq  = Parameters->neq;
	nnodes  = Parameters->nnodes;
	lm = *lm_out;
	lmaux = *lmaux_out;

	CSR_List = (NodeListType **)mycalloc("CSR_List of 'CSR_Initialization'",neq,sizeof(NodeListType *));

	for (K = 0; K<neq ; K++)
		CSR_List[K] = NULL;	

	nnzero = 0;
	for (K = 0; K < nel ; K++)
	{
		I1 = lm[K][0];
		I2 = lm[K][1];
		I3 = lm[K][2];
		I4 = lm[K][3];
		
		if (I1 != neq){
			csr_List_insertA (CSR_List, I1, I1, &nnzero);
			if (I2 != neq)
				csr_List_insertA (CSR_List, I1, I2, &nnzero);
			if (I3 != neq)
				csr_List_insertA (CSR_List, I1, I3, &nnzero);
			if (I4 != neq)
				csr_List_insertA (CSR_List, I1, I4, &nnzero);
		}

		if (I2 != neq){
			csr_List_insertA (CSR_List, I2, I2, &nnzero);
			if (I1 != neq)
				csr_List_insertA (CSR_List, I2, I1, &nnzero);
			if (I3 != neq)
				csr_List_insertA (CSR_List, I2, I3, &nnzero);
			if (I4 != neq)
				csr_List_insertA (CSR_List, I2, I4, &nnzero);
		}
		
		if (I3 != neq){
			csr_List_insertA (CSR_List, I3, I3, &nnzero);
			if (I1 != neq)
				csr_List_insertA (CSR_List, I3, I1, &nnzero);
			if (I2 != neq)
				csr_List_insertA (CSR_List, I3, I2, &nnzero);
			if (I4 != neq)
				csr_List_insertA (CSR_List, I3, I4, &nnzero);
		}
		
		if (I4 != neq){
			csr_List_insertA (CSR_List, I4, I4, &nnzero);
			if (I1 != neq)
				csr_List_insertA (CSR_List, I4, I1, &nnzero);
			if (I2 != neq)
				csr_List_insertA (CSR_List, I4, I2, &nnzero);
			if (I3 != neq)
				csr_List_insertA (CSR_List, I4, I3, &nnzero);
		}

	}

	Parameters->nnzero = nnzero;
	JA = (int *)mycalloc("JA of 'CSR_Initialization'",nnzero, sizeof(int));
	IA = (int *)mycalloc("IA of 'CSR_Initialization'",neq+1, sizeof(int));
	count = 0;
	for (K = 0 ; K < neq ; K++){
		current = CSR_List[K];
		while (current != NULL){
			JA[count++] = current->J;	
			IA[K+1]++;
			current = current->next;
		}	
	}	 

	for (K = 1; K < neq; K++)
		IA[K+1] += IA[K];	
	
	CSR_by_Element = (int**) mycalloc("CSR_by_Element of 'CSR_Initialization'", nel, sizeof(int*));
	for (K = 0; K < nel; K++)
		CSR_by_Element[K] = (int*) mycalloc("CSR_by_Element of 'CSR_Initialization'", 16, sizeof(int));
		

	for (K = 0; K < nel ; K++)
		for (I = 0; I < 16; I++)
			CSR_by_Element[K][I] = nnzero;

	for (K = 0; K < nel ; K++){
		I1 = lm[K][0];
		I2 = lm[K][1];
		I3 = lm[K][2];
		I4 = lm[K][3];
		
		//printf("I1 = %d \t I2 = %d \t I3 = %d \t I4 = %d\n", I1, I2, I3, I4); getchar();
		
		if (I1 != neq){
			pos1 = csr_search( I1, I1, CSR_List[I1]);
			CSR_by_Element[K][0] = IA[I1] + pos1;
			if (I2 != neq){
				pos2 = csr_search( I1, I2, CSR_List[I1]);
				CSR_by_Element[K][1] = IA[I1] + pos2;
			}
			if (I3 != neq){
				pos3 = csr_search( I1, I3, CSR_List[I1]);
				CSR_by_Element[K][2] = IA[I1] + pos3;
				
			}
			if (I4 != neq){
				pos4 = csr_search( I1, I4, CSR_List[I1]);
				CSR_by_Element[K][3] = IA[I1] + pos4;
				
			}
		}
		
		if (I2 != neq){
			pos2 = csr_search( I2, I2, CSR_List[I2]);
			CSR_by_Element[K][5] = IA[I2] + pos2;
			if (I1 != neq){
				pos1 = csr_search( I2, I1, CSR_List[I2]);
				CSR_by_Element[K][4] = IA[I2] + pos1;
			}
			if (I3 != neq){
				pos3 = csr_search( I2, I3, CSR_List[I2]);
				CSR_by_Element[K][6] = IA[I2] + pos3;
				
			}
			if (I4 != neq){
				pos4 = csr_search( I2, I4, CSR_List[I2]);
				CSR_by_Element[K][7] = IA[I2] + pos4;
				
			}
		}
		
		if (I3 != neq){
			pos3 = csr_search( I3, I3, CSR_List[I3]);
			CSR_by_Element[K][10] = IA[I3] + pos3;
			if (I1 != neq){
				pos1 = csr_search( I3, I1, CSR_List[I3]);
				CSR_by_Element[K][8] = IA[I3] + pos1;
			}
			if (I2 != neq){
				pos2 = csr_search( I3, I2, CSR_List[I3]);
				CSR_by_Element[K][9] = IA[I3] + pos2;
				
			}
			if (I4 != neq){
				pos4 = csr_search( I3, I4, CSR_List[I3]);
				CSR_by_Element[K][11] = IA[I3] + pos4;
				
			}
		}

		if (I4 != neq){
			pos4 = csr_search( I4, I4, CSR_List[I4]);
			//printf("pos4 = %d\n",pos4);
			CSR_by_Element[K][15] = IA[I4] + pos4;
			if (I1 != neq){
				pos1 = csr_search( I4, I1, CSR_List[I4]);
				CSR_by_Element[K][12] = IA[I4] + pos1;
			}
			if (I2 != neq){
				pos2 = csr_search( I4, I2, CSR_List[I4]);
				CSR_by_Element[K][13] = IA[I4] + pos2;
				
			}
			if (I3 != neq){
				pos3 = csr_search( I4, I3, CSR_List[I4]);
				CSR_by_Element[K][14] = IA[I4] + pos3;
				
			}
		}
		
		//printf("valores do CSR by element do elemento %d\n", K);
		//printf("%d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \t %d \n", CSR_by_Element[K][0],
		  //      CSR_by_Element[K][1], CSR_by_Element[K][2], CSR_by_Element[K][3], CSR_by_Element[K][4], CSR_by_Element[K][5], CSR_by_Element[K][6],
		  //      CSR_by_Element[K][7], CSR_by_Element[K][8], CSR_by_Element[K][9], CSR_by_Element[K][10], CSR_by_Element[K][11], CSR_by_Element[K][12],
		  //      CSR_by_Element[K][13], CSR_by_Element[K][14], CSR_by_Element[K][15]);
		//getchar();
		
	}
//	free(lm);
//	free(lmaux);

//	lm = NULL;
//	lmaux = NULL;
 	
	
	for (K = 0; K < neq; K++){
		current = CSR_List[K];
		while (current != NULL){
			temp = current;		
			current	= current->next;
			free(temp);
		}
	}

	free(CSR_List);

	int *PermCSR = mycalloc("PermCSR of 'csr_Inititalization'",nnzero,sizeof(int));
	int *perm = mycalloc("perm of 'csr_Initialization'",neq,sizeof(int));

	for (I=0; I<nnzero; I++)
		PermCSR[I] = I;
	for (I=0; I<neq; I++)
		perm[I]=I;

	reordering (Parameters,JA,IA,perm,PermCSR);

	int *invPermCSR = mycalloc("invPermCSR of 'csr_Inititalization'",nnzero,sizeof(int));
	int *invperm = mycalloc("invPermCSR of 'csr_Inititalization'",neq,sizeof(int));

	for (I=0;I<nnzero; I++){
		invPermCSR[PermCSR[I]] = I;
	}


	for (I=0;I<neq;I++)
		invperm[perm[I]]=I;

	int aux[16];
	for (K=0; K<nel; K++){
		for (I=0;I<16;I++)
			if (CSR_by_Element[K][I]<nnzero)
				aux[I] = invPermCSR[CSR_by_Element[K][I]];
			else
				aux[I] = nnzero;	
		for (I=0;I<16;I++)
			CSR_by_Element[K][I] = aux[I];
	}

	for (K=0; K<nnodes; K++){
		if (Node[K].id!=-1)	
			Node[K].id = invperm[Node[K].id];
	}

	free(perm);	
	free(invperm);	
	free(PermCSR);
	free(invPermCSR);
	*lm_out = lm;
	*lmaux_out = lmaux;
	*CSR_by_Element_out = CSR_by_Element;
	*JA_out = JA; 
	*IA_out = IA;
	*perm_out = perm;
	*invperm_out = invperm;
}


