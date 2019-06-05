#include "NavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/Reordering/reordering.h"

void csr_Initialization(ParametersType *Parameters, NodeType *Node, int **JA_out, int **IA_out, int **perm_out, int  **invperm_out,
			int ***lm_out, int **lmaux_out, int ***CSR_by_Element_out)
{
	int neq, nel, nnodes, nnzero, count, pos,  K, I, J, II, JJ;
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
		for (I=0;I<9;I++){
			II = lm[K][I];
			if (II != neq){
				csr_List_insertA( CSR_List, II,II, &nnzero);
				for (J=0; J<9; J++){
					JJ = lm[K][J];
					if (JJ != II && JJ != neq)
						csr_List_insertA( CSR_List, II,JJ, &nnzero);
				}
			}
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
	
	
	CSR_by_Element = (int**) mycalloc("CSR_by_Element of 'CSR_Initialization'", nel,sizeof(int*));
	for (K = 0; K < nel; K++)
		CSR_by_Element[K] = (int*) mycalloc("CSR_by_Element of 'CSR_Initialization'",81,sizeof(int));
		

	for (K = 0; K < nel ; K++)
		for (I = 0; I < 81; I++)
			CSR_by_Element[K][I] = nnzero;

	for (K = 0; K < nel ; K++){
		for (I=0; I<9; I++){
			II = lm[K][I];
			if (II != neq){
				for (J=0; J<9; J++){	
					JJ = lm[K][J];					
					if (JJ != neq){
						pos = csr_search(II, JJ, CSR_List[II]);
						CSR_by_Element[K][9*I+J] = IA[II] + pos;
					}

				}
			}
		}

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

	int aux[81];
	for (K=0; K<nel; K++){
		for (I=0;I<81;I++)
			if (CSR_by_Element[K][I]<nnzero)
				aux[I] = invPermCSR[CSR_by_Element[K][I]];
			else
				aux[I] = nnzero;	
		for (I=0;I<81;I++)
			CSR_by_Element[K][I] = aux[I];
	}

	for (K=0; K<nnodes; K++){
		for (I=0; I<3; I++){
			if (Node[K].id[I]!=-1)	
				Node[K].id[I] = invperm[Node[K].id[I]];
		}
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


void csr_List_insertA(NodeListType **CSR_List, int I, int J, int *nnzero_out)
{
	int nnzero;
	NodeListType *current, *previous, *new;

	nnzero = *nnzero_out;
	nnzero++;
	new = mycalloc("new of 'CSR_Initialization'", 1, sizeof(NodeListType));
	new->J = J;
	new->next = NULL;
	//only to initialization
	previous = new;
	previous->next = NULL;

	if (CSR_List[I]==NULL)
		CSR_List[I] = new;
	else if (CSR_List[I]->J > J){
		new->next = CSR_List[I];
		CSR_List[I] = new;
	}
	else {
		current = CSR_List[I];		
		while (current != NULL){
			if (J <= current->J)			
				break;
			previous = current;	
			current = current->next;
		}
		if (current == NULL)
			previous->next = new;
		else if (current->J == J){
			nnzero--;
			free(new);
		}
		else {
			new->next = current;
			previous->next = new;
		}
	}
	*nnzero_out = nnzero;
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



