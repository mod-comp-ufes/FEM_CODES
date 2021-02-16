#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

void ede_Initialization(ParametersType *Parameters, int **order, int ***lm_out, int **lmaux_out, int ***EDGE_by_Element_out)
{
	int neq, nel, nedge, K, J, I1, I2, I3, count;
	int **EDGE_by_Element, **lm, *lmaux, *Indicator;
	NodeListType *current, *temp, **EDGE_List;

	nel = Parameters->nel;
	neq = Parameters->neq;

	lm = *lm_out;
	lmaux = *lmaux_out;

	EDGE_List = (NodeListType **)mycalloc("EDGE_List of 'EDE_Inicialization'",neq+1,sizeof(NodeListType *));

	for (K = 0; K < neq+1 ; K++)
		EDGE_List[K] = NULL;	

	nedge = 0;
	for (K = 0; K < nel ; K++)
	{
		if (lm[K][0]<lm[K][1] && lm[K][0]<lm[K][2]){
			I1 = lm[K][0];
			if (lm[K][1] < lm[K][2]){
				I2 = lm[K][1];	
				I3 = lm[K][2];
				order[K][0] = 0;
				order[K][1] = 1;
				order[K][2] = 2;

			}else{
				I2 = lm[K][2];	
				I3 = lm[K][1];
				order[K][0] = 0;
				order[K][1] = 2;
				order[K][2] = 1;
			}
		}
		else if (lm[K][1] < lm[K][2]){
			I1 = lm[K][1];
			if (lm[K][0] < lm[K][2]){
				I2 = lm[K][0];	
				I3 = lm[K][2];
				order[K][0] = 1;
				order[K][1] = 0;
				order[K][2] = 2;
		
			}else{ 
				I2 = lm[K][2];	
				I3 = lm[K][0];
				order[K][0] = 1;
				order[K][1] = 2;
				order[K][2] = 0;

			}
		}
		else{
			I1 = lm[K][2];
			if (lm[K][0] < lm[K][1]){
				I2 = lm[K][0];
				I3 = lm[K][1];
				order[K][0] = 2;
				order[K][1] = 0;
				order[K][2] = 1;

			}else{

				I2 = lm[K][1];
				I3 = lm[K][0];
				order[K][0] = 2;
				order[K][1] = 1;
				order[K][2] = 0;

			}

		}			
		ede_List_insertA( EDGE_List, I1,I2, &nedge);
		ede_List_insertA( EDGE_List, I1,I3, &nedge);
		ede_List_insertA( EDGE_List, I2,I3, &nedge);

	}
 	
	Indicator = (int *)mycalloc("Incator of 'EDE_Inicialization'",neq+1, sizeof(int));

	for (K = 0 ; K < neq ; K++){
		current = EDGE_List[K];
		while (current != NULL){
			Indicator[K+1]++;	
			current = current->next;
		}	
	}	 

	for (K = 1; K < neq; K++)
		Indicator[K+1] += Indicator[K];	

	EDGE_by_Element = (int**) mycalloc("EDGE_by_Element of 'EDE_Inicialization'",nel,sizeof(int*));
	for (K = 0; K < nel; K++)
		EDGE_by_Element[K] = (int*) mycalloc("EDGE_by_Element of 'EDE_Inicialization'",3,sizeof(int));
		

	nedge--;
	for (K = 0; K < nel ; K++){
		I1 = lm[K][order[K][0]];
		I2 = lm[K][order[K][1]];
		I3 = lm[K][order[K][2]];

		current = EDGE_List[I1];
		count = 0;
		while (current->J != I2){
			count++;
			current = current->next;
		}
		EDGE_by_Element[K][0] = Indicator[I1] + count;
		
		current = EDGE_List[I1];
		count = 0;
		while (current->J != I3){
			count++;
			current = current->next;
		}	
		EDGE_by_Element[K][1] = Indicator[I1] + count;

		current = EDGE_List[I2];
		count = 0;
		while (current->J != I3){
			count++;
			current = current->next;
		}	
		EDGE_by_Element[K][2] = Indicator[I2] + count;

	}
	
	lm = (int**) mycalloc("lm of 'EDE_Inicialization'",nedge, sizeof(int*));
	lmaux = (int*) mycalloc("lmaux of 'EDE_Inicialization'",2*nedge, sizeof(int));

	for (K = 0; K < nedge; K++)
		lm[K] = &lmaux[2*K];

	for (K = 0; K < neq; K++){
		current = EDGE_List[K];
		count = Indicator[K];
		while (current != NULL){
			J = current->J;
			lm[count][0] = K;
			lm[count][1] = J;
			count++;
			current = current->next;
		}
	}
		
	free(Indicator);
	
	for (K = 0; K < neq + 1; K++){
		current = EDGE_List[K];
		while (current != NULL){
			temp = current;		
			current	= current->next;
			free(temp);
		}

	}

	free(EDGE_List);
	
	Parameters->nedge = nedge;
	*lm_out = lm;
	*lmaux_out = lmaux;
	*EDGE_by_Element_out = EDGE_by_Element;
}


void ede_List_insertA(NodeListType **EDGE_List, int I, int J, int *nedge_out)
{
	int nedge;
	NodeListType *current, *previous, *new;

	nedge = *nedge_out;
	nedge++;
	new = mycalloc("new of 'EDE_Inicialization'", 1,sizeof(NodeListType));
	new->J = J;
	new->next = NULL;
	//only to initialization
	previous = new;
	previous->next = NULL;

	if (EDGE_List[I]==NULL)
		EDGE_List[I] = new;
	else if (EDGE_List[I]->J > J){
		new->next = EDGE_List[I];
		EDGE_List[I] = new;
	}
	else {
		current = EDGE_List[I];		
		while (current != NULL){
			if (J <= current->J)			
				break;
			previous = current;	
			current = current->next;
		}
		if (current == NULL)
			previous->next = new;
		else if (current->J == J){
			nedge--;
			free(new);
		}
		else {
			new->next = current;
			previous->next = new;
		}
	}
	*nedge_out = nedge;
}


