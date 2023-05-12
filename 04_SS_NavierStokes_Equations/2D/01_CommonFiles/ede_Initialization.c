#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

void ede_Initialization(int nnodes, int neq, int nel, NodeType *Node, ElementType *Element, int *nedge_out, int **order, int ***lm_out, int **lmaux_out, int ***EDGE_by_Element_out)
{
	int nedge, K, I, I1, I2, I3, J1, J2, J3;
	int **EDGE_by_Element, **lm, *lmaux, Io[3];;
	NodeListType *current, *temp, **EDGE_List;


	lm = *lm_out;
	lmaux = *lmaux_out;

	EDGE_List = (NodeListType **)mycalloc("EDGE_List of 'EDE_Inicialization'",nnodes,sizeof(NodeListType *));

	for (K = 0; K < nnodes; K++)
		EDGE_List[K] = NULL;	

	nedge = 0;
	for (K = 0; K < nel ; K++)
	{
		J1 = Element[K].Vertex[0];
		J2 = Element[K].Vertex[1];
		J3 = Element[K].Vertex[2];

		if (J1<J2 && J1<J3){
			I1 = J1;
			if (J2<J3){
				I2 = J2; 
				I3 = J3;
				order[K][0] = 0;
				order[K][1] = 1;
				order[K][2] = 2;
			}
			else{
				I2 = J3;
				I3 = J2;
				order[K][0] = 0;
				order[K][1] = 2;
				order[K][2] = 1;
			}

		}
		else if (J2<J3){
			I1 = J2;
			if (J1<J3){
				I2 = J1;
				I3 = J3;
				order[K][0] = 1;
				order[K][1] = 0;
				order[K][2] = 2;
			}
			else{
				I2 = J3;
				I3 = J1;
				order[K][0] = 1;
				order[K][1] = 2;
				order[K][2] = 0;
			}

		}
		else{
			I1 = J3;
			if (J1<J2){
				I2 = J1;
				I3 = J2;
				order[K][0] = 2;
				order[K][1] = 0;
				order[K][2] = 1;
			}
			else{
				I2 = J2;
				I3 = J1;
				order[K][0] = 2;
				order[K][1] = 1;
				order[K][2] = 0;
			}

		}

		ede_List_insertA( EDGE_List, I1,I2, &nedge);
		ede_List_insertA( EDGE_List, I1,I3, &nedge);
		ede_List_insertA( EDGE_List, I2,I3, &nedge);

	}

	EDGE_by_Element = (int**) mycalloc("EDGE_by_Element of 'EDE_Inicialization'",nel,sizeof(int*));
	for (K = 0; K < nel; K++)
		EDGE_by_Element[K] = (int*) mycalloc("EDGE_by_Element of 'EDE_Inicialization'", NNOEL,sizeof(int));
		
	int Reverse[3];

	for (K = 0; K < nel ; K++){
	
		I1 = Element[K].Vertex[0];
		I2 = Element[K].Vertex[1];
		I3 = Element[K].Vertex[2];

		Reverse[order[K][0]]=0;
		Reverse[order[K][1]]=1;
		Reverse[order[K][2]]=2;
		
		Io[Reverse[0]] = I1;
		Io[Reverse[1]] = I2;
		Io[Reverse[2]] = I3;
		
		current = EDGE_List[Io[0]];
		while (current->J != Io[1]){
			current = current->next;
		}
		EDGE_by_Element[K][0] = current->count;
		
		current = EDGE_List[Io[0]];
		while (current->J != Io[2]){
			current = current->next;
		}
		EDGE_by_Element[K][1] = current->count;

		current = EDGE_List[Io[1]];
		while (current->J != Io[2]){
			current = current->next;
		}
		EDGE_by_Element[K][2] = current->count;


	}
	
	myfree(lm);
	myfree(lmaux);

	lm = (int**) mycalloc("lm of 'EDE_Inicialization'",nedge, sizeof(int*));
	lmaux = (int*) mycalloc("lmaux of 'EDE_Inicialization'",2*NDOF*nedge, sizeof(int));

	for (K = 0; K < nedge; K++)
		lm[K] = &lmaux[2*NDOF*K];

	for (K = 0; K < nel; K++){

		I1 = Element[K].Vertex[0];
		I2 = Element[K].Vertex[1];
		I3 = Element[K].Vertex[2];

		Reverse[order[K][0]]=0;
		Reverse[order[K][1]]=1;
		Reverse[order[K][2]]=2;
		
		Io[Reverse[0]] = I1;
		Io[Reverse[1]] = I2;
		Io[Reverse[2]] = I3;
		
		for (I=0 ; I<NDOF; I++){

			if (Node[Io[0]].id[I]==-1){
				lm[EDGE_by_Element[K][0]][I] = neq; 
				lm[EDGE_by_Element[K][1]][I] = neq; 
			}
			else{
				lm[EDGE_by_Element[K][0]][I] = Node[Io[0]].id[I]; 
				lm[EDGE_by_Element[K][1]][I] = Node[Io[0]].id[I]; 

			}
			
			if (Node[Io[1]].id[I]==-1){
				lm[EDGE_by_Element[K][0]][I+NDOF] = neq; 
				lm[EDGE_by_Element[K][2]][I] = neq; 
			}
			else{
				lm[EDGE_by_Element[K][0]][I+NDOF] = Node[Io[1]].id[I]; 
				lm[EDGE_by_Element[K][2]][I] = Node[Io[1]].id[I]; 
			}
			
			if (Node[Io[2]].id[I]==-1){
				lm[EDGE_by_Element[K][1]][I+NDOF] = neq; 
				lm[EDGE_by_Element[K][2]][I+NDOF] = neq; 
			}
			else{
				lm[EDGE_by_Element[K][1]][I+NDOF] = Node[Io[2]].id[I]; 
				lm[EDGE_by_Element[K][2]][I+NDOF] = Node[Io[2]].id[I]; 
			}
		}
	

	}
		
	
	for (K = 0; K < nnodes; K++){
		current = EDGE_List[K];
		while (current != NULL){
			temp = current;		
			current	= current->next;
			myfree(temp);
		}

	}

	myfree(EDGE_List);
	
	*nedge_out = nedge;
	*lm_out = lm;
	*lmaux_out = lmaux;
	*EDGE_by_Element_out = EDGE_by_Element;
}


void ede_List_insertA(NodeListType **EDGE_List, int I, int J, int *nedge_out)
{
	int nedge;
	NodeListType *current=NULL, *previous=NULL, *new=NULL;

	nedge = *nedge_out;
	new = mycalloc("new of 'EDE_Inicialization'", 1,sizeof(NodeListType));
	new->count = nedge++;
	new->J = J;
	new->next = NULL;

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
			myfree(new);
		}
		else {
			new->next = current;
			previous->next = new;
		}
	}
	*nedge_out = nedge;
}


