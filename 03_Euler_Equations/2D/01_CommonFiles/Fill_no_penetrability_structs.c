#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocations/allocations.h"

int Fill_no_penetrability_structs(ParametersType *Parameters, FemStructsType *FemStructs)
{
	int I, J, J1, J2, J3;
	double x1, x2, x3;
	int **neighbors, *neighborsAux;
	int nonpnodes = Parameters->nonpnodes;
	int nel = Parameters->nel;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStruts->Element;

	if (nonpnodes==0)
		return 0;
	
	NeighborsAux = mycalloc("neighborsAux of 'Fill_no_penetrability_structs'",2*nonpnodes,sizeof(int));
	Neighbors = mycalloc("neighbors of 'Fill_no_penetrability_structs'",nonpnodes,sizeof(int*));		//Neighbors[I][0] is the neighbor before
	npnodeInNode =mycalloc("npnodeInNode of 'Fill_no_penetrability_structs'",nonpnodes,size(int));		//Neighbors[I][1] is the neighbor after
	for (I=0, I<nonpnodes;I++){
		Neighbors[I] = &NeighborsAux[2*I];
		Neighbors[I][0] = -1;
		Neighbors[I][1] = -1;
	}

	count = 0;
	for (I=0;I<nel;I++){
		J1 = Element[I].Vertex[0];
		J2 = Element[I].Vertex[1];
		J3 = Element[I].Vertex[2];
		
		x1 = Node[J1].x;

		x2 = Node[J2].x;
		
		x3 = Node[J3].x;

		if (Node[J1].v1Type==-1){
			if (Node[J2].v1Type<=0){
				if (x1>x2){
					if (Neighbors[count][0]==-1&&Neigbors[count][1]==-1){
						count++;
						Neighbors[count-1][0] = J2;
					}
					else
						Neighbors[count][0] = J2;
				}
				else{
					if (Neighbors[count][0]==-1&&Neigbors[count][1]==-1){
						count++;
						Neighbors[count-1][1] = J2;
					}
					else
						Neighbors[count][1] = J2;
								
				}

			}
			if (Node[J3].v1Type<=0){
				if (x1>x3){
					if (Neighbors[count][0]==-1&&Neigbors[count][1]==-1){
						count++;
						Neighbors[count-1][0] = J3;
					}
					else
						Neighbors[count][0] = J3;
				}
				else{
					if (Neighbors[count][0]==-1&&Neigbors[count][1]==-1){
						count++;
						Neighbors[count-1][1] = J3;
					}
					else
						Neighbors[count][1] = J3;
								
				}

			}
		}

		if (Node[J2].v1Type==-1){
			if (Node[J1].v1Type<=0){
				if (x2>x1){
					if (Neighbors[count][0]==-1&&Neigbors[count][1]==-1){
						count++;
						Neighbors[count-1][0] = J1;
					}
					else
						Neighbors[count][0] = J1;
				}
				else{
					if (Neighbors[count][0]==-1&&Neigbors[count][1]==-1){
						count++;
						Neighbors[count-1][1] = J1;
					}
					else
						Neighbors[count][1] = J1;
								
				}

			}
			if (Node[J3].v1Type<=0){
				if (x2>x3){
					if (Neighbors[count][0]==-1&&Neigbors[count][1]==-1){
						count++;
						Neighbors[count-1][0] = J3;
					}
					else
						Neighbors[count][0] = J3;
				}
				else{
					if (Neighbors[count][0]==-1&&Neigbors[count][1]==-1){
						count++;
						Neighbors[count-1][1] = J3;
					}
					else
						Neighbors[count][1] = J3;
								
				}

			}
		}



			


	}	
	
	


}
