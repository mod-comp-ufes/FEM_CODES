#include "preconditioners.h"
#include "../Allocation_Operations/allocations.h"

int NOBlockDiag2_precond_EDE_setup(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I, J, E, I1, I2, I3;
	double **A = MatrixData->A;
	int **order = MatrixData->order;
	int **EDGE_by_Element = MatrixData->Scheme_by_Element;
	int neq = Parameters-> neq;
	int nel = Parameters-> nel;
	int nedge = Parameters-> nedge;
	int nnodes = Parameters->nnodes;
	int **Id, *IdAux;
	double **BlockDiag, *BlockDiagAux;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	int I_diag[2][32] = {{0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27},{36,37,38,39,44,45,46,47,52,53,54,55,60,61,62,63}};

	if (tag==1){
		IdAux = (int*) mycalloc("IdAux of 'BlockBlockDiag_precond_EDE'",4*nnodes,sizeof(int));
		Id = (int**) mycalloc("Id of 'BlockBlockDiag_precond_EDE'", nnodes,sizeof(int*));
		BlockDiagAux = (double*) mycalloc("BlockDiagAux of 'BlockBlockDiag_precond_EDE'",16*nnodes,sizeof(double));
		BlockDiag = (double**) mycalloc("BlockDiag of 'BlockBlockDiag_precond_EDE'",nnodes,sizeof(double*));

		for (I=0;I<nnodes;I++){
			BlockDiag[I] = &BlockDiagAux[16*I]; 
			Id[I] = &IdAux[4*I];
		}
	
		for (I=0;I<nnodes;I++){

			if (Node[I].id[0] < 0) 
				Id[I][0] = neq;
			else
				Id[I][0] = Node[I].id[0];

			if (Node[I].id[1] < 0) 
				Id[I][1] = neq; 
			else
				Id[I][1] = Node[I].id[1];

			if (Node[I].id[2] < 0) 
				Id[I][2] = neq; 
			else
				Id[I][2] = Node[I].id[2];

			if (Node[I].id[3] < 0)
				Id[I][3] = neq; 
			else
				Id[I][3] = Node[I].id[3];
		}
		MatrixData->BlockDiag = BlockDiag;
		MatrixData->BlockDiagAux = BlockDiagAux;
		MatrixData->Id = Id;
	}
	else{
		Id = MatrixData->Id;
		BlockDiag = MatrixData->BlockDiag;
		BlockDiagAux = MatrixData->BlockDiagAux;
		memset(BlockDiagAux,0,16*nnodes*sizeof(double));
	}	


	int Io[3], Reverse[3];
	int **EDGE, *EDGEAux;

	EDGEAux = (int*) mycalloc("EDGEAux of 'BlocBlockDiag_precond_EDE'", 2*nedge, sizeof(int));
	EDGE = (int**) mycalloc("EDGEAux of 'BlocBlockDiag_precond_EDE'", nedge, sizeof(int*));

	for (I=0;I<nedge; I++)
		EDGE[I] = &EDGEAux[2*I];
		
	for (E = 0; E < nel ; E++){
	
		I1 = Element[E].Vertex[0];
		I2 = Element[E].Vertex[1];
		I3 = Element[E].Vertex[2];

		Reverse[order[E][0]]=0;
		Reverse[order[E][1]]=1;
		Reverse[order[E][2]]=2;
		
		Io[Reverse[0]] = I1;
		Io[Reverse[1]] = I2;
		Io[Reverse[2]] = I3;

		EDGE[EDGE_by_Element[E][0]][0] = Io[0];
		EDGE[EDGE_by_Element[E][0]][1] = Io[1];
		EDGE[EDGE_by_Element[E][1]][0] = Io[0];
		EDGE[EDGE_by_Element[E][1]][1] = Io[2];
		EDGE[EDGE_by_Element[E][2]][0] = Io[1];
		EDGE[EDGE_by_Element[E][2]][1] = Io[2];

	}

	for (E=0;E<nedge;E++){

		I1 = EDGE[E][0];	
		I2 = EDGE[E][1];
	
		for (I=0;I<4;I++){
			for(J=0;J<4;J++){
				BlockDiag[I1][4*I+J] += A[E][I_diag[0][4*I+J]];
				BlockDiag[I2][4*I+J] += A[E][I_diag[1][4*I+J]];
			}
		}
	}

	myfree(EDGEAux);
	myfree(EDGE);

	return 0;
}





