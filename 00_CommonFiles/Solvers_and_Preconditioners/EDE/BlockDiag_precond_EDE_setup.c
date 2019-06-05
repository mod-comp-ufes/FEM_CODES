#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int EDE_matrix_Inverse(double *, double *);

int BlockDiag_precond_EDE_setup(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
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
	double **BlockDiag, *BlockDiagAux, **invBlockDiag, *invBlockDiagAux;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	int I_diag[2][32] = {{0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27},{36,37,38,39,44,45,46,47,52,53,54,55,60,61,62,63}};

	if (tag==1){
		IdAux = (int*) mycalloc("IdAux of 'BlockBlockDiag_precond_EDE'",4*nnodes,sizeof(int));
		Id = (int**) mycalloc("Id of 'BlockBlockDiag_precond_EDE'", nnodes,sizeof(int*));
		BlockDiagAux = (double*) mycalloc("BlockDiagAux of 'BlockBlockDiag_precond_EDE'",16*nnodes,sizeof(double));
		BlockDiag = (double**) mycalloc("BlockDiag of 'BlockBlockDiag_precond_EDE'",nnodes,sizeof(double*));
		invBlockDiagAux = (double*) mycalloc("invBlockDiagAux of 'BlockBlockDiag_precond_EDE'",16*nnodes,sizeof(double));
		invBlockDiag = (double**) mycalloc("invBlockDiag of 'BlockBlockDiag_precond_EDE'",nnodes,sizeof(double*));

		for (I=0;I<nnodes;I++){
			BlockDiag[I] = &BlockDiagAux[16*I]; 
			invBlockDiag[I] = &invBlockDiagAux[16*I]; 
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
		MatrixData->invBlockDiag = invBlockDiag;
	}
	else{
		Id = MatrixData->Id;
		BlockDiag = MatrixData->BlockDiag;
		BlockDiagAux = MatrixData->BlockDiagAux;
		memset(BlockDiagAux,0,16*nnodes*sizeof(double));
		invBlockDiag = MatrixData->invBlockDiag;
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

	for (I=0;I<nnodes;I++)
		EDE_matrix_Inverse(BlockDiag[I],invBlockDiag[I]);

	/* Preconditioning of F */
	BlockDiag_precond(Parameters, MatrixData, FemStructs, F, F);

	free(EDGEAux);
	free(EDGE);

	return 0;
}

int EDE_matrix_Inverse(double *BlockDiag, double *invBlockDiag)
{
	double a, b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,det;

	/*This function calculates the inverse of a 4 x 4 matrix 

	   |a b c d|	
	A =|e f g h| 	
	   |i j k l|	
	   |m n o p|	*/

	a = BlockDiag[0];
	b = BlockDiag[1];
	c = BlockDiag[2];
	d = BlockDiag[3];
	e = BlockDiag[4];
	f = BlockDiag[5];
	g = BlockDiag[6];
	h = BlockDiag[7];
	i = BlockDiag[8];
	j = BlockDiag[9];
	k = BlockDiag[10];
	l = BlockDiag[11];
	m = BlockDiag[12];
	n = BlockDiag[13];
	o = BlockDiag[14];
	p = BlockDiag[15];

	/* Determinant calculations */
	det = (a*f*k*p - a*f*l*o - a*g*j*p + a*g*l*n + a*h*j*o - a*h*k*n
			-b*e*k*p  + b*e*l*o + b*g*i*p - b*g*l*m - b*h*i*o + b*h*k*m +
			c*(e*j*p - e*l*n - f*i*p + f*l*m + h*i*n - h*j*m) +
			d*(-e*j*o + e*k*n + f*i*o - f*k*m - g*i*n + g*j*m));

	double invdet = 1.0/det;

	/* Inverse of matrix calculations */
	invBlockDiag[0] = invdet*(-h*k*n + g*l*n + h*j*o - f*l*o - g*j*p + f*k*p);		 	
	invBlockDiag[1] = invdet*(d*k*n - c*l*n - d*j*o + b*l*o + c*j*p - b*k*p);		 	
	invBlockDiag[2] = invdet*(-d*g*n + c*h*n + d*f*o - b*h*o - c*f*p + b*g*p);		 	
	invBlockDiag[3] = invdet*(d*g*j - c*h*j - d*f*k + b*h*k + c*f*l - b*g*l);		 	
	invBlockDiag[4] = invdet*(h*k*m - g*l*m - h*i*o + e*l*o + g*i*p - e*k*p);		 	
	invBlockDiag[5] = invdet*(-d*k*m + c*l*m + d*i*o - a*l*o - c*i*p + a*k*p);
	invBlockDiag[6] = invdet*(d*g*m - c*h*m - d*e*o + a*h*o + c*e*p - a*g*p);		 	
	invBlockDiag[7] = invdet*(-d*g*i + c*h*i + d*e*k - a*h*k - c*e*l + a*g*l);
	invBlockDiag[8] = invdet*(-h*j*m + f*l*m + h*i*n - e*l*n - f*i*p + e*j*p);
	invBlockDiag[9] = invdet*(d*j*m - b*l*m - d*i*n + a*l*n + b*i*p - a*j*p);		 	
	invBlockDiag[10] = invdet*(-d*f*m + b*h*m + d*e*n - a*h*n - b*e*p + a*f*p);
	invBlockDiag[11] = invdet*(d*f*i - b*h*i - d*e*j + a*h*j + b*e*l - a*f*l);		 	
	invBlockDiag[12] = invdet*(g*j*m - f*k*m - g*i*n + e*k*n + f*i*o - e*j*o);		 	
	invBlockDiag[13] = invdet*(-c*j*m + b*k*m + c*i*n - a*k*n - b*i*o + a*j*o);
	invBlockDiag[14] = invdet*(c*j*m - b*g*m - c*e*n + a*g*n + b*e*o - a*f*o);		 	
	invBlockDiag[15] = invdet*(-c*f*i + b*g*i + c*e*j - a*g*j - b*e*k + a*f*k);

	return 0;
}		







