#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int fill_Block_Ae(double [NDOF][NDOF], int, int, double *, int [NNOEL*NDOF][NOEL*NDOF]);
int matrix_Bondary_Conditions(double [NDOF][NDOF], double [NDOF][NDOF], double [NDOF][NDOF], int, int);
int matrix_Scaling(double *, int, double *, double [NDOF][NDOF], double [NDOF][NDOF], double [NDOF][NDOF], int [NNOEL*NDOF][NNOEL*NDOF]);

int BlockLU_precond_EBE_setup(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	int I, J, E, I1, I2, I3;
	double **A = MatrixData->A;
	int neq = Parameters-> neq;
	int nel = Parameters-> nel;
	int nnodes = Parameters->nnodes;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	int **Id, *Idaux;
	double **invDiag;

	BlockDiag_precond_EBE_setup(Parameters, MatrixData, FemStructs, tag, F);
	
	Id = MatrixData->Id;
	invDiag = MatrixData->invBlockDiag;

	double A11[NDOF][NDOF], A12[NDOF][NDOF], A13[NDOF][NDOF], A21[NDOF][NDOF], A22[NDOF][NDOF], A23[NDOF][NDOF], A31[NDOF][NDOF], A32[NDOF][NDOF], A33[NDOF][NDOF];
	
	for (E=0;E<nel;E++){
		I1 = Element[E].Vertex[0];	
		I2 = Element[E].Vertex[1];	
		I3 = Element[E].Vertex[2];	

		fill_Block_Ae(A11,0,0,A[E],I_general);
		fill_Block_Ae(A12,0,1,A[E],I_general);
		fill_Block_Ae(A13,0,2,A[E],I_general);
		fill_Block_Ae(A21,1,0,A[E],I_general);
		fill_Block_Ae(A22,1,1,A[E],I_general);
		fill_Block_Ae(A23,1,2,A[E],I_general);
		fill_Block_Ae(A31,2,0,A[E],I_general);
		fill_Block_Ae(A32,2,1,A[E],I_general);
		fill_Block_Ae(A33,2,2,A[E],I_general);

		/* Adjusting boundary conditions on matrices M of the elements*/
		for (I=0; I<NDOF; I++){
			matrix_Bondary_Conditions(A11, A12, A13, Node[I1].id[I], I);
			matrix_Bondary_Conditions(A22, A21, A23, Node[I2].id[I], I);
			matrix_Bondary_Conditions(A33, A31, A32, Node[I3].id[I], I);
		}

		/* Preconditioning of the Ae = Inv*Me */	
		matrix_Preconditioning(A[E],0, invDiag[I1], A11, A12, A13, I_general);
		matrix_Preconditioning(A[E],1, invDiag[I2], A21, A22, A23, I_general);
		matrix_Preconditioning(A[E],2, invDiag[I3], A31, A32, A33, I_general);
	}

	return 0;
}

int fill_Block_Ae(double M[NDOF][NDOF],int I, int J, double *A, int I_g[NNOEL*NDOF][NNOEL*NDOF])
{
	int i, j;

	for (i=0; i <NDOF; i++)
		for (j=0; j<NDOF; j++)
			M[i][j] = A[I_g[NDOF*I+i][NDOF*J+j]];

	return 0;
}


int matrix_Bondary_Conditions(double M1[NDOF][NDOF], double M2[NDOF][NDOF], double M3[NDOF][NDOF], int NodeId, int tag)
{
	int i;

	if (NodeId < 0){
		
		for (i=0; i<NDOF; i++){
			M1[tag][i] = 0;
			M1[i][tag] = 0;
			M2[tag][i] = 0;
			M3[tag][i] = 0;
		}
		M1[tag][tag] = 1;
	}

	return 0;	
}


int matrix_Scaling(double *A, int tag, double *invDiag, double M1[NDOF][NDOF], double M2[NDOF][NDOF], double M3[NDOF][NDOF], int I_g[NNOEL*NDOF][NNOEL*NDOF])
{
	int i, j, k;

	for (i=0; i<NDOF; i++){
		for (k=0; k<NDOF; k++){
			m1 = k;
			m2 = k + NDOF;
			m3 = k + 2*NDOF;
			n = NDOF*tag + i;
			A[I_g[n][m1]] = 0;
			A[I_g[n][m2]] = 0;
			A[I_g[n][m3]] = 0;
			for (j=0; j<NDOF; j++){
				A[I_g[n][m1] += invDiag[j]*M1[j][k]];
				A[I_g[n][m2] += invDiag[j]*M2[j][k]];
				A[I_g[n][m3] += invDiag[j]*M3[j][k]];
			}
		}
	}
/*
	A[I_g[4*tag][0]]     = invDiag[0]*M1[0][0]  + invDiag[1]*M1[1][0]  + invDiag[2]*M1[2][0]  + invDiag[3]*M1[3][0];
	A[I_g[4*tag][1]]     = invDiag[0]*M1[0][1]  + invDiag[1]*M1[1][1]  + invDiag[2]*M1[2][1]  + invDiag[3]*M1[3][1];
	A[I_g[4*tag][2]]     = invDiag[0]*M1[0][2]  + invDiag[1]*M1[1][2]  + invDiag[2]*M1[2][2]  + invDiag[3]*M1[3][2];
	A[I_g[4*tag][3]]     = invDiag[0]*M1[0][3]  + invDiag[1]*M1[1][3]  + invDiag[2]*M1[2][3]  + invDiag[3]*M1[3][3];

	A[I_g[4*tag][4]]     = invDiag[0]*M2[0][0]  + invDiag[1]*M2[1][0]  + invDiag[2]*M2[2][0]  + invDiag[3]*M2[3][0];
	A[I_g[4*tag][5]]     = invDiag[0]*M2[0][1]  + invDiag[1]*M2[1][1]  + invDiag[2]*M2[2][1]  + invDiag[3]*M2[3][1];
	A[I_g[4*tag][6]]     = invDiag[0]*M2[0][2]  + invDiag[1]*M2[1][2]  + invDiag[2]*M2[2][2]  + invDiag[3]*M2[3][2];
	A[I_g[4*tag][7]]     = invDiag[0]*M2[0][3]  + invDiag[1]*M2[1][3]  + invDiag[2]*M2[2][3]  + invDiag[3]*M2[3][3];

	A[I_g[4*tag][8]]     = invDiag[0]*M3[0][0]  + invDiag[1]*M3[1][0]  + invDiag[2]*M3[2][0]  + invDiag[3]*M3[3][0];
	A[I_g[4*tag][9]]     = invDiag[0]*M3[0][1]  + invDiag[1]*M3[1][1]  + invDiag[2]*M3[2][1]  + invDiag[3]*M3[3][1];
	A[I_g[4*tag][10]]    = invDiag[0]*M3[0][2]  + invDiag[1]*M3[1][2]  + invDiag[2]*M3[2][2]  + invDiag[3]*M3[3][2];
	A[I_g[4*tag][11]]    = invDiag[0]*M3[0][3]  + invDiag[1]*M3[1][3]  + invDiag[2]*M3[2][3]  + invDiag[3]*M3[3][3];

	A[I_g[4*tag+1][0]]     = invDiag[4]*M1[0][0]  + invDiag[5]*M1[1][0]  + invDiag[6]*M1[2][0]  + invDiag[7]*M1[3][0];
	A[I_g[4*tag+1][1]]     = invDiag[4]*M1[0][1]  + invDiag[5]*M1[1][1]  + invDiag[6]*M1[2][1]  + invDiag[7]*M1[3][1];
	A[I_g[4*tag+1][2]]     = invDiag[4]*M1[0][2]  + invDiag[5]*M1[1][2]  + invDiag[6]*M1[2][2]  + invDiag[7]*M1[3][2];
	A[I_g[4*tag+1][3]]     = invDiag[4]*M1[0][3]  + invDiag[5]*M1[1][3]  + invDiag[6]*M1[2][3]  + invDiag[7]*M1[3][3];

	A[I_g[4*tag+1][4]]     = invDiag[4]*M2[0][0]  + invDiag[5]*M2[1][0]  + invDiag[6]*M2[2][0]  + invDiag[7]*M2[3][0];
	A[I_g[4*tag+1][5]]     = invDiag[4]*M2[0][1]  + invDiag[5]*M2[1][1]  + invDiag[6]*M2[2][1]  + invDiag[7]*M2[3][1];
	A[I_g[4*tag+1][6]]     = invDiag[4]*M2[0][2]  + invDiag[5]*M2[1][2]  + invDiag[6]*M2[2][2]  + invDiag[7]*M2[3][2];
	A[I_g[4*tag+1][7]]     = invDiag[4]*M2[0][3]  + invDiag[5]*M2[1][3]  + invDiag[6]*M2[2][3]  + invDiag[7]*M2[3][3];

	A[I_g[4*tag+1][8]]     = invDiag[4]*M3[0][0]  + invDiag[5]*M3[1][0]  + invDiag[6]*M3[2][0]  + invDiag[7]*M3[3][0];
	A[I_g[4*tag+1][9]]     = invDiag[4]*M3[0][1]  + invDiag[5]*M3[1][1]  + invDiag[6]*M3[2][1]  + invDiag[7]*M3[3][1];
	A[I_g[4*tag+1][10]]    = invDiag[4]*M3[0][2]  + invDiag[5]*M3[1][2]  + invDiag[6]*M3[2][2]  + invDiag[7]*M3[3][2];
	A[I_g[4*tag+1][11]]    = invDiag[4]*M3[0][3]  + invDiag[5]*M3[1][3]  + invDiag[6]*M3[2][3]  + invDiag[7]*M3[3][3];

	A[I_g[4*tag+2][0]]     = invDiag[8]*M1[0][0]  + invDiag[9]*M1[1][0]  + invDiag[10]*M1[2][0]  + invDiag[11]*M1[3][0];
	A[I_g[4*tag+2][1]]     = invDiag[8]*M1[0][1]  + invDiag[9]*M1[1][1]  + invDiag[10]*M1[2][1]  + invDiag[11]*M1[3][1];
	A[I_g[4*tag+2][2]]     = invDiag[8]*M1[0][2]  + invDiag[9]*M1[1][2]  + invDiag[10]*M1[2][2]  + invDiag[11]*M1[3][2];
	A[I_g[4*tag+2][3]]     = invDiag[8]*M1[0][3]  + invDiag[9]*M1[1][3]  + invDiag[10]*M1[2][3]  + invDiag[11]*M1[3][3];

	A[I_g[4*tag+2][4]]     = invDiag[8]*M2[0][0]  + invDiag[9]*M2[1][0]  + invDiag[10]*M2[2][0]  + invDiag[11]*M2[3][0];
	A[I_g[4*tag+2][5]]     = invDiag[8]*M2[0][1]  + invDiag[9]*M2[1][1]  + invDiag[10]*M2[2][1]  + invDiag[11]*M2[3][1];
	A[I_g[4*tag+2][6]]     = invDiag[8]*M2[0][2]  + invDiag[9]*M2[1][2]  + invDiag[10]*M2[2][2]  + invDiag[11]*M2[3][2];
	A[I_g[4*tag+2][7]]     = invDiag[8]*M2[0][3]  + invDiag[9]*M2[1][3]  + invDiag[10]*M2[2][3]  + invDiag[11]*M2[3][3];

	A[I_g[4*tag+2][8]]     = invDiag[8]*M3[0][0]  + invDiag[9]*M3[1][0]  + invDiag[10]*M3[2][0]  + invDiag[11]*M3[3][0];
	A[I_g[4*tag+2][9]]     = invDiag[8]*M3[0][1]  + invDiag[9]*M3[1][1]  + invDiag[10]*M3[2][1]  + invDiag[11]*M3[3][1];
	A[I_g[4*tag+2][10]]    = invDiag[8]*M3[0][2]  + invDiag[9]*M3[1][2]  + invDiag[10]*M3[2][2]  + invDiag[11]*M3[3][2];
	A[I_g[4*tag+2][11]]    = invDiag[8]*M3[0][3]  + invDiag[9]*M3[1][3]  + invDiag[10]*M3[2][3]  + invDiag[11]*M3[3][3];

	A[I_g[4*tag+3][0]]     = invDiag[12]*M1[0][0]  + invDiag[13]*M1[1][0]  + invDiag[14]*M1[2][0]  + invDiag[15]*M1[3][0];
	A[I_g[4*tag+3][1]]     = invDiag[12]*M1[0][1]  + invDiag[13]*M1[1][1]  + invDiag[14]*M1[2][1]  + invDiag[15]*M1[3][1];
	A[I_g[4*tag+3][2]]     = invDiag[12]*M1[0][2]  + invDiag[13]*M1[1][2]  + invDiag[14]*M1[2][2]  + invDiag[15]*M1[3][2];
	A[I_g[4*tag+3][3]]     = invDiag[12]*M1[0][3]  + invDiag[13]*M1[1][3]  + invDiag[14]*M1[2][3]  + invDiag[15]*M1[3][3];

	A[I_g[4*tag+3][4]]     = invDiag[12]*M2[0][0]  + invDiag[13]*M2[1][0]  + invDiag[14]*M2[2][0]  + invDiag[15]*M2[3][0];
	A[I_g[4*tag+3][5]]     = invDiag[12]*M2[0][1]  + invDiag[13]*M2[1][1]  + invDiag[14]*M2[2][1]  + invDiag[15]*M2[3][1];
	A[I_g[4*tag+3][6]]     = invDiag[12]*M2[0][2]  + invDiag[13]*M2[1][2]  + invDiag[14]*M2[2][2]  + invDiag[15]*M2[3][2];
	A[I_g[4*tag+3][7]]     = invDiag[12]*M2[0][3]  + invDiag[13]*M2[1][3]  + invDiag[14]*M2[2][3]  + invDiag[15]*M2[3][3];

	A[I_g[4*tag+3][8]]     = invDiag[12]*M3[0][0]  + invDiag[13]*M3[1][0]  + invDiag[14]*M3[2][0]  + invDiag[15]*M3[3][0];
	A[I_g[4*tag+3][9]]     = invDiag[12]*M3[0][1]  + invDiag[13]*M3[1][1]  + invDiag[14]*M3[2][1]  + invDiag[15]*M3[3][1];
	A[I_g[4*tag+3][10]]    = invDiag[12]*M3[0][2]  + invDiag[13]*M3[1][2]  + invDiag[14]*M3[2][2]  + invDiag[15]*M3[3][2];
	A[I_g[4*tag+3][11]]    = invDiag[12]*M3[0][3]  + invDiag[13]*M3[1][3]  + invDiag[14]*M3[2][3]  + invDiag[15]*M3[3][3];
*/
	return 0;

}
