#include "EulerEquations.h"

void ede_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*Me)[12])
{
	double **M;
	int I1, I2, I3, J1, J2, J3, *order, *EDGE_by_Element;

	M = MatrixData->A;
	order = MatrixData->order[E];
	I1 = order[0]*NDOF;
	I2 = order[1]*NDOF;
	I3 = order[2]*NDOF;
	EDGE_by_Element = MatrixData->Scheme_by_Element[E];
	J1 = EDGE_by_Element[0];
	J2 = EDGE_by_Element[1];
	J3 = EDGE_by_Element[2];

	//Edge I1_I2 with 64 contributions
	edge_JI_assembly(J1,I1,I2, Me, M);
	

	//Edge I1_I3 with 64 contributions
	edge_JI_assembly(J2,I1,I3, Me, M);


	//Edge I2_I3 with 64 contributions
	edge_JI_assembly(J3,I2,I3, Me, M);

}

void edge_JI_assembly(int J, int I1, int I2, double (*Me)[12], double **M)
{
	M[J][0] += 0.5*Me[I1][I1];
	M[J][1] += 0.5*Me[I1][I1+1];
	M[J][2] += 0.5*Me[I1][I1+2];
	M[J][3] += 0.5*Me[I1][I1+3];

	M[J][4] += Me[I1][I2];
	M[J][5] += Me[I1][I2+1];
	M[J][6] += Me[I1][I2+2];
	M[J][7] += Me[I1][I2+3];
	
	M[J][8] += 0.5*Me[I1+1][I1];
	M[J][9] += 0.5*Me[I1+1][I1+1];
	M[J][10] += 0.5*Me[I1+1][I1+2];
	M[J][11] += 0.5*Me[I1+1][I1+3];

	M[J][12] += Me[I1+1][I2];
	M[J][13] += Me[I1+1][I2+1];
	M[J][14] += Me[I1+1][I2+2];
	M[J][15] += Me[I1+1][I2+3];

	M[J][16] += 0.5*Me[I1+2][I1];
	M[J][17] += 0.5*Me[I1+2][I1+1];
	M[J][18] += 0.5*Me[I1+2][I1+2];
	M[J][19] += 0.5*Me[I1+2][I1+3];

	M[J][20] += Me[I1+2][I2];
	M[J][21] += Me[I1+2][I2+1];
	M[J][22] += Me[I1+2][I2+2];
	M[J][23] += Me[I1+2][I2+3];

	M[J][24] += 0.5*Me[I1+3][I1];
	M[J][25] += 0.5*Me[I1+3][I1+1];
	M[J][26] += 0.5*Me[I1+3][I1+2];
	M[J][27] += 0.5*Me[I1+3][I1+3];

	M[J][28] += Me[I1+3][I2];
	M[J][29] += Me[I1+3][I2+1];
	M[J][30] += Me[I1+3][I2+2];
	M[J][31] += Me[I1+3][I2+3];

	M[J][32] += Me[I2][I1];
	M[J][33] += Me[I2][I1+1];
	M[J][34] += Me[I2][I1+2];
	M[J][35] += Me[I2][I1+3];

	M[J][36] += 0.5*Me[I2][I2];
	M[J][37] += 0.5*Me[I2][I2+1];
	M[J][38] += 0.5*Me[I2][I2+2];
	M[J][39] += 0.5*Me[I2][I2+3];

	M[J][40] += Me[I2+1][I1];
	M[J][41] += Me[I2+1][I1+1];
	M[J][42] += Me[I2+1][I1+2];
	M[J][43] += Me[I2+1][I1+3];

	M[J][44] += 0.5*Me[I2+1][I2];
	M[J][45] += 0.5*Me[I2+1][I2+1];
	M[J][46] += 0.5*Me[I2+1][I2+2];
	M[J][47] += 0.5*Me[I2+1][I2+3];

	M[J][48] += Me[I2+2][I1];
	M[J][49] += Me[I2+2][I1+1];
	M[J][50] += Me[I2+2][I1+2];
	M[J][51] += Me[I2+2][I1+3];

	M[J][52] += 0.5*Me[I2+2][I2];
	M[J][53] += 0.5*Me[I2+2][I2+1];
	M[J][54] += 0.5*Me[I2+2][I2+2];
	M[J][55] += 0.5*Me[I2+2][I2+3];

	M[J][56] += Me[I2+3][I1];
	M[J][57] += Me[I2+3][I1+1];
	M[J][58] += Me[I2+3][I1+2];
	M[J][59] += Me[I2+3][I1+3];

	M[J][60] += 0.5*Me[I2+3][I2];
	M[J][61] += 0.5*Me[I2+3][I2+1];
	M[J][62] += 0.5*Me[I2+3][I2+2];
	M[J][63] += 0.5*Me[I2+3][I2+3];
}

