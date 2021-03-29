#include "SSTranspEquation.h"

void ede_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*ke)[3])
{
	double **K;
	int I1, I2, I3, J1, J2, J3, *order, *EDGE_by_Element;
	int nedge = Parameters->nedge;

	K = MatrixData->A;
	order = MatrixData->order[E];
	I1 = order[0];
	I2 = order[1];
	I3 = order[2];
	EDGE_by_Element = MatrixData->Scheme_by_Element[E];
	J1 = EDGE_by_Element[0];
	J2 = EDGE_by_Element[1];
	J3 = EDGE_by_Element[2];
	K[J1][0] += 0.5*ke[I1][I1];
	K[J1][1] += ke[I1][I2];
	K[J1][2] += ke[I2][I1];
	K[J1][3] += 0.5*ke[I2][I2];
	K[J2][0] += 0.5*ke[I1][I1];
	K[J2][1] += ke[I1][I3];
	K[J2][2] += ke[I3][I1];
	K[J2][3] += 0.5*ke[I3][I3];
	K[J3][0] += 0.5*ke[I2][I2];
	K[J3][1] += ke[I2][I3];
	K[J3][2] += ke[I3][I2];
	K[J3][3] += 0.5*ke[I3][I3];
	K[nedge][0] = 0;
	K[nedge][1] = 0;
	K[nedge][2] = 0;
	K[nedge][3] = 0;

}




