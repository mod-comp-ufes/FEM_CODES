#include "TranspEquation.h"

int uB_InitialSolution(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *u, double *uB)
{
	int nel, I, J1, J2, J3, I1, I2, I3;
	double X[3],Y[3];
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	
	nel = Parameters->nel;

	for (I=0;I<nel;I++){

		uB[I] = 0;

		J1 = Element[I].Vertex[0];
		J2 = Element[I].Vertex[1];
		J3 = Element[I].Vertex[2];

		I1 = Node[J1].id;
		I2 = Node[J2].id;
		I3 = Node[J3].id;
	
		X[0] = Node[J1].x;
		X[1] = Node[J2].x;
		X[2] = Node[J3].x;
		Y[0] = Node[J1].y;
		Y[1] = Node[J2].y;
		Y[2] = Node[J3].y;

		if (I1>=0)
			uB[I] += u[I1];
		else
			uB[I] += FemFunctions->upresc(X[0],Y[0]);
		
		if (I2>=0)
			uB[I] += u[I2];
		else
			uB[I] += FemFunctions->upresc(X[1],Y[1]);
	
		if (I3>=0)
			uB[I] += u[I3];
		else
			uB[I] += FemFunctions->upresc(X[2],Y[2]);
	
		uB[I]/=3;
	
	}
	return 0;
}
