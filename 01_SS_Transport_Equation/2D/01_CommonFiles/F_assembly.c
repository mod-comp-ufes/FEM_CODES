#include "SSTranspEquation.h"

void F_assembly(NodeType *Node, int J1, int J2, int J3, double Fe1, double Fe2, double Fe3, double *F, double ke[3][3],double X[3], double Y[3], FemFunctionsType *FemFunctions)
{
	if (Node[J1].Type >= 1){ 
		F[Node[J1].id] += Fe1;	
		if (Node[J2].Type == 0)
			F[Node[J1].id] -= ke[0][1]*FemFunctions->upresc(X[1],Y[1]);
		if (Node[J3].Type == 0)
			F[Node[J1].id] -= ke[0][2]*FemFunctions->upresc(X[2],Y[2]);
	}

	if (Node[J2].Type >= 1){ 
		F[Node[J2].id] += Fe2;	
		if (Node[J1].Type == 0)
			F[Node[J2].id] -= ke[1][0]*FemFunctions->upresc(X[0],Y[0]);
		if (Node[J3].Type == 0)
			F[Node[J2].id] -= ke[1][2]*FemFunctions->upresc(X[2],Y[2]);
	}

	if (Node[J3].Type >= 1){ 
		F[Node[J3].id] += Fe3;	
		if (Node[J1].Type == 0)
			F[Node[J3].id] -= ke[2][0]*FemFunctions->upresc(X[0],Y[0]);
		if (Node[J2].Type == 0)
			F[Node[J3].id] -= ke[2][1]*FemFunctions->upresc(X[1],Y[1]);
	}

	// Fluxo
	double h1, h2, h3, l;
	double sixth = 1/6;
	h1 = FemFunctions->hflux(X[0], Y[0]);
	h2 = FemFunctions->hflux(X[1], Y[1]);
	h3 = FemFunctions->hflux(X[2], Y[2]);

	if (Node[J1].Type == 2 && Node[J2].Type == 2){ 
		l = sqrt((X[0] - X[1])*(X[0] - X[1]) + (Y[0] - Y[1])*(Y[0] - Y[1]));
		F[Node[J1].id] += sixth*l*(2*h1 + h2);
		F[Node[J2].id] += sixth*l*(h1 + 2*h2);
	}
	else if (Node[J1].Type == 2 && Node[J3].Type == 2){ 
		l = sqrt((X[0] - X[2])*(X[0] - X[2]) + (Y[0] - Y[2])*(Y[0] - Y[2]));
		F[Node[J1].id] += sixth*l*(2*h1 + h3);
		F[Node[J3].id] += sixth*l*(h1 + 2*h3);
	}
	else if (Node[J1].Type == 2 && Node[J2].Type == 2){ 
		l = sqrt((X[1] - X[2])*(X[1] - X[2]) + (Y[1] - Y[2])*(Y[1] - Y[2]));
		F[Node[J2].id] += sixth*l*(2*h2 + h3);
		F[Node[J3].id] += sixth*l*(h2 + 2*h3);
	}
}


