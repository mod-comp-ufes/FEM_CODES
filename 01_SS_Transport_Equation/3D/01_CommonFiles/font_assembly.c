#include "SSTransportEquation3D.h"

void font_assembly(int e, double Fe[4], double (*Ke)[4], FemFunctionsType *FemFunctions, FemStructsType *FemStructs, int neq, double A, double Lx)
{
	int **lm, J1, J2, J3, J4;
	double *F;
	double g[4];
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	
	F = FemStructs->F;
	lm = FemStructs->lm;
	
	// Global node that composes the element
	J1 = Element[e].Vertex[0];
	J2 = Element[e].Vertex[1];
	J3 = Element[e].Vertex[2];
	J4 = Element[e].Vertex[3];
	
	// Fill auxiliar vector g (type == 1 incognita, type == 0 prescribed)
	if(Node[J1].uType == 1){
		g[0] = 0.0;
	}else{
		g[0] = FemFunctions->upresc(Node[J1].x, Node[J1].y, Node[J1].z, A, Lx);
	}
	
	if(Node[J2].uType == 1){ 
		g[1] = 0.0;
	}else{
		g[1] = FemFunctions->upresc(Node[J2].x, Node[J2].y, Node[J2].z, A, Lx);
	}
	
	if(Node[J3].uType == 1){ 
		g[2] = 0.0;
	}else{
		g[2] = FemFunctions->upresc(Node[J3].x, Node[J3].y, Node[J3].z, A, Lx);
	}
	
	if(Node[J4].uType == 1){ 
		g[3] = 0.0;
	}else{
		g[3] = FemFunctions->upresc(Node[J4].x, Node[J4].y, Node[J4].z, A, Lx);
	}
	
	Fe[0] = Fe[0] - (Ke[0][0]*g[0] + Ke[0][1]*g[1] + Ke[0][2]*g[2] + Ke[0][3]*g[3]);
	Fe[1] = Fe[1] - (Ke[1][0]*g[0] + Ke[1][1]*g[1] + Ke[1][2]*g[2] + Ke[1][3]*g[3]);
	Fe[2] = Fe[2] - (Ke[2][0]*g[0] + Ke[2][1]*g[1] + Ke[2][2]*g[2] + Ke[2][3]*g[3]);
	Fe[3] = Fe[3] - (Ke[3][0]*g[0] + Ke[3][1]*g[1] + Ke[3][2]*g[2] + Ke[3][3]*g[3]);
	
	// preenchendo vetor fonte F
	F[lm[e][0]]  = F[lm[e][0]] + Fe[0];
	F[lm[e][1]]  = F[lm[e][1]] + Fe[1];
	F[lm[e][2]]  = F[lm[e][2]] + Fe[2];
	F[lm[e][3]]  = F[lm[e][3]] + Fe[3];
	
	F[neq] = 0.0;
}
