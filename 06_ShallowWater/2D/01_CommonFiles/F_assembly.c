#include "ShalowWater.h"


void F_assembly(int e, double Fe[9], double De[9][9], FemFunctionsType *FemFunctions, FemStructsType *FemStructs, int neq)
{
	int **lm, J1, J2, J3;
	double *F;
	double g[9];
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	
	F = FemStructs->F;
	lm = FemStructs->lm;
	
	// Global node that composes the element
	J1 = Element[e].Vertex[0];
	J2 = Element[e].Vertex[1];
	J3 = Element[e].Vertex[2];
	
	// Fill auxiliar vector g (type == 1 incognita, type == 0 prescribed)
	// height
	if(Node[J1].hType == 1){
		g[0] = 0.0;
	}else{
		g[0] = FemFunctions->hpresc(Node[J1].x, Node[J1].y);
	}
	
	if(Node[J2].hType == 1){ 
		g[3] = 0.0;
	}else{
		g[3] = FemFunctions->hpresc(Node[J2].x, Node[J2].y);
	}
	
	if(Node[J3].hType == 1){ 
		g[6] = 0.0;
	}else{
		g[6] = FemFunctions->hpresc(Node[J3].x, Node[J3].y);
	}

	// discharge in X
	if(Node[J1].qxType == 1){
		g[1] = 0.0;
	}else{
		g[1] = FemFunctions->qxpresc(Node[J1].x, Node[J1].y);
	}
	
	if(Node[J2].qxType == 1){ 
		g[4] = 0.0;
	}else{
		g[4] = FemFunctions->qxpresc(Node[J2].x, Node[J2].y);
	}
	
	if(Node[J3].qxType == 1){ 
		g[7] = 0.0;
	}else{
		g[7] = FemFunctions->qxpresc(Node[J3].x, Node[J3].y);
	}

	// discharge in Y
	if(Node[J1].qyType == 1){
		g[2] = 0.0;
	}else{
		g[2] = FemFunctions->qypresc(Node[J1].x, Node[J1].y);
	}
	
	if(Node[J2].qyType == 1){ 
		g[5] = 0.0;
	}else{
		g[5] = FemFunctions->qypresc(Node[J2].x, Node[J2].y);
	}
	
	if(Node[J3].qyType == 1){ 
		g[8] = 0.0;
	}else{
		g[8] = FemFunctions->qypresc(Node[J3].x, Node[J3].y);
	}


	Fe[0] -= De[0][0]*g[0] + De[0][1]*g[1] + De[0][2]*g[2] + De[0][3]*g[3] + De[0][4]*g[4] + De[0][5]*g[5] + De[0][6]*g[6] + De[0][7]*g[7] + De[0][8]*g[8];
	Fe[1] -= De[1][0]*g[0] + De[1][1]*g[1] + De[1][2]*g[2] + De[1][3]*g[3] + De[1][4]*g[4] + De[1][5]*g[5] + De[1][6]*g[6] + De[1][7]*g[7] + De[1][8]*g[8];
	Fe[2] -= De[2][0]*g[0] + De[2][1]*g[1] + De[2][2]*g[2] + De[2][3]*g[3] + De[2][4]*g[4] + De[2][5]*g[5] + De[2][6]*g[6] + De[2][7]*g[7] + De[2][8]*g[8];
	Fe[3] -= De[3][0]*g[0] + De[3][1]*g[1] + De[3][2]*g[2] + De[3][3]*g[3] + De[3][4]*g[4] + De[3][5]*g[5] + De[3][6]*g[6] + De[3][7]*g[7] + De[3][8]*g[8];
	Fe[4] -= De[4][0]*g[0] + De[4][1]*g[1] + De[4][2]*g[2] + De[4][3]*g[3] + De[4][4]*g[4] + De[4][5]*g[5] + De[4][6]*g[6] + De[4][7]*g[7] + De[4][8]*g[8];
	Fe[5] -= De[5][0]*g[0] + De[5][1]*g[1] + De[5][2]*g[2] + De[5][3]*g[3] + De[5][4]*g[4] + De[5][5]*g[5] + De[5][6]*g[6] + De[5][7]*g[7] + De[5][8]*g[8];
	Fe[6] -= De[6][0]*g[0] + De[6][1]*g[1] + De[6][2]*g[2] + De[6][3]*g[3] + De[6][4]*g[4] + De[6][5]*g[5] + De[6][6]*g[6] + De[6][7]*g[7] + De[6][8]*g[8];
	Fe[7] -= De[7][0]*g[0] + De[7][1]*g[1] + De[7][2]*g[2] + De[7][3]*g[3] + De[7][4]*g[4] + De[7][5]*g[5] + De[7][6]*g[6] + De[7][7]*g[7] + De[7][8]*g[8];
	Fe[8] -= De[8][0]*g[0] + De[8][1]*g[1] + De[8][2]*g[2] + De[8][3]*g[3] + De[8][4]*g[4] + De[8][5]*g[5] + De[8][6]*g[6] + De[8][7]*g[7] + De[8][8]*g[8];


	// preenchendo vetor fonte F
	F[lm[e][0]]  = F[lm[e][0]] + Fe[0];
	F[lm[e][1]]  = F[lm[e][1]] + Fe[1];
	F[lm[e][2]]  = F[lm[e][2]] + Fe[2];
	F[lm[e][3]]  = F[lm[e][3]] + Fe[3];
	F[lm[e][4]]  = F[lm[e][4]] + Fe[4];
	F[lm[e][5]]  = F[lm[e][5]] + Fe[5];
	F[lm[e][6]]  = F[lm[e][6]] + Fe[6];
	F[lm[e][7]]  = F[lm[e][7]] + Fe[7];
	F[lm[e][8]]  = F[lm[e][8]] + Fe[8];

	F[neq] = 0.0;
}
