#include "ShalowWater.h"


void F_assembly(int e, double *Fe, double (*De)[9], FemFunctionsType *FemFunctions, FemStructsType *FemStructs, int neq)
{
	int i, **lm, J1, J2, J3;
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
	if(Node[J1].hType == 1)
		g[0] = 0.0;
	else
		g[0] = FemFunctions->hpresc(Node[J1].x, Node[J1].y);
	
	if(Node[J2].hType == 1)
		g[3] = 0.0;
	else
		g[3] = FemFunctions->hpresc(Node[J2].x, Node[J2].y);
	
	if(Node[J3].hType == 1)
		g[6] = 0.0;
	else
		g[6] = FemFunctions->hpresc(Node[J3].x, Node[J3].y);

	// discharge in X
	if(Node[J1].qxType == 1)
		g[1] = 0.0;
	else
		g[1] = FemFunctions->qxpresc(Node[J1].x, Node[J1].y);
	
	if(Node[J2].qxType == 1)
		g[4] = 0.0;
	else
		g[4] = FemFunctions->qxpresc(Node[J2].x, Node[J2].y);
	
	if(Node[J3].qxType == 1)
		g[7] = 0.0;
	else
		g[7] = FemFunctions->qxpresc(Node[J3].x, Node[J3].y);

	// discharge in Y
	if(Node[J1].qyType == 1)
		g[2] = 0.0;
	else
		g[2] = FemFunctions->qypresc(Node[J1].x, Node[J1].y);
	
	if(Node[J2].qyType == 1)
		g[5] = 0.0;
	else
		g[5] = FemFunctions->qypresc(Node[J2].x, Node[J2].y);
	
	if(Node[J3].qyType == 1)
		g[8] = 0.0;
	else
		g[8] = FemFunctions->qypresc(Node[J3].x, Node[J3].y);
	
	for(i=0; i<9; i++)
	{
		Fe[i] -= (De[i][0]*g[0] + De[i][1]*g[1] + De[i][2]*g[2] + De[i][3]*g[3] + De[i][4]*g[4] + De[i][5]*g[5] + De[i][6]*g[6] + De[i][7]*g[7] + De[i][8]*g[8]);
		
		// preenchendo vetor fonte F
		F[lm[e][i]] += Fe[i];
	}

	F[neq] = 0.0;
}
