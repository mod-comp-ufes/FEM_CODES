#include "SSTransportEquation3D.h"

void eval_U_Space(ParametersType *Parameters,FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *Us)
{
	int I;
	double x, y, z;
	int nnodes = Parameters->nnodes;
	double Cst = Parameters->ConstApli;
	double xSup = Parameters->xSup;
	double *us = FemStructs->u;
	NodeType *Node = FemStructs->Node;

	for (I = 0; I < nnodes; I++){

		x = Node[I].x;
		y = Node[I].y;
		z = Node[I].z;
				
		if (Node[I].id >= 0){
			Us[I] = us[Node[I].id];
		}
		else{
			Us[I] = FemFunctions->upresc(x, y, z, Cst, xSup);
		}
	
	}
}



