#include "SSNavierStokesEquations.h"

void eval_U_dU(ParametersType *Parameters,FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *U,double *dU)
{
	int I;
	double x, y;
	int nnodes = Parameters->nnodes;
	double *u = FemStructs->u;
	double *du = FemStructs->du;
	NodeType *Node = FemStructs->Node;

	for (I=0; I<nnodes; I++){

		x = Node[I].x;
		y = Node[I].y;

		// Velocity x
		if (Node[I].id[0] >= 0){
			U[3*I] = u[Node[I].id[0]];
			dU[3*I] = du[Node[I].id[0]];
		}
		else{
			U[3*I] = FemFunctions->v1presc(x,y);
			dU[3*I] = 0.0;
		}

		// Velocity y
		if (Node[I].id[1] >= 0){
			U[3*I+1] = u[Node[I].id[1]];
			dU[3*I+1] = du[Node[I].id[1]];
		}
		else{
			U[3*I+1] = FemFunctions->v2presc(x, y);
			dU[3*I+1] = 0.0;
		}

		// Press
		if (Node[I].id[2] >= 0){
			U[3*I+2] = u[Node[I].id[2]];
		}
		else{
			U[3*I+2] = FemFunctions->ppresc(x, y);
		}

	}
}



