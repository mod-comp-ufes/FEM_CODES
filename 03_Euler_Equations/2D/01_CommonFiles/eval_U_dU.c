#include "EulerEquations.h"

void eval_U_dU(ParametersType *Parameters,FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *U,double *dU)
{
	int I;
	double x, y;
	int nnodes = Parameters->nnodes;
	double *u = FemStructs->u;
	double *du = FemStructs->du;
	//double theta;
	NodeType *Node = FemStructs->Node;

	for (I=0; I<nnodes; I++){

		x = Node[I].x;
		y = Node[I].y;
				
		//if(Node[I].v1Type < 0){
		//	theta = FemFunctions->BC_theta(x,y);
		//}

		//density
		if (Node[I].id[0] >= 0){
			U[4*I] = u[Node[I].id[0]];
			dU[4*I] = du[Node[I].id[0]];
		}
		else{
			U[4*I] = FemFunctions->rhopresc(x,y);
			dU[4*I] = 0.0;
		}
		/*
		// Velocity x
		if (Node[I].id[1] >= 0){
			if(Node[I].v1Type < 0){
				U[4*I+1] = cos( theta ) * u[Node[I].id[1]];
	   			dU[4*I+1] = cos( theta ) * du[Node[I].id[1]];
			}
			else{
	   			U[4*I+1] = u[Node[I].id[1]];
				dU[4*I+1] = du[Node[I].id[1]];
			}
		}
		else{
			U[4*I+1] = FemFunctions->v1presc(x, y);
			dU[4*I+1] = 0.0;
		}

		// Velocity y
		if (Node[I].id[2] >= 0){
			U[4*I+2] = u[Node[I].id[2]];
			dU[4*I+2] = du[Node[I].id[2]];
		}
		else if(Node[I].v1Type < 0){
			U[4*I+2] = sin( theta ) * u[Node[I].id[1]];
			dU[4*I+2] = sin( theta ) * du[Node[I].id[1]];
		}
		else{
			U[4*I+2] = FemFunctions->v2presc(x, y);
			dU[4*I+2] = 0.0;
		}
		
		*/
		// Velocity x
		if (Node[I].id[1] >= 0){
			U[4*I+1] = u[Node[I].id[1]];
			dU[4*I+1] = du[Node[I].id[1]];
		}
		else{
			U[4*I+1] = FemFunctions->v1presc(x, y);
			dU[4*I+1] = 0.0;
		}

		// Velocity y
		if (Node[I].id[2] >= 0){
			U[4*I+2] = u[Node[I].id[2]];
			dU[4*I+2] = du[Node[I].id[2]];
		}
		else{
			U[4*I+2] = FemFunctions->v2presc(x, y);
			dU[4*I+2] = 0.0;
		}
		
		// Energy
		if (Node[I].id[3] >= 0){
			U[4*I+3] = u[Node[I].id[3]];
			dU[4*I+3] = du[Node[I].id[3]];
		}
		else{
			U[4*I+3] = FemFunctions->epresc(x, y);
			dU[4*I+3] = 0.0;
		}

	}
}



