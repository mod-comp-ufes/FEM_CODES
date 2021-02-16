#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_Galerkin(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int J, E, J1, J2, J3;
	int neq, nel;
	double kx, ky, Be_x, Be_y, gamma, TwoA, invArea, Area;
	double third=1.0/3.0, sixth = 1.0/6.0, twelfth= 1.0/12.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Beta[2], ke[3][3];
	double fe1, fe2, fe3, Fe1, Fe2, Fe3;
	double D12, D13, D23, C11, C12, C13, R11, R12, XB, YB;
	double *U;
	double *F = FemStructs->F;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	neq = Parameters->neq;
	nel = Parameters->nel;
	kx = FemFunctions->Condutivity();
	ky = FemFunctions->Condutivity();
	gamma = FemFunctions->Reaction();

	setzeros(Parameters,MatrixData);
	U = mycalloc("U of 'Build_Galerkin'", Parameters->nnodes, sizeof(double));
	eval_U(Parameters,FemStructs,FemFunctions,U);
	
	for (J=0;J<neq;J++)
		F[J] = 0;

	for (E=0; E<nel; E++)
	{
		J1 = Element[E].Vertex[0];
		J2 = Element[E].Vertex[1];
		J3 = Element[E].Vertex[2];

		X[0] = Node[J1].x;
		X[1] = Node[J2].x;
		X[2] = Node[J3].x;
		Y[0] = Node[J1].y;
		Y[1] = Node[J2].y;
		Y[2] = Node[J3].y;

		y23 = Y[1] - Y[2];
		y31 = Y[2] - Y[0];
		y12 = Y[0] - Y[1];

		x32 = X[2] - X[1];
		x13 = X[0] - X[2];
		x21 = X[1] - X[0];

		TwoA = x21*y31 - x13*y12;
		Area = 0.5*TwoA;
		invArea = 1.0/Area;

		XB = third*(X[0]+X[1]+X[2]);
		YB = third*(Y[0]+Y[1]+Y[2]);

		//Beta velocity calculation
		FemFunctions->Velocity(XB, YB, Beta);

		Be_x = Beta[0];
		Be_y = Beta[1];

		//Calculation of Font
		fe1 = FemFunctions->Font(X[0], Y[0], kx, gamma, Be_x, Be_y);
		fe2 = FemFunctions->Font(X[1], Y[1], kx, gamma, Be_x, Be_y);
		fe3 = FemFunctions->Font(X[2], Y[2], kx, gamma, Be_x, Be_y);

		// Galerkin diffusion matrix
		D12 = 0.25*invArea*((kx * y23 * y31) + (ky * x32 * x13));
		D13 = 0.25*invArea*((kx * y23 * y12) + (ky * x32 * x21));
		D23 = 0.25*invArea*((kx * y31 * y12) + (ky * x13 * x21));

		// Galerkin Convection matrix
		C11 = sixth*(Be_x*y23 + Be_y*x32);
		C12 = sixth*(Be_x*y31 + Be_y*x13);
		C13 = sixth*(Be_x*y12 + Be_y*x21);

		// Galerkin reaction matrix
		R11 = sixth*(gamma*Area);
		R12 = twelfth*(gamma * Area);

		/********************************************************************************************************/

	    // Matrix ke
		ke[0][0] = -(D12 + D13) + C11 + R11;
		ke[0][1] =      D12     + C12 + R12;
		ke[0][2] = 		D13     + C13 + R12;
		ke[1][0] = 		D12     + C11 + R12;
		ke[1][1] = -(D12 + D23) + C12 + R11;
		ke[1][2] = 		D23     + C13 + R12;
		ke[2][0] = 		D13     + C11 + R12;
		ke[2][1] = 		D23     + C12 + R12;
		ke[2][2] = -(D13 + D23) + C13 + R11;

        // PAREI AQUI
		// Boundary conditions on vector F and Fill Matrix Diagonal 
		Fe1 = Area*twelfth*(2*fe1 +  fe2  +  fe3 );
		Fe2 = Area*twelfth*(  fe1 + 2*fe2 +  fe3 );
		Fe3 = Area*twelfth*(  fe1 +  fe2  + 2*fe3);

		F_assembly(Node,J1,J2,J3,Fe1,Fe2,Fe3, FemStructs->F, ke, X, Y, FemFunctions);

		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, E, ke);

	}

	free(U);

	return 0;
}
