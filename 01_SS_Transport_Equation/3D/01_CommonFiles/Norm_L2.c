#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"


void Norm_L2 (double *U, ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions){
	
	int e, i, nel, J1, J2, J3, J4;
	double x[4], y[4], z[4], X, Y, Z, Ue = 0.0, Uh;
	double w[5], Xi[5][3];
	double Volume, SomaInt, Error;//, normUUh, detJ;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	CoefFormFuncType *CFF = FemStructs->CFF;
	
/*	double **J;
	J = (double**) mycalloc("J line of 'Norm_L2'", 3, sizeof(double*));
	for (i = 0; i < 3; i++){
		J[i] = (double*) mycalloc("J colum of 'Norm_L2'", 3, sizeof(double));
	}*/
	
	nel = Parameters->nel;
	
	// Fills in with the gaussian quadrature weights and points
	Xi[0][0] = Xi[0][1] = Xi[0][2] = 1.0/6.0;  //ponto 1
	Xi[1][0] = Xi[1][1] = 1.0/6.0; Xi[1][2] = 1.0/2.0;  //ponto 2
	Xi[2][0] = 1.0/6.0;  Xi[2][1] = 1.0/2.0; Xi[2][2] = 1.0/6.0; //ponto 3
	Xi[3][0] = 1.0/2.0;  Xi[3][1] = Xi[3][2] = 1.0/6.0;  //ponto 4
	Xi[4][0] = Xi[4][1] = Xi[4][2] = 1.0/4.0;  //ponto 5
	
	w[0] = w[1] = w[2] = w[3] = 9.0/20.0;  
	w[4] = -16.0/20.0;  
	
	SomaInt = 0.0;

	for(e = 0; e < nel; e++){
		Volume = CFF[e].volume;
		
		J1 = Element[e].Vertex[0];
		J2 = Element[e].Vertex[1];
		J3 = Element[e].Vertex[2];
		J4 = Element[e].Vertex[3];
		
		// Nodal coordinates
		x[0] = Node[J1].x;
		x[1] = Node[J2].x;
		x[2] = Node[J3].x;
		x[3] = Node[J4].x;
		
		y[0] = Node[J1].y;
		y[1] = Node[J2].y;
		y[2] = Node[J3].y;
		y[3] = Node[J4].y;
		
		z[0] = Node[J1].z;
		z[1] = Node[J2].z;
		z[2] = Node[J3].z;
		z[3] = Node[J4].z;
		
		// Jacobian Matriz
	/*	J[0][0] = x[1] - x[0]; J[0][1] = y[1] - y[0]; J[0][2] = z[1] - z[0];
		J[1][0] = x[2] - x[0]; J[1][1] = y[2] - y[0]; J[1][2] = z[2] - z[0];
		J[2][0] = x[3] - x[0]; J[2][1] = y[3] - y[0]; J[2][2] = z[3] - z[0];*/
		
		//detJ = determinant(J, 3);
		
		for(i = 0; i < 5; i++){
			// Transformed coordinates
			X = (1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2])*x[0] + Xi[i][0]*x[1] + Xi[i][1]*x[2] + Xi[i][2]*x[3];
			Y = (1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2])*y[0] + Xi[i][0]*y[1] + Xi[i][1]*y[2] + Xi[i][2]*y[3];
			Z = (1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2])*z[0] + Xi[i][0]*z[1] + Xi[i][1]*z[2] + Xi[i][2]*z[3];
			
			// Exact solucion
			Ue = FemFunctions->ExactSolution(X,Y,Z,Parameters->ConstApli);
			
			// Interpolated solucion
			Uh = U[J1]*(1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2]) + U[J2]*Xi[i][0] + U[J3]*Xi[i][1] + U[J4]*Xi[i][2];
			
			// Local Error
			Error = Ue - Uh;
			
			SomaInt = SomaInt + Volume*pow(Error,2)*w[i];
		}
		
	}
	
	//printf("\nNORM L2 = %lf \n", normUUh);
	
	Parameters->NormL2 = sqrt(SomaInt);
	
/*	for(i = 0; i < 3; i++){
		free(J[i]);
	}
	free(J);*/
}
