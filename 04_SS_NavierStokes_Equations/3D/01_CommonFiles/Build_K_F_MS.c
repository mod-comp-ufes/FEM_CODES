#include "NavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_K_F_MS(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int e;
	int nel;
	int J1, J2, J3, J4;
	double X[4], Y[4], Z[4], a[4], b[4], c[4], d[4], sixVolume, Volume, invVolume;
	double E11, E12, E13, E14, E22, E23, E24, E33, E34, E44, V1, V2, V3, V4;
	
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	
	
	
	for(e = 0; e < nel; e++){
		// Global node that composes the element
		J1 = Element[e].Vertex[0];
		J2 = Element[e].Vertex[1];
		J3 = Element[e].Vertex[2];
		J4 = Element[e].Vertex[3];

		// Nodal coordinates
		X[0] = Node[J1].x;
		X[1] = Node[J2].x;
		X[2] = Node[J3].x;
		X[3] = Node[J4].x;
		
		Y[0] = Node[J1].y;
		Y[1] = Node[J2].y;
		Y[2] = Node[J3].y;
		Y[4] = Node[J4].y;
		
		Z[0] = Node[J1].z;
		Z[1] = Node[J2].z;
		Z[2] = Node[J3].z;
		Z[3] = Node[J4].z;
		
		// Component of the form function
		a[0] = X[1]*Y[2]*Z[3] + Y[1]*Z[2]*X[3] + Z[1]*X[2]*Y[3] - Z[1]*Y[2]*X[3] - X[1]*Z[2]*Y[3] - Y[1]*X[2]*Z[3]; // a1
		a[1] = X[2]*Y[3]*Z[0] + Y[2]*Z[3]*X[0] + Z[2]*X[3]*Y[0] - Z[2]*Y[3]*X[0] - X[2]*Z[3]*Y[0] - Y[2]*X[3]*Z[0]; // a2
		a[2] = X[3]*Y[0]*Z[1] + Y[3]*Z[0]*X[1] + Z[3]*X[0]*Y[1] - Z[3]*Y[0]*X[1] - X[3]*Z[0]*Y[1] - Y[3]*X[0]*Z[1]; // a3
		a[3] = X[0]*Y[1]*Z[2] + Y[0]*Z[1]*X[2] + Z[0]*X[1]*Y[2] - Z[0]*Y[1]*X[2] - X[0]*Z[1]*Y[2] - Y[0]*X[1]*Z[2]; // a4
		
		b[0] =   Y[1]*(Z[3]-Z[2]) - Y[2]*(Z[3]-Z[1]) + Y[3]*(Z[2]-Z[1]); // b1
		b[1] = - Y[0]*(Z[3]-Z[2]) + Y[2]*(Z[3]-Z[0]) - Y[3]*(Z[2]-Z[0]); // b2
		b[2] =   Y[0]*(Z[3]-Z[1]) - Y[1]*(Z[3]-Z[0]) + Y[3]*(Z[1]-Z[0]); // b3
		b[3] = - Y[0]*(Z[2]-Z[1]) + Y[1]*(Z[2]-Z[0]) - Y[2]*(Z[1]-Z[0]); // b4
		
		c[0] = - X[1]*(Z[3]-Z[2]) + X[2]*(Z[3]-Z[1]) - X[3]*(Z[2]-Z[1]); // c1
		c[1] =   X[0]*(Z[3]-Z[2]) - X[2]*(Z[3]-Z[0]) + X[3]*(Z[2]-Z[0]); // c2
		c[2] = - X[0]*(Z[3]-Z[1]) + X[1]*(Z[3]-Z[0]) - X[3]*(Z[1]-Z[0]); // c3
		c[3] =   X[0]*(Z[2]-Z[1]) - X[1]*(Z[2]-Z[0]) + X[2]*(Z[1]-Z[0]); // c4
		
		d[0] =   X[1]*(Z[3]-Z[2]) - X[2]*(Z[3]-Z[1]) + X[3]*(Z[2]-Z[1]); // d1
		d[1] = - X[0]*(Z[3]-Z[2]) + X[2]*(Z[3]-Z[0]) - X[3]*(Z[2]-Z[0]); // d2
		d[2] =   X[0]*(Z[3]-Z[1]) - X[1]*(Z[3]-Z[0]) + X[3]*(Z[1]-Z[0]); // d3
		d[3] = - X[0]*(Z[2]-Z[1]) + X[1]*(Z[2]-Z[0]) - X[2]*(Z[1]-Z[0]); // d4
		
		// Volume
		sixVolume = (X[1]-X[0])*((Y[2]-Y[0])*(Z[3]-Z[0]) - (Y[3]-Y[0])*(Z[2]-Z[0])) + (Y[1]-Y[0])*((X[3]-X[0])*(Z[2]-Z[0]) - (X[2]-X[0])*(Z[3]-Z[0])) + (Z[1]-Z[0])*((X[2]-X[0])*(Y[3]-Y[0]) - (X[3]-X[0])*(Y[2]-Y[0]));
		Volume = sixVolume/6.0;
		invVolume = 1.0/Volume;
		
		// Stiffness matrix viscous term R_hh^1
		E11 = b[0]^2 + c[0]^2 + d[0]^2;
		E12 = b[0]*b[1] + c[0]*c[1] + d[0]*d[1];
		E13 = b[0]*b[2] + c[0]*c[2] + d[0]*d[2];
		E14 = b[0]*b[3] + c[0]*c[3] + d[0]*d[3];
		E22 = b[1]^2 + c[1]^2 + d[1]^2;
		E23 = b[1]*b[2] + c[1]*c[2] + d[1]*d[2];
		E24 = b[1]*b[3] + c[1]*c[3] + d[1]*d[3];
		E33 = b[2]^2 + c[2]^2 + d[2]^2;
		E34 = b[2]*b[3] + c[2]*c[3] + d[2]*d[3];
		E44 = b[3]^2 + c[3]^2 + d[3]^2;
		
		// Stiffness matrix viscous term R_hb^1
		V1 = 
		
	}
	
	
}
