#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

double Compute_Residual(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{	
	double normR = DBL_MAX;
	
	if (strcasecmp(Parameters->ComputeResidual,"Sistem")==0){
		int i;
		double *R, *KU, *U = FemStructs->u, *F = FemStructs->F;
		double neq = Parameters->neq;
		R = (double*) mycalloc("R of 'Compute_Residual'", neq+1, sizeof(double));
		KU = (double*) mycalloc("KU of 'Compute_Residual'", neq+1, sizeof(double));
	
		// multiplica K por U
		FemFunctions->ProductMatrixVector(Parameters, MatrixData, FemStructs, U, KU);
	
		for(i = 0; i < neq; i++){
			R[i] =  F[i]- KU[i];
		}
	
		normR = sqrt(ddot(neq, R, R));
	
		myfree(R);
		myfree(KU);
	
	}else if (strcasecmp(Parameters->ComputeResidual,"Equation")==0){
		int i, e, J1, J2, J3, J4, nel, nnodes;
		double betaX, betaY, betaZ, sigma;
		double X[4], Y[4], Z[4], b[4], c[4], d[4], xb, yb, zb, Fbari, Ubari;
		double *U, Ue[4], sixV, ConstGradU, gradU[3], *R;
		NodeType *Node = FemStructs->Node;
		ElementType *Element = FemStructs->Element;
		CoefFormFuncType *CFF = FemStructs->CFF;	
	
		nel = Parameters->nel;
		nnodes = Parameters->nnodes;
	
		R = (double*) mycalloc("R of 'Compute_Residual'", nel, sizeof(double));
		U = (double*) mycalloc("U of 'Compute_Residual'", nnodes, sizeof(double));
		eval_U_Space(Parameters, FemStructs, FemFunctions, U);
	
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
			Y[3] = Node[J4].y;
		
			Z[0] = Node[J1].z;
			Z[1] = Node[J2].z;
			Z[2] = Node[J3].z;
			Z[3] = Node[J4].z;
		
			// baricentric coordinate
			xb = (X[0] + X[1] + X[2] + X[3])/4.0;
			yb = (Y[0] + Y[1] + Y[2] + Y[3])/4.0;
			zb = (Z[0] + Z[1] + Z[2] + Z[3])/4.0;
			
			for(i = 0; i < 4; i++){
				b[i] = CFF[e].b[i];
				c[i] = CFF[e].c[i];
				d[i] = CFF[e].d[i];
			}	
			sixV = 6.0*CFF[e].volume;
			
			// calculates the delayed u in space
			Ue[0] = U[J1];
			Ue[1] = U[J2];
			Ue[2] = U[J3];
			Ue[3] = U[J4];
	
			// U and F in baricentric of tetrahedron
			Ubari = (Ue[0] + Ue[1] + Ue[2] + Ue[3])/4.0;
			Fbari = FemFunctions->f(xb, yb, zb, Parameters->ConstApli);
			
			ConstGradU = 1.0/sixV;
			
			gradU[0] = ConstGradU*(Ue[0]*b[0] + Ue[1]*b[1] + Ue[2]*b[2] + Ue[3]*b[3]);
			gradU[1] = ConstGradU*(Ue[0]*c[0] + Ue[1]*c[1] + Ue[2]*c[2] + Ue[3]*c[3]);
			gradU[2] = ConstGradU*(Ue[0]*d[0] + Ue[1]*d[1] + Ue[2]*d[2] + Ue[3]*d[3]);
		
			sigma = FemFunctions->sigma(xb, yb, zb);
			FemFunctions->beta(&betaX, &betaY, &betaZ, xb, yb, zb, Parameters->ConstApli, Parameters->Vmax, Parameters->xSup, Parameters->ySup);
	
			R[e] = betaX*gradU[0] + betaY*gradU[1] + betaZ*gradU[2] + sigma*Ubari - Fbari;
		
		}//end for elemento
	
		normR = sqrt(ddot(nel, R, R));
	
		myfree(U);
		myfree(R);
	
	}else if (strcasecmp(Parameters->ComputeResidual,"EquationL2")==0){
		int e, i, nel, nnodes, J1, J2, J3, J4;
		double x[4], y[4], z[4], X, Y, Z, Ue[4], Uh;
		double w[5], Xi[5][3];
		double Volume, SomaInt;
		double *U, b[4], c[4], d[4], gradU[3], sixV, ConstGradU;
		double sigma, betaX, betaY, betaZ, R, F;
		NodeType *Node = FemStructs->Node;
		ElementType *Element = FemStructs->Element;
		CoefFormFuncType *CFF = FemStructs->CFF;

		nel = Parameters->nel;
		nnodes = Parameters->nnodes;

		U = (double*) mycalloc("U of 'Compute_Residual'", nnodes, sizeof(double));
		eval_U_Space(Parameters, FemStructs, FemFunctions, U);

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
			sixV = 6*Volume;

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

			for(i = 0; i < 4; i++){
				b[i] = CFF[e].b[i];
				c[i] = CFF[e].c[i];
				d[i] = CFF[e].d[i];
			}
			
			// calculates the delayed u in space
			Ue[0] = U[J1];
			Ue[1] = U[J2];
			Ue[2] = U[J3];
			Ue[3] = U[J4];

			ConstGradU = 1.0/sixV;

			gradU[0] = ConstGradU*(Ue[0]*b[0] + Ue[1]*b[1] + Ue[2]*b[2] + Ue[3]*b[3]);
			gradU[1] = ConstGradU*(Ue[0]*c[0] + Ue[1]*c[1] + Ue[2]*c[2] + Ue[3]*c[3]);
			gradU[2] = ConstGradU*(Ue[0]*d[0] + Ue[1]*d[1] + Ue[2]*d[2] + Ue[3]*d[3]);
		
			for(i = 0; i < 5; i++){
				// Transformed coordinates
				X = (1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2])*x[0] + Xi[i][0]*x[1] + Xi[i][1]*x[2] + Xi[i][2]*x[3];
				Y = (1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2])*y[0] + Xi[i][0]*y[1] + Xi[i][1]*y[2] + Xi[i][2]*y[3];
				Z = (1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2])*z[0] + Xi[i][0]*z[1] + Xi[i][1]*z[2] + Xi[i][2]*z[3];

				sigma = FemFunctions->sigma(X, Y, Z);
				FemFunctions->beta(&betaX, &betaY, &betaZ, X, Y, Z, Parameters->ConstApli, Parameters->Vmax, Parameters->xSup, Parameters->ySup);

				// Source term
				F = FemFunctions->f(X, Y, Z, Parameters->ConstApli);

				// Interpolated solucion
				Uh = Ue[0]*(1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2]) + Ue[1]*Xi[i][0] + Ue[2]*Xi[i][1] + Ue[3]*Xi[i][2];

				// Local resídue
				R = betaX*gradU[0] + betaY*gradU[1] + betaZ*gradU[2] + sigma*Uh - F;

				SomaInt = SomaInt + Volume*pow(R,2)*w[i];
			}

		}

		normR = sqrt(SomaInt);
		
	}else{
		printf("Compute residue is not defined correctly!\n");
		exit(1);
	}
	
	return normR;
}
