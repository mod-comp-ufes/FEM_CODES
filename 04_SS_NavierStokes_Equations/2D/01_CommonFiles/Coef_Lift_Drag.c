#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Coef_Lift_Drag(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
//
// Esta funccao calcula os coeficientes C_Lift e C_Drag do problema do Cylindro para Re = 20
//
	int E, J1, J2, J3, nel;
	double Re;
	double twoArea, invArea, Area, visc, rho; 
	double third=1.0/3.0, sixth = 1.0/6.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], ux1, ux2, ux3, uy1, uy2, uy3, ux, uy;
	double C1, C2, C3;
	//double Nv[6], Ndeltav[6], Kv[6], Ksv[6], Gv[6], Gdeltav[6], GTv[3], Cv[3], Nphiv[3], Gphiv[3], Res[9];
	//double M11, M13, Mdelta11, Mdelta31, Mdelta51;
	double K11, K13, K15, K22, K24, K26, K33, K35, K44, K46, K55, K66; 
	double G11, G21, G31, G41, G51, G61, N11, N13, N15;
	double p1, p2, p3;
	double Clift, Cdrag, P1, P2, DP;
	
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	nel = Parameters->nel;

	double *U = (double*) mycalloc("U of 'Build_K_F_SUPG_PSPG'", 3*Parameters->nnodes, sizeof(double));
	eval_U(Parameters,FemStructs,FemFunctions,U);

	setzeros(Parameters,MatrixData);
	
	Clift = 0.0;
	Cdrag = 0.0;
	P1 = 0.0;
	P2 = 0.0;
	DP = 0.0;

	for (E=0; E<nel; E++){
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

		twoArea =  fabs(x21*y31 - x13*y12);
		Area = 0.5*twoArea;
		invArea = 1.0/Area;		
				
	 	//x Velocity  
		ux1 = U[3*J1];
		ux2 = U[3*J2];
		ux3 = U[3*J3];

		//y Velocity
		uy1 = U[3*J1+1];
		uy2 = U[3*J2+1];
		uy3 = U[3*J3+1];
		
		//pression 
		p1 = U[3*J1+2];
		p2 = U[3*J2+2];
		p3 = U[3*J3+2];

		//*****************u at barycentre of the element***************
		ux = third*( ux1 + ux2 + ux3 );
		uy = third*( uy1 + uy2 + uy3 );
		
		//*****************Coefficients C_i******************************
		C1 = 0.5*invArea*( ux*y23 + uy*x32 );
		C2 = 0.5*invArea*( ux*y31 + uy*x13 );
		C3 = 0.5*invArea*( ux*y12 + uy*x21 );

		//*************Calculation of tau_SUPG=tau_PSPG=tau*********
		Re = Parameters->ReynoldsNumber;		
		visc =1./Re;
		rho = 1.;		

		//*************************Matrix N****************************
		N11 = rho*third*Area*C1;		
		N13 = rho*third*Area*C2;		
		N15 = rho*third*Area*C3;

		//********* Calculo de \nu(\nabla u,\nabla w) *****************
		K11 = visc*0.25*invArea*( y23*y23 + x32*x32 );
		K13 = visc*0.25*invArea*( y23*y31 + x32*x13 );
		K15 = visc*0.25*invArea*( y23*y12 + x32*x21 );
		
		K22 = visc*0.25*invArea*( y23*y23 + x32*x32 );
		K24 = visc*0.25*invArea*( y23*y31 + x32*x13 );
		K26 = visc*0.25*invArea*( y23*y12 + x32*x21 );

		K33 = visc*0.25*invArea*( y31*y31 + x13*x13 );
		K35 = visc*0.25*invArea*( y31*y12 + x13*x21 );
		
		K44 = visc*0.25*invArea*( y31*y31 + x13*x13 );
		K46 = visc*0.25*invArea*( y31*y12 + x13*x21 );
	
		K55 = visc*0.25*invArea*( y12*y12 + x21*x21 );		
			
		K66 = visc*0.25*invArea*( y12*y12 + x21*x21);	

		//************************Matrix G****************************
		G11 = sixth*y23;
		G21 = sixth*x32;
		G31 = sixth*y31;
		G41 = sixth*x13;
		G51 = sixth*y12;
		G61 = sixth*x21;

		//( 3*N11*ux1 + 3*N13*ux2 + 3*N15*uy3 + 3*N11*uy1 + 3*N13*uy2 + 3*N15*uy3); 
		
		//( (K11+K13+K15)*ux1 + (K13+K33+K35)*ux2 + (K15+K35+K55)*ux3 + (K22+K24+K26)*uy1 + (K24+K44+K46)*uy2 + (K26+K46+K66)*uy3); 
			
		//(G11+G21+G31+G41+G51+G61)*(p1+p2+p3);
					
		// Apenas nó 1 sobre o cylindro.  
		if (Node[J1].v1Type == -2 && Node[J2].v1Type == 1 && Node[J3].v1Type == 1){
			Cdrag = Cdrag - (N11*ux1 + N13*ux2 + N15*ux3) - (K11*ux1 + K13*ux2 + K15*ux3) + G11*(p1+p2+p3);
			Clift = Clift - (N11*uy1 + N13*uy2 + N15*uy3) - (K22*uy1 + K24*uy2 + K66*uy3) + G21*(p1+p2+p3);
//			printf(" \n Nó 1 ux1 = %f, ux2 = %f, ux3 = %f, uy1 = %f, uy2 = %f, uy3 = %f ", ux1, ux2, ux3, uy1, uy2, uy3);
				
		}
		// Apenas nó 2 sobre o cylindro
		else if (Node[J1].v1Type == 1 && Node[J2].v1Type == -2 && Node[J3].v1Type == 1){
			Cdrag = Cdrag - (N11*ux1 + N13*ux2 + N15*ux3) - (K13*ux1 + K33*ux2 + K35*ux3) + G31*(p1+p2+p3);
			Clift = Clift - (N11*uy1 + N13*uy2 + N15*uy3) - (K24*uy1 + K44*uy2 + K46*uy3) + G41*(p1+p2+p3);
//			printf(" \n Nó 2 ux1 = %f, ux2 = %f, ux3 = %f, uy1 = %f, uy2 = %f, uy3 = %f ", ux1, ux2, ux3, uy1, uy2, uy3);
		}
		// Apenas nó 3 sobre o cylindro		
		else if (Node[J1].v1Type == 1 && Node[J2].v1Type == 1 && Node[J3].v1Type == -2){
			Cdrag = Cdrag - (N11*ux1 + N13*ux2 + N15*ux3) - (K15*ux1 + K35*ux2 + K55*ux3) + G51*(p1+p2+p3);
			Clift = Clift - (N11*uy1 + N13*uy2 + N15*uy3) - (K26*uy1 + K46*uy2 + K66*uy3) + G61*(p1+p2+p3);
//			printf(" \n Nó 3 ux1 = %f, ux2 = %f, ux3 = %f, uy1 = %f, uy2 = %f, uy3 = %f ", ux1, ux2, ux3, uy1, uy2, uy3);
		}
		// Apenas nó 1 e 2 sobre o cylindro
		else if (Node[J1].v1Type == -2 && Node[J2].v1Type == -2 && Node[J3].v1Type == 1){
			Cdrag = Cdrag - (N11*ux1 + N13*ux2 + N15*ux3) - (K11*ux1 + K13*ux2 + K15*ux3) + G11*(p1+p2+p3) - (N11*ux1 + N13*ux2 + N15*ux3) - (K13*ux1 + K33*ux2 + K35*ux3) + G31*(p1+p2+p3);
			Clift = Clift - (N11*uy1 + N13*uy2 + N15*uy3) - (K22*uy1 + K24*uy2 + K66*uy3) + G21*(p1+p2+p3) - (N11*uy1 + N13*uy2 + N15*uy3) - (K24*uy1 + K44*uy2 + K46*uy3) + G41*(p1+p2+p3);
		}
		// Apenas nó 1 e 3 sobre o cylindro
		else if (Node[J1].v1Type == -2 && Node[J2].v1Type == 1 && Node[J3].v1Type == -2){
			Cdrag = Cdrag - (N11*ux1 + N13*ux2 + N15*ux3) - (K11*ux1 + K13*ux2 + K15*ux3) + G11*(p1+p2+p3) - (N11*ux1 + N13*ux2 + N15*ux3) - (K15*ux1 + K35*ux2 + K55*ux3) + G51*(p1+p2+p3);
			Clift = Clift - (N11*uy1 + N13*uy2 + N15*uy3) - (K22*uy1 + K24*uy2 + K66*uy3) + G21*(p1+p2+p3) - (N11*uy1 + N13*uy2 + N15*uy3) - (K26*uy1 + K46*uy2 + K66*uy3) + G61*(p1+p2+p3);
		}
		// Apenas nó 2 e 3 sobre o cylindro
		else if (Node[J1].v1Type == 1 && Node[J2].v1Type == -2 && Node[J3].v1Type == -2){
			Cdrag = Cdrag - (N11*ux1 + N13*ux2 + N15*ux3) - (K13*ux1 + K33*ux2 + K35*ux3) + G31*(p1+p2+p3) - (N11*ux1 + N13*ux2 + N15*ux3) - (K15*ux1 + K35*ux2 + K55*ux3) + G51*(p1+p2+p3);
			Clift = Clift - (N11*uy1 + N13*uy2 + N15*uy3) - (K24*uy1 + K44*uy2 + K46*uy3) + G41*(p1+p2+p3) - (N11*uy1 + N13*uy2 + N15*uy3) - (K26*uy1 + K46*uy2 + K66*uy3) + G61*(p1+p2+p3);
		}
		
		if(fabs(Y[0]-0.2)<=1e-15 || fabs(Y[1]-0.2)<=1e-15 || fabs(Y[2]-0.2)<=1e-15){
			if(fabs(X[0]-0.15)<=1e-15){
				P1 = p1;
			//	printf("\n Passei aqui 01, P1 = %f \n", P1);				
			}
			else if(fabs(X[1]-0.15)<=1e-15){
				P1 = p2;
			//	printf("\n Passei aqui 02, P1 = %f \n", P1);				
			}
			     else if(fabs(Y[2]-0.15)<=1e-15){
				P1 = p3;
			//	printf("\n Passei aqui 03, P1 = %f \n", P1);				
			}
			if(fabs(X[0]-0.25)<=1e-15){
				P2 = p1;
			//	printf("\n Passei aqui 04, P2 = %f \n", P2);				
			}
			else if(fabs(X[1]-0.25)<=1e-15){
				     P2 = p2;
			//	     printf("\n Passei aqui 05, P2 = %f \n", P2);				
			}
			     else if(fabs(Y[2]-0.25)<=1e-15){
				     P2 = p3;
			//	     printf("\n Passei aqui 06, P2 = %f \n", P2);				
			}
			
			DP = P1-P2;		
		}
	
	}
	printf("\n Cdrag = %f, Clift = %f, Delta P = %f", 500*Cdrag, 500*Clift, DP);
	printf("\n\n");
	return 0;

}

















