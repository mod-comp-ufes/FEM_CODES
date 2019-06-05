#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

static int varglobal = 0;

int Build_K_F_SUPG_PSPG(ParametersType *Parameters, MatrixDataType *MatrixData, NodeType *Node, ElementType *Element, double *F, double *u, int **lm,
void (*assembly)(ParametersType *, int , int,  double (*)[9], MatrixDataType *, int **), void (*end_assembly)(ParametersType *, int, MatrixDataType *), double (*ppresc)(double, double), double (*v1presc)(double, double), double (*v2presc)(double, double))
//, double (*Viscosity)(),double (*Rho)() , double (*U_ref)(), double (*L_ref)())	
{
	int i, j, neq, nel;//, op;
	int J, E, J1, J2, J3, lm_aux[9];
	//double Uref, Lref;
	double Re, x, y;
	double TwoA, invArea, Area, tau, h, visc, rho, invrho, nu, unorm, auxtau; 
	double soma=0.0, third=1.0/3.0, sixth = 1.0/6.0, twelfth = 1.0/12.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Ke[9][9], Keb[9][9], Fe[9], ux1, ux2, ux3, uy1, uy2, uy3, ux, uy, duxdx, duxdy, duydx, duydy;
	double fx1, fx2, fx3, fy1, fy2, fy3, fxB, fyB, f1, f2, fdelta1, fdelta2, fdelta3, fdelta4, fdelta5, fdelta6, fphi1, fphi2, fphi3;
	double C1, C2, C3;
	double Nv[6], Ndeltav[6], Kv[6], Gv[6], Gdeltav[6], GTv[3], Nphiv[3], Gphiv[3], RQres[9];
	//double M11, M13, Mdelta11, Mdelta31, Mdelta51;
	double K11, K12, K13, K14, K15, K16, K22, K23, K24, K25, K26, K33, K34, K35, K36, K44, K45, K46, K55, K56, K66; 
	double GT11, GT12, GT13, GT14, GT15, GT16, N11, N13, N15;
	double Ndelta11, Ndelta13, Ndelta15, Ndelta22, Ndelta24, Ndelta26, Ndelta33, Ndelta35, Ndelta44, Ndelta46, Ndelta55, Ndelta66;
	double Gdelta11, Gdelta12, Gdelta13, Gdelta21, Gdelta22, Gdelta23, Gdelta31, Gdelta32, Gdelta33, Gdelta41, Gdelta42, Gdelta43;
	double Gdelta51, Gdelta52, Gdelta53, Gdelta61, Gdelta62, Gdelta63;
	//double Mphi11, Mphi12, Mphi21, Mphi22, Mphi31, Mphi32;
	double Gphi11, Gphi12, Gphi13, Gphi22, Gphi23, Gphi33;
	double Np11, Np12, Np21, Np22, Np31, Np32, Np41, Np42;
	double Npdelta11, Npdelta12, Npdelta21, Npdelta22, Npdelta31, Npdelta32, Npdelta41, Npdelta42, Npdelta51, Npdelta52, Npdelta61, Npdelta62; 
	double Nppdelta11, Nppdelta12, Nppdelta21, Nppdelta22, Nppdelta31, Nppdelta32, Nppdelta41, Nppdelta42, Nppdelta51, Nppdelta52, Nppdelta61, Nppdelta62;
	double Npphi11, Npphi12, Npphi21, Npphi22, Npphi31, Npphi32;
	double p1, p2, p3, p, dpdx, dpdy;
	double *vet_bound;
	
	vet_bound = (double*) mycalloc("vet_bound of 'Build_M_K_F_DD_Estatico'", 12, sizeof(double)); 

	nel = Parameters->nel;
	neq = Parameters->neq;
	
	for (J=0; J<neq+1; J++)
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

		TwoA =  fabs(x21*y31 - x13*y12);
		Area = 0.5*TwoA;
		invArea = 1.0/Area;		
				
	 

		//(u_x,u_y) e p of the node 0
		if (Node[J1].id[0]>=0){
	   		ux1 = u[Node[J1].id[0]];
	   	}
	   	else{
	   		ux1 = (*v1presc)(X[0],Y[0]);
	   	}
		if (Node[J1].id[1]>=0){
	   		uy1 = u[Node[J1].id[1]];
	   	}
	   	else{
	   		uy1 = (*v2presc)(X[0],Y[0]);
	   	}
		if (Node[J1].id[2]>=0){
	   		p1 = u[Node[J1].id[2]];
	   	}
	   	else{
	   		p1 = (*ppresc)(X[0],Y[0]);
	   	}
		//(u_x,u_y) e p of the node 1
		if (Node[J2].id[0]>=0){
	   		ux2 = u[Node[J2].id[0]];
	   	}
	   	else{
	   		ux2 = (*v1presc)(X[1],Y[1]);
	   	}
		if (Node[J2].id[1]>=0){
	   		uy2 = u[Node[J2].id[1]];
	   	}
	   	else{
	   		uy2 = (*v2presc)(X[1],Y[1]);
	   	}
		if (Node[J2].id[2]>=0){
	   		p2 = u[Node[J2].id[2]];
	   	}
	   	else{
	   		p2 = (*ppresc)(X[1],Y[1]);
	   	}
		//(u_x,u_y) e p of the node 2
		if (Node[J3].id[0]>=0){
	   		ux3 = u[Node[J3].id[0]];
	   	}
	   	else{
	   		ux3 = (*v1presc)(X[2],Y[2]);
	   	}
		if (Node[J3].id[1]>=0){
	   		uy3 = u[Node[J3].id[1]];
	   	}
	   	else{
	   		uy3 = (*v2presc)(X[2],Y[2]);
	   	}
		if (Node[J3].id[2]>=0){
	   		p3 = u[Node[J3].id[2]];
	   	}
	   	else{
	   		p3 = (*ppresc)(X[2],Y[2]);
	   	}
		
		//*****************u at barycentre of the element***************
		ux = third*( ux1 + ux2 + ux3 );
		uy = third*( uy1 + uy2 + uy3 );
		
		//*****************************Gradient of u********************
		duxdx = 0.5*invArea*( ux1*y23 + ux2*y31 + ux3*y12 ); 
		duxdy = 0.5*invArea*( ux1*x32 + ux2*x13 + ux3*x21 );
		duydx = 0.5*invArea*( uy1*y23 + uy2*y31 + uy3*y12 );
		duydy = 0.5*invArea*( uy1*x32 + uy2*x13 + uy3*x21 );
			
		//*****************p at barycentre of the element***************
		p = third*( p1 + p2 + p3 );

		//**************************Gradient of p**********************
		dpdx = 0.5*invArea*( p1*y23 + p2*y31 + p3*y12 );
		dpdy = 0.5*invArea*( p1*x32 + p2*x13 + p3*x21 );

		//*****************Coefficients C_i******************************
		C1 = 0.5*invArea*( ux*y23 + uy*x32 );
		C2 = 0.5*invArea*( ux*y31 + uy*x13 );
		C3 = 0.5*invArea*( ux*y12 + uy*x21 );

		//*************Calculation of tau_SUPG=tau_PSPG=tau*********
		Re = Parameters->ReynoldsNumber;		
		visc = 1./Re;
		rho = 1.;		
		//visc = 1./sqrt(Re);
		//rho = sqrt(Re);		

		//visc = (*Viscosity)(); 
		//rho = (*Rho)();
		//Uref = (*U_ref)();
		//Lref = (*L_ref)();
		invrho = 1.0/rho;
		nu = visc/rho;
		h = sqrt( 4*Area/PI );		
		unorm = sqrt( ux*ux + uy*uy );
		auxtau = ( 2*unorm/h )*( 2*unorm/h ) + 9*( 4*nu/(h*h) )*( 4*nu/(h*h) );
		auxtau = sqrt( auxtau);
		tau = 1.0/auxtau;
		//tau=tau_SUPG=tau_PSPG
		
	
		//************Matrices of the Galerkin formulation*************

		//*************************Matrix M****************************
		//M11 = rho*2*Area*twelfth;
		//M13 = rho*Area*twelfth; 
		
		//*************************Matrix N****************************
		N11 = rho*third*Area*C1;		
		N13 = rho*third*Area*C2;		
		N15 = rho*third*Area*C3;		

		//************************Matrix K******************************
		K11 = visc*0.25*invArea*( 2*y23*y23 + x32*x32 );
		K12 = visc*0.25*invArea*( x32*y23 );
		K13 = visc*0.25*invArea*( 2*y23*y31 + x32*x13 );
		K14 = visc*0.25*invArea*( x32*y31 );
		K15 = visc*0.25*invArea*( 2*y23*y12 + x32*x21 );
		K16 = visc*0.25*invArea*( x32*y12 );

		K22 = visc*0.25*invArea*( y23*y23 + 2*x32*x32 );
		K23 = visc*0.25*invArea*( y23*x13 );
		K24 = visc*0.25*invArea*( y23*y31 + 2*x32*x13 );
		K25 = visc*0.25*invArea*( y23*x21 );
		K26 = visc*0.25*invArea*( y23*y12 + 2*x32*x21 );

		K33 = visc*0.25*invArea*( 2*y31*y31 + x13*x13 );
		K34 = visc*0.25*invArea*( x13*y31 );
		K35 = visc*0.25*invArea*( 2*y31*y12 + x13*x21 );
		K36 = visc*0.25*invArea*( x13*y12 );

		K44 = visc*0.25*invArea*( y31*y31 + 2*x13*x13 );
		K45 = visc*0.25*invArea*( y31*x21 );
		K46 = visc*0.25*invArea*( y31*y12 + 2*x13*x21 );

		K55 = visc*0.25*invArea*( 2*y12*y12 + x21*x21 );		
		K56 = visc*0.25*invArea*( x21*y12 );

		K66 = visc*0.25*invArea*( y12*y12 + 2*x21*x21);

		//************************Matrix G^T****************************
		GT11 = sixth*y23;
		GT12 = sixth*x32;
		GT13 = sixth*y31;
		GT14 = sixth*x13;
		GT15 = sixth*y12;
		GT16 = sixth*x21;
		
		//***************Matrices of the SUPG formulation****************

		//*******************Matrix Mdelta=tau*N^T***********************
		//Mdelta11 = tau*N11;
		//Mdelta31 = tau*N13;
		//Mdelta51 = tau*N15;

		//**********************Matrix Ndelta****************************
		Ndelta11 = rho*tau*Area*C1*C1;
		Ndelta13 = rho*tau*Area*C1*C2;
		Ndelta15 = rho*tau*Area*C1*C3;
	
		Ndelta22 = Ndelta11;
		Ndelta24 = Ndelta13;
		Ndelta26 = Ndelta15;
	
		Ndelta33 = rho*tau*Area*C2*C2;
		Ndelta35 = rho*tau*Area*C2*C3;
	
		Ndelta44 = Ndelta33;
		Ndelta46 = Ndelta35;
	
		Ndelta55 = rho*tau*Area*C3*C3;
	
		Ndelta66 = Ndelta55;
		
		//**********************Matrix Gdelta****************************
		Gdelta11 = 0.5*tau*C1*y23;
		Gdelta12 = 0.5*tau*C1*y31;
		Gdelta13 = 0.5*tau*C1*y12;

		Gdelta21 = 0.5*tau*C1*x32;
		Gdelta22 = 0.5*tau*C1*x13;
		Gdelta23 = 0.5*tau*C1*x21;

		Gdelta31 = 0.5*tau*C2*y23;
		Gdelta32 = 0.5*tau*C2*y31;
		Gdelta33 = 0.5*tau*C2*y12;

		Gdelta41 = 0.5*tau*C2*x32;
		Gdelta42 = 0.5*tau*C2*x13;
		Gdelta43 = 0.5*tau*C2*x21;

		Gdelta51 = 0.5*tau*C3*y23;
		Gdelta52 = 0.5*tau*C3*y31;
		Gdelta53 = 0.5*tau*C3*y12;

		Gdelta61 = 0.5*tau*C3*x32;
		Gdelta62 = 0.5*tau*C3*x13;
		Gdelta63 = 0.5*tau*C3*x21;

		//***************Matrices of the PSPG formulation****************

		//**********************Matrix Mphi******************************
		//Mphi11 = tau*sixth*y23;
		//Mphi12 = tau*sixth*x32;
	
		//Mphi21 = tau*sixth*y31;
		//Mphi22 = tau*sixth*x13;
	
		//Mphi31 = tau*sixth*y12;
		//Mphi32 = tau*sixth*x21;

		//******************Matrix Nphi=Gdelta^T**************************

		//**********************Matrix Gphi*******************************
		Gphi11 = invrho*0.25*tau*invArea*( y23*y23 + x32*x32 );
		Gphi12 = invrho*0.25*tau*invArea*( y23*y31 + x32*x13 );
		Gphi13 = invrho*0.25*tau*invArea*( y23*y12 + x32*x21 );
		
		Gphi22 = invrho*0.25*tau*invArea*( y31*y31 + x13*x13 );
		Gphi23 = invrho*0.25*tau*invArea*( y31*y12 + x13*x21 );
		
		Gphi33 = invrho*0.25*tau*invArea*( y12*y12 + x21*x21 );	
		
		//*********************Incremental matrices***********************
		
		//*************************Matrix Np*******************************
		Np11 = rho*Area*sixth*duxdx;
		Np12 = rho*Area*sixth*duxdy;
		Np21 = rho*Area*sixth*duydx;
		Np22 = rho*Area*sixth*duydy;

		Np31 = rho*Area*twelfth*duxdx;
		Np32 = rho*Area*twelfth*duxdy;
		Np41 = rho*Area*twelfth*duydx;
		Np42 = rho*Area*twelfth*duydy;
		
		//**********************Matrix Npdelta*******************************
		Npdelta11 = rho*tau*Area*third*C1*duxdx;
		Npdelta12 = rho*tau*Area*third*C1*duxdy;
		Npdelta21 = rho*tau*Area*third*C1*duydx;
		Npdelta22 = rho*tau*Area*third*C1*duydy;
	
		Npdelta31 = rho*tau*Area*third*C2*duxdx;
		Npdelta32 = rho*tau*Area*third*C2*duxdy;
		Npdelta41 = rho*tau*Area*third*C2*duydx;
		Npdelta42 = rho*tau*Area*third*C2*duydy;
	
		Npdelta51 = rho*tau*Area*third*C3*duxdx;
		Npdelta52 = rho*tau*Area*third*C3*duxdy;
		Npdelta61 = rho*tau*Area*third*C3*duydx;
		Npdelta62 = rho*tau*Area*third*C3*duydy;

		//*********************Matrix Nppdelta*******************************
		Nppdelta11 = rho*tau*sixth*y23*( duxdx*ux + duxdy*uy );
		Nppdelta12 = rho*tau*sixth*x32*( duxdx*ux + duxdy*uy );
		Nppdelta21 = rho*tau*sixth*y23*( duydx*ux + duydy*uy );
		Nppdelta22 = rho*tau*sixth*x32*( duydx*ux + duydy*uy );

		Nppdelta31 = rho*tau*sixth*y31*( duxdx*ux + duxdy*uy );
		Nppdelta32 = rho*tau*sixth*x13*( duxdx*ux + duxdy*uy );
		Nppdelta41 = rho*tau*sixth*y31*( duydx*ux + duydy*uy );
		Nppdelta42 = rho*tau*sixth*x13*( duydx*ux + duydy*uy );

		Nppdelta51 = rho*tau*sixth*y12*( duxdx*ux + duxdy*uy );
		Nppdelta52 = rho*tau*sixth*x21*( duxdx*ux + duxdy*uy );
		Nppdelta61 = rho*tau*sixth*y12*( duydx*ux + duydy*uy );
		Nppdelta62 = rho*tau*sixth*x21*( duydx*ux + duydy*uy );

		//************************Matrix Npphi*******************************
		Npphi11 = tau*sixth*( duxdx*y23 + duydx*x32 );
		Npphi12 = tau*sixth*( duxdy*y23 + duydy*x32 );

		Npphi21 = tau*sixth*( duxdx*y31 + duydx*x13 );
		Npphi22 = tau*sixth*( duxdy*y31 + duydy*x13 );

		Npphi31 = tau*sixth*( duxdx*y12 + duydx*x21 );
		Npphi32 = tau*sixth*( duxdy*y12 + duydy*x21 );
	
		//****************************Font***********************************
/*		fx1 = (*Font)(X[0],Y[0],1);*/
/*		fx2 = (*Font)(X[1],Y[1],1);*/
/*		fx3 = (*Font)(X[2],Y[2],1);*/
/*		*/
/*		fy1 = (*Font)(X[0],Y[0],2);*/
/*		fy2 = (*Font)(X[1],Y[1],2);*/
/*		fy3 = (*Font)(X[2],Y[2],2);*/

		fx1 = 0.0;
		fx2 = 0.0;
		fx3 = 0.0;
		
		fy1 = 0.0;
		fy2 = 0.0;
		fy3 = 0.0;


		//*********************Fonte Sol Exata Conhecida********************
		x = X[0];
		y = Y[0];
		
		fx1 = 2*x + pow(x,4)*(1.2 - 2.4*y) + pow(x,3)*(-2.4 + 4.8*y) + 
 y*(-0.4 + 1.2*y - 0.8*pow(y,2)) + 8*pow((-1 + x),3)*pow(x,3)*(-1 + 2*x)*pow(y,2)*pow((1 - 3*y + 2*pow(y,2)),2) + x*y*(2.4 - 7.2*y + 4.8*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,3)*(1 - 3*x + 2*pow(x,2))*pow((-1 + y),2)*pow(y,2)*(1 - 6*y + 6*pow(y,2)) + pow(x,2)*(1.2 - 4.8*y + 7.2*pow(y,2) - 4.8*pow(y,3));
		
		fy1 =	-2*y - 1.2*pow((1 - y),2)*pow(y,2) + 
 8*pow(x,2)*pow((1 - 3*x + 2*pow(x,2)),2)*pow((-1 + y),3)*pow(y,3)*(-1 + 2*y) + 
 pow(x,2)*(-1.2 + 7.2*y - 7.2*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,2)*(1 - 6*x + 6*pow(x,2))*pow((-1 + y),2)*pow(y,3)*(1 - 3*y + 2*pow(y,2)) + pow(x,3)*(0.8 - 4.8*y + 4.8*pow(y,2)) + x*(0.4 - 2.4*y + 4.8*pow(y,2) - 4.8*pow(y,3) + 2.4*pow(y,4));	

		x = X[1];
		y = Y[1];
		
		fx2 = 2*x + pow(x,4)*(1.2 - 2.4*y) + pow(x,3)*(-2.4 + 4.8*y) + 
 y*(-0.4 + 1.2*y - 0.8*pow(y,2)) + 8*pow((-1 + x),3)*pow(x,3)*(-1 + 2*x)*pow(y,2)*pow((1 - 3*y + 2*pow(y,2)),2) + x*y*(2.4 - 7.2*y + 4.8*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,3)*(1 - 3*x + 2*pow(x,2))*pow((-1 + y),2)*pow(y,2)*(1 - 6*y + 6*pow(y,2)) + pow(x,2)*(1.2 - 4.8*y + 7.2*pow(y,2) - 4.8*pow(y,3));
		
		fy2 =	-2*y - 1.2*pow((1 - y),2)*pow(y,2) + 
 8*pow(x,2)*pow((1 - 3*x + 2*pow(x,2)),2)*pow((-1 + y),3)*pow(y,3)*(-1 + 2*y) + 
 pow(x,2)*(-1.2 + 7.2*y - 7.2*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,2)*(1 - 6*x + 6*pow(x,2))*pow((-1 + y),2)*pow(y,3)*(1 - 3*y + 2*pow(y,2)) + pow(x,3)*(0.8 - 4.8*y + 4.8*pow(y,2)) + x*(0.4 - 2.4*y + 4.8*pow(y,2) - 4.8*pow(y,3) + 2.4*pow(y,4));
		
		x = X[2];
		y = Y[2];
		
		fx3 = 2*x + pow(x,4)*(1.2 - 2.4*y) + pow(x,3)*(-2.4 + 4.8*y) + 
 y*(-0.4 + 1.2*y - 0.8*pow(y,2)) + 8*pow((-1 + x),3)*pow(x,3)*(-1 + 2*x)*pow(y,2)*pow((1 - 3*y + 2*pow(y,2)),2) + x*y*(2.4 - 7.2*y + 4.8*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,3)*(1 - 3*x + 2*pow(x,2))*pow((-1 + y),2)*pow(y,2)*(1 - 6*y + 6*pow(y,2)) + pow(x,2)*(1.2 - 4.8*y + 7.2*pow(y,2) - 4.8*pow(y,3));
		
		fy3 =	-2*y - 1.2*pow((1 - y),2)*pow(y,2) + 
 8*pow(x,2)*pow((1 - 3*x + 2*pow(x,2)),2)*pow((-1 + y),3)*pow(y,3)*(-1 + 2*y) + 
 pow(x,2)*(-1.2 + 7.2*y - 7.2*pow(y,2)) - 
 4*pow((-1 + x),2)*pow(x,2)*(1 - 6*x + 6*pow(x,2))*pow((-1 + y),2)*pow(y,3)*(1 - 3*y + 2*pow(y,2)) + pow(x,3)*(0.8 - 4.8*y + 4.8*pow(y,2)) + x*(0.4 - 2.4*y + 4.8*pow(y,2) - 4.8*pow(y,3) + 2.4*pow(y,4));

		fxB = third*( fx1 + fx2 + fx3 );
		fyB = third*( fy1 + fy2 + fy3 );
		
		//***************Font of the Galerkin formulation********************
		f1 = rho*third*Area*fxB;
		f2 = rho*third*Area*fyB;

		//f1 = f2 = fxB = fyB = 0.0;
		
		//*****************Font of the SUPG formulation**********************
		fdelta1 = rho*tau*Area*C1*fxB;
		fdelta2 = rho*tau*Area*C1*fyB;
		fdelta3 = rho*tau*Area*C2*fxB;
		fdelta4 = rho*tau*Area*C2*fyB;
		fdelta5 = rho*tau*Area*C3*fxB;
		fdelta6 = rho*tau*Area*C3*fyB;

		//*****************Font of the PSPG formulation**********************
		fphi1 = tau*Area*( y23*fxB + x32*fyB );
		fphi2 = tau*Area*( y31*fxB + x13*fyB );
		fphi3 = tau*Area*( y12*fxB + x21*fyB );

		//*******************************vector Fe***************************
		Fe[0] = f1 + fdelta1;
		Fe[1] = f2 + fdelta2;
		Fe[2] = f1 + fdelta3;
		Fe[3] = f2 + fdelta4;
		Fe[4] = f1 + fdelta5;
		Fe[5] = f2 + fdelta6;
		Fe[6] = fphi1;
		Fe[7] = fphi2;
		Fe[8] = fphi3;

		//**********************Calculation of the residue*******************		
		Nv[0] = rho*third*Area*( duxdx*ux + duxdy*uy );
		Nv[1] = rho*third*Area*( duydx*ux + duydy*uy );
		Nv[2] = rho*third*Area*( duxdx*ux + duxdy*uy );
		Nv[3] = rho*third*Area*( duydx*ux + duydy*uy );
		Nv[4] = rho*third*Area*( duxdx*ux + duxdy*uy );
		Nv[5] = rho*third*Area*( duydx*ux + duydy*uy );


		Ndeltav[0] = rho*tau*Area*C1*( duxdx*ux + duxdy*uy );
		Ndeltav[1] = rho*tau*Area*C1*( duydx*ux + duydy*uy );
		Ndeltav[2] = rho*tau*Area*C2*( duxdx*ux + duxdy*uy );
		Ndeltav[3] = rho*tau*Area*C2*( duydx*ux + duydy*uy );
		Ndeltav[4] = rho*tau*Area*C3*( duxdx*ux + duxdy*uy );
		Ndeltav[5] = rho*tau*Area*C3*( duydx*ux + duydy*uy );


		Kv[0] = 0.5*visc*( 2*y23*duxdx + x32*( duxdy + duydx ) );
		Kv[1] = 0.5*visc*( 2*x32*duydy + y23*( duxdy + duydx ) );
		Kv[2] = 0.5*visc*( 2*y31*duxdx + x13*( duxdy + duydx ) );
		Kv[3] = 0.5*visc*( 2*x13*duydy + y31*( duxdy + duydx ) );
		Kv[4] = 0.5*visc*( 2*y12*duxdx + x21*( duxdy + duydx ) );
		Kv[5] = 0.5*visc*( 2*x21*duydy + y12*( duxdy + duydx ) );


		Gv[0] = 0.5*p*y23; 
		Gv[1] = 0.5*p*x32;
		Gv[2] = 0.5*p*y31;
		Gv[3] = 0.5*p*x13;
		Gv[4] = 0.5*p*y12;
		Gv[5] = 0.5*p*x21;

		Gdeltav[0] = tau*Area*C1*dpdx;
		Gdeltav[1] = tau*Area*C1*dpdy;
		Gdeltav[2] = tau*Area*C2*dpdx;
		Gdeltav[3] = tau*Area*C2*dpdy;
		Gdeltav[4] = tau*Area*C3*dpdx;
		Gdeltav[5] = tau*Area*C3*dpdy;
		
		//Residue R
		RQres[0] = Fe[0] - (Nv[0] + Ndeltav[0] + Kv[0] - ( Gv[0] + Gdeltav[0] ) );
		RQres[1] = Fe[1] - (Nv[1] + Ndeltav[1] + Kv[1] - ( Gv[1] + Gdeltav[1] ) );
		RQres[3] = Fe[2] - (Nv[2] + Ndeltav[2] + Kv[2] - ( Gv[2] + Gdeltav[2] ) );
		RQres[4] = Fe[3] - (Nv[3] + Ndeltav[3] + Kv[3] - ( Gv[3] + Gdeltav[3] ) );
		RQres[6] = Fe[4] - (Nv[4] + Ndeltav[4] + Kv[4] - ( Gv[4] + Gdeltav[4] ) );
		RQres[7] = Fe[5] - (Nv[5] + Ndeltav[5] + Kv[5] - ( Gv[5] + Gdeltav[5] ) );
				
		GTv[0] = third*Area*( duxdx + duydy);
		GTv[1] = third*Area*( duxdx + duydy);
		GTv[2] = third*Area*( duxdx + duydy);
		
		Nphiv[0] = 0.5*tau*(( duxdx*ux + duxdy*uy )*y23 + ( duydx*ux + duydy*uy )*x32);
		Nphiv[1] = 0.5*tau*(( duxdx*ux + duxdy*uy )*y31 + ( duydx*ux + duydy*uy )*x13);
		Nphiv[2] = 0.5*tau*(( duxdx*ux + duxdy*uy )*y12 + ( duydx*ux + duydy*uy )*x21);

		Gphiv[0] = 0.5*tau*invrho*( y23*dpdx + x32*dpdy); 
		Gphiv[1] = 0.5*tau*invrho*( y31*dpdx + x13*dpdy);
		Gphiv[2] = 0.5*tau*invrho*( y12*dpdx + x21*dpdy);

		//Residue Q

		RQres[2] = Fe[6] - ( GTv[0] + Nphiv[0] + Gphiv[0] );
		RQres[5] = Fe[7] - ( GTv[1] + Nphiv[1] + Gphiv[1] );
		RQres[8] = Fe[8] - ( GTv[2] + Nphiv[2] + Gphiv[2] );
				
		//**********************Tangent matrix*******************************
		//op = 1; // 0 se ISI e 1 se NI
		//if(op==0){		
		if(varglobal < 5){ //varglobal (<5 ISI) (>=5 NI)

		Ke[0][0] = N11 + Ndelta11 + K11;
		Ke[0][1] =       K12;
		Ke[0][3] = N13 + Ndelta13 + K13;
		Ke[0][4] =		  K14;
		Ke[0][6] = N15 + Ndelta15 + K15;
		Ke[0][7] = 		  K16;
		Ke[0][2] = -( GT11 + Gdelta11 );
		Ke[0][5] = -( GT11 + Gdelta12 );
		Ke[0][8] = -( GT11 + Gdelta13 );
	
		Ke[1][0] = 		  K12;
		Ke[1][1] = N11 + Ndelta22 + K22;
		Ke[1][3] = 		  K23;
		Ke[1][4] = N13 + Ndelta24 + K24;
		Ke[1][6] = 		  K25;
		Ke[1][7] = N15 + Ndelta26 + K26;
		Ke[1][2] = -( GT12 + Gdelta21 );
		Ke[1][5] = -( GT12 + Gdelta22 );
		Ke[1][8] = -( GT12 + Gdelta23 );
		

		Ke[3][0] = N11 + Ndelta13 + K13;
		Ke[3][1] =        K23;
		Ke[3][3] = N13 + Ndelta33 + K33;
		Ke[3][4] =		   K34;
		Ke[3][6] = N15 + Ndelta35 + K35;
		Ke[3][7] = 		   K36;
		Ke[3][2] = -( GT13 + Gdelta31 );
		Ke[3][5] = -( GT13 + Gdelta32 );
		Ke[3][8] = -( GT13 + Gdelta33 );


		Ke[4][0] = 		  K14;
		Ke[4][1] = N11 + Ndelta24 + K24;
		Ke[4][3] = 		  K34;
		Ke[4][4] = N13 + Ndelta44 + K44;
		Ke[4][6] = 		  K45;
		Ke[4][7] = N15 + Ndelta46 + K46;
		Ke[4][2] = -( GT14 + Gdelta41 );
		Ke[4][5] = -( GT14 + Gdelta42 );
		Ke[4][8] = -( GT14 + Gdelta43 );


		Ke[6][0] = N11 + Ndelta15 + K15;
		Ke[6][1] =       K25;
		Ke[6][3] = N13 + Ndelta35 + K35;
		Ke[6][4] =		 K45;
		Ke[6][6] = N15 + Ndelta55 + K55;
		Ke[6][7] = 		  K56;
		Ke[6][2] = -( GT15 + Gdelta51 );
		Ke[6][5] = -( GT15 + Gdelta52 );
		Ke[6][8] = -( GT15 + Gdelta53 );


		Ke[7][0] = 		  K16;
		Ke[7][1] = N11 + Ndelta26 + K26;
		Ke[7][3] = 		  K36;
		Ke[7][4] = N13 + Ndelta46 + K46;
		Ke[7][6] = 		  K56;
		Ke[7][7] = N15 + Ndelta66 + K66;
		Ke[7][2] = -( GT16 + Gdelta61 );
		Ke[7][5] = -( GT16 + Gdelta62 );
		Ke[7][8] = -( GT16 + Gdelta63 );


		Ke[2][0] = GT11 + Gdelta11;
		Ke[2][1] = GT12 + Gdelta21;
		Ke[2][3] = GT13 + Gdelta31;
		Ke[2][4] = GT14 + Gdelta41;
		Ke[2][6] = GT15 + Gdelta51;
		Ke[2][7] = GT16 + Gdelta61;
		Ke[2][2] = Gphi11;
		Ke[2][5] = Gphi12;
		Ke[2][8] = Gphi13;

	
		Ke[5][0] = GT11 + Gdelta12;
		Ke[5][1] = GT12 + Gdelta22;
		Ke[5][3] = GT13 + Gdelta32;
		Ke[5][4] = GT14 + Gdelta42;
		Ke[5][6] = GT15 + Gdelta52;
		Ke[5][7] = GT16 + Gdelta62;
		Ke[5][2] = Gphi12;
		Ke[5][5] = Gphi22;
		Ke[5][8] = Gphi23;


		Ke[8][0] = GT11 + Gdelta13;
		Ke[8][1] = GT12 + Gdelta23;
		Ke[8][3] = GT13 + Gdelta33;
		Ke[8][4] = GT14 + Gdelta43;
		Ke[8][6] = GT15 + Gdelta53;
		Ke[8][7] = GT16 + Gdelta63;
		Ke[8][2] = Gphi13;
		Ke[8][5] = Gphi23;
		Ke[8][8] = Gphi33;
		
		}else{

		Ke[0][0] = N11 + Np11 + Ndelta11 + Npdelta11 + Nppdelta11 + K11;
		Ke[0][1] =       Np12 			 + Npdelta12 + Nppdelta12 + K12;
		Ke[0][3] = N13 + Np31 + Ndelta13 + Npdelta11 + Nppdelta11 + K13;
		Ke[0][4] =		 Np32 	 		 + Npdelta12 + Nppdelta12 + K14;
		Ke[0][6] = N15 + Np31 + Ndelta15 + Npdelta11 + Nppdelta11 + K15;
		Ke[0][7] = 		 Np32 	 		 + Npdelta12 + Nppdelta12 + K16;
		Ke[0][2] = -( GT11 + Gdelta11 );
		Ke[0][5] = -( GT11 + Gdelta12 );
		Ke[0][8] = -( GT11 + Gdelta13 );
	
		Ke[1][0] = 		 Np21            + Npdelta21 + Nppdelta21 + K12;
		Ke[1][1] = N11 + Np22 + Ndelta22 + Npdelta22 + Nppdelta22 + K22;
		Ke[1][3] = 		 Np41			 + Npdelta21 + Nppdelta21 + K23;
		Ke[1][4] = N13 + Np42 + Ndelta24 + Npdelta22 + Nppdelta22 + K24;
		Ke[1][6] = 		 Np41 			 + Npdelta21 + Nppdelta21 + K25;
		Ke[1][7] = N15 + Np42 + Ndelta26 + Npdelta22 + Nppdelta22 + K26;
		Ke[1][2] = -( GT12 + Gdelta21 );
		Ke[1][5] = -( GT12 + Gdelta22 );
		Ke[1][8] = -( GT12 + Gdelta23 );
		

		Ke[3][0] = N11 + Np31 + Ndelta13 + Npdelta31 + Nppdelta31 + K13;
		Ke[3][1] =       Np32 			 + Npdelta32 + Nppdelta32 + K23;
		Ke[3][3] = N13 + Np11 + Ndelta33 + Npdelta31 + Nppdelta31 + K33;
		Ke[3][4] =		 Np12 	 		 + Npdelta32 + Nppdelta32 + K34;
		Ke[3][6] = N15 + Np31 + Ndelta35 + Npdelta31 + Nppdelta31 + K35;
		Ke[3][7] = 		 Np32 	 		 + Npdelta32 + Nppdelta32 + K36;
		Ke[3][2] = -( GT13 + Gdelta31 );
		Ke[3][5] = -( GT13 + Gdelta32 );
		Ke[3][8] = -( GT13 + Gdelta33 );


		Ke[4][0] = 		 Np41            + Npdelta41 + Nppdelta41 + K14;
		Ke[4][1] = N11 + Np42 + Ndelta24 + Npdelta42 + Nppdelta42 + K24;
		Ke[4][3] = 		 Np21			 + Npdelta41 + Nppdelta41 + K34;
		Ke[4][4] = N13 + Np22 + Ndelta44 + Npdelta42 + Nppdelta42 + K44;
		Ke[4][6] = 		 Np41 			 + Npdelta41 + Nppdelta41 + K45;
		Ke[4][7] = N15 + Np42 + Ndelta46 + Npdelta42 + Nppdelta42 + K46;
		Ke[4][2] = -( GT14 + Gdelta41 );
		Ke[4][5] = -( GT14 + Gdelta42 );
		Ke[4][8] = -( GT14 + Gdelta43 );


		Ke[6][0] = N11 + Np31 + Ndelta15 + Npdelta51 + Nppdelta51 + K15;
		Ke[6][1] =       Np32 			 + Npdelta52 + Nppdelta52 + K25;
		Ke[6][3] = N13 + Np31 + Ndelta35 + Npdelta51 + Nppdelta51 + K35;
		Ke[6][4] =		 Np32 	 		 + Npdelta52 + Nppdelta52 + K45;
		Ke[6][6] = N15 + Np11 + Ndelta55 + Npdelta51 + Nppdelta51 + K55;
		Ke[6][7] = 		 Np12 	 		 + Npdelta52 + Nppdelta52 + K56;
		Ke[6][2] = -( GT15 + Gdelta51 );
		Ke[6][5] = -( GT15 + Gdelta52 );
		Ke[6][8] = -( GT15 + Gdelta53 );


		Ke[7][0] = 		 Np41            + Npdelta61 + Nppdelta61 + K16;
		Ke[7][1] = N11 + Np42 + Ndelta26 + Npdelta62 + Nppdelta62 + K26;
		Ke[7][3] = 		 Np41			 + Npdelta61 + Nppdelta61 + K36;
		Ke[7][4] = N13 + Np42 + Ndelta46 + Npdelta62 + Nppdelta62 + K46;
		Ke[7][6] = 		 Np21 			 + Npdelta61 + Nppdelta61 + K56;
		Ke[7][7] = N15 + Np22 + Ndelta66 + Npdelta62 + Nppdelta62 + K66;
		Ke[7][2] = -( GT16 + Gdelta61 );
		Ke[7][5] = -( GT16 + Gdelta62 );
		Ke[7][8] = -( GT16 + Gdelta63 );


		Ke[2][0] = GT11 + Gdelta11 + Npphi11;
		Ke[2][1] = GT12 + Gdelta21 + Npphi12;
		Ke[2][3] = GT13 + Gdelta31 + Npphi11;
		Ke[2][4] = GT14 + Gdelta41 + Npphi12;
		Ke[2][6] = GT15 + Gdelta51 + Npphi11;
		Ke[2][7] = GT16 + Gdelta61 + Npphi12;
		Ke[2][2] = Gphi11;
		Ke[2][5] = Gphi12;
		Ke[2][8] = Gphi13;

	
		Ke[5][0] = GT11 + Gdelta12 + Npphi21;
		Ke[5][1] = GT12 + Gdelta22 + Npphi22;
		Ke[5][3] = GT13 + Gdelta32 + Npphi21;
		Ke[5][4] = GT14 + Gdelta42 + Npphi22;
		Ke[5][6] = GT15 + Gdelta52 + Npphi21;
		Ke[5][7] = GT16 + Gdelta62 + Npphi22;
		Ke[5][2] = Gphi12;
		Ke[5][5] = Gphi22;
		Ke[5][8] = Gphi23;


		Ke[8][0] = GT11 + Gdelta13 + Npphi31;
		Ke[8][1] = GT12 + Gdelta23 + Npphi32;
		Ke[8][3] = GT13 + Gdelta33 + Npphi31;
		Ke[8][4] = GT14 + Gdelta43 + Npphi32;
		Ke[8][6] = GT15 + Gdelta53 + Npphi31;
		Ke[8][7] = GT16 + Gdelta63 + Npphi32;
		Ke[8][2] = Gphi13;
		Ke[8][5] = Gphi23;
		Ke[8][8] = Gphi33;

		}

		//================
		//== Cond Contorno
		//================
		Keb[0][0] = N11 + Ndelta11 + K11;
		Keb[0][1] =       K12;
		Keb[0][3] = N13 + Ndelta13 + K13;
		Keb[0][4] =		  K14;
		Keb[0][6] = N15 + Ndelta15 + K15;
		Keb[0][7] = 		  K16;
		Keb[0][2] = -( GT11 + Gdelta11 );
		Keb[0][5] = -( GT11 + Gdelta12 );
		Keb[0][8] = -( GT11 + Gdelta13 );
	
		Keb[1][0] = 		  K12;
		Keb[1][1] = N11 + Ndelta22 + K22;
		Keb[1][3] = 		  K23;
		Keb[1][4] = N13 + Ndelta24 + K24;
		Keb[1][6] = 		  K25;
		Keb[1][7] = N15 + Ndelta26 + K26;
		Keb[1][2] = -( GT12 + Gdelta21 );
		Keb[1][5] = -( GT12 + Gdelta22 );
		Keb[1][8] = -( GT12 + Gdelta23 );

		Keb[3][0] = N11 + Ndelta13 + K13;
		Keb[3][1] =        K23;
		Keb[3][3] = N13 + Ndelta33 + K33;
		Keb[3][4] =		   K34;
		Keb[3][6] = N15 + Ndelta35 + K35;
		Keb[3][7] = 		   K36;
		Keb[3][2] = -( GT13 + Gdelta31 );
		Keb[3][5] = -( GT13 + Gdelta32 );
		Keb[3][8] = -( GT13 + Gdelta33 );

		Keb[4][0] = 		  K14;
		Keb[4][1] = N11 + Ndelta24 + K24;
		Keb[4][3] = 		  K34;
		Keb[4][4] = N13 + Ndelta44 + K44;
		Keb[4][6] = 		  K45;
		Keb[4][7] = N15 + Ndelta46 + K46;
		Keb[4][2] = -( GT14 + Gdelta41 );
		Keb[4][5] = -( GT14 + Gdelta42 );
		Keb[4][8] = -( GT14 + Gdelta43 );

		Keb[6][0] = N11 + Ndelta15 + K15;
		Keb[6][1] =       K25;
		Keb[6][3] = N13 + Ndelta35 + K35;
		Keb[6][4] =		 K45;
		Keb[6][6] = N15 + Ndelta55 + K55;
		Keb[6][7] = 		  K56;
		Keb[6][2] = -( GT15 + Gdelta51 );
		Keb[6][5] = -( GT15 + Gdelta52 );
		Keb[6][8] = -( GT15 + Gdelta53 );

		Keb[7][0] = 		  K16;
		Keb[7][1] = N11 + Ndelta26 + K26;
		Keb[7][3] = 		  K36;
		Keb[7][4] = N13 + Ndelta46 + K46;
		Keb[7][6] = 		  K56;
		Keb[7][7] = N15 + Ndelta66 + K66;
		Keb[7][2] = -( GT16 + Gdelta61 );
		Keb[7][5] = -( GT16 + Gdelta62 );
		Keb[7][8] = -( GT16 + Gdelta63 );

		Keb[2][0] = GT11 + Gdelta11;
		Keb[2][1] = GT12 + Gdelta21;
		Keb[2][3] = GT13 + Gdelta31;
		Keb[2][4] = GT14 + Gdelta41;
		Keb[2][6] = GT15 + Gdelta51;
		Keb[2][7] = GT16 + Gdelta61;
		Keb[2][2] = Gphi11;
		Keb[2][5] = Gphi12;
		Keb[2][8] = Gphi13;
	
		Keb[5][0] = GT11 + Gdelta12;
		Keb[5][1] = GT12 + Gdelta22;
		Keb[5][3] = GT13 + Gdelta32;
		Keb[5][4] = GT14 + Gdelta42;
		Keb[5][6] = GT15 + Gdelta52;
		Keb[5][7] = GT16 + Gdelta62;
		Keb[5][2] = Gphi12;
		Keb[5][5] = Gphi22;
		Keb[5][8] = Gphi23;

		Keb[8][0] = GT11 + Gdelta13;
		Keb[8][1] = GT12 + Gdelta23;
		Keb[8][3] = GT13 + Gdelta33;
		Keb[8][4] = GT14 + Gdelta43;
		Keb[8][6] = GT15 + Gdelta53;
		Keb[8][7] = GT16 + Gdelta63;
		Keb[8][2] = Gphi13;
		Keb[8][5] = Gphi23;
		Keb[8][8] = Gphi33;

		// Vetor Fe local do Elemento
		//Node 1
		if (Node[J1].id[0]==-1){ // nó prescrito
			vet_bound[0] = v1presc(X[0],Y[0]);
			
		}else{ // nó incógnita
			vet_bound[0] = 0.0;
			} 
		if (Node[J1].id[1]==-1){ // nó prescrito
			vet_bound[1] = v2presc(X[0],Y[0]);
		}else{ // nó incógnita
			vet_bound[1] = 0.0;
		}
		if (Node[J1].id[2]==-1){ // nó prescrito
			vet_bound[2] = ppresc(X[0],Y[0]);
		}else{ // nó incógnita
			vet_bound[2] = 0.0;
		} 
				
		//Node 2
		if (Node[J2].id[0]==-1){ // nó prescrito
			vet_bound[3] = v1presc(X[1],Y[1]);
		}else{ // nó incógnita
			vet_bound[3] = 0.0;
			} 
		if (Node[J2].id[1]==-1){ // nó prescrito
			vet_bound[4] = v2presc(X[1],Y[1]);
		}else{ // nó incógnita
			vet_bound[4] = 0.0;
		}
		if (Node[J2].id[2]==-1){ // nó prescrito
			vet_bound[5] = ppresc(X[1],Y[1]);
		}else{ // nó incógnita
			vet_bound[5] = 0.0;
		} 
	
		//Node 3
		if (Node[J3].id[0]==-1){ // nó prescrito
			vet_bound[6] = v1presc(X[2],Y[2]);
		}else{ // nó incógnita
			vet_bound[6] = 0.0;
			} 
		if (Node[J3].id[1]==-1){ // nó prescrito
			vet_bound[7] = v2presc(X[2],Y[2]);
		}else{ // nó incógnita
			vet_bound[7] = 0.0;
		}
		if (Node[J3].id[2]==-1){ // nó prescrito
			vet_bound[8] = ppresc(X[2],Y[2]);
		}else{ // nó incógnita
			vet_bound[8] = 0.0;
		} 

		/* direciona para um lixo durante o produto matriz vetor os elementos da matriz
		 que deveriam ser nulos (os que foram para o vetor força devido a c.c.). 
		Nos da a posicao de assemble global do vetor forca
		Aqui neq == numequacao + 1. */
		for (i = 0; i<NDOF; i++){
			if (Node[J1].id[i]==-1)
				lm_aux[i] = neq;
			else
				lm_aux[i] = Node[J1].id[i];

			if (Node[J2].id[i]==-1)
				lm_aux[i + NDOF] = neq;
			else
				lm_aux[i + NDOF] = Node[J2].id[i];
		
			if (Node[J3].id[i]==-1)
				lm_aux[i + 2*NDOF] = neq;
			else
				lm_aux[i + 2*NDOF] = Node[J3].id[i];
		}

		// Condicoes de fronteira no vetor Fe
		for (i = 0; i < 9; i++)
		{
			soma = 0.0;             
			for (j = 0; j < 9; j++)
				soma = soma + Keb[i][j]*vet_bound[j];
		//	F[lm_aux[i]] -= soma;
		}

		// Assemble global do vetor independente F de Au=F 
		for (i = 0; i < 9; i++)
			F[lm_aux[i]] +=  RQres[i];
		
		F[neq] = 0;


		
		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		assembly(Parameters, 0, E, Ke, MatrixData, lm);
		
	}//for elemento
		
	printf(" \n\n VARGLOBAL = %d \n\n", varglobal);
	varglobal++;

/*	if(varglobal==1){*/
/*	int j;*/
/*	printf("==PARA OCTAVE==\n\n");*/
/*	printf("nel=%d; \n",nel);*/
/*	printf("neq=%d; \n",neq);*/
/*	printf("nodes=%d; \n",Parameters->nnodes);	*/


/*	printf("Elements=[\n");*/
/*	for (i=0;i<nel;i++)*/
/*		printf("%d %d %d;\n",Element[i].Vertex[0]+1,Element[i].Vertex[1]+1,Element[i].Vertex[2]+1);*/
/*	printf("];\n");*/

/*	printf("Nodes=[\n");*/
/*	for (i=0;i<Parameters->nnodes;i++){*/
/*		if (Node[i].id[0] == -1)*/
/*			printf("%d ",neq+1);*/
/*		else*/
/*			printf("%d ",Node[i].id[0]+1);*/
/*		if (Node[i].id[1] == -1)*/
/*			printf("%d ",neq+1);*/
/*		else*/
/*			printf("%d ",Node[i].id[1]+1);*/
/*		if (Node[i].id[2] == -1)*/
/*			printf("%d ",neq+1);*/
/*		else*/
/*			printf("%d ",Node[i].id[2]+1);*/
/*		printf(";\n");*/
/*	}	*/
/*	printf("];\n");*/

/*	printf("lm=[\n");*/
/*	for (i=0;i<nel;i++){*/
/*		for (j=0;j<9;j++)*/
/*			printf("%d ",lm[i][j]+1);*/
/*		printf(";\n");*/
/*	}*/
/*	printf("];\n");*/

/*	printf("A=[\n");*/
/*	for(i=0;i<nel;i++){*/
/*		for (j=0;j<81;j++){*/
/*			printf("%lf ",MatrixData->A[0][i][j]);*/
/*			//printf("\n passei aqui \n");*/
/*		}		*/
/*		printf(";\n");*/
/*	}*/
/*	printf("];\n");*/
/*	*/
/*	printf("\n M = monta_Matriz_Global(A,lm, neq, nel);\n");*/
/*	printf("cond(M)\n");*/
/*	printf("[Diag,invDiag,PA] = monta_DiagSemFonte(A,Elements,Nodes,lm,nodes,neq,nel);\n");*/
/*	printf("\n MP = monta_Matriz_Global(PA,lm, neq, nel);\n");*/
/*	printf("cond(MP)\n");*/

/*	printf("\n==OCTAVE ATE AQUI==\n");*/
/*	//exit(1);*/
/*	}*/
/*
	printf("F=[\n");
	for (i=0;i<neq;i++)
		printf("%lf;\n",F[i]);
	printf("];\n");

	exit(1);		
*/	

	end_assembly(Parameters, 0, MatrixData);

/*	int j; 
	FILE *Out;
	Out = fopen("/home/lmuniz/Cavity1.mtx","w");
	fprintf(Out,"%d\t%d\t%d\n",Parameters->neq,Parameters->neq,Parameters->nnzero);
	for (i=0; i<Parameters->neq;i++)
		for (j=MatrixData->IA[0][i];j<MatrixData->IA[0][i+1];j++)
			fprintf(Out,"%d\t%d\t%.14lf\n",i+1,MatrixData->JA[0][j]+1,MatrixData->AA[0][j]);
	fclose(Out);
	exit(1);
*/	



	// Liberacao dos espacos alocados
//	free(vet_bound);

	return 0;

}
