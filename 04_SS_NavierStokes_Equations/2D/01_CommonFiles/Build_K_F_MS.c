#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

static int varglobal = 0;

int Build_K_F_MS(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int i, neq, nel;//, op;
	int J, E, J1, J2, J3;
	//double Uref, Lref;
	double Re, VelMax, ReLoc, ReLocMax, ReLocMin, ReLocMed, epsilon, coefnorm;
	double normres, normgradU;//, normresF;
	double hRGN;//, GNUx, GNUy, NGNU, r1, r2;
	double twoArea, invArea, Area, h, visc, rho, nu, unorm, tau_l, delta_aux, delta, deltam; //, invrho
	double third=1.0/3.0, sixth = 1.0/6.0,  ninefortieth = 9.0/40.0, thirdtwentieth = 3.0/20.0, ninetwentieth = 9.0/20.0, twelfth = 1.0/12.0;
	double y23, y31, y12, x32, x13, x21, X[3], Y[3], Ke[9][9], Fe[9], ux1, ux2, ux3, uy1, uy2, uy3, ux, uy, uxB, uyB, duxdx, duxdy, duydx, duydy;
	double fx1, fx2, fx3, fy1, fy2, fy3, fxB, fyB, f1, f2, fdelta1, fdelta2, fdelta3, fdelta4, fdelta5, fdelta6, fphi1, fphi2, fphi3;
	double f1B, f2B;	
	double C1, C2, C3;
	double Nv[6], Ndeltav[6], Kv[6], Ksv[6], Kddv[6], Gv[6], Gdeltav[6], GTv[3], Cv[3], Nphiv[3], Gphiv[3], Res[9];
	//double M11, M13, Mdelta11, Mdelta31, Mdelta51;
	double K11, K12, K13, K14, K15, K16, K22, K23, K24, K25, K26, K33, K34, K35, K36, K44, K45, K46, K55, K56, K66; 
	double Ks11, Ks12, Ks13, Ks14, Ks15, Ks16, Ks21, Ks22, Ks23, Ks24, Ks25, Ks26, Ks31, Ks32, Ks33, Ks34, Ks35, Ks36;
	double Ks41, Ks42, Ks43, Ks44, Ks45, Ks46, Ks51, Ks52, Ks53, Ks54, Ks55, Ks56, Ks61, Ks62, Ks63, Ks64, Ks65, Ks66;
	double Kdd11, Kdd12, Kdd13, Kdd22, Kdd23, Kdd33; 
	double GT11, GT12, GT13, GT14, GT15, GT16, N11, N13, N15, C11, C12;
	double Ndelta11, Ndelta12, Ndelta13, Ndelta14, Ndelta15, Ndelta16, Ndelta21, Ndelta22, Ndelta23, Ndelta24, Ndelta25, Ndelta26;
	double Ndelta31, Ndelta32, Ndelta33, Ndelta34, Ndelta35, Ndelta36, Ndelta41, Ndelta42, Ndelta43, Ndelta44, Ndelta45, Ndelta46;
	double Ndelta51, Ndelta52, Ndelta53, Ndelta54, Ndelta55, Ndelta56, Ndelta61, Ndelta62, Ndelta63, Ndelta64, Ndelta65, Ndelta66;
	double Gdelta11, Gdelta12, Gdelta13, Gdelta21, Gdelta22, Gdelta23, Gdelta31, Gdelta32, Gdelta33, Gdelta41, Gdelta42, Gdelta43;
	double Gdelta51, Gdelta52, Gdelta53, Gdelta61, Gdelta62, Gdelta63;
	//double Mphi11, Mphi12, Mphi21, Mphi22, Mphi31, Mphi32;
	double Nphi11, Nphi12, Nphi13, Nphi14, Nphi15, Nphi16, Nphi21, Nphi22, Nphi23, Nphi24, Nphi25, Nphi26, Nphi31, Nphi32, Nphi33;
	double Nphi34, Nphi35, Nphi36;	
	double Gphi11, Gphi12, Gphi13, Gphi21, Gphi22, Gphi23, Gphi31, Gphi32, Gphi33;
	double Np11, Np12, Np21, Np22, Np31, Np32, Np41, Np42;
	double Ndeltap11, Ndeltap12, Ndeltap21, Ndeltap22, Ndeltap31, Ndeltap32, Ndeltap41, Ndeltap42, Ndeltap51, Ndeltap52, Ndeltap61, Ndeltap62; 
	//double Nppdelta11, Nppdelta12, Nppdelta21, Nppdelta22, Nppdelta31, Nppdelta32, Nppdelta41, Nppdelta42, Nppdelta51, Nppdelta52, Nppdelta61, Nppdelta62;
	double Nphip11, Nphip12, Nphip21, Nphip22, Nphip31, Nphip32;
	double p1, p2, p3, p, dpdx, dpdy, duxduy;
	double ne1, ne2, cont, auxh;		
		
	
	double *delta_old = FemStructs->delta_old; 
	double *normres_old = FemStructs->normres_old; 
	//double *deltam_old = FemStructs->deltam_old; 
	double *F = FemStructs->F;
	int **lm = FemStructs->lm;
	
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	nel = Parameters->nel;
	neq = Parameters->neq;
	VelMax = Parameters->VelMax;
	
	double *U = (double*) mycalloc("U of 'Build_K_F_VMS_DCDD'", 3*Parameters->nnodes, sizeof(double));
	eval_U(Parameters,FemStructs,FemFunctions,U);

	setzeros(Parameters,MatrixData);

	for (J=0; J<neq+1; J++)
		F[J] = 0;
	
	ReLoc = 0.0;
	ReLocMax = 0.0;
	ReLocMin = 1000.0;
	ReLocMed = 0.0;
	
	epsilon = 1e-7;	//Usar 10-7
	coefnorm = 0.8;	

	double coe01, coe02, delta_aux01, delta_aux02; 
	// coe01 = coe02 = 0.0 .=. metodo galerkim enriquecido - VMS
	// coe01 = 1.0 e coe02 = 0.0 .=. metodo NSGS1-NS
	// coe01 = 1.0 e coe02 = 1.0 .=. metodo NSGS2-NS
	coe01 = 1.0;  
	coe02 = 1.0;

	int conta=0, conta2=0;

	for (E=0; E<nel; E++)
	{
		conta++;

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
		ux = third*(ux1 + ux2 + ux3);
		uy = third*(uy1 + uy2 + uy3);
		
		//*****************************Gradient of u********************
		duxdx = 0.5*invArea*( ux1*y23 + ux2*y31 + ux3*y12 ); 
		duxdy = 0.5*invArea*( ux1*x32 + ux2*x13 + ux3*x21 );
		duydx = 0.5*invArea*( uy1*y23 + uy2*y31 + uy3*y12 );
		duydy = 0.5*invArea*( uy1*x32 + uy2*x13 + uy3*x21 );
			
		duxduy = ux1*y23 + uy1*x32 + ux2*y31 + uy2*x13 + ux3*y12 + uy3*x21;


		//*****************p at barycentre of the element***************
		p = third*( p1 + p2 + p3 );

		//**************************Gradient of p**********************
		dpdx = 0.5*invArea*( p1*y23 + p2*y31 + p3*y12 );
		dpdy = 0.5*invArea*( p1*x32 + p2*x13 + p3*x21 );

		//***Coefficients C_i em funcao da vel macro (u_h) apenas *********
		C1 = 0.5*invArea*( ux*y23 + uy*x32 );
		C2 = 0.5*invArea*( ux*y31 + uy*x13 );
		C3 = 0.5*invArea*( ux*y12 + uy*x21 );

		//************* Dados de entrada *********
		Re = Parameters->ReynoldsNumber;
		//h = sqrt( 4*Area/PI );		//definido por Renato Elias em sua dissertacao  
		h = sqrt(twoArea);		
		visc =1./Re;
		rho = 1.;	
		nu = visc/rho;	
		unorm = sqrt( ux*ux + uy*uy );

		//****************calculo do Re local************************
		if(varglobal==0){		
			ReLoc = VelMax*h/nu;
			ReLocMin = (ReLoc < ReLocMin)? ReLoc: ReLocMin;
			ReLocMax = (ReLoc > ReLocMax)? ReLoc: ReLocMax;
			ReLocMed += ReLoc;
			if(E==nel-1)
				ReLocMed = ReLocMed/nel;
		}
		
		//*************Calculation of tau_LSIC stabilization*********
		tau_l = 0.333*h*unorm;		
		//tau_l = 0.5*h*unorm;		
		//tau_l = sqrt(nu*nu + (0.5*unorm*h)*(0.5*unorm*h));
		//tau_l = 0.0;
	
		//==============================================
		//*********************Fonte********************
		//==============================================
		fx1 = FemFunctions->f1ext(X[0],Y[0]);
		fy1 = FemFunctions->f2ext(X[0],Y[0]);

		fx2 = FemFunctions->f1ext(X[1],Y[1]);
		fy2 = FemFunctions->f2ext(X[1],Y[1]);

		fx3 = FemFunctions->f1ext(X[2],Y[2]);
		fy3 = FemFunctions->f2ext(X[2],Y[2]);
		
		fxB = third*( fx1 + fx2 + fx3 );
		fyB = third*( fy1 + fy2 + fy3 );
		
		//***************Font of the Galerkin formulation********************
		f1 = rho*third*Area*fxB;
		f2 = rho*third*Area*fyB;

		//f1 = f2 = fxB = fyB = 0.0;
		
		//***************Font of the Bobble formulation********************
		f1B = rho*ninetwentieth*Area*fxB;
		f2B = rho*ninetwentieth*Area*fyB;


		//*******************************vector Fe***************************
		Fe[0] = f1;
		Fe[1] = f2;
		Fe[2] = 0.0;
		Fe[3] = f1;
		Fe[4] = f2;
		Fe[5] = 0.0;
		Fe[6] = f1;
		Fe[7] = f2;
		Fe[8] = 0.0;
		
		//============================================================
		//************Calculo da velocidade na micro escala **********
		//============================================================
		ne1 = ux*duxdx + uy*duydx + dpdx - f1;
		ne2 = ux*duxdy + uy*duydy + dpdy - f2;
		cont = duxdx + duydy;
		
		//***************Norma do residuo ****************************
		// Com Div v = 0 
		//normres = sqrt(ne1*ne1 + ne2*ne2 + cont*cont); 
		// Sem Div v = 0 
		normres = sqrt(ne1*ne1 + ne2*ne2);
					
		//*****************Norma do grad de U = (Velocidade, pressao)
		//normgradU = sqrt(duxdx*duxdx + duxdy*duxdy + duydx*duydx + duydy*duydy + dpdx*dpdx + dpdy*dpdy);
		//*****************Norma do grad de U = (Velocidade)
		normgradU = sqrt(duxdx*duxdx + duxdy*duxdy + duydx*duydx + duydy*duydy);
	
		//******* Parametro de malha h  *************
		//
		double tolnormgradU = 1e-5;
		double w = 0.5;		
		/*
		if(ux>toluxy || uy>toluxy){
			auxh = fabs(ux*y23 + uy*x32) + fabs(ux*y31 + uy*x13) + fabs(ux*y12 + uy*x21);
			hRGN = 4*Area*unorm/auxh;
		*/
		auxh = fabs(ux*y23 + uy*x32) + fabs(ux*y31 + uy*x13) + fabs(ux*y12 + uy*x21);
		if(fabs(auxh) >= 1e-6){
			hRGN = 4*Area*unorm/auxh;
		}else{
			hRGN = 0.0;
		}
		//
		//******* Parametro de difusao artificial := delta ******
		if(normgradU > tolnormgradU){
			delta_aux01 = coe01*hRGN*normres/normgradU;
			delta_aux02 = coe02*hRGN*hRGN*sqrt(cont*cont);
			//delta_aux = (delta_aux01 > delta_aux02)? delta_aux01 : delta_aux02;			
			//delta_aux = sqrt(delta_aux01*delta_aux01 + delta_aux02*delta_aux02); 
			delta_aux = (delta_aux01 + delta_aux02)/2.0;
			//delta_aux = delta_aux01;
		}else{
			delta_aux = 0.0;
		}

		if(varglobal==0 || normres < coefnorm*normres_old[E]){
			delta =  w*delta_aux + (1-w)*delta_old[E];
		}else{
			delta = delta_old[E];
		}
		deltam = delta;
		//delta = 0.0;
				
		//****** Matrix InvKBBd = (KBB + KBBdelta)^{-1}**********
		double A, B, InvKBBd11, InvKBBd12, InvKBBd22, C, D, KBB11, KBB12, KBB22, NBB;

		A = y23*y23 + y31*y31 + y12*y12 + y23*y31 + y23*y12 + y31*y12;
		B = x32*x32 + x13*x13 + x21*x21 + x32*x13 + x32*x21 + x13*x21;
		C = y23*x32 + y31*x13 + y12*x21;
		
		//  termo convectivo da micro escala com xb = 0.34 e yb = 0.34 *********
		double cg1=0.311998, cg2=0.311998, cg3=0.594661;
		
		NBB = cg1*(ux*y23 + uy*x32) + cg2*(ux*y31 + uy*x13) + cg3*(ux*y12 + uy*x21);
		NBB = 0.0;		
		//printf("NBB_%d = %f \n", E, NBB);
		
		KBB11 = 81.0*((2*visc+deltam)*A + (visc+deltam)*B)/(40.0*Area) + NBB;
		KBB12 = 81.0*visc*C/(80.0*Area);
		KBB22 = 81.0*((2*visc+deltam)*B + (visc+deltam)*A)/(40.0*Area) + NBB;		

		//******** D eh o determinande da matriz (KBB + KBBdelta) 
		D = KBB11*KBB22 - KBB12*KBB12;

		//****** Matrix InvKBBd **********
		InvKBBd11 = KBB22/D;
		InvKBBd12 = - KBB12/D;
		InvKBBd22 = KBB11/D;

		//*************** Matrix NBh ***************************
		double NBh11, NBh13, NBh15;
		NBh11 = rho*ninefortieth*(ux*y23 + uy*x32);
		NBh13 = rho*ninefortieth*(ux*y31 + uy*x13);
		NBh15 = rho*ninefortieth*(ux*y12 + uy*x21);

		//************** Matrix GBh ****************************
		double GBh11, GBh12, GBh13, GBh21, GBh22, GBh23;		
		GBh11 = - ninefortieth*y23;		
		GBh12 = - ninefortieth*y31;		
		GBh13 = - ninefortieth*y12;		
		GBh21 = - ninefortieth*x32;		
		GBh22 = - ninefortieth*x13;		
		GBh23 = - ninefortieth*x21;		
				
		uxB = InvKBBd11*(f1B - (NBh11*ux1 + NBh13*ux2 + NBh15*ux3) + (GBh11*p1 + GBh12*p2 + GBh13*p3)) + InvKBBd12*(f2B - (NBh11*uy1 + NBh13*uy2 + NBh15*uy3) + (GBh21*p1 + GBh22*p2 + GBh23*p3));
		uyB = InvKBBd12*(f1B - (NBh11*ux1 + NBh13*ux2 + NBh15*ux3) + (GBh11*p1 + GBh12*p2 + GBh13*p3)) + InvKBBd22*(f2B - (NBh11*uy1 + NBh13*uy2 + NBh15*uy3) + (GBh21*p1 + GBh22*p2 + GBh23*p3));

		//==== Fim do Claculo da velocidade na micro escala=======================
		//if(E % 100 ==0)		
		//	printf("\nE = %d Velocidades ux = %lf, uxB = %lf, uy = %lf, uyB = %lf",E, ux, uxB, uy, uyB);
		 
		ux = ux + uxB;
		uy = uy + uyB;		

		//==== Fim do Claculo da velocidade total = micro + macro escalas =========
		//=========================================================================

		//***Coefficients C_i em funcao da vel macro (u_h) e micro (u_B) *********
		C1 = 0.5*invArea*( ux*y23 + uy*x32 );
		C2 = 0.5*invArea*( ux*y31 + uy*x13 );
		C3 = 0.5*invArea*( ux*y12 + uy*x21 );

		//************Matrices of the Galerkin formulation*************

		//*************************Matrix M****************************
		//M11 = rho*2*Area*twelfth;
		//M13 = rho*Area*twelfth; 
		
		//*************************Matrix N****************************
		N11 = rho*third*Area*C1;		
		N13 = rho*third*Area*C2;		
		N15 = rho*third*Area*C3;	
		
		//**************Incremental Matrix Np**************************
		Np11 = rho*Area*sixth*duxdx;
		Np12 = rho*Area*sixth*duxdy;
		Np21 = rho*Area*sixth*duydx;
		Np22 = rho*Area*sixth*duydy;

		Np31 = rho*Area*twelfth*duxdx;
		Np32 = rho*Area*twelfth*duxdy;
		Np41 = rho*Area*twelfth*duydx;
		Np42 = rho*Area*twelfth*duydy;	

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
		
		//************************Matrix Cij - media nula*************************
		C11 = epsilon*Area*sixth;
		C12 = epsilon*Area*twelfth;
		//***************Matrices of the LSIC formulation****************
		
		//***************Matrices K_s ***********************************
		
		Ks11 = tau_l*rho*0.25*invArea*y23*y23;
		Ks12 = tau_l*rho*0.25*invArea*x32*y23;
		Ks13 = tau_l*rho*0.25*invArea*y31*y23;
		Ks14 = tau_l*rho*0.25*invArea*x13*y23;
		Ks15 = tau_l*rho*0.25*invArea*y12*y23;
		Ks16 = tau_l*rho*0.25*invArea*x21*y23;

		Ks21 = tau_l*rho*0.25*invArea*y23*x32;
		Ks22 = tau_l*rho*0.25*invArea*x32*x32;
		Ks23 = tau_l*rho*0.25*invArea*y31*x32;
		Ks24 = tau_l*rho*0.25*invArea*x13*x32;
		Ks25 = tau_l*rho*0.25*invArea*y12*x32;
		Ks26 = tau_l*rho*0.25*invArea*x21*x32;

		Ks31 = tau_l*rho*0.25*invArea*y23*y31;
		Ks32 = tau_l*rho*0.25*invArea*x32*y31;
		Ks33 = tau_l*rho*0.25*invArea*y31*y31;
		Ks34 = tau_l*rho*0.25*invArea*x13*y31;
		Ks35 = tau_l*rho*0.25*invArea*y12*y31;
		Ks36 = tau_l*rho*0.25*invArea*x21*y31;

		Ks41 = tau_l*rho*0.25*invArea*y23*x13;
		Ks42 = tau_l*rho*0.25*invArea*x32*x13;
		Ks43 = tau_l*rho*0.25*invArea*y31*x13;
		Ks44 = tau_l*rho*0.25*invArea*x13*x13;
		Ks45 = tau_l*rho*0.25*invArea*y12*x13;
		Ks46 = tau_l*rho*0.25*invArea*x21*x13;
		
		Ks51 = tau_l*rho*0.25*invArea*y23*y12;
		Ks52 = tau_l*rho*0.25*invArea*x32*y12;
		Ks53 = tau_l*rho*0.25*invArea*y31*y12;
		Ks54 = tau_l*rho*0.25*invArea*x13*y12;
		Ks55 = tau_l*rho*0.25*invArea*y12*y12;
		Ks56 = tau_l*rho*0.25*invArea*x21*y12;

		Ks61 = tau_l*rho*0.25*invArea*y23*x21;
		Ks62 = tau_l*rho*0.25*invArea*x32*x21;
		Ks63 = tau_l*rho*0.25*invArea*y31*x21;
		Ks64 = tau_l*rho*0.25*invArea*x13*x21;
		Ks65 = tau_l*rho*0.25*invArea*y12*x21;
		Ks66 = tau_l*rho*0.25*invArea*x21*x21;
		

		
		//**********************Calculation of the residue*******************		
		Nv[0] = rho*third*Area*( duxdx*ux + duxdy*uy );
		Nv[1] = rho*third*Area*( duydx*ux + duydy*uy );
		Nv[2] = rho*third*Area*( duxdx*ux + duxdy*uy );
		Nv[3] = rho*third*Area*( duydx*ux + duydy*uy );
		Nv[4] = rho*third*Area*( duxdx*ux + duxdy*uy );
		Nv[5] = rho*third*Area*( duydx*ux + duydy*uy );

		Kv[0] = 0.5*visc*( 2*y23*duxdx + x32*( duxdy + duydx ) );
		Kv[1] = 0.5*visc*( 2*x32*duydy + y23*( duxdy + duydx ) );
		Kv[2] = 0.5*visc*( 2*y31*duxdx + x13*( duxdy + duydx ) );
		Kv[3] = 0.5*visc*( 2*x13*duydy + y31*( duxdy + duydx ) );
		Kv[4] = 0.5*visc*( 2*y12*duxdx + x21*( duxdy + duydx ) );
		Kv[5] = 0.5*visc*( 2*x21*duydy + y12*( duxdy + duydx ) );
		
		Ksv[0] = tau_l*rho*0.25*invArea*y23*duxduy;
		Ksv[1] = tau_l*rho*0.25*invArea*x32*duxduy;
		Ksv[2] = tau_l*rho*0.25*invArea*y31*duxduy;
		Ksv[3] = tau_l*rho*0.25*invArea*x13*duxduy;
		Ksv[4] = tau_l*rho*0.25*invArea*y12*duxduy;
		Ksv[5] = tau_l*rho*0.25*invArea*x21*duxduy;

		Gv[0] = 0.5*p*y23; 
		Gv[1] = 0.5*p*x32;
		Gv[2] = 0.5*p*y31;
		Gv[3] = 0.5*p*x13;
		Gv[4] = 0.5*p*y12;
		Gv[5] = 0.5*p*x21;

		GTv[0] = third*Area*( duxdx + duydy);
		GTv[1] = third*Area*( duxdx + duydy);
		GTv[2] = third*Area*( duxdx + duydy);

		Cv[0] = C11*p1 + C12*(p2+p3);
		Cv[1] = C11*p2 + C12*(p1+p3);
		Cv[2] = C11*p3 + C12*(p1+p2);

		//*****************Residuo para o calculo do coeficiente da difusao artificial 
		// *****************Norma euclidiana do residuo da formulacao variacional ||R(U,P)||
/*	
		Res[0] = Fe[0] - (Nv[0] + Kv[0] - Gv[0]);
		Res[1] = Fe[1] - (Nv[1] + Kv[1] - Gv[1]);
		Res[2] = Fe[2] - GTv[0] ;
		Res[3] = Fe[3] - (Nv[2] + Kv[2] - Gv[2] );
		Res[4] = Fe[4] - (Nv[3] + Kv[3] - Gv[3] );
		Res[5] = Fe[5] - GTv[1] ;
		Res[6] = Fe[6] - (Nv[4] + Kv[4] - Gv[4]);
		Res[7] = Fe[7] - (Nv[5] + Kv[5] - Gv[5]);
		Res[8] = Fe[8] - GTv[2];
		normres = sqrt(Res[0]*Res[0] + Res[1]*Res[1] + Res[2]*Res[2] + Res[3]*Res[3] + Res[4]*Res[4] + Res[5]*Res[5] + Res[6]*Res[6] + Res[7]*Res[7] + Res[8]*Res[8]);
		normres = sqrt(Res[0]*Res[0] + Res[1]*Res[1] + Res[3]*Res[3] + Res[4]*Res[4] + Res[6]*Res[6] + Res[7]*Res[7]);
*/	
		//*****************Norma euclidiana do residuo ||R(U,P)|| da equacao forte		
		ne1 = ux*duxdx + uy*duydx + dpdx - f1;
		ne2 = ux*duxdy + uy*duydy + dpdy - f2;
		//cont = duxdx + duydy;
		//
		//normres_oldaux = normres_old[E];
		//
		// Com Div v = 0 
		//normres = sqrt(ne1*ne1 + ne2*ne2 + cont*cont); 
		// Sem Div v = 0 
		normres = sqrt(ne1*ne1 + ne2*ne2); 
		//		
				
		//*****************Norma do grad de U = (Velocidade, pressao)
		//normgradU = sqrt(duxdx*duxdx + duxdy*duxdy + duydx*duydx + duydy*duydy + dpdx*dpdx + dpdy*dpdy);
		//*****************Norma do grad de U = (Velocidade)
		normgradU = sqrt(duxdx*duxdx + duydy*duydy + duxdy*duxdy+ duydx*duydx);
		
		//********** norma de u = u_h + u_B *********************
		unorm = sqrt( ux*ux + uy*uy );

		//*** Vetor: Grad da Norma de U = < GNUx, GNUy >*(unorm)********
//		GNUx = ux*duxdx + uy*duydx;
//		GNUy = ux*duxdy + uy*duydy;
		
		//*** Escalar: Norma do Grad da Norma de U = (somente Velocidade) ***********
//		NGNU = sqrt(GNUx*GNUx + GNUy*GNUy);

		//*** Vetor r: Grad da Norma de U(Velocidade) normalizado Ref. Tezduyar DCDD 2001**********
/*		if(NGNU>toluxy){
			r1 = GNUx/NGNU;
			r2 = GNUy/NGNU;
		}else{
			r1 = 0.0;
			r2 = 0.0;
		}
*/
		//******* Parametro de malha h_RGN Ref. Tezduyar DCDD 2001*************
/*		if(r1>toluxy || r2>toluxy){
			//auxh = fabs(r1*y23 + r2*x32) + fabs(r1*y31 + r2*x13) + fabs(r1*y12 + r2*x21);
			//hRGN = 4*Area/auxh;			
		}else{
			hRGN = 0.0;
		}
*/		
		//*******Parametro de malha tradicional
		/*		
		if(ux>toluxy || uy>toluxy){
			auxh = fabs(ux*y23 + uy*x32) + fabs(ux*y31 + uy*x13) + fabs(ux*y12 + uy*x21);
			hRGN = 4*Area*unorm/auxh;
		*/
		auxh = fabs(ux*y23 + uy*x32) + fabs(ux*y31 + uy*x13) + fabs(ux*y12 + uy*x21);
		if(fabs(auxh) >= 1e-6){
			hRGN = 4*Area*unorm/auxh;
		}else{
			hRGN = 0.0;
		}
		
		//******* Parametro de difusao artificial := delta ******
		if(normgradU > tolnormgradU){
			delta_aux01 = coe01*hRGN*normres/normgradU;
			delta_aux02 = coe02*hRGN*hRGN*sqrt(cont*cont);
			//delta_aux = (delta_aux01 > delta_aux02)? delta_aux01 : delta_aux02;			
			//delta_aux = sqrt(delta_aux01*delta_aux01 + delta_aux02*delta_aux02); 
			delta_aux = (delta_aux01 + delta_aux02)/2.0;
			//delta_aux = delta_aux01;
		}else{
			delta_aux = 0.0;
		}

		//if(varglobal == 2 || varglobal == 6)
		/*	if((conta % 3000) == 0)	{		
				//printf("Elem = %d, x = %f, y = %f, |Res| = %f, |GradU| = %f, deltam = %f, nu = %f  \n\n", E, X[0], Y[0], normres, normgradU, deltam, nu); 
				printf("\n |Res|/|ResOld| = %f", normres/normres_old[E]); 
		}
		*/
		if(varglobal==0 || normres < coefnorm*normres_old[E]){
			delta =  w*delta_aux + (1-w)*delta_old[E];
			delta_old[E] = delta;
			conta2++;
		}else{
			delta = delta_old[E];
		}
		//
		normres_old[E] = normres;
		//
		deltam = delta;
		delta = 0.0;
		//
		//if(varglobal == 2 || varglobal == 6)
		/*	if((conta % 3000) == 0)	{		
				//printf("Elem = %d, x = %f, y = %f, |Res| = %f, |GradU| = %f, deltam = %f, nu = %f  \n\n", E, X[0], Y[0], normres, normgradU, deltam, nu); 
				printf("\n |Res|= %f, |GradU|= %f, del_m= %.1e, del_aux01= %.1e, del_aux02= %.1e, nu= %.1e \n", normres, normgradU, deltam, delta_aux01, delta_aux02, nu); 
		}		
		*/		
		//****** Matrix InvKBBd = (KBB + KBBdelta)^{-1}**********
		A = y23*y23 + y31*y31 + y12*y12 + y23*y31 + y23*y12 + y31*y12;
		B = x32*x32 + x13*x13 + x21*x21 + x32*x13 + x32*x21 + x13*x21;
		C = y23*x32 + y31*x13 + y12*x21;
		//
		KBB11 = 81.0*((2*visc+deltam)*A + (visc+deltam)*B)/(40.0*Area) + NBB;
		KBB12 = 81.0*visc*C/(80.0*Area);
		KBB22 = 81.0*((2*visc+deltam)*B + (visc+deltam)*A)/(40.0*Area) + NBB;
		//
		//******** D eh o determinande da matriz (KBB + KBBdelta) 
		D = KBB11*KBB22 - KBB12*KBB12;
		//
		//****** Matrix InvKBBd **********
		InvKBBd11 = KBB22/D;
		InvKBBd12 = - KBB12/D;
		InvKBBd22 = KBB11/D;

		// *** Matriz Difusao Dinamica Macro - DD 
		double delta4area = delta / (4.0 * Area);
		
		Kdd11 = delta4area * (y23*y23 + x32*x32);
		Kdd12 = delta4area * (y23*y31 + x32*x13);
		Kdd13 = delta4area * (y23*y12 + x32*x21);
		Kdd22 = delta4area * (y31*y31 + x13*x13);
		Kdd23 = delta4area * (y31*y12 + x13*x21);
		Kdd33 = delta4area * (y12*y12 + x21*x21);

		// ======= Matrizes de estabilização advinda da micro escala ====
		// ******Analogous matrices of the SUPG formulation**************

		//*********Matrix Ndelta = NhB(KBB + KBBdelta)^{-1}NBh**********
	
		//*************** Matrix NBh ***************************
		NBh11 = rho*ninefortieth*(ux*y23 + uy*x32);
		NBh13 = rho*ninefortieth*(ux*y31 + uy*x13);
		NBh15 = rho*ninefortieth*(ux*y12 + uy*x21);
		

		//************** Matrix NhB = - NBh^T *******************		
		double NhB11, NhB31, NhB51;
		NhB11 = -NBh11;
		NhB31 = -NBh13;
		NhB51 = -NBh15;

		//***************Incremental Matrix NBh+ ****************
		double NBhp11, NBhp12, NBhp21, NBhp22; 		
		NBhp11 = rho*thirdtwentieth*Area*duxdx;
		NBhp12 = rho*thirdtwentieth*Area*duxdy;
		NBhp21 = rho*thirdtwentieth*Area*duydx;
		NBhp22 = rho*thirdtwentieth*Area*duydy;

		
		//******** D eh o determinande da matriz (KBB + KBBdelta) 
/*		D = ((2.0*visc+deltam)*(2.0*visc+deltam) + (visc+deltam)*(visc+deltam))*A*B + (2.0*visc+deltam)*(visc+deltam)*(A*A + B*B) - visc*visc*C*C/4.0;
		
		double coe = 1.;
		InvKBBd11 = coe*40.0*Area*((2.0*visc+deltam)*B + (visc+deltam)*A)/(81.0*D);
		InvKBBd12 = - coe*40.0*Area*visc*C/(81.0*D*2.0);
		InvKBBd22 = coe*40.0*Area*((2*visc+deltam)*A + (visc+deltam)*B)/(81.0*D);
*/
		//if(varglobal == 3 || varglobal == 14)
		//	printf("InvK11 = %f, InvK12 = %f, InvK22 = %f  \n", InvKBBd11, InvKBBd12, InvKBBd22);
		
		//****** Matriz Ndelta = NhB*InvKBB*NBh ******************
		Ndelta11 = NhB11*InvKBBd11*NBh11;
		Ndelta12 = NhB11*InvKBBd12*NBh11;
		Ndelta13 = NhB11*InvKBBd11*NBh13;
		Ndelta14 = NhB11*InvKBBd12*NBh13;
		Ndelta15 = NhB11*InvKBBd11*NBh15;
		Ndelta16 = NhB11*InvKBBd12*NBh15;
		
		//if(varglobal == 3 || varglobal == 14)
		//	printf("Ndelta11 = %e, Ndelta12 = %e, Ndelta13 = %e, Ndelta14 = %e, Ndelta15 = %e, Ndelta16 = %e  \n", Ndelta11, Ndelta12, Ndelta13, Ndelta14, Ndelta15, Ndelta16);

		Ndelta21 = Ndelta12;
		Ndelta22 = NhB11*InvKBBd22*NBh11;
		Ndelta23 = Ndelta14; 		
		Ndelta24 = NhB11*InvKBBd22*NBh13;
		Ndelta25 = Ndelta16; 		
		Ndelta26 = NhB11*InvKBBd22*NBh15;

		Ndelta31 = NhB31*InvKBBd11*NBh11;
		Ndelta32 = NhB31*InvKBBd12*NBh11;
		Ndelta33 = NhB31*InvKBBd11*NBh13;
		Ndelta34 = NhB31*InvKBBd12*NBh13;
		Ndelta35 = NhB31*InvKBBd11*NBh15;
		Ndelta36 = NhB31*InvKBBd12*NBh15;
		
		Ndelta41 = Ndelta32;
		Ndelta42 = NhB31*InvKBBd22*NBh11;
		Ndelta43 = Ndelta34;
		Ndelta44 = NhB31*InvKBBd22*NBh13;
		Ndelta45 = Ndelta36;		
		Ndelta46 = NhB31*InvKBBd22*NBh15;

		Ndelta51 = NhB51*InvKBBd11*NBh11;
		Ndelta52 = NhB51*InvKBBd12*NBh11;
		Ndelta53 = NhB51*InvKBBd11*NBh13;
		Ndelta54 = NhB51*InvKBBd12*NBh13;
		Ndelta55 = NhB51*InvKBBd11*NBh15;
		Ndelta56 = NhB51*InvKBBd12*NBh15;
	
		Ndelta61 = Ndelta52;
		Ndelta62 = NhB51*InvKBBd22*NBh11;
		Ndelta63 = Ndelta54;
		Ndelta64 = NhB51*InvKBBd22*NBh13;
		Ndelta65 = Ndelta56;
		Ndelta66 = NhB51*InvKBBd22*NBh15;

		//****** Incremental Matriz Ndelta+ = NhB*InvKBB*NBhp **********
		Ndeltap11 = NhB11*InvKBBd11*NBhp11 + NhB11*InvKBBd12*NBhp21;
		Ndeltap12 = NhB11*InvKBBd11*NBhp12 + NhB11*InvKBBd12*NBhp22;
		Ndeltap21 = NhB11*InvKBBd12*NBhp11 + NhB11*InvKBBd22*NBhp21;
		Ndeltap22 = NhB11*InvKBBd12*NBhp12 + NhB11*InvKBBd22*NBhp22;

		Ndeltap31 = NhB31*InvKBBd11*NBhp11 + NhB31*InvKBBd12*NBhp21;
		Ndeltap32 = NhB31*InvKBBd11*NBhp12 + NhB31*InvKBBd12*NBhp22;
		Ndeltap41 = NhB31*InvKBBd12*NBhp11 + NhB31*InvKBBd22*NBhp21;
		Ndeltap42 = NhB31*InvKBBd12*NBhp12 + NhB31*InvKBBd22*NBhp22;

		Ndeltap51 = NhB51*InvKBBd11*NBhp11 + NhB51*InvKBBd12*NBhp21;
		Ndeltap52 = NhB51*InvKBBd11*NBhp12 + NhB51*InvKBBd12*NBhp22;
		Ndeltap61 = NhB51*InvKBBd12*NBhp11 + NhB51*InvKBBd22*NBhp21;
		Ndeltap62 = NhB51*InvKBBd12*NBhp12 + NhB51*InvKBBd22*NBhp22;

		//*********Matrix Gdelta = NhB(KBB + KBBdelta)^{-1}GBh**********
		//		
		//************** Matrix GBh ****************************
		// Estah sendo calculada acima, junto com a vel na micro escala		
		
		//*************Matrix Gdelta****************************
		Gdelta11 = NhB11*InvKBBd11*GBh11 + NhB11*InvKBBd12*GBh21;
		Gdelta12 = NhB11*InvKBBd11*GBh12 + NhB11*InvKBBd12*GBh22;
		Gdelta13 = NhB11*InvKBBd11*GBh13 + NhB11*InvKBBd12*GBh23;
		
		Gdelta21 = NhB11*InvKBBd12*GBh11 + NhB11*InvKBBd22*GBh21;
		Gdelta22 = NhB11*InvKBBd12*GBh12 + NhB11*InvKBBd22*GBh22;
		Gdelta23 = NhB11*InvKBBd12*GBh13 + NhB11*InvKBBd22*GBh23;
			
		Gdelta31 = NhB31*InvKBBd11*GBh11 + NhB31*InvKBBd12*GBh21;
		Gdelta32 = NhB31*InvKBBd11*GBh12 + NhB31*InvKBBd12*GBh22;
		Gdelta33 = NhB31*InvKBBd11*GBh13 + NhB31*InvKBBd12*GBh23;
		
		Gdelta41 = NhB31*InvKBBd12*GBh11 + NhB31*InvKBBd22*GBh21;
		Gdelta42 = NhB31*InvKBBd12*GBh12 + NhB31*InvKBBd22*GBh22;
		Gdelta43 = NhB31*InvKBBd12*GBh13 + NhB31*InvKBBd22*GBh23;
		
		Gdelta51 = NhB51*InvKBBd11*GBh11 + NhB51*InvKBBd12*GBh21;
		Gdelta52 = NhB51*InvKBBd11*GBh12 + NhB51*InvKBBd12*GBh22;
		Gdelta53 = NhB51*InvKBBd11*GBh13 + NhB51*InvKBBd12*GBh23;
		
		Gdelta61 = NhB51*InvKBBd12*GBh11 + NhB51*InvKBBd22*GBh21;
		Gdelta62 = NhB51*InvKBBd12*GBh12 + NhB51*InvKBBd22*GBh22;
		Gdelta63 = NhB51*InvKBBd12*GBh13 + NhB51*InvKBBd22*GBh23;
		
		//*******Analogous matrices of the PSPG formulation**************

		//**********************Matrix Mphi******************************
		//Mphi11 = tau*sixth*y23;
		//Mphi12 = tau*sixth*x32;
	
		//Mphi21 = tau*sixth*y31;
		//Mphi22 = tau*sixth*x13;
	
		//Mphi31 = tau*sixth*y12;
		//Mphi32 = tau*sixth*x21;

		//** Matrix Nphi = GhB(KBB+KBBdelta)^{-1}NBh = - Gdelta^T*************
		Nphi11 = - Gdelta11;		
		Nphi12 = - Gdelta21;		
		Nphi13 = - Gdelta31;		
		Nphi14 = - Gdelta41;		
		Nphi15 = - Gdelta51;		
		Nphi16 = - Gdelta61;		
		
		Nphi21 = - Gdelta12;		
		Nphi22 = - Gdelta22;		
		Nphi23 = - Gdelta32;		
		Nphi24 = - Gdelta42;		
		Nphi25 = - Gdelta52;		
		Nphi26 = - Gdelta62;		
		
		Nphi31 = - Gdelta13;		
		Nphi32 = - Gdelta23;		
		Nphi33 = - Gdelta33;		
		Nphi34 = - Gdelta43;		
		Nphi35 = - Gdelta53;		
		Nphi36 = - Gdelta63;	

		//******** Incremental Matrix Nphi+ ****************************
		Nphip11 = (GBh11*InvKBBd11 + GBh21*InvKBBd12)*NBhp11 + (GBh11*InvKBBd12 + GBh21*InvKBBd22)*NBhp21;
		Nphip12 = (GBh11*InvKBBd11 + GBh21*InvKBBd12)*NBhp12 + (GBh11*InvKBBd12 + GBh21*InvKBBd22)*NBhp22;

		Nphip21 = (GBh12*InvKBBd11 + GBh22*InvKBBd12)*NBhp11 + (GBh12*InvKBBd12 + GBh22*InvKBBd22)*NBhp21;
		Nphip22 = (GBh12*InvKBBd11 + GBh22*InvKBBd12)*NBhp12 + (GBh12*InvKBBd12 + GBh22*InvKBBd22)*NBhp22;

		Nphip31 = (GBh13*InvKBBd11 + GBh23*InvKBBd12)*NBhp11 + (GBh13*InvKBBd12 + GBh23*InvKBBd22)*NBhp21;
		Nphip32 = (GBh13*InvKBBd11 + GBh23*InvKBBd12)*NBhp12 + (GBh13*InvKBBd12 + GBh23*InvKBBd22)*NBhp22;
		
		//** Matrix Gphi = GhB(KBB+KBBdelta)^{-1}GBh, with GhB = GBh^T*************
 	
		Gphi11 = GBh11*(GBh11*InvKBBd11 + GBh21*InvKBBd12) + GBh21*(GBh11*InvKBBd12 + GBh21*InvKBBd22);
		Gphi12 = GBh12*(GBh11*InvKBBd11 + GBh21*InvKBBd12) + GBh22*(GBh11*InvKBBd12 + GBh21*InvKBBd22);
		Gphi13 = GBh13*(GBh11*InvKBBd11 + GBh21*InvKBBd12) + GBh23*(GBh11*InvKBBd12 + GBh21*InvKBBd22);
		
		Gphi21 = GBh11*(GBh12*InvKBBd11 + GBh22*InvKBBd12) + GBh21*(GBh12*InvKBBd12 + GBh22*InvKBBd22);
		Gphi22 = GBh12*(GBh12*InvKBBd11 + GBh22*InvKBBd12) + GBh22*(GBh12*InvKBBd12 + GBh22*InvKBBd22);
		Gphi23 = GBh13*(GBh12*InvKBBd11 + GBh22*InvKBBd12) + GBh23*(GBh12*InvKBBd12 + GBh22*InvKBBd22);
		
		Gphi31 = GBh11*(GBh13*InvKBBd11 + GBh23*InvKBBd12) + GBh21*(GBh13*InvKBBd12 + GBh23*InvKBBd22);
		Gphi32 = GBh12*(GBh13*InvKBBd11 + GBh23*InvKBBd12) + GBh22*(GBh13*InvKBBd12 + GBh23*InvKBBd22);
		Gphi33 = GBh13*(GBh13*InvKBBd11 + GBh23*InvKBBd12) + GBh23*(GBh13*InvKBBd12 + GBh23*InvKBBd22);

		//*****************Font of the Multiscale formulation analougus to SUPG****
		fdelta1 = - (NhB11*InvKBBd11*f1B + NhB11*InvKBBd12*f2B);
		fdelta2 = - (NhB11*InvKBBd12*f1B + NhB11*InvKBBd22*f2B);
		fdelta3 = - (NhB31*InvKBBd11*f1B + NhB31*InvKBBd12*f2B);
		fdelta4 = - (NhB31*InvKBBd12*f1B + NhB31*InvKBBd22*f2B);
		fdelta5 = - (NhB51*InvKBBd11*f1B + NhB51*InvKBBd12*f2B);
		fdelta6 = - (NhB51*InvKBBd12*f1B + NhB51*InvKBBd22*f2B);

		//*****************Font of the Multiscale formulation analougus to PSPG ****
		fphi1 = - ((GBh11*InvKBBd11 + GBh21*InvKBBd12)*f1B + (GBh11*InvKBBd12 + GBh21*InvKBBd22)*f2B) ;
		fphi2 = - ((GBh12*InvKBBd11 + GBh22*InvKBBd12)*f1B + (GBh12*InvKBBd12 + GBh22*InvKBBd22)*f2B);
		fphi3 = - ((GBh13*InvKBBd11 + GBh23*InvKBBd12)*f1B + (GBh13*InvKBBd12 + GBh23*InvKBBd22)*f2B);

		//printf("Fonte fphi_1 = %E \n", fdelta1);
		//*******************************vector Fe***************************
		Fe[0] = f1 + fdelta1;
		Fe[1] = f2 + fdelta2;
		Fe[2] = fphi1;
		Fe[3] = f1 + fdelta3;
		Fe[4] = f2 + fdelta4;
		Fe[5] = fphi2;
		Fe[6] = f1 + fdelta5;
		Fe[7] = f2 + fdelta6;
		Fe[8] = fphi3;
		
		//*************Calculation of the residue continue*******************
		Kddv[0] = 0.5*delta*(y23*duxdx + x32*duxdy);
		Kddv[1] = 0.5*delta*(y23*duydx + x32*duydy);
		Kddv[2] = 0.5*delta*(y31*duxdx + x13*duxdy);
		Kddv[3] = 0.5*delta*(y31*duydx + x13*duydy);
		Kddv[4] = 0.5*delta*(y12*duxdx + x21*duxdy);
		Kddv[5] = 0.5*delta*(y12*duydx + x21*duydy);  

		Ndeltav[0] = Ndelta11*ux1 +  Ndelta12*uy1 +  Ndelta13*ux2 +  Ndelta14*uy2 +  Ndelta15*ux3 +  Ndelta16*uy3;
		Ndeltav[1] = Ndelta21*ux1 +  Ndelta22*uy1 +  Ndelta23*ux2 +  Ndelta24*uy2 +  Ndelta25*ux3 +  Ndelta26*uy3;
		Ndeltav[2] = Ndelta31*ux1 +  Ndelta32*uy1 +  Ndelta33*ux2 +  Ndelta34*uy2 +  Ndelta35*ux3 +  Ndelta36*uy3;
		Ndeltav[3] = Ndelta41*ux1 +  Ndelta42*uy1 +  Ndelta43*ux2 +  Ndelta44*uy2 +  Ndelta45*ux3 +  Ndelta46*uy3;
		Ndeltav[4] = Ndelta51*ux1 +  Ndelta52*uy1 +  Ndelta53*ux2 +  Ndelta54*uy2 +  Ndelta55*ux3 +  Ndelta56*uy3;
		Ndeltav[5] = Ndelta61*ux1 +  Ndelta62*uy1 +  Ndelta63*ux2 +  Ndelta64*uy2 +  Ndelta65*ux3 +  Ndelta66*uy3;

		Gdeltav[0] = Gdelta11*p1 + Gdelta12*p2 + Gdelta13*p3;
		Gdeltav[1] = Gdelta21*p1 + Gdelta22*p2 + Gdelta23*p3;
		Gdeltav[2] = Gdelta31*p1 + Gdelta32*p2 + Gdelta33*p3;
		Gdeltav[3] = Gdelta41*p1 + Gdelta42*p2 + Gdelta43*p3;
		Gdeltav[4] = Gdelta51*p1 + Gdelta52*p2 + Gdelta53*p3;
		Gdeltav[5] = Gdelta61*p1 + Gdelta62*p2 + Gdelta63*p3;

		Nphiv[0] = Nphi11*ux1 + Nphi12*uy1 + Nphi13*ux2 + Nphi14*uy2 + Nphi15*ux3 + Nphi16*uy3;
		Nphiv[1] = Nphi21*ux1 + Nphi22*uy1 + Nphi23*ux2 + Nphi24*uy2 + Nphi25*ux3 + Nphi26*uy3;
		Nphiv[2] = Nphi31*ux1 + Nphi32*uy1 + Nphi33*ux2 + Nphi34*uy2 + Nphi35*ux3 + Nphi36*uy3;

		Gphiv[0] = Gphi11*p1 + Gphi12*p2 + Gphi13*p3; 
		Gphiv[1] = Gphi21*p1 + Gphi22*p2 + Gphi23*p3;
		Gphiv[2] = Gphi31*p1 + Gphi32*p2 + Gphi33*p3;

		//*****************Residuo do metodo de Newton*****************
		
		Res[0] = Fe[0] - (Nv[0] - Ndeltav[0] + Kv[0] + Ksv[0] + Kddv[0] - ( Gv[0] - Gdeltav[0] ) );
		Res[1] = Fe[1] - (Nv[1] - Ndeltav[1] + Kv[1] + Ksv[1] + Kddv[1] - ( Gv[1] - Gdeltav[1] ) );
		Res[2] = Fe[2] - ( GTv[0] + Cv[0] - Nphiv[0] + Gphiv[0] );
		Res[3] = Fe[3] - (Nv[2] - Ndeltav[2] + Kv[2] + Ksv[2] + Kddv[2] - ( Gv[2] - Gdeltav[2] ) );
		Res[4] = Fe[4] - (Nv[3] - Ndeltav[3] + Kv[3] + Ksv[3] + Kddv[3] - ( Gv[3] - Gdeltav[3] ) );
		Res[5] = Fe[5] - ( GTv[1] + Cv[1] - Nphiv[1] + Gphiv[1] );
		Res[6] = Fe[6] - (Nv[4] - Ndeltav[4] + Kv[4] + Ksv[4] + Kddv[4] - ( Gv[4] - Gdeltav[4] ) );
		Res[7] = Fe[7] - (Nv[5] - Ndeltav[5] + Kv[5] + Ksv[5] + Kddv[5] - ( Gv[5] - Gdeltav[5] ) );
		Res[8] = Fe[8] - ( GTv[2] + Cv[2] - Nphiv[2] + Gphiv[2] );
/*
		Res[0] = Fe[0] + Nv[0];
		Res[1] = Fe[1] + Nv[1];
		Res[2] = Fe[2];
		Res[3] = Fe[3] + Nv[2];
		Res[4] = Fe[4] + Nv[3];
		Res[5] = Fe[5];
		Res[6] = Fe[6] + Nv[4];
		Res[7] = Fe[7] + Nv[5];
		Res[8] = Fe[8];
*/				
		//**********************Tangent matrix*******************************
				
		if(varglobal < 5e+40){ //varglobal (<5 ISI, Ponto Fixo) (>=5 NI. Met. de Newton)

		Ke[0][0] = N11 - Ndelta11 + K11 + Ks11 + Kdd11;
		Ke[0][1] =       K12 + Ks12;
		Ke[0][2] = -( GT11 - Gdelta11 );
		Ke[0][3] = N13 - Ndelta13 + K13 + Ks13 + Kdd12;
		Ke[0][4] =		  K14 + Ks14;
		Ke[0][5] = -( GT11 - Gdelta12 );
		Ke[0][6] = N15 - Ndelta15 + K15 + Ks15 + Kdd13;
		Ke[0][7] = 		  K16 + Ks16;
		Ke[0][8] = -( GT11 - Gdelta13 );
	
		Ke[1][0] = 		  K12 + Ks21;
		Ke[1][1] = N11 - Ndelta22 + K22 + Ks22 + Kdd11;
		Ke[1][2] = -( GT12 - Gdelta21 );
		Ke[1][3] = 		  K23 + Ks23;
		Ke[1][4] = N13 - Ndelta24 + K24 + Ks24 + Kdd12;
		Ke[1][5] = -( GT12 - Gdelta22 );
		Ke[1][6] = 		  K25 + Ks25;
		Ke[1][7] = N15 - Ndelta26 + K26 + Ks26 + Kdd13;
		Ke[1][8] = -( GT12 - Gdelta23 );
		
		Ke[2][0] = GT11 - Nphi11;
		Ke[2][1] = GT12 - Nphi12;
		Ke[2][2] = Gphi11 + C11;
		Ke[2][3] = GT13 - Nphi13;
		Ke[2][4] = GT14 - Nphi14;
		Ke[2][5] = Gphi12 + C12;
		Ke[2][6] = GT15 - Nphi15;
		Ke[2][7] = GT16 - Nphi16;
		Ke[2][8] = Gphi13 + C12;

		Ke[3][0] = N11 - Ndelta13 + K13 + Ks31 + Kdd12;
		Ke[3][1] =        K23 + Ks32;
		Ke[3][2] = -( GT13 - Gdelta31 );
		Ke[3][3] = N13 - Ndelta33 + K33 + Ks33 + Kdd22;
		Ke[3][4] =	  K34 + Ks34;
		Ke[3][5] = -( GT13 - Gdelta32 );
		Ke[3][6] = N15 - Ndelta35 + K35 + Ks35 + Kdd23;
		Ke[3][7] = 	  K36 + Ks36;
		Ke[3][8] = -( GT13 - Gdelta33 );

		Ke[4][0] = 	 K14 + Ks41;
		Ke[4][1] = N11 - Ndelta24 + K24 + Ks42 + Kdd12;
		Ke[4][2] = -( GT14 - Gdelta41 );
		Ke[4][3] = 	 K34 + Ks43;
		Ke[4][4] = N13 - Ndelta44 + K44 + Ks44 + Kdd22;
		Ke[4][5] = -( GT14 - Gdelta42 );
		Ke[4][6] = 	 K45 + Ks45;
		Ke[4][7] = N15 - Ndelta46 + K46 + Ks46 + Kdd23;
		Ke[4][8] = -( GT14 - Gdelta43 );

		Ke[5][0] = GT11 - Nphi21;
		Ke[5][1] = GT12 - Nphi22;
		Ke[5][2] = Gphi12 + C12;
		Ke[5][3] = GT13 - Nphi23;
		Ke[5][4] = GT14 - Nphi24;
		Ke[5][5] = Gphi22 + C11;
		Ke[5][6] = GT15 - Nphi25;
		Ke[5][7] = GT16 - Nphi26;
		Ke[5][8] = Gphi23 + C12;

		Ke[6][0] = N11 - Ndelta15 + K15 + Ks51 + Kdd13;
		Ke[6][1] =       K25 + Ks52;
		Ke[6][2] = -( GT15 - Gdelta51 );
		Ke[6][3] = N13 + Ndelta35 + K35 + Ks53 + Kdd23;
		Ke[6][4] =	 K45 + Ks54;
		Ke[6][5] = -( GT15 - Gdelta52 );
		Ke[6][6] = N15 - Ndelta55 + K55 + Ks55 + Kdd33;
		Ke[6][7] = 	 K56 + Ks56;
		Ke[6][8] = -( GT15 - Gdelta53 );

		Ke[7][0] = 	 K16 + Ks61;
		Ke[7][1] = N11 - Ndelta26 + K26 + Ks62 + Kdd13;
		Ke[7][2] = -( GT16 - Gdelta61 );
		Ke[7][3] = 	 K36 + Ks63;
		Ke[7][4] = N13 - Ndelta46 + K46 + Ks64 + Kdd23;
		Ke[7][5] = -( GT16 - Gdelta62 );
		Ke[7][6] = 	 K56 + Ks65;
		Ke[7][7] = N15 - Ndelta66 + K66 + Ks66 + Kdd33;
		Ke[7][8] = -( GT16 - Gdelta63 );

		Ke[8][0] = GT11 - Nphi31;
		Ke[8][1] = GT12 - Nphi32;
		Ke[8][2] = Gphi13 + C12;
		Ke[8][3] = GT13 - Nphi33;
		Ke[8][4] = GT14 - Nphi34;
		Ke[8][5] = Gphi23 + C12;
		Ke[8][6] = GT15 - Nphi35;
		Ke[8][7] = GT16 - Nphi36;
		Ke[8][8] = Gphi33 + C11;
		
		}else{

		Ke[0][0] = N11 + Np11 - Ndelta11 - Ndeltap11 + K11 + Ks11 + Kdd11;
		Ke[0][1] =	Np12 - Ndeltap12 + K12 + Ks12;
		Ke[0][2] = -( GT11 - Gdelta11 );
		Ke[0][3] = N13 + Np31 - Ndelta13 - Ndeltap11 + K13 + Ks13 + Kdd12;
		Ke[0][4] =	Np32 - Ndeltap12 + K14 + Ks14;
		Ke[0][5] = -( GT11 - Gdelta12 );
		Ke[0][6] = N15 + Np31 - Ndelta15 - Ndeltap11 + K15 + Ks15 + Kdd13;
		Ke[0][7] = 	Np32 - Ndeltap12 + K16 + Ks16;
		Ke[0][8] = -( GT11 - Gdelta13 );
	
		Ke[1][0] =	Np21 - Ndeltap21 + K12 + Ks21;
		Ke[1][1] = N11 + Np22 - Ndelta22 - Ndeltap22 + K22 + Ks22 + Kdd11;
		Ke[1][2] = -( GT12 - Gdelta21 );
		Ke[1][3] =	Np41 - Ndeltap21 + K23 + Ks23;
		Ke[1][4] = N13 + Np42 - Ndelta24 - Ndeltap22 + K24 + Ks24 + Kdd12;
		Ke[1][5] = -( GT12 - Gdelta22 );
		Ke[1][6] = 	Np41 - Ndeltap21 + K25 + Ks25;
		Ke[1][7] = N15 + Np42 - Ndelta26 - Ndeltap22 + K26 + Ks26 + Kdd13;
		Ke[1][8] = -( GT12 - Gdelta23 );
		
		Ke[2][0] = GT11 - Nphi11 - Nphip11;
		Ke[2][1] = GT12 - Nphi12 - Nphip12;
		Ke[2][2] = Gphi11 + C11;
		Ke[2][3] = GT13 - Nphi13 - Nphip11;
		Ke[2][4] = GT14 - Nphi14 - Nphip12;
		Ke[2][5] = Gphi12 + C12;
		Ke[2][6] = GT15 - Nphi15 - Nphip11;
		Ke[2][7] = GT16 - Nphi16 - Nphip12;
		Ke[2][8] = Gphi13 + C12;

		Ke[3][0] = N11 + Np31- Ndelta13 - Ndeltap31 + K13 + Ks31 + Kdd12;
		Ke[3][1] =	Np32 - Ndeltap32 + K23 + Ks32;
		Ke[3][2] = -( GT13 - Gdelta31 );
		Ke[3][3] = N13 + Np11 - Ndelta33 - Ndeltap31 + K33 + Ks33 + Kdd22;
		Ke[3][4] =	Np12 - Ndeltap32 + K34 + Ks34;
		Ke[3][5] = -( GT13 - Gdelta32 );
		Ke[3][6] = N15 + Np31 - Ndelta35 - Ndeltap31 + K35 + Ks35 + Kdd23;
		Ke[3][7] =	Np32 - Ndeltap32 + K36 + Ks36;
		Ke[3][8] = -( GT13 - Gdelta33 );

		Ke[4][0] =	Np41 - Ndeltap41 + K14 + Ks41;
		Ke[4][1] = N11 + Np42 - Ndelta24 - Ndeltap42 + K24 + Ks42 + Kdd12;
		Ke[4][2] = -( GT14 - Gdelta41 );
		Ke[4][3] =	Np21 - Ndeltap41 + K34 + Ks43;
		Ke[4][4] = N13 + Np22 - Ndelta44 - Ndeltap42 + K44 + Ks44 + Kdd22;
		Ke[4][5] = -( GT14 - Gdelta42 );
		Ke[4][6] =	Np41 - Ndeltap41 + K45 + Ks45;
		Ke[4][7] = N15 + Np42 - Ndelta46 - Ndeltap42 + K46 + Ks46 + Kdd23;
		Ke[4][8] = -( GT14 - Gdelta43 );

		Ke[5][0] = GT11 - Nphi21 - Nphip21;
		Ke[5][1] = GT12 - Nphi22 - Nphip22;
		Ke[5][2] = Gphi12 + C12;
		Ke[5][3] = GT13 - Nphi23 - Nphip21;
		Ke[5][4] = GT14 - Nphi24 - Nphip22;
		Ke[5][5] = Gphi22 + C11;
		Ke[5][6] = GT15 - Nphi25 - Nphip21;
		Ke[5][7] = GT16 - Nphi26 - Nphip22;
		Ke[5][8] = Gphi23 + C12;

		Ke[6][0] = N11 + Np31 - Ndelta15 - Ndeltap51 + K15 + Ks51 + Kdd13;
		Ke[6][1] =	Np32 - Ndeltap52 + K25 + Ks52;
		Ke[6][2] = -( GT15 - Gdelta51 );
		Ke[6][3] = N13 + Np31 + Ndelta35 - Ndeltap51 + K35 + Ks53 + Kdd23;
		Ke[6][4] =	Np32 - Ndeltap52 + K45 + Ks54;
		Ke[6][5] = -( GT15 - Gdelta52 );
		Ke[6][6] = N15 + Np11 - Ndelta55 - Ndeltap51 + K55 + Ks55 + Kdd33;
		Ke[6][7] =	Np12 - Ndeltap52 + K56 + Ks56;
		Ke[6][8] = -( GT15 - Gdelta53 );

		Ke[7][0] =	Np41 - Ndeltap61 + K16 + Ks61;
		Ke[7][1] = N11 + Np42 - Ndelta26 - Ndeltap62 + K26 + Ks62 + Kdd13;
		Ke[7][2] = -( GT16 - Gdelta61 );
		Ke[7][3] =	Np41 - Ndeltap61 + K36 + Ks63;
		Ke[7][4] = N13 + Np42 - Ndelta46 - Ndeltap62 + K46 + Ks64 + Kdd23;
		Ke[7][5] = -( GT16 - Gdelta62 );
		Ke[7][6] =	Np21 - Ndeltap61 + K56 + Ks65;
		Ke[7][7] = N15 + Np22 - Ndelta66 - Ndeltap62  + K66 + Ks66 + Kdd33;
		Ke[7][8] = -( GT16 - Gdelta63 );

		Ke[8][0] = GT11 - Nphi31 - Nphip31;
		Ke[8][1] = GT12 - Nphi32 - Nphip32;
		Ke[8][2] = Gphi13 + C12;
		Ke[8][3] = GT13 - Nphi33 - Nphip31;
		Ke[8][4] = GT14 - Nphi34 - Nphip32;
		Ke[8][5] = Gphi23 + C12;
		Ke[8][6] = GT15 - Nphi35 - Nphip31;
		Ke[8][7] = GT16 - Nphi36 - Nphip32;
		Ke[8][8] = Gphi33 + C11;
		
		}



		/* direciona para um lixo durante o produto matriz vetor os elementos da matriz
		 que deveriam ser nulos (os que foram para o vetor força devido a c.c.). 
		Nos da a posicao de assemble global do vetor forca
		Aqui neq == numequacao + 1. */
/*		for (i = 0; i<NDOF; i++){
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
*/
		// Assemble global do vetor independente F de Au=F 
		for (i = 0; i < 9; i++)
			F[lm[E][i]] +=  Res[i];
		
		F[neq] = 0;


		
		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, E, Ke);
		
	}//for elemento
	
	printf("\n Qt delta_m atualizado = %d \n", conta2);
	if(varglobal==0)
		printf("\n Re Mínimo = %f, Máximo = %f, Médio = %f \n", ReLocMin, ReLocMax, ReLocMed);		
		
	myfree(U);		
	//printf(" \n\n VARGLOBAL = %d \n\n", varglobal);
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

// nao tem mais	end_assembly(Parameters, 0, MatrixData);

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




	return 0;

}
