#include "SSTransportEquation3D.h"
#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_CAU(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int e, i;
	int nel, neq;
	double X[4], Y[4], Z[4], a[4], b[4], c[4], d[4], xb, yb, zb;
	double sixV, kappaX, kappaY, kappaZ, normK, betaX, betaY, betaZ, normBeta, sigma;
	double peclet, xi, tau;
	double Volume, MaxVolume, Ke[4][4], Fe[4];
	double h, hc;
	//double lixo1[4], lixo2[Parameters->nnodes], lixo3; // Variáveis usadas pra descartar valores na função Element_Configuration
	double *F = FemStructs->F;
	double Cst = Parameters->ConstApli;
	double Vmax = Parameters->Vmax;
	double xSup = Parameters->xSup;
	double ySup = Parameters->ySup;
	double *uSpace, Ue[4], uB;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	CoefFormFuncType *CFF = FemStructs->CFF;
	
	nel = Parameters->nel;
	neq = Parameters->neq;
	MaxVolume = 0.0;
	
	double *PeL, MedPeL = 0.0, MaxPeL = 0.0, MinPeL = DBL_MAX; // convection-diffusion realtion
	double *RDRL, MedRDRL = 0.0, MaxRDRL = 0.0, MinRDRL = DBL_MAX; // reaction-diffusion relation
	
	PeL = (double*) mycalloc("PeL of 'Build_DD'", nel, sizeof(double));
	dzero(nel, PeL);

	RDRL = (double*) mycalloc("RDRL of 'Build_DD'", nel, sizeof(double));
	dzero(nel, RDRL);

	dzero(neq+1, F);

	setZeros(Parameters,MatrixData);

	uSpace = (double*) mycalloc("uSpace of 'Build_CAU'", Parameters->nnodes, sizeof(double));
	eval_U_Space(Parameters, FemStructs, FemFunctions, uSpace);

	for(e = 0; e < nel; e++){
		
		//printf("---- Elemento %d ----\n", e);
		
		Volume = Element_Configuration(e, Element, Node, X, Y, Z, a, b, c, d, Ue, uSpace, &uB);
		CFF[e].volume = Volume;
		sixV = 6.0*Volume;
		if(Volume > MaxVolume){
			MaxVolume = Volume;
		}
		
		// baricentric coordinate
		xb = (X[0] + X[1] + X[2] + X[3])/4.0;
		yb = (Y[0] + Y[1] + Y[2] + Y[3])/4.0;
		zb = (Z[0] + Z[1] + Z[2] + Z[3])/4.0;
		
		//printf("\n\n VOLUME = %lf\n",Volume);
		
		for(i = 0; i < 4; i++){
			CFF[e].a[i] = a[i];
			CFF[e].b[i] = b[i];
			CFF[e].c[i] = c[i];
			CFF[e].d[i] = d[i];
		}	
		
		double b1b1, b2b2, b3b3, b4b4, c1c1, c2c2, c3c3, c4c4, d1d1, d2d2, d3d3, d4d4;
		double b1b2, b1b3, b1b4, b2b3, b2b4, b3b4, c1c2, c1c3, c1c4, c2c3, c2c4, c3c4, d1d2, d1d3, d1d4, d2d3, d2d4, d3d4;
		
		b1b1 = b[0]*b[0]; c1c1 = c[0]*c[0]; d1d1 = d[0]*d[0];
		b2b2 = b[1]*b[1]; c2c2 = c[1]*c[1]; d2d2 = d[1]*d[1];
		b3b3 = b[2]*b[2]; c3c3 = c[2]*c[2]; d3d3 = d[2]*d[2];
		b4b4 = b[3]*b[3]; c4c4 = c[3]*c[3]; d4d4 = d[3]*d[3];
		
		b1b2 = b[0]*b[1]; c1c2 = c[0]*c[1]; d1d2 = d[0]*d[1];
		b1b3 = b[0]*b[2]; c1c3 = c[0]*c[2]; d1d3 = d[0]*d[2];
		b1b4 = b[0]*b[3]; c1c4 = c[0]*c[3]; d1d4 = d[0]*d[3];
		b2b3 = b[1]*b[2]; c2c3 = c[1]*c[2]; d2d3 = d[1]*d[2];
		b2b4 = b[1]*b[3]; c2c4 = c[1]*c[3]; d2d4 = d[1]*d[3];
		b3b4 = b[2]*b[3]; c3c4 = c[2]*c[3]; d3d4 = d[2]*d[3];
		
		// Coefficients of the equation
		// diffusion
		FemFunctions->kappa(&kappaX, &kappaY, &kappaZ, xb, yb, zb, Cst);
		normK = sqrt(kappaX*kappaX + kappaY*kappaY + kappaZ*kappaZ);
		
		// convection
		FemFunctions->beta(&betaX, &betaY, &betaZ, xb, yb, zb, Cst, Vmax, xSup, ySup);
		normBeta = sqrt(betaX*betaX + betaY*betaY + betaZ*betaZ);
		
		// reaction
		sigma = FemFunctions->sigma(xb, yb, zb);
		
		// D = int grad(v). K grad(u) dOmega
		double E11, E12, E13, E14, E22, E23, E24, E33, E34, E44;
		double ConstD = 1.0/(36.0*Volume);
		
		E11 = kappaX*b1b1 + kappaY*c1c1 + kappaZ*d1d1;
		E12 = kappaX*b1b2 + kappaY*c1c2 + kappaZ*d1d2;
		E13 = kappaX*b1b3 + kappaY*c1c3 + kappaZ*d1d3;
		E14 = kappaX*b1b4 + kappaY*c1c4 + kappaZ*d1d4;
		E22 = kappaX*b2b2 + kappaY*c2c2 + kappaZ*d2d2;
		E23 = kappaX*b2b3 + kappaY*c2c3 + kappaZ*d2d3;
		E24 = kappaX*b2b4 + kappaY*c2c4 + kappaZ*d2d4;
		E33 = kappaX*b3b3 + kappaY*c3c3 + kappaZ*d3d3;
		E34 = kappaX*b3b4 + kappaY*c3c4 + kappaZ*d3d4;
		E44 = kappaX*b4b4 + kappaY*c4c4 + kappaZ*d4d4;
		
		// C = int v beta . grad(u) dOmega
		double C1, C2, C3, C4;
		double ConstC = 1.0/24.0;
		
		C1 = betaX*b[0] + betaY*c[0] + betaZ*d[0];
		C2 = betaX*b[1] + betaY*c[1] + betaZ*d[1];
		C3 = betaX*b[2] + betaY*c[2] + betaZ*d[2];
		C4 = betaX*b[3] + betaY*c[3] + betaZ*d[3];
		
		// R = int sigma v u dOmega
		double ConstR = (sigma*Volume)/20.0;
	
		double dXidx, dEtadx, dZetadx, dXidy, dEtady, dZetady, dXidz, dEtadz, dZetadz;
			
		dXidx = b[1]/sixV;
		dXidy = c[1]/sixV;
		dXidz = d[1]/sixV;
		dEtadx = b[2]/sixV;
		dEtady = c[2]/sixV;
		dEtadz = d[2]/sixV;
		dZetadx = b[3]/sixV;
		dZetady = c[3]/sixV;
		dZetadz = d[3]/sixV;
		
		// Paramêtro de malha
		double Be[3], normBe;
			
		// Be_i = SUM (dXi_i/dx_j beta_j)
		Be[0] = dXidx*betaX + dXidy*betaY + dXidz*betaZ;
		Be[1] = dEtadx*betaX + dEtady*betaY + dEtadz*betaZ;
		Be[2] = dZetadx*betaX + dZetady*betaY + dZetadz*betaZ;
			
		normBe = sqrt(Be[0]*Be[0] + Be[1]*Be[1] + Be[2]*Be[2]);
		
		/*if(e == 571){
			printf("dXidx = %E\n dXidy = %E\n dXidz = %E\n dEtadx = %E\n dEtady = %E\n dEtadz = %E\n dZetadx = %E\n dZetady = %E\n dZetadz = %E\n", dXidx, dXidy, dXidz, dEtadx, dEtady, dEtadz, dZetadx, dZetady, dZetadz);
			printf("betaX = %lf\n betaY = %lf\n betaZ = %lf\n Be[0] = %E\n Be[1] = %E\n Be[2] = %E\n", betaX, betaY, betaZ, Be[0], Be[1], Be[2]);
			getchar();
		}*/
		
		if (normBe < 1e-5){
			h = 0.0;
		}else{	
			h = (2*normBeta)/normBe;
		}
		
		//Tau CAU
		peclet = (normBeta*h)/(2.0*normK); 
		
		xi = fmax(0, 1-(1/peclet));  //cosh(peclet)/sinh(peclet) - 1.0/peclet; 
		tau = (h*xi)/(2.0*normBeta);
		
		// Csupg = int beta . grad(v) tau beta . grad(u) dOmega
		double Cs11, Cs12, Cs13, Cs14, Cs22, Cs23, Cs24, Cs33, Cs34, Cs44;
		double ConstCsupg = tau/(36.0*Volume);
		
		Cs12 = betaX*betaX*b[0]*b[1] + betaX*betaY*c[0]*b[1] + betaX*betaZ*d[0]*b[1] + betaX*betaY*b[0]*c[1] + betaY*betaY*c[0]*c[1] + betaY*betaZ*d[0]*c[1] + betaX*betaZ*b[0]*d[1] + betaY*betaZ*c[0]*d[1] + betaZ*betaZ*d[0]*d[1];
		Cs13 = betaX*betaX*b[0]*b[2] + betaX*betaY*c[0]*b[2] + betaX*betaZ*d[0]*b[2] + betaX*betaY*b[0]*c[2] + betaY*betaY*c[0]*c[2] + betaY*betaZ*d[0]*c[2] + betaX*betaZ*b[0]*d[2] + betaY*betaZ*c[0]*d[2] + betaZ*betaZ*d[0]*d[2];
		Cs14 = betaX*betaX*b[0]*b[3] + betaX*betaY*c[0]*b[3] + betaX*betaZ*d[0]*b[3] + betaX*betaY*b[0]*c[3] + betaY*betaY*c[0]*c[3] + betaY*betaZ*d[0]*c[3] + betaX*betaZ*b[0]*d[3] + betaY*betaZ*c[0]*d[3] + betaZ*betaZ*d[0]*d[3];
		Cs11 = -(Cs12 + Cs13 + Cs14);
		Cs23 = betaX*betaX*b[1]*b[2] + betaX*betaY*c[1]*b[2] + betaX*betaZ*d[1]*b[2] + betaX*betaY*b[1]*c[2] + betaY*betaY*c[1]*c[2] + betaY*betaZ*d[1]*c[2] + betaX*betaZ*b[1]*d[2] + betaY*betaZ*c[1]*d[2] + betaZ*betaZ*d[1]*d[2];
		Cs24 = betaX*betaX*b[1]*b[3] + betaX*betaY*c[1]*b[3] + betaX*betaZ*d[1]*b[3] + betaX*betaY*b[1]*c[3] + betaY*betaY*c[1]*c[3] + betaY*betaZ*d[1]*c[3] + betaX*betaZ*b[1]*d[3] + betaY*betaZ*c[1]*d[3] + betaZ*betaZ*d[1]*d[3];
		Cs22 = -(Cs12 + Cs23 + Cs24);
		Cs34 = betaX*betaX*b[2]*b[3] + betaX*betaY*c[2]*b[3] + betaX*betaZ*d[2]*b[3] + betaX*betaY*b[2]*c[3] + betaY*betaY*c[2]*c[3] + betaY*betaZ*d[2]*c[3] + betaX*betaZ*b[2]*d[3] + betaY*betaZ*c[2]*d[3] + betaZ*betaZ*d[2]*d[3];
		Cs33 = -(Cs13 + Cs23 + Cs34);
		Cs44 = -(Cs14 + Cs24 + Cs34);
		
		// Rsupg = int beta . grad(v) tau sigma u dOmega
		double ConstRsupg = (sigma*tau)/24.0;
		
		// delta CAU
		double betaH[3], pecletc, xic, ErrBeta[3], normErrBeta;
		double delta, gradU[3], normGradU, dotGradU, res, Fbari, modR;
		double ConstGradU = 1.0/sixV;
		
		gradU[0] = ConstGradU*(Ue[0]*b[0] + Ue[1]*b[1] + Ue[2]*b[2] + Ue[3]*b[3]);
		gradU[1] = ConstGradU*(Ue[0]*c[0] + Ue[1]*c[1] + Ue[2]*c[2] + Ue[3]*c[3]);
		gradU[2] = ConstGradU*(Ue[0]*d[0] + Ue[1]*d[1] + Ue[2]*d[2] + Ue[3]*d[3]);
		
		dotGradU = ddot(3, gradU, gradU); //gradU[0]*gradU[0] + gradU[1]*gradU[1] + gradU[2]*gradU[2];
		normGradU = sqrt(dotGradU); //sqrt(gradU[0]*gradU[0] + gradU[1]*gradU[1] + gradU[2]*gradU[2]);
		
		Fbari = FemFunctions->f(xb, yb, zb, Cst);
		
		res = betaX*gradU[0] + betaY*gradU[1] + betaZ*gradU[2] + sigma*uB - Fbari;
		modR = fabs(res);
		
		if(normGradU > 1e-5){
			betaH[0] = betaX - (modR/dotGradU)*gradU[0];
			betaH[1] = betaY - (modR/dotGradU)*gradU[1];
			betaH[2] = betaZ - (modR/dotGradU)*gradU[2];
		}else{
			betaH[0] = betaX;
			betaH[1] = betaY;
			betaH[2] = betaZ;
		}
		
		ErrBeta[0] = betaX - betaH[0];
		ErrBeta[1] = betaY - betaH[1];
		ErrBeta[2] = betaZ - betaH[2];
		
		normErrBeta = sqrt(ddot(3, ErrBeta, ErrBeta));

		double Bec[3], normBec;
		
		// Be_i = SUM (dXi_i/dx_j ErrBeta_j)
		Bec[0] = dXidx*ErrBeta[0] + dXidy*ErrBeta[1] + dXidz*ErrBeta[2];
		Bec[1] = dEtadx*ErrBeta[0] + dEtady*ErrBeta[1] + dEtadz*ErrBeta[2];
		Bec[2] = dZetadx*ErrBeta[0] + dZetady*ErrBeta[1] + dZetadz*ErrBeta[2];
		
		normBec = sqrt(ddot(3, Bec, Bec)); //sqrt(Be[0]*Be[0] + Be[1]*Be[1] + Be[2]*Be[2]);
		
		if (normBec < 1e-5){ // igual a zero
			hc = 0.0;
		}else{
			hc = (2*normErrBeta)/normBec;
		}
		
		pecletc = (hc*normErrBeta)/(2.0*normK);
		
		if (pecletc != 0.0){
			xic = fmax(0, 1-(1/pecletc));   //cosh(pecletc)/sinh(pecletc) - 1.0/pecletc;
		}else{
			xic = 0.0;
		}
		
		double TOL;
		if (xi < 1e-5){ // igual a zero
			TOL = 0.0;
		}else{
			TOL = (xic*hc)/(xi*h);;
		}
		
		double NUM;
		
		if(normGradU < 1e-5){
			NUM = 0.0;
		}else{
			NUM = modR/(normBeta*normGradU);
		}
		
		/*if(e == 571){
			printf("normBeta = %E\n normBe = %E\n h = %E \n normK = %E\n peclet = %E\n conta = %E\n",normBeta, normBe, h, normK, peclet, 1-(1/peclet));
			printf("xi = %E\n tau = %E\n", xi, tau);
			printf("normErrBeta = %E\n normBec = %E\n hc = %E \n pecletc = %E\n conta = %E\n",normErrBeta, normBec, hc, pecletc, 1-(1/pecletc));
			printf("xic = %E\n hc = %E\n xi = %E\n h = %E\n modR = %E\n normBeta = %E\n normGradU = %E\n", xic, hc, xi, h, modR, normBeta, normGradU);
			printf("TOL = %E\n NUM = %E\n", TOL, NUM);
			getchar();
		}*/
		
		if(NUM >= TOL){
			delta = 0.0;
			//printf(" if delta = %E\n", delta); getchar();
		}else{
			delta = (xi*h/2.0)*(TOL-NUM)*(modR/normGradU);
			//printf("else delta = %E\n", delta); getchar();
		}
		
		// Dcau = int delta grad(v) . grad(u) dOmega
		double Ec11, Ec12, Ec13, Ec14, Ec22, Ec23, Ec24, Ec33, Ec34, Ec44;
		double ConstDcau = delta/(36.0*Volume);

		Ec11 = b1b1 + c1c1 + d1d1;
		Ec12 = b1b2 + c1c2 + d1d2;
		Ec13 = b1b3 + c1c3 + d1d3;
		Ec14 = b1b4 + c1c4 + d1d4;
		Ec22 = b2b2 + c2c2 + d2d2;
		Ec23 = b2b3 + c2c3 + d2d3;
		Ec24 = b2b4 + c2c4 + d2d4;
		Ec33 = b3b3 + c3c3 + d3d3;
		Ec34 = b3b4 + c3c4 + d3d4;
		Ec44 = b4b4 + c4c4 + d4d4;
		
		// Elementary Matrix
		Ke[0][0] = ConstD*E11 + ConstDcau*Ec11 + ConstC*C1 + ConstCsupg*Cs11 + ConstRsupg*C1 + 2*ConstR; 
		Ke[1][0] = ConstD*E12 + ConstDcau*Ec12 + ConstC*C1 + ConstCsupg*Cs12 + ConstRsupg*C2 + ConstR; 
		Ke[2][0] = ConstD*E13 + ConstDcau*Ec13 + ConstC*C1 + ConstCsupg*Cs13 + ConstRsupg*C3 + ConstR; 
		Ke[3][0] = ConstD*E14 + ConstDcau*Ec14 + ConstC*C1 + ConstCsupg*Cs14 + ConstRsupg*C4 + ConstR; 
		                                                               
		Ke[0][1] = ConstD*E12 + ConstDcau*Ec12 + ConstC*C2 + ConstCsupg*Cs12 + ConstRsupg*C1 + ConstR; 
		Ke[1][1] = ConstD*E22 + ConstDcau*Ec22 + ConstC*C2 + ConstCsupg*Cs22 + ConstRsupg*C2 + 2*ConstR;            
		Ke[2][1] = ConstD*E23 + ConstDcau*Ec23 + ConstC*C2 + ConstCsupg*Cs23 + ConstRsupg*C3 + ConstR;            
		Ke[3][1] = ConstD*E24 + ConstDcau*Ec24 + ConstC*C2 + ConstCsupg*Cs24 + ConstRsupg*C4 + ConstR;            
		                                                                
		Ke[0][2] = ConstD*E13 + ConstDcau*Ec13 + ConstC*C3 + ConstCsupg*Cs13 + ConstRsupg*C1 + ConstR; 
		Ke[1][2] = ConstD*E23 + ConstDcau*Ec23 + ConstC*C3 + ConstCsupg*Cs23 + ConstRsupg*C2 + ConstR;              
		Ke[2][2] = ConstD*E33 + ConstDcau*Ec33 + ConstC*C3 + ConstCsupg*Cs33 + ConstRsupg*C3 + 2*ConstR;              
		Ke[3][2] = ConstD*E34 + ConstDcau*Ec34 + ConstC*C3 + ConstCsupg*Cs34 + ConstRsupg*C4 + ConstR;              
		                                                                
		Ke[0][3] = ConstD*E14 + ConstDcau*Ec14 + ConstC*C4 + ConstCsupg*Cs14 + ConstRsupg*C1 + ConstR;            
		Ke[1][3] = ConstD*E24 + ConstDcau*Ec24 + ConstC*C4 + ConstCsupg*Cs24 + ConstRsupg*C2 + ConstR;            
		Ke[2][3] = ConstD*E34 + ConstDcau*Ec34 + ConstC*C4 + ConstCsupg*Cs34 + ConstRsupg*C3 + ConstR;            
		Ke[3][3] = ConstD*E44 + ConstDcau*Ec44 + ConstC*C4 + ConstCsupg*Cs44 + ConstRsupg*C4 + 2*ConstR;
		
/*		printf("\nMatriz Ke E: %d\n", e);
		for(i = 0; i < 4; i++){
			printf("%E\t%E\t%E\t%E\n", Ke[i][0], Ke[i][1], Ke[i][2], Ke[i][3]);
		}
		getchar();*/
		                                                      
		// F = int v f dOmega 
		double f[4], Faux;
		double ConstF = Volume/20.0;
		
		for(i = 0; i < 4; i++){
			f[i] = FemFunctions->f(X[i], Y[i], Z[i], Cst);
		}

		Faux = f[0] + f[1] + f[2] + f[3];
		
		// Fsupg = int beta . grad(v) tau f dOmega
		double ConstFsupg = tau/24.0;
		
		Fe[0]  = ConstF*(Faux + f[0]) + ConstFsupg*C1*Faux;
		Fe[1]  = ConstF*(Faux + f[1]) + ConstFsupg*C2*Faux;
		Fe[2]  = ConstF*(Faux + f[2]) + ConstFsupg*C3*Faux;
		Fe[3]  = ConstF*(Faux + f[3]) + ConstFsupg*C4*Faux;			
		
/*		printf("\nVetor Fe E: %d\n", e);
		for(i = 0; i < 4; i++){
			printf("%E\n",Fe[i]);
		}
		getchar();*/
		
		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->MatrixAssembly(Parameters, MatrixData, FemStructs, e, Ke);

//		AssemblyGlobalMatrix(MatrixData, FemStructs, e, Ke, neq);
		
		// Font Assembly 
		font_assembly(e, Fe, Ke, FemFunctions, FemStructs, neq, Cst, xSup);
		
		PeL[e] = (sqrt(betaX*betaX + betaY*betaY + betaZ*betaZ) * pow(Volume,1.0/3.0)) / (2.0 * sqrt(kappaX*kappaX + kappaY*kappaY + kappaZ*kappaZ));
		RDRL[e] = (sigma * pow(Volume,2.0/3.0)) / sqrt(kappaX*kappaX + kappaY*kappaY + kappaZ*kappaZ);
		
	}//end for element
	
	// Keeping de maximum and minimun Peclet local
	for(e=0; e<nel; e++){
		MedPeL = MedPeL + PeL[e];

		if(PeL[e] > MaxPeL)
			MaxPeL = PeL[e];

		if(PeL[e] < MinPeL)
			MinPeL = PeL[e];
	}
	Parameters->MaxPecletLocal = MaxPeL;
	Parameters->MinPecletLocal = MinPeL;
	Parameters->MedPecletLocal = MedPeL/nel;

	// Keeping de maximum and minimun RDR local
	for(e=0; e<nel; e++){
		MedRDRL = MedRDRL + RDRL[e];

		if(RDRL[e] > MaxRDRL)
			MaxRDRL = RDRL[e];

		if(RDRL[e] < MinRDRL)
			MinRDRL = RDRL[e];
	}
	Parameters->MaxReaDifRelationLocal = MaxRDRL;
	Parameters->MinReaDifRelationLocal = MinRDRL;
	Parameters->MedReaDifRelationLocal = MedRDRL/nel;
	
	Parameters->MaxVolume = MaxVolume;
	
	myfree(uSpace);
	myfree(PeL);

	return 0;
		
}
