#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

double TauMacro(int e, ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *Ue)
{	
	
	int J1, J2, J3, J4, i;
	double tau=0.0, X[4], Y[4], Z[4], xb, yb, zb, sixV, b[4], c[4], d[4], uB;
	double kappaX, kappaY, kappaZ, betaX, betaY, betaZ, sigma, normKappa, normBeta;
	double ConstGradU, gradU[3], dotGradU, normGradU, res, modR, Fbari;
	double Cst = Parameters->ConstApli;
	double Vmax = Parameters->Vmax;
	double xSup = Parameters->xSup;
	double ySup = Parameters->ySup;
	double tolGradU = Parameters->tolGradU;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	CoefFormFuncType *CFF = FemStructs->CFF;
	
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
	
	sixV = 6.0*CFF[e].volume;

	// baricentric coordinate
	xb = (X[0] + X[1] + X[2] + X[3])/4.0;
	yb = (Y[0] + Y[1] + Y[2] + Y[3])/4.0;
	zb = (Z[0] + Z[1] + Z[2] + Z[3])/4.0;
	
	uB = (Ue[0] + Ue[1] + Ue[2] + Ue[3])/4.0;
	
	for(i = 0; i < 4; i++){
		b[i] = CFF[e].b[i];
		c[i] = CFF[e].c[i];
		d[i] = CFF[e].d[i];
	}
	
	// diffusion
	FemFunctions->kappa(&kappaX, &kappaY, &kappaZ, xb, yb, zb, Cst);
	normKappa = Cst; //sqrt(kappaX*kappaX + kappaY*kappaY + kappaZ*kappaZ);
		
	// convection
	FemFunctions->beta(&betaX, &betaY, &betaZ, xb, yb, zb, Cst, Vmax, xSup, ySup);
	normBeta = sqrt(betaX*betaX + betaY*betaY + betaZ*betaZ);
		
	// reaction
	sigma = FemFunctions->sigma(xb, yb, zb);
	
	ConstGradU = 1.0/sixV;
		
	gradU[0] = ConstGradU*(Ue[0]*b[0] + Ue[1]*b[1] + Ue[2]*b[2] + Ue[3]*b[3]);
	gradU[1] = ConstGradU*(Ue[0]*c[0] + Ue[1]*c[1] + Ue[2]*c[2] + Ue[3]*c[3]);
	gradU[2] = ConstGradU*(Ue[0]*d[0] + Ue[1]*d[1] + Ue[2]*d[2] + Ue[3]*d[3]);
		
	dotGradU = ddot(3, gradU, gradU);
	normGradU = sqrt(dotGradU);
//printf("\n MACRO normGradU = %e\n", normGradU);getchar();	
	Fbari = FemFunctions->f(xb, yb, zb, Cst);
	
	res = betaX*gradU[0] + betaY*gradU[1] + betaZ*gradU[2] + sigma*uB - Fbari;
	modR = fabs(res);
	
	if (strcasecmp(Parameters->TauMacro,"PROP1")==0){
		
		double w = Parameters->wMacro;
		double h, Cbtil, Cb;
		double *CbOld = FemStructs->CbOldMacro;
		double Pe, FuncPe;
		
		// Paramêtro de malha
		h = LengthMesh(e, Parameters, Element, betaX, betaY, betaZ, CFF[e].volume, X, Y, Z, modR, gradU);
		
		Pe = (normBeta*h)/(2.0*normKappa);
		FuncPe = max(0, 1 - 1/Pe);
		 
		if(normGradU < tolGradU){
			Cbtil = 0.0;
		}else{
			Cbtil = (h*FuncPe*modR)/(2.0*normGradU);
		}
		
		Cb = w*Cbtil + (1-w)*CbOld[e];
		CbOld[e] = Cbtil;
	
		tau = Cb;
		
	}else if (strcasecmp(Parameters->TauMacro,"PROP2")==0){
		
		double w = Parameters->wMacro;
		double h, Cbtil, Cb;
		double *CbOld = FemStructs->CbOldMacro;
		double Pe, FuncPe;
		
		// Paramêtro de malha
		h = LengthMesh(e, Parameters, Element, betaX, betaY, betaZ, CFF[e].volume, X, Y, Z, modR, gradU);
		
		Pe = (modR*h)/(normGradU*2.0*normKappa); // modificação da normBeta = modR/normGradU
		FuncPe = max(0, 1 - 1/Pe);
		 
		if(normGradU < tolGradU){
			Cbtil = 0.0;
		}else{
			Cbtil = (h*FuncPe*modR)/(2.0*normGradU);
		}
		
		Cb = w*Cbtil + (1-w)*CbOld[e];
		CbOld[e] = Cbtil;
	
		tau = Cb;
		
	}else if (strcasecmp(Parameters->TauMacro,"CAU")==0) {

		double dXidx, dEtadx, dZetadx, dXidy, dEtady, dZetady, dXidz, dEtadz, dZetadz;
		double Be[3], normBe, h, peclet, xi, TOL, NUM;
		double betaH[3], hc, pecletc, xic, ErrBeta[3], normErrBeta;
		double Bec[3], normBec;
	
		dXidx = b[1]/sixV;
		dXidy = c[1]/sixV;
		dXidz = d[1]/sixV;
		dEtadx = b[2]/sixV;
		dEtady = c[2]/sixV;
		dEtadz = d[2]/sixV;
		dZetadx = b[3]/sixV;
		dZetady = c[3]/sixV;
		dZetadz = d[3]/sixV;
	
		// Be_i = SUM (dXi_i/dx_j beta_j)
		Be[0] = dXidx*betaX + dXidy*betaY + dXidz*betaZ;
		Be[1] = dEtadx*betaX + dEtady*betaY + dEtadz*betaZ;
		Be[2] = dZetadx*betaX + dZetady*betaY + dZetadz*betaZ;
			
		normBe = sqrt(Be[0]*Be[0] + Be[1]*Be[1] + Be[2]*Be[2]);
		normBeta = sqrt(betaX*betaX + betaY*betaY + betaZ*betaZ);
			
		h = (2*normBeta)/normBe;
		
		peclet = (normBeta*h)/(2.0*normKappa); 
		
		xi = fmax(0, 1-(1/peclet)); 
	
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
			
		pecletc = (hc*normErrBeta)/(2.0*normKappa);
		
		xic = fmax(0, 1-(1/pecletc));   //cosh(pecletc)/sinh(pecletc) - 1.0/pecletc;
		
		TOL = (xic*hc)/(xi*h);
		NUM = modR/(normBeta*normGradU);
		
		if(NUM >= TOL){
			tau = 0.0;
		}else{
			tau = (xi*h/2.0)*(TOL-NUM)*(modR/normGradU);
		}
		
	}else if (strcasecmp(Parameters->TauMacro,"NSGS")==0){
		
		double w = Parameters->wMacro;
		double h, Cbtil, Cb;
		double *CbOld = FemStructs->CbOldMacro;
		
		// Paramêtro de malha
		h = LengthMesh(e, Parameters, Element, betaX, betaY, betaZ, CFF[e].volume, X, Y, Z, modR, gradU);
		
		//double peclet = (normBeta*h)/(2.0*normKappa);
		 
		//if((normGradU < 1e-6)||(peclet < 1)){
		if(normGradU < 1e-6){
			Cbtil = 0.0;
		}else{
			Cbtil = modR/(2.0*normGradU);
		}
		
		Cb = w*Cbtil + (1-w)*CbOld[e];
		CbOld[e] = Cbtil;
	
		tau = h*Cb;
		
	}else if (strcasecmp(Parameters->TauMacro,"NSGS2")==0){
		
		double w = Parameters->wMacro;
		double tetha = Parameters->tetha;
		double h, Cbtil, Cb;
		double *CbOld = FemStructs->CbOldMacro;
		double peclet, RDR;
		
		// Paramêtro de malha
		//h = LengthMesh(e, Parameters, Element, betaX, betaY, betaZ, CFF[e].volume, X, Y, Z, modR, gradU);
		h = LengthMesh(e, Parameters, Element, betaX, betaY, betaZ, CFF[e].volume, X, Y, Z, res, gradU);
//printf("hMacro = %lf\n", h);		
		
		if (strcasecmp(Parameters->UsePeclet,"YES")==0){
			peclet = (normBeta*h)/(2.0*normKappa);
			RDR = (sigma*h*h)/normKappa;
			
			/*if(peclet > 1.0){
				if (normGradU < tolGradU){
					Cbtil = 0.0;
				}else{
					Cbtil = 0.5*( modR/(normGradU + tetha) );
				}
			}*/
			
			if( (normGradU*(peclet - 1) > tolGradU) || (normGradU*(RDR - 1) > tolGradU) ){
				Cbtil = (h/2.0)*( modR/(normGradU + tetha) );
			}else{
				Cbtil = 0.0;
			}
			
		}else{
			if(normGradU < tolGradU){
				Cbtil = 0.0;
			}else{
				Cbtil = (h/2.0)*( modR/(normGradU + tetha) );
			}
		}
		
		Cb = w*Cbtil + (1-w)*CbOld[e];
		CbOld[e] = Cb;
	
		tau = Cb;
		
	}else if (strcasecmp(Parameters->TauMacro,"NEW")==0){
		
		double w = Parameters->wMacro;
		double h, Cbtil, Cb;
		double *CbOld = FemStructs->CbOldMacro;
		
		// Paramêtro de malha
		h = LengthMesh(e, Parameters, Element, betaX, betaY, betaZ, CFF[e].volume, X, Y, Z, modR, gradU);
		
		if(normGradU < 1e-6){
			Cbtil = 0.0;
		}else{
			Cbtil = pow(modR/normGradU,2.0);
		}
		Cb = w*Cbtil + (1-w)*CbOld[e];
		CbOld[e] = Cbtil;
		
		tau = 0.5*h*Cb;
		
	}else{
		printf("Tau Macro is not defined correctly!\n");
		exit(1);
	}
	
	return tau;
}
