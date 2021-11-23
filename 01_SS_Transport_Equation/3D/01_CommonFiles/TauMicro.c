#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

double TauMicro(int e, ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *Ue)
{	
	
	int J1, J2, J3, J4, i;
	double X[4], Y[4], Z[4], xb, yb, zb, sixV, b[4], c[4], d[4], uB;
	double tau = 0.0, h, ConstGradU, gradU[3], normGradU, res, modR, Fbari, Cbtil;
	double kappaX, kappaY, kappaZ, normKappa, normBeta;
	double betaX, betaY, betaZ, sigma;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	CoefFormFuncType *CFF = FemStructs->CFF;
	double w = Parameters->wMicro;
	double Cb = 0.0;
	double *CbOld = FemStructs->CbOldMicro;
	double Cst = Parameters->ConstApli;
	double Vmax = Parameters->Vmax;
	double xSup = Parameters->xSup;
	double ySup = Parameters->ySup;
	double peclet, RDR;
	double tolGradU = Parameters->tolGradU;
	
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
		
	normGradU = sqrt(gradU[0]*gradU[0] + gradU[1]*gradU[1] + gradU[2]*gradU[2]);
//printf("\n MICRO normGradU = %e\n", normGradU);getchar();
	Fbari = FemFunctions->f(xb, yb, zb, Cst);
		
	res = betaX*gradU[0] + betaY*gradU[1] + betaZ*gradU[2] + sigma*uB - Fbari;
	modR = fabs(res);
	
	if (strcasecmp(Parameters->TauMicro,"PROP1")==0){
		
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
		
	}else if (strcasecmp(Parameters->TauMicro,"PROP2")==0){
		
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
		
	}else if (strcasecmp(Parameters->TauMicro,"NSGS")==0){
		// Paramêtro de malha
		h = LengthMesh(e, Parameters, Element, betaX, betaY, betaZ, CFF[e].volume, X, Y, Z, modR, gradU);
		
		peclet = (normBeta*h)/(2.0*normKappa); 
	
		//if((normGradU < 1e-6)||(peclet < 1)){
		if(normGradU < 1e-6){
			Cbtil = 0.0;
		}else{
			Cbtil = modR/(2.0*normGradU);
		}
		
		Cb = w*Cbtil + (1-w)*CbOld[e];
		CbOld[e] = Cbtil;
	
		tau = h*Cb;
	
	}else if (strcasecmp(Parameters->TauMicro,"NSGS2")==0){
		
		double tetha = Parameters->tetha;	
		
		// Paramêtro de malha
		//h = LengthMesh(e, Parameters, Element, betaX, betaY, betaZ, CFF[e].volume, X, Y, Z, modR, gradU);
		h = LengthMesh(e, Parameters, Element, betaX, betaY, betaZ, CFF[e].volume, X, Y, Z, res, gradU);
//printf("hMicro = %lf\n", h);		
		
		if (strcasecmp(Parameters->UsePeclet,"YES")==0){
			peclet = (normBeta*h)/(2.0*normKappa);
			RDR = (sigma*h*h)/normKappa;
			
			/*if(peclet > 1.0){
				if (normGradU < tolGradU){
					Cbtil = 0.0;
				}else{
					Cbtil = 0.5*( modR/(normGradU + tetha) );
				}
			} */
			
			if( (normGradU*(peclet - 1) > tolGradU)||(normGradU*(RDR - 1) > tolGradU) ){
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
		
	}else{
		
		printf("Tau Micro is not defined correctly!\n");
		exit(1);
		
	}
	
	
	
	return tau;	
}
