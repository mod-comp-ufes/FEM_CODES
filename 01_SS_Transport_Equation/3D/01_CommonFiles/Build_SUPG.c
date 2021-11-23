#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

/* Dada a equação
 * - div(kappa grad(u)) + beta . grad(u) + sigma u = f
 *  
 * Obtemos a formulação
 * int (grad(v) . kappa grad(u)) + int (v beta . grad(u)) + int (sigma v u) = int (v f)
 * 
 * Que dá origem ao sistema
 * Ku = F
 *
 * Ke*Ue = Fe
 * Ke_{4x4} e Fe_{4x1}
 * 
 * Ke: int (grad(v).grad(u))
 * Fe: int (v f)
 *  */

int Build_SUPG(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int e, i;
	int nel, neq;
	double X[4], Y[4], Z[4], a[4], b[4], c[4], d[4], xb, yb, zb;
	double kappaX, kappaY, kappaZ, normK, betaX, betaY, betaZ, normB, sigma;
	double h, peclet, xi, tau;
	double Volume, MaxVolume, Ke[4][4], Fe[4];
	double lixo1[4], lixo2[Parameters->nnodes], lixo3; // Variáveis usadas pra descartar valores na função Element_Configuration
	double *F = FemStructs->F;
	double Cst = Parameters->ConstApli;
	double Vmax = Parameters->Vmax;
	double xSup = Parameters->xSup;
	double ySup = Parameters->ySup;
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

	for(e = 0; e < nel; e++){
		
		Volume = Element_Configuration(e, Element, Node, X, Y, Z, a, b, c, d, lixo1, lixo2, &lixo3);
		CFF[e].volume = Volume;
		h = pow(Volume,1.0/3.0); // mesh parameter
		if(Volume > MaxVolume)
			MaxVolume = Volume;
		
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
		
		// D = int grad(v). K grad(u) dOmega
		double E11, E12, E13, E14, E22, E23, E24, E33, E34, E44;
		double ConstD = 1.0/(36.0*Volume);
		
		FemFunctions->kappa(&kappaX, &kappaY, &kappaZ, xb, yb, zb, Cst);
		normK = sqrt(kappaX*kappaX + kappaY*kappaY + kappaZ*kappaZ);
		
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
		
		FemFunctions->beta(&betaX, &betaY, &betaZ, xb, yb, zb, Cst, Vmax, xSup, ySup);
		normB = sqrt(betaX*betaX + betaY*betaY + betaZ*betaZ);
		
		C1 = betaX*b[0] + betaY*c[0] + betaZ*d[0];
		C2 = betaX*b[1] + betaY*c[1] + betaZ*d[1];
		C3 = betaX*b[2] + betaY*c[2] + betaZ*d[2];
		C4 = betaX*b[3] + betaY*c[3] + betaZ*d[3];
		
		// R = int sigma v u dOmega
		sigma = FemFunctions->sigma(xb, yb, zb);
		
		double ConstR = (sigma*Volume)/20.0;
		
		//Tau SUPG
		peclet = (normB*h)/(2.0*normK);
		xi = fmax(0, 1-1/peclet); //coth(peclet) - 1.0/peclet; 
		tau = (h*xi)/(2.0*normB);
		
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
		
		// Elementary Matrix
		Ke[0][0] = ConstD*E11 + ConstC*C1 + ConstCsupg*Cs11 + ConstRsupg*C1 + 2*ConstR; 
		Ke[1][0] = ConstD*E12 + ConstC*C1 + ConstCsupg*Cs12 + ConstRsupg*C2 + ConstR; 
		Ke[2][0] = ConstD*E13 + ConstC*C1 + ConstCsupg*Cs13 + ConstRsupg*C3 + ConstR; 
		Ke[3][0] = ConstD*E14 + ConstC*C1 + ConstCsupg*Cs14 + ConstRsupg*C4 + ConstR; 
		                                                    
		Ke[0][1] = ConstD*E12 + ConstC*C2 + ConstCsupg*Cs12 + ConstRsupg*C1 + ConstR; 
		Ke[1][1] = ConstD*E22 + ConstC*C2 + ConstCsupg*Cs22 + ConstRsupg*C2 + 2*ConstR;            
		Ke[2][1] = ConstD*E23 + ConstC*C2 + ConstCsupg*Cs23 + ConstRsupg*C3 + ConstR;            
		Ke[3][1] = ConstD*E24 + ConstC*C2 + ConstCsupg*Cs24 + ConstRsupg*C4 + ConstR;            
		                                                    
		Ke[0][2] = ConstD*E13 + ConstC*C3 + ConstCsupg*Cs13 + ConstRsupg*C1 + ConstR; 
		Ke[1][2] = ConstD*E23 + ConstC*C3 + ConstCsupg*Cs23 + ConstRsupg*C2 + ConstR;              
		Ke[2][2] = ConstD*E33 + ConstC*C3 + ConstCsupg*Cs33 + ConstRsupg*C3 + 2*ConstR;              
		Ke[3][2] = ConstD*E34 + ConstC*C3 + ConstCsupg*Cs34 + ConstRsupg*C4 + ConstR;              
		                                                    
		Ke[0][3] = ConstD*E14 + ConstC*C4 + ConstCsupg*Cs14 + ConstRsupg*C1 + ConstR;            
		Ke[1][3] = ConstD*E24 + ConstC*C4 + ConstCsupg*Cs24 + ConstRsupg*C2 + ConstR;            
		Ke[2][3] = ConstD*E34 + ConstC*C4 + ConstCsupg*Cs34 + ConstRsupg*C3 + ConstR;            
		Ke[3][3] = ConstD*E44 + ConstC*C4 + ConstCsupg*Cs44 + ConstRsupg*C4 + 2*ConstR;
		
		
		/*printf("\nMatriz Ke E: %d\n", e);
		for(i = 0; i < 4; i++){
			printf("%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", Ke[i][0], Ke[i][1], Ke[i][2], Ke[i][3]);
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
		
	/*	printf("\nVetor Fe E: %d\n", e);
		for(i = 0; i < 16; i++){
			printf("%lf\n",Fe[i]);
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

	myfree(PeL);
	
	return 0;
		
}
