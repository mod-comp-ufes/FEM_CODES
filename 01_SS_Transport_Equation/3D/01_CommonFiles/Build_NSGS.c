#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_NSGS(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int e, i;
	int nel, neq;
	double X[4], Y[4], Z[4], a[4], b[4], c[4], d[4], xb, yb, zb;
	double kappaX, kappaY, kappaZ, betaX, betaY, betaZ, sigma;
	double Volume, MaxVolume, Ke[4][4], Fe[4];
	double *uSpace, Ue[4], uB;
	double Cb = 0.0;
	double *CbOld = FemStructs->CbOldMicro;
	double w;
	double *F = FemStructs->F;
	double Cst = Parameters->ConstApli;
	double Vmax = Parameters->Vmax;
	double xSup = Parameters->xSup;
	double ySup = Parameters->ySup;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	CoefFormFuncType *CFF = FemStructs->CFF;
	
	w = Parameters->wMicro;
	nel = Parameters->nel;
	neq = Parameters->neq;
	MaxVolume = 0.0;
	
	dzero(neq+1, F);
	setZeros(Parameters,MatrixData);
	
	uSpace = (double*) mycalloc("uSpace of 'Build_NSGS'", Parameters->nnodes, sizeof(double));
	eval_U_Space(Parameters, FemStructs, FemFunctions, uSpace);
	
	for(e = 0; e < nel; e++){
		
		Volume = Element_Configuration(e, Element, Node, X, Y, Z, a, b, c, d, Ue, uSpace, &uB);		
		CFF[e].volume = Volume;
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
		
		// Dhh = int grad(v). K grad(u) dOmega
		double Dhh11, Dhh12, Dhh13, Dhh14, Dhh22, Dhh23, Dhh24, Dhh33, Dhh34, Dhh44;
		double ConstDhh = 1.0/(36.0*Volume);
		
		FemFunctions->kappa(&kappaX, &kappaY, &kappaZ, xb, yb, zb, Cst);
		
		Dhh11 = kappaX*b1b1 + kappaY*c1c1 + kappaZ*d1d1;
		Dhh12 = kappaX*b1b2 + kappaY*c1c2 + kappaZ*d1d2;
		Dhh13 = kappaX*b1b3 + kappaY*c1c3 + kappaZ*d1d3;
		Dhh14 = kappaX*b1b4 + kappaY*c1c4 + kappaZ*d1d4;
		Dhh22 = kappaX*b2b2 + kappaY*c2c2 + kappaZ*d2d2;
		Dhh23 = kappaX*b2b3 + kappaY*c2c3 + kappaZ*d2d3;
		Dhh24 = kappaX*b2b4 + kappaY*c2c4 + kappaZ*d2d4;
		Dhh33 = kappaX*b3b3 + kappaY*c3c3 + kappaZ*d3d3;
		Dhh34 = kappaX*b3b4 + kappaY*c3c4 + kappaZ*d3d4;
		Dhh44 = kappaX*b4b4 + kappaY*c4c4 + kappaZ*d4d4;
		
		// Chh = int v beta . grad(u) dOmega
		double Chh1, Chh2, Chh3, Chh4;
		double ConstChh = 1.0/24.0;
		
		FemFunctions->beta(&betaX, &betaY, &betaZ, xb, yb, zb, Cst, Vmax, xSup, ySup);
		
		Chh1 = betaX*b[0] + betaY*c[0] + betaZ*d[0];
		Chh2 = betaX*b[1] + betaY*c[1] + betaZ*d[1];
		Chh3 = betaX*b[2] + betaY*c[2] + betaZ*d[2];
		Chh4 = betaX*b[3] + betaY*c[3] + betaZ*d[3];
		
		// Rhh = int sigma v u dOmega
		sigma = FemFunctions->sigma(xb, yb, zb);
		
		double ConstRhh = (sigma*Volume)/20.0;
		
		// Ahh = Dhh + Chh + Rhh
		double Ahh[4][4];
		
		Ahh[0][0] = ConstDhh*Dhh11 + ConstChh*Chh1 + 2*ConstRhh; 
		Ahh[1][0] = ConstDhh*Dhh12 + ConstChh*Chh1 + ConstRhh; 
		Ahh[2][0] = ConstDhh*Dhh13 + ConstChh*Chh1 + ConstRhh; 
		Ahh[3][0] = ConstDhh*Dhh14 + ConstChh*Chh1 + ConstRhh; 
		
		Ahh[0][1] = ConstDhh*Dhh12 + ConstChh*Chh2 + ConstRhh; 
		Ahh[1][1] = ConstDhh*Dhh22 + ConstChh*Chh2 + 2*ConstRhh;            
		Ahh[2][1] = ConstDhh*Dhh23 + ConstChh*Chh2 + ConstRhh;            
		Ahh[3][1] = ConstDhh*Dhh24 + ConstChh*Chh2 + ConstRhh;            
		
		Ahh[0][2] = ConstDhh*Dhh13 + ConstChh*Chh3 + ConstRhh; 
		Ahh[1][2] = ConstDhh*Dhh23 + ConstChh*Chh3 + ConstRhh;              
		Ahh[2][2] = ConstDhh*Dhh33 + ConstChh*Chh3 + 2*ConstRhh;              
		Ahh[3][2] = ConstDhh*Dhh34 + ConstChh*Chh3 + ConstRhh;              
		          
		Ahh[0][3] = ConstDhh*Dhh14 + ConstChh*Chh4 + ConstRhh;            
		Ahh[1][3] = ConstDhh*Dhh24 + ConstChh*Chh4 + ConstRhh;            
		Ahh[2][3] = ConstDhh*Dhh34 + ConstChh*Chh4 + ConstRhh;            
		Ahh[3][3] = ConstDhh*Dhh44 + ConstChh*Chh4 + 2*ConstRhh;    
		
		// Chb = int v beta . grad(u') dOmega
		double ConstChb = 9.0/560;
		double Chb1, Chb2, Chb3, Chb4;
		
		Chb1 = betaX*(b[1] + b[2] + b[3]) + betaY*(c[1] + c[2] + c[3]) + betaZ*(d[1] + d[2] + d[3]);
		Chb2 = betaX*(b[0] + b[2] + b[3]) + betaY*(c[0] + c[2] + c[3]) + betaZ*(d[0] + d[2] + d[3]);
		Chb3 = betaX*(b[0] + b[1] + b[3]) + betaY*(c[0] + c[1] + c[3]) + betaZ*(d[0] + d[1] + d[3]);
		Chb4 = betaX*(b[0] + b[1] + b[2]) + betaY*(c[0] + c[1] + c[2]) + betaZ*(d[0] + d[1] + d[2]);
		
		// Rhb = int sigma v u' dOmega
		double ConstRhb = (sigma*Volume*27.0)/1120.0;
		
		// Ahb = Chb + Rhb (Dhb = 0)
		double Ahb[4];
		
		Ahb[0] = ConstChb*Chb1 + ConstRhb; 
		Ahb[1] = ConstChb*Chb2 + ConstRhb;
		Ahb[2] = ConstChb*Chb3 + ConstRhb;
		Ahb[3] = ConstChb*Chb4 + ConstRhb;
		
		// Cbh = int v' beta . grad(u) dOmega
		double ConstCbh = 9.0/560;
		double Cbh1, Cbh2, Cbh3, Cbh4;
		
		Cbh1 = betaX*b[0] + betaY*c[0] + betaZ*d[0];
		Cbh2 = betaX*b[1] + betaY*c[1] + betaZ*d[1];
		Cbh3 = betaX*b[2] + betaY*c[2] + betaZ*d[2];
		Cbh4 = betaX*b[3] + betaY*c[3] + betaZ*d[3];
		
		// Rbh = int sigma v' u dOmega
		double ConstRbh = (sigma*Volume*27.0)/1120.0;
		
		// Abh = Cbh + Rbh (Dbh = 0)
		double Abh[4];
		
		Abh[0] = ConstCbh*Cbh1 + ConstRbh; 
		Abh[1] = ConstCbh*Cbh2 + ConstRbh;
		Abh[2] = ConstCbh*Cbh3 + ConstRbh;
		Abh[3] = ConstCbh*Cbh4 + ConstRbh;
		
		// Rbb = int sigma v' u' dOmega
		double ConstRbb = (sigma*Volume*243.0)/15400.0;
		
		// tau_nsgs
		double Tnsgs, gradU[3], normGradU, res, Fbari, modR, modB;
		double h = pow(Volume,1.0/3.0);
		double ConstGradU = 1.0/(6*Volume);
		
		gradU[0] = ConstGradU*(Ue[0]*b[0] + Ue[1]*b[1] + Ue[2]*b[2] + Ue[3]*b[3]);
		gradU[1] = ConstGradU*(Ue[0]*c[0] + Ue[1]*c[1] + Ue[2]*c[2] + Ue[3]*c[3]);
		gradU[2] = ConstGradU*(Ue[0]*d[0] + Ue[1]*d[1] + Ue[2]*d[2] + Ue[3]*d[3]);
		
		normGradU = sqrt(gradU[0]*gradU[0] + gradU[1]*gradU[1] + gradU[2]*gradU[2]);
		Fbari = FemFunctions->f(xb, yb, zb, Cst);
		
		res = betaX*gradU[0] + betaY*gradU[1] + betaZ*gradU[2] + sigma*uB - Fbari;
		modR = fabs(res);
		
		if(normGradU != 0.0){
			modB = modR/normGradU;
		}else{
			modB = 0.0;
		}
		
		Cb = w*(0.5*modB) + (1-w)*CbOld[e];
		CbOld[e] = Cb;
		Tnsgs = h*Cb;
		
		// Dbb = int grad(v'). K grad(u') dOmega
		double ConstDbb = 27.0/(1120.0*Volume);
		double D1, D2, D3;
		
		D1 = kappaX*(b1b1 + b2b2 + b3b3 + b4b4 + b1b2 + b1b3 + b1b4 + b2b3 + b2b4 + b3b4); 
		D2 = kappaY*(c1c1 + c2c2 + c3c3 + c4c4 + c1c2 + c1c3 + c1c4 + c2c3 + c2c4 + c3c4);
		D3 = kappaZ*(d1d1 + d2d2 + d3d3 + d4d4 + d1d2 + d1d3 + d1d4 + d2d3 + d2d4 + d3d4);
		
		// Dbb_nsgs = int tau_nsgs grad(v').grad(u') dOmega
		double ConstDbbn = (Tnsgs*27.0)/(1120.0*Volume);
		double Dn1, Dn2, Dn3;
		
		Dn1 = (b1b1 + b2b2 + b3b3 + b4b4 + b1b2 + b1b3 + b1b4 + b2b3 + b2b4 + b3b4); 
		Dn2 = (c1c1 + c2c2 + c3c3 + c4c4 + c1c2 + c1c3 + c1c4 + c2c3 + c2c4 + c3c4);
		Dn3 = (d1d1 + d2d2 + d3d3 + d4d4 + d1d2 + d1d3 + d1d4 + d2d3 + d2d4 + d3d4);
		
		// Abb = Dbb + Dbb_nsgs + Rbb (Cbb = 0)
		double Abb = ConstDbb*(D1 + D2 + D3) + ConstDbbn*(Dn1 + Dn2 + Dn3) + ConstRbb;
		double invAbb = 1.0/Abb;
		
		// Elementary Matrix (Ahh - AhbAbb^{-1}Abh)
		Ke[0][0] = Ahh[0][0] - Ahb[0]*invAbb*Abh[0];
		Ke[1][0] = Ahh[1][0] - Ahb[1]*invAbb*Abh[0];
		Ke[2][0] = Ahh[2][0] - Ahb[2]*invAbb*Abh[0];
		Ke[3][0] = Ahh[3][0] - Ahb[3]*invAbb*Abh[0];
		                              
		Ke[0][1] = Ahh[0][1] - Ahb[0]*invAbb*Abh[1];
		Ke[1][1] = Ahh[1][1] - Ahb[1]*invAbb*Abh[1];            
		Ke[2][1] = Ahh[2][1] - Ahb[2]*invAbb*Abh[1];          
		Ke[3][1] = Ahh[3][1] - Ahb[3]*invAbb*Abh[1];          
		                              
		Ke[0][2] = Ahh[0][2] - Ahb[0]*invAbb*Abh[2];
		Ke[1][2] = Ahh[1][2] - Ahb[1]*invAbb*Abh[2];            
		Ke[2][2] = Ahh[2][2] - Ahb[2]*invAbb*Abh[2];              
		Ke[3][2] = Ahh[3][2] - Ahb[3]*invAbb*Abh[2];            
		                              
		Ke[0][3] = Ahh[0][3] - Ahb[0]*invAbb*Abh[3];          
		Ke[1][3] = Ahh[1][3] - Ahb[1]*invAbb*Abh[3];          
		Ke[2][3] = Ahh[2][3] - Ahb[2]*invAbb*Abh[3];          
		Ke[3][3] = Ahh[3][3] - Ahb[3]*invAbb*Abh[3];            
		
	/*	printf("\nMatriz Ke E: %d\n", e);
		for(i = 0; i < 4; i++){
			printf("%E\t%E\t%E\t%E\n", Ke[i][0], Ke[i][1], Ke[i][2], Ke[i][3]);
		}
		getchar();*/
		
		// Fh = int v f dOmega
		double f[4], Faux;
		double ConstFh = Volume/20.0;
		
		for(i = 0; i < 4; i++){
			f[i] = FemFunctions->f(X[i], Y[i], Z[i], Cst);
		}

		Faux = f[0] + f[1] + f[2] + f[3];
		
		// Fb = int v' f dOmega
		double ConstFb = (27.0*Volume*Faux)/1120.0;

		// Fe = Fh - AhbAbb^{-1}Fb
		Fe[0]  = ConstFh*(Faux + f[0])- Ahb[0]*invAbb*ConstFb;
		Fe[1]  = ConstFh*(Faux + f[1])- Ahb[1]*invAbb*ConstFb;
		Fe[2]  = ConstFh*(Faux + f[2])- Ahb[2]*invAbb*ConstFb;
		Fe[3]  = ConstFh*(Faux + f[3])- Ahb[3]*invAbb*ConstFb;			
		
	/*	printf("\nVetor Fe E: %d\n", e);
		for(i = 0; i < 4; i++){
			printf("%E\n",Fe[i]);
		}
		getchar();*/
		
		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->MatrixAssembly(Parameters, MatrixData, FemStructs, e, Ke);

//		AssemblyGlobalMatrix(MatrixData, FemStructs, e, Ke, neq);
		
		// Font Assembly 
		font_assembly(e, Fe, Ke, FemFunctions, FemStructs, neq, Cst, xSup);
		
	}//end for element
	
	Parameters->MaxVolume = MaxVolume;
	
	myfree(uSpace);
	
	return 0;
		
}
