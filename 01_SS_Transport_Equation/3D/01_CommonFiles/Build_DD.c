#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Build_DD(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int e, i;
	int nel, neq;
	double X[4], Y[4], Z[4], a[4], b[4], c[4], d[4], xb, yb, zb;
	double kappaX, kappaY, kappaZ, betaX, betaY, betaZ, sigma;
	double Volume, MaxVolume, Ke[4][4], Fe[4];
	double *uSpace, Ue[4], uB;
	double *F = FemStructs->F;
	double Cst = Parameters->ConstApli;
	double Vmax = Parameters->Vmax;
	double xSup = Parameters->xSup;
	double ySup = Parameters->ySup;
	double Tmacro=0.0, Tmicro=0.0;
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

	dzero(neq+1, F); // zerando termo de fonte
	setZeros(Parameters,MatrixData); // zerando matriz elementar

	uSpace = (double*) mycalloc("uSpace of 'Build_DD'", Parameters->nnodes, sizeof(double));
	eval_U_Space(Parameters, FemStructs, FemFunctions, uSpace);
	
	for(e = 0; e < nel; e++){
		//printf("Elemento = %d\n", e);
		Volume = Element_Configuration(e, Element, Node, X, Y, Z, a, b, c, d, Ue, uSpace, &uB);		
		CFF[e].volume = Volume;
		if(Volume > MaxVolume)
			MaxVolume = Volume;
		
		// baricentric coordinate
		xb = (X[0] + X[1] + X[2] + X[3])/4.0;
		yb = (Y[0] + Y[1] + Y[2] + Y[3])/4.0;
		zb = (Z[0] + Z[1] + Z[2] + Z[3])/4.0;
		
		for(i = 0; i < 4; i++){
			CFF[e].a[i] = a[i];
			CFF[e].b[i] = b[i];
			CFF[e].c[i] = c[i];
			CFF[e].d[i] = d[i];
			//printf("a[%d] = %lf\t b[%d] = %lf\t c[%d] = %lf\t d[%d] = %lf\n", i, a[i], i, b[i], i, c[i], i, d[i]);
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
			// convection
		FemFunctions->beta(&betaX, &betaY, &betaZ, xb, yb, zb, Cst, Vmax, xSup, ySup);
			// reaction
		sigma = FemFunctions->sigma(xb, yb, zb);
	
		// Artificial diffusion
		Tmacro = TauMacro(e, Parameters, FemStructs, FemFunctions, Ue);
		Tmicro = TauMicro(e, Parameters, FemStructs, FemFunctions, Ue);
	//printf("TauMacro = %lf \t TauMicro = %lf\n", Tmacro, Tmicro);
		// Dhh = int grad(v). K grad(u) dOmega
		double Dhh11, Dhh12, Dhh13, Dhh14, Dhh22, Dhh23, Dhh24, Dhh33, Dhh34, Dhh44;
		double ConstDhh = 1.0/(36.0*Volume);
		
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
		
		// Dhhd = int tau_dd grad(v).grad(u) dOmega
		double Dhhd11, Dhhd12, Dhhd13, Dhhd14, Dhhd22, Dhhd23, Dhhd24, Dhhd33, Dhhd34, Dhhd44;
		double ConstDhhd = Tmacro/(36.0*Volume);
		
		Dhhd11 = b1b1 + c1c1 + d1d1;
		Dhhd12 = b1b2 + c1c2 + d1d2;
		Dhhd13 = b1b3 + c1c3 + d1d3;
		Dhhd14 = b1b4 + c1c4 + d1d4;
		Dhhd22 = b2b2 + c2c2 + d2d2;
		Dhhd23 = b2b3 + c2c3 + d2d3;
		Dhhd24 = b2b4 + c2c4 + d2d4;
		Dhhd33 = b3b3 + c3c3 + d3d3;
		Dhhd34 = b3b4 + c3c4 + d3d4;
		Dhhd44 = b4b4 + c4c4 + d4d4;
		
		// Chh = int v beta . grad(u) dOmega
		double Chh1, Chh2, Chh3, Chh4;
		double ConstChh = 1.0/24.0;
		
		Chh1 = betaX*b[0] + betaY*c[0] + betaZ*d[0];
		Chh2 = betaX*b[1] + betaY*c[1] + betaZ*d[1];
		Chh3 = betaX*b[2] + betaY*c[2] + betaZ*d[2];
		Chh4 = betaX*b[3] + betaY*c[3] + betaZ*d[3];
		
		// Rhh = int sigma v u dOmega		
		double ConstRhh = (sigma*Volume)/20.0;
		
		// Ahh = Dhh + Dhhd + Chh + Rhh
		double Ahh[4][4];
		
		Ahh[0][0] = ConstDhh*Dhh11 + ConstDhhd*Dhhd11 + ConstChh*Chh1 + 2*ConstRhh; 
		Ahh[1][0] = ConstDhh*Dhh12 + ConstDhhd*Dhhd12 + ConstChh*Chh1 + ConstRhh; 
		Ahh[2][0] = ConstDhh*Dhh13 + ConstDhhd*Dhhd13 + ConstChh*Chh1 + ConstRhh; 
		Ahh[3][0] = ConstDhh*Dhh14 + ConstDhhd*Dhhd14 + ConstChh*Chh1 + ConstRhh; 
		                                     
		Ahh[0][1] = ConstDhh*Dhh12 + ConstDhhd*Dhhd12 + ConstChh*Chh2 + ConstRhh; 
		Ahh[1][1] = ConstDhh*Dhh22 + ConstDhhd*Dhhd22 + ConstChh*Chh2 + 2*ConstRhh;            
		Ahh[2][1] = ConstDhh*Dhh23 + ConstDhhd*Dhhd23 + ConstChh*Chh2 + ConstRhh;            
		Ahh[3][1] = ConstDhh*Dhh24 + ConstDhhd*Dhhd24 + ConstChh*Chh2 + ConstRhh;            
		                                    
		Ahh[0][2] = ConstDhh*Dhh13 + ConstDhhd*Dhhd13 + ConstChh*Chh3 + ConstRhh; 
		Ahh[1][2] = ConstDhh*Dhh23 + ConstDhhd*Dhhd23 + ConstChh*Chh3 + ConstRhh;              
		Ahh[2][2] = ConstDhh*Dhh33 + ConstDhhd*Dhhd33 + ConstChh*Chh3 + 2*ConstRhh;              
		Ahh[3][2] = ConstDhh*Dhh34 + ConstDhhd*Dhhd34 + ConstChh*Chh3 + ConstRhh;              
		                                    
		Ahh[0][3] = ConstDhh*Dhh14 + ConstDhhd*Dhhd14 + ConstChh*Chh4 + ConstRhh;            
		Ahh[1][3] = ConstDhh*Dhh24 + ConstDhhd*Dhhd24 + ConstChh*Chh4 + ConstRhh;            
		Ahh[2][3] = ConstDhh*Dhh34 + ConstDhhd*Dhhd34 + ConstChh*Chh4 + ConstRhh;            
		Ahh[3][3] = ConstDhh*Dhh44 + ConstDhhd*Dhhd44 + ConstChh*Chh4 + 2*ConstRhh;    

	/*	printf("\n");
		for (i = 0; i < 4; i++){
			for(j = 0; j < 4; j++){
				printf("Ahh[%d][%d] = %lf\t", i,j,Ahh[i][j]);
			}
			printf("\n");
		} */
		
		// Chb = int v beta . grad(u') dOmega
		double ConstChb = 16.0/315.0;
		double Chb1, Chb2, Chb3, Chb4;
		
		Chb1 = betaX*(b[1] + b[2] + b[3]) + betaY*(c[1] + c[2] + c[3]) + betaZ*(d[1] + d[2] + d[3]);
		Chb2 = betaX*(b[0] + b[2] + b[3]) + betaY*(c[0] + c[2] + c[3]) + betaZ*(d[0] + d[2] + d[3]);
		Chb3 = betaX*(b[0] + b[1] + b[3]) + betaY*(c[0] + c[1] + c[3]) + betaZ*(d[0] + d[1] + d[3]);
		Chb4 = betaX*(b[0] + b[1] + b[2]) + betaY*(c[0] + c[1] + c[2]) + betaZ*(d[0] + d[1] + d[2]);
		
		// Rhb = int sigma v u' dOmega
		double ConstRhb = (sigma*Volume*8.0)/105.0;
		
		// Ahb = Chb + Rhb (Dhb = 0)
		double Ahb[4];
		
		Ahb[0] = ConstChb*Chb1 + ConstRhb; 
		Ahb[1] = ConstChb*Chb2 + ConstRhb;
		Ahb[2] = ConstChb*Chb3 + ConstRhb;
		Ahb[3] = ConstChb*Chb4 + ConstRhb;
		
	/*	printf("\n");
		for (i = 0; i < 4; i++){
			printf("Ahb[%d] = %lf\n", i, Ahb[i]);
		} */

		// Cbh = int v' beta . grad(u) dOmega
		double ConstCbh = 16.0/315.0;
		double Cbh1, Cbh2, Cbh3, Cbh4;
		
		Cbh1 = betaX*b[0] + betaY*c[0] + betaZ*d[0];
		Cbh2 = betaX*b[1] + betaY*c[1] + betaZ*d[1];
		Cbh3 = betaX*b[2] + betaY*c[2] + betaZ*d[2];
		Cbh4 = betaX*b[3] + betaY*c[3] + betaZ*d[3];
		
		// Rbh = int sigma v' u dOmega
		double ConstRbh = (sigma*Volume*8.0)/105.0;
		
		// Abh = Cbh + Rbh (Dbh = 0)
		double Abh[4];
		
		Abh[0] = ConstCbh*Cbh1 + ConstRbh; 
		Abh[1] = ConstCbh*Cbh2 + ConstRbh;
		Abh[2] = ConstCbh*Cbh3 + ConstRbh;
		Abh[3] = ConstCbh*Cbh4 + ConstRbh;
		
	/*	printf("\n");
		for (i = 0; i < 4; i++){
			printf("Abh[%d] = %lf\n", i, Abh[i]);
		} */

		// Dbb = int grad(v'). K grad(u') dOmega
		double ConstDbb = 2048.0/(8505.0*Volume);
		double D1, D2, D3;

		//D1 = kappaX*(b1b1 + b2b2 + b3b3 + b4b4 + b1b2 + b1b3 + b1b4 + b2b3 + b2b4 + b3b4); 
		//D2 = kappaY*(c1c1 + c2c2 + c3c3 + c4c4 + c1c2 + c1c3 + c1c4 + c2c3 + c2c4 + c3c4);
		//D3 = kappaZ*(d1d1 + d2d2 + d3d3 + d4d4 + d1d2 + d1d3 + d1d4 + d2d3 + d2d4 + d3d4);

		D1 = kappaX*(b2b2 + b3b3 + b4b4 + b2b3 + b2b4 + b3b4); 
		D2 = kappaY*(c2c2 + c3c3 + c4c4 + c2c3 + c2c4 + c3c4);
		D3 = kappaZ*(d2d2 + d3d3 + d4d4 + d2d3 + d2d4 + d3d4);
		
		// Dbb_dd = int tau_dd grad(v').grad(u') dOmega
		double ConstDbbd = (Tmicro*2048.0)/(8505.0*Volume);
		double Dd1, Dd2, Dd3;
		
		//Dd1 = (b1b1 + b2b2 + b3b3 + b4b4 + b1b2 + b1b3 + b1b4 + b2b3 + b2b4 + b3b4); 
		//Dd2 = (c1c1 + c2c2 + c3c3 + c4c4 + c1c2 + c1c3 + c1c4 + c2c3 + c2c4 + c3c4);
		//Dd3 = (d1d1 + d2d2 + d3d3 + d4d4 + d1d2 + d1d3 + d1d4 + d2d3 + d2d4 + d3d4);

		Dd1 = (b2b2 + b3b3 + b4b4 + b2b3 + b2b4 + b3b4); 
		Dd2 = (c2c2 + c3c3 + c4c4 + c2c3 + c2c4 + c3c4);
		Dd3 = (d2d2 + d3d3 + d4d4 + d2d3 + d2d4 + d3d4);

				//printf("Dentro da Build\n");
		//printf("b1+b2+b3+b4 = %e\n", b[1]+b[2]+b[3]+b[0]);
		//printf("c1+c2+c3+c4 = %e\n", c[1]+c[2]+c[3]+c[0]);
		//printf("d1+d2+d3+d4 = %e\n", d[1]+d[2]+d[3]+d[0]);

		//printf("b1b1+b1b2+b1b3+b1b4 = %e\n", b1b1+b1b2+b1b3+b1b4);
		//printf("c1c1+c1c2+c1c3+c1c4 = %e\n", c1c1+c1c2+c1c3+c1c4);
		//printf("d1d1+d1d2+d1d3+d1d4 = %e\n", d1d1+d1d2+d1d3+d1d4);

		//getchar();


		// Rbb = int sigma v' u' dOmega
		double ConstRbb = (sigma*Volume*8192.0)/4725.0;
		
		// Abb = Dbb + Dbb_dd + Rbb (Cbb = 0)
		double Abb = ConstDbb*(D1 + D2 + D3) + ConstDbbd*(Dd1 + Dd2 + Dd3) + ConstRbb;
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
		
	/*	printf("\n");
		for (i = 0; i < 4; i++){
			for(j = 0; j < 4; j++){
				printf("Ke[%d][%d] = %lf\t", i,j,Ke[i][j]);
			}
			printf("\n");
		} */

		// Fh = int v f dOmega
		double f[4], Faux;
		double ConstFh = Volume/20.0;
		
		for(i = 0; i < 4; i++){
			f[i] = FemFunctions->f(X[i], Y[i], Z[i], Cst);
		}

		Faux = f[0] + f[1] + f[2] + f[3];
		
		// Fb = int v' f dOmega
		double ConstFb = (8.0*Volume*Faux)/105.0;

		// Fe = Fh - AhbAbb^{-1}Fb
		Fe[0]  = ConstFh*(Faux + f[0])- Ahb[0]*invAbb*ConstFb;
		Fe[1]  = ConstFh*(Faux + f[1])- Ahb[1]*invAbb*ConstFb;
		Fe[2]  = ConstFh*(Faux + f[2])- Ahb[2]*invAbb*ConstFb;
		Fe[3]  = ConstFh*(Faux + f[3])- Ahb[3]*invAbb*ConstFb;			
		
	/*	printf("\n");
		for (i = 0; i < 4; i++){
			printf("Fe[%d] = %lf\n", i, Fe[i]);
		} */

		// Matrix assembly according to chosen storage scheme (EBE, EDE or CSR)
		FemFunctions->MatrixAssembly(Parameters, MatrixData, FemStructs, e, Ke);
		
		// Font Assembly 
		font_assembly(e, Fe, Ke, FemFunctions, FemStructs, neq, Cst, xSup);

	/*	printf("\nVetor F\n");
		for(i = 0; i <= neq; i++){
			printf("%d: %E\n", i, FemStructs->F[i]);
			if(i % 20 == 0){
				getchar();
			}
		}

		printf("\nMatriz Global\n");
		for(i = 0; i < neq; i++){
			if(MatrixData->IA[i+1]-MatrixData->IA[i] != 0){
				for (j = MatrixData->IA[i]; j < MatrixData->IA[i+1]; j++){
					printf("M[%d][%d] = %E\n", i+1, MatrixData->JA[j]+1, MatrixData->AA[j]);
					if(MatrixData->JA[j] % 20 == 0){
						getchar();
					}
				}
			}
		} */

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

	// Keeping de maximum volume
	Parameters->MaxVolume = MaxVolume;
	
	// Freeinf memory
	myfree(uSpace);
	myfree(PeL);
	
	return 0;
		
}
