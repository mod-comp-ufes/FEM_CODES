#include "SSTransportEquation3D.h"

double LengthMesh(int e, ParametersType *Parameters, ElementType *Element, double betaX, double betaY, double betaZ, double Volume, double *X, double *Y, double *Z, double modR, double *gradU)
{
	double h, h2, h3 ,h4 ,h5, hold;
	double R3Vol = pow(Volume,1.0/3.0);
	double betaB[3], normBetaB, B2[3], normB2, B3[3], normB3, B4[3], normB4, B5[3], normB5;
	double sixV, x21, x31, x41, y21, y31, y41, z21, z31, z41;
	double normGradU, dotGradU;
		
	sixV = 6.0*Volume;
		
	x21 = X[1] - X[0];
	x31 = X[2] - X[0];
	x41 = X[3] - X[0];
	y21 = Y[1] - Y[0];
	y31 = Y[2] - Y[0];
	y41 = Y[3] - Y[0];
	z21 = Z[1] - Z[0];
	z31 = Z[2] - Z[0];
	z41 = Z[3] - Z[0];
	
	normGradU = sqrt(gradU[0]*gradU[0] + gradU[1]*gradU[1] + gradU[2]*gradU[2]);
	dotGradU = gradU[0]*gradU[0] + gradU[1]*gradU[1] + gradU[2]*gradU[2];
//printf("dotGradU: %lf\n", dotGradU);
	
	betaB[0] = (modR/dotGradU)*gradU[0];
	betaB[1] = (modR/dotGradU)*gradU[1];
	betaB[2] = (modR/dotGradU)*gradU[2];
	
	normBetaB = sqrt(betaB[0]*betaB[0] + betaB[1]*betaB[1] + betaB[2]*betaB[2]);
//printf("normBetaB Lenght: %lf\n", normBetaB);

	// transpose of matrix dXi/dx 		
	double dXidxT[3][3];
	// Coluna 1
	dXidxT[0][0] = ((y31*z41) - (z31*y41))/sixV;   // b2/6V
	dXidxT[1][0] = (- (x31*z41) + (z31*x41))/sixV; // c2/6V
	dXidxT[2][0] = ((x31*y41) - (y31*x41))/sixV;   // d2/6V
	// Coluna 2
	dXidxT[0][1] = (- (y21*z41) + (z21*y41))/sixV; // b3/6V 
	dXidxT[1][1] = ((x21*z41) - (z21*x41))/sixV;   // c3/6V
	dXidxT[2][1] = (- (x21*y41) + (y21*x41))/sixV; // d3/6V
	// Coluna 3
	dXidxT[0][2] = ((y21*z31) - (z21*y31))/sixV;   // b4/6V
	dXidxT[1][2] = (- (x21*z31) + (z21*x31))/sixV; // c4/6V
	dXidxT[2][2] = ((x21*y31) - (x31*y21))/sixV;   // d4/6V
		
	B2[0] = betaB[0]*dXidxT[0][0] + betaB[1]*dXidxT[1][0] + betaB[2]*dXidxT[2][0];
	B2[1] = betaB[0]*dXidxT[0][1] + betaB[1]*dXidxT[1][1] + betaB[2]*dXidxT[2][1];
	B2[2] = betaB[0]*dXidxT[0][2] + betaB[1]*dXidxT[1][2] + betaB[2]*dXidxT[2][2];
		
	normB2 = sqrt(B2[0]*B2[0] + B2[1]*B2[1] + B2[2]*B2[2]);
//printf("normB2 Lenght: %lf\n", normB2);
		
	// matrix dXi/dx 		
	double dXidx[3][3];
	// Linha 1
	dXidx[0][0] = ((y31*z41) - (z31*y41))/sixV;   // b2/6V
	dXidx[0][1] = (- (x31*z41) + (z31*x41))/sixV; // c2/6V
	dXidx[0][2] = ((x31*y41) - (y31*x41))/sixV;   // d2/6V
	// Linha 2
	dXidx[1][0] = (- (y21*z41) + (z21*y41))/sixV; // b3/6V 
	dXidx[1][1] = ((x21*z41) - (z21*x41))/sixV;   // c3/6V
	dXidx[1][2] = (- (x21*y41) + (y21*x41))/sixV; // d3/6V
	// Linha 3
	dXidx[2][0] = ((y21*z31) - (z21*y31))/sixV;   // b4/6V
	dXidx[2][1] = (- (x21*z31) + (z21*x31))/sixV; // c4/6V
	dXidx[2][2] = ((x21*y31) - (x31*y21))/sixV;   // d4/6V
		
	B3[0] = betaB[0]*dXidx[0][0] + betaB[1]*dXidx[1][0] + betaB[2]*dXidx[2][0];
	B3[1] = betaB[0]*dXidx[0][1] + betaB[1]*dXidx[1][1] + betaB[2]*dXidx[2][1];
	B3[2] = betaB[0]*dXidx[0][2] + betaB[1]*dXidx[1][2] + betaB[2]*dXidx[2][2];
			
	normB3 = sqrt(B3[0]*B3[0] + B3[1]*B3[1] + B3[2]*B3[2]);
//printf("normB3 Lenght: %lf\n", normB3);

	// transpose of matrix dx/dXi 		
	double dxdXiT[3][3];
	// Coluna 1
	dxdXiT[0][0] = x21;
	dxdXiT[1][0] = x31;
	dxdXiT[2][0] = x41;
	// Coluna 2
	dxdXiT[0][1] = y21;
	dxdXiT[1][1] = y31;
	dxdXiT[2][1] = y41;
	// Coluna 3
	dxdXiT[0][2] = z21; 
	dxdXiT[1][2] = z31; 
	dxdXiT[2][2] = z41; 
			
	B4[0] = betaB[0]*dxdXiT[0][0] + betaB[1]*dxdXiT[1][0] + betaB[2]*dxdXiT[2][0];
	B4[1] = betaB[0]*dxdXiT[0][1] + betaB[1]*dxdXiT[1][1] + betaB[2]*dxdXiT[2][1];
	B4[2] = betaB[0]*dxdXiT[0][2] + betaB[1]*dxdXiT[1][2] + betaB[2]*dxdXiT[2][2];
			
	normB4 = sqrt(B4[0]*B4[0] + B4[1]*B4[1] + B4[2]*B4[2]);
//printf("normB4 Lenght: %lf\n", normB4);

	// matrix dx/dXi 		
	double dxdXi[3][3];
	// Coluna 1
	dxdXi[0][0] = x21;
	dxdXi[1][0] = y21;
	dxdXi[2][0] = z21;
	// Coluna 2
	dxdXi[0][1] = x31;
	dxdXi[1][1] = y31;
	dxdXi[2][1] = z31;
	// Coluna 3
	dxdXi[0][2] = x41; 
	dxdXi[1][2] = y41; 
	dxdXi[2][2] = z41; 	
			
	B5[0] = betaB[0]*dxdXi[0][0] + betaB[1]*dxdXi[1][0] + betaB[2]*dxdXi[2][0];
	B5[1] = betaB[0]*dxdXi[0][1] + betaB[1]*dxdXi[1][1] + betaB[2]*dxdXi[2][1];
	B5[2] = betaB[0]*dxdXi[0][2] + betaB[1]*dxdXi[1][2] + betaB[2]*dxdXi[2][2];
			
	normB5 = sqrt(B5[0]*B5[0] + B5[1]*B5[1] + B5[2]*B5[2]);
//printf("normB5 Lenght: %lf\n", normB5);
					
	if((normGradU < Parameters->tolGradU)){ //||(modR < 1e-6)){
		h2 = 0.0;
		h3 = 0.0;
		h4 = 0.0;
		h5 = 0.0;
	}else{
		h2 = (2*normBetaB)/normB2;
		h3 = (normBetaB)/normB3;
		h4 = (2*normBetaB)/normB4;
		h5 = (2*normBetaB)/normB5;
	}
	
	if(Parameters->TipoH == 1){
		h = R3Vol;
	}else if(Parameters->TipoH == 2){
		h = h2;
	}else if(Parameters->TipoH == 3){
		h = h3;
		//printf("## res = %e \t Norm GRadU = %e \t h = %lf\n", modR, normGradU, h); getchar();
	}else if(Parameters->TipoH == 4){
		h = h4;
	}else if(Parameters->TipoH == 5){
		h = h5;
	}
	
	double nx0[3], nx1[3], ny0[3], ny1[3], nz0[3], nz1[3];
	double betaN1, betaN2, betaN3, betaN4, betaN5, betaN6;
	
	// Vetores normais
	nx0[0] = -1.0; nx0[1] = 0.0;  nx0[2] = 0.0;  // plano x = 0
	nx1[0] = 1.0;  nx1[1] = 0.0;  nx1[2] = 0.0;  // plano x = 1
	ny0[0] = 0.0;  ny0[1] = -1.0; ny0[2] = 0.0;  // plano y = 0
	ny1[0] = 0.0;  ny1[1] = 1.0;  ny1[2] = 0.0;  // plano y = 1
	nz0[0] = 0.0;  nz0[1] = 0.0;  nz0[2] = -1.0; // plano z = 0
	nz1[0] = 0.0;  nz1[1] = 0.0;  nz1[2] = 1.0;  // plano z = 1
	
	betaN1 = betaX*ny0[0] + betaY*ny0[1] + betaZ*ny0[2];
	betaN2 = betaX*nx1[0] + betaY*nx1[1] + betaZ*nx1[2];
	betaN3 = betaX*ny1[0] + betaY*ny1[1] + betaZ*ny1[2];
	betaN4 = betaX*nx0[0] + betaY*nx0[1] + betaZ*nx0[2];
	betaN5 = betaX*nz1[0] + betaY*nz1[1] + betaZ*nz1[2];
	betaN6 = betaX*nz0[0] + betaY*nz0[1] + betaZ*nz0[2];
	
	
	// Analise se elemento é de saida de fluxo
	if (strncmp(Parameters->OutputFlow,"YES",3)==0){
		
		int OutFlow = atoi(&(Parameters->OutputFlow[3]));
		
		if (OutFlow == 1){ // na saída h recebe 3*(raíz cúbica de volume)
			hold = 3*R3Vol;
		}else if (OutFlow == 2){ // na saída h recebe 3*h3
			hold = 3*h3;
		}else if (OutFlow == 3){ // na saída h recebe a média entre a raíz cúbica e a norma calculada em h3
			hold = 1.5*(R3Vol + h3);
		}
		
	//if(strcasecmp(Parameters->OutputFlow,"YES")==0){
		if(Element[e].Type == 1){
			if(betaN1 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 2){
			if(betaN2 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 3){
			if(betaN3 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 4){
			if(betaN4 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 5){
			if(betaN5 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 6){
			if(betaN6 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 7){
			if(betaN1 > 0.0 || betaN4 > 0.0 || betaN6 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 8){
			if(betaN1 > 0.0 || betaN2 > 0.0 || betaN6 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 9){
			if(betaN3 > 0.0 || betaN2 > 0.0 || betaN6 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 10){
			if(betaN3 > 0.0 || betaN4 > 0.0 || betaN6 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 11){
			if(betaN1 > 0.0 || betaN4 > 0.0 || betaN5 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 12){
			if(betaN1 > 0.0 || betaN2 > 0.0 || betaN5 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 13){
			if(betaN3 > 0.0 || betaN2 > 0.0 || betaN5 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 14){
			if(betaN3 > 0.0 || betaN4 > 0.0 || betaN5 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 15){
			if(betaN1 > 0.0 || betaN5 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 16){
			if(betaN1 > 0.0 || betaN2 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 17){
			if(betaN1 > 0.0 || betaN6 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 18){
			if(betaN1 > 0.0 || betaN4 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 19){
			if(betaN2 > 0.0 || betaN6 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 20){
			if(betaN2 > 0.0 || betaN3 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 21){
			if(betaN2 > 0.0 || betaN5 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 22){
			if(betaN3 > 0.0 || betaN6 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 23){
			if(betaN3 > 0.0 || betaN4 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 24){
			if(betaN3 > 0.0 || betaN5 > 0.0){
				h = hold;
			}
		}else if(Element[e].Type == 25){
			if(betaN4 > 0.0 || betaN5 > 0.0){
				h = hold;
			}
		}					
	}
	
	return h;
}
