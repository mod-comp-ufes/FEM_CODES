#include "SSTranspEquation.h"

double h_shock_sqrtArea(ParametersType *Parameters, ElementType *Element, int e, double Be_x, double Be_y, double u1, double u2, double u3, double y23, double y31, double y12, double x32, double x13, double x21, double Area)
{
	double h, nx0[2], nx1[2], ny0[2], ny1[2], betaN1, betaN2, betaN3, betaN4;
	
	//printf("Elemento = %d Tipo = %d\n", e, Element[e].Type); getchar();
	
	h = 2*sqrt(fabs(Area));
	
	// Vetores normais
	nx0[0] = -1.0; nx0[1] = 0.0;  // reta x = 0
	nx1[0] = 1.0;  nx1[1] = 0.0;  // reta x = 1
	ny0[0] = 0.0;  ny0[1] = -1.0; // reta y = 0
	ny1[0] = 0.0;  ny1[1] = 1.0;  // reta y = 1
	
	betaN1 = Be_x*nx0[0] + Be_y*nx0[1];
	betaN2 = Be_x*ny0[0] + Be_y*ny0[1];
	betaN3 = Be_x*nx1[0] + Be_y*nx1[1];
	betaN4 = Be_x*ny1[0] + Be_y*ny1[1];
	
	// Analise se elemento é de saida de fluxo
	if(strcasecmp(Parameters->OutputFlow,"YES")==0){
		if(Element[e].Type == 1){ // reta x = 0
			if(betaN1 > 0.0){
				h = 2*sqrt(fabs(Area));
			}
		}else if(Element[e].Type == 2){ // reta y = 0
			if(betaN2 > 0.0){
				h = 2*sqrt(fabs(Area));
			}
		}else if(Element[e].Type == 3){ // reta x = 1
			if(betaN3 > 0.0){
				h = 2*sqrt(fabs(Area));
			}
		}else if(Element[e].Type == 4){ // reta y = 1
			if(betaN4 > 0.0){
				h = 2*sqrt(fabs(Area));
			}
		}else if(Element[e].Type == 5){ // esquina entre x = 0 e y = 0
			if(betaN1 > 0.0 || betaN2 > 0.0){
				h = 2*sqrt(fabs(Area));
			}
		}else if(Element[e].Type == 6){ // esquina entre y = 0 e x = 1
			if(betaN2 > 0.0 || betaN3 > 0.0){
				h = 2*sqrt(fabs(Area));
			}
		}else if(Element[e].Type == 7){ // esquina entre x = 1 e y = 1
			if(betaN3 > 0.0 || betaN4 > 0.0){
				h = 2*sqrt(fabs(Area));
			}
		}else if(Element[e].Type == 8){ // esquina entre y = 1 e x = 0
			if(betaN4 > 0.0 || betaN1 > 0.0){
				h = 2*sqrt(fabs(Area));
			}
		}		
	}
		
	return h;

}


