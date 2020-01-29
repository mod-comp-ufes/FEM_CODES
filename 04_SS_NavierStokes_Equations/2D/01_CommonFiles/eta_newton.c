#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include  <math.h>

double eta_newton(double *Fold, double *F, double etaold, int pricek, int it, double talnl, double eta0, ParametersType *Parameters)
{
	double alpha, gamma, beta, etaa, etaaux, eta=0;
	double delta, delta1, delta2, old, ak, bk, thetak;
	int caso, neq;	

	caso = Parameters->SolverToleranceCase; 
	neq = Parameters->neq;
	delta1 = sqrt(ddot(neq, F, F));
	delta2 = sqrt(ddot(neq, Fold, Fold));
	
	delta = delta1/delta2;
	
	switch(caso){
		//case 0:
		//	eta = eta0;
		//	break;
		case 1:			// etapp
			eta = (eta0 < delta)? eta0: delta*delta;
			break;
		case 2:			// etaewk
			alpha = 2.;
			gamma = 0.9;
			beta = 0.1;
			etaa = gamma*pow(delta,alpha);
			old = gamma*pow(etaold,alpha);
			if(old < beta){
				eta = (eta0 < etaa)? eta0: etaa;	
			}else{
				etaaux = (etaa > old)? etaa: old;
				eta = (eta0 < etaaux)? eta0: etaaux;
			}
			break;			 
		case 3:			// etaewc
			alpha = (1.+sqrt(5))*0.5;
			gamma = 1.0;
			beta = 0.;
			etaa = gamma*pow(delta,alpha);			 
			old = gamma*pow(etaold,alpha);
			if(old < beta){
				eta = (eta0 < etaa)? eta0: etaa;	
			}else{
				etaaux = (etaa > old)? etaa: old;
				eta = (eta0 < etaaux)? eta0: etaaux;
			}break;
		case 4:			// etaglt
			ak = log10(delta1) - log10(delta2);
			bk = log10(pricek);
			thetak = bk/sqrt(ak*ak + bk*bk);
			gamma = pow(1./(it+1),1.1)*cos(thetak)*cos(thetak);
					
			etaa = gamma*delta;			 
			eta = (eta0 < etaa)? eta0: etaa;	
			break;
	}
	old = 0.5*talnl/delta1;
	etaaux = (eta > old)? eta: old;
	eta = (eta0 < etaaux)? eta0: etaaux; 
	if(caso == 0)
		eta = eta0;
	printf("\n Eta: %3.2E, Iteracao Newton: %d, Pricek: %d \n", eta, it+1, pricek);
				
	return eta;
}

