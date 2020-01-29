#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int i, it, itmax, neq, nel;
	int iter, iterold, pricek;	
	double norm_F, norm_Fold, delta, delta2, epsilon, tol, norms, normu, etaold, eta0;	
	double *s, *Fold, *u, *F, *delta_old, *normres_old;//, *deltam_old;
	

	char FileName[300], FileName2[300];
	FILE *OutFile, *OutFile2;
	
	setProblem(Parameters, FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters,FemOtherFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions);

	u = FemStructs->u;
	F = FemStructs->F;


	//==== print residuo ======	
	sprintf(FileName,"../03_output/Residuo_%s_%s_%s_%s_%s_N%d_E%d.dat", Parameters->ProblemTitle, Parameters->Experiments, Parameters->StabilizationForm, Parameters->MatrixVectorProductScheme, Parameters->Preconditioner,Parameters->nnodes, Parameters->nel);
	OutFile = myfopen(FileName,"w");

	//=====print correcao do fator peso da nova solucao
	sprintf(FileName2,"../03_output/FatorCorr_%s_%s_%s_%s_%s_N%d_E%d.dat", Parameters->ProblemTitle, Parameters->Experiments, Parameters->StabilizationForm, Parameters->MatrixVectorProductScheme, Parameters->Preconditioner,Parameters->nnodes, Parameters->nel);
	OutFile2 = myfopen(FileName2,"w");

	//========= Inicio Iteração de Newton===========
	neq = Parameters->neq;
	nel = Parameters->nel;
	
	s = (double*) mycalloc("s of 'Process'", neq+1, sizeof(double));
	Fold = (double*) mycalloc("Fold of 'Process'", neq+1, sizeof(double));
	delta_old = (double*) mycalloc("delta_old of 'Process'", nel, sizeof(double));
	normres_old = (double*) mycalloc("normres_old of 'Process'", nel, sizeof(double));
	//deltam_old = (double*) mycalloc("delta_old of 'Process'", nel, sizeof(double));
	
	FemStructs->delta_old = delta_old;
	FemStructs->normres_old = normres_old;
	//FemStructs->deltam_old = deltam_old;

	for(i=0; i<neq+1; i++){
		u[i] = 0.0;
		s[i] = 0.0;
		Fold[i] = 0.0;
	}
	for(i=0; i<nel; i++){
		delta_old[i] = 0.0;
		normres_old[i] = 0.0;
	}
	

	//====== Constroi matriz e vetor força (Residuo) ======	
	FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
	norm_F = sqrt(ddot(neq, F, F));	
	printf("\n Norma de F_0. |F_0| = %3.2E \n", norm_F);	

	//====== Precondiona sistema ======
	it = 1;
	FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, it, F);
//	printf("\n AQUI \n");
	norm_F = sqrt(ddot(neq, F, F));
	delta2 = 1;
	delta = 1;	
	tol = Parameters->SolverToleranceNonLin;	
	epsilon = tol*norm_F;	
	printf("\n Norma de F_0. |F_0| = %3.2E === Saida de |F| = %3.2E  \n", norm_F, epsilon);
	itmax = 2000;
	iter = 0;

	//== Parametros para globalizacao ======
	double w_min=0.05, w_max=1.00, c1=1.001, c2=1.1, c3=1.001, c4=0.9, c5=0.5; 
	// teste double w_min=0.05, w_max=1.00, c1=1.001, c2=1.1, c3=1.01, c4=0.9, c5=0.9; 
	double w = w_max, w_old;
	int fd, cont, qtcorr, cont2;

	//================Parameters->SolverTolerance = eta_0 (0.1, 0.5, 0.9) 
	eta0 = Parameters->SolverTolerance;
	pricek = 1;
	dcopy(neq, F, Fold);   //Fould = F
//	int I, J;
	//while((delta1 > epsilon || delta2 > tol) && it < itmax){
	//while(delta2 > tol && it < itmax){
	while(norm_F > epsilon && it < itmax){
		it++;		
		iterold = iter;		
		norm_Fold = sqrt(ddot(neq, Fold, Fold));
		//======== decaimento do residuo=======
		fprintf(OutFile,"\t   %d   %.12lf\n", it-1, log10(norm_F));

		//========================================	
//		FILE *Out;

/*		Out = fopen("octave_solution.m","w");
		fprintf(Out,"A=sparse(%d,%d);\n",neq,neq);
		for (I=0;I<neq;I++)
			for (J=MatrixData->IA[I]; J<MatrixData->IA[I+1]; J++)
				fprintf(Out,"A(%d,%d)=%.15lf;\n",I+1,MatrixData->JA[J]+1,MatrixData->AA[J]);
		for (I=0;I<neq;I++)
			fprintf(Out,"F(%d)=%.15lf;\n",I+1,F[I]);
*/				
		//====== Resolve sitema linear ======		
		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, s);
		
/*		for (I=0;I<neq;I++)
			fprintf(Out,"s(%d)=%.15lf;\n",I+1,s[I]);
		fclose(Out);		
		break;
*/
		iter = Parameters->iterations;
		pricek = iter - iterold + 1; // tenho duvidas!!!!!
		
		fd = 1;
		cont = 1;
		qtcorr = 0;
		w_old = 0.0;
		printf("=====================???????==============================");			
		//dcopy(neq, F, Fold);   //Fould = F
		do{
			w_old = w;		
			daxpy(neq, w, s, u);		//u = uold + w*s
			//dcopy(neq, F, Fold);   //Fould = F
			//====== Constroi matriz e vetor força ======		
			FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
			//====== Precondiona sistema ======
			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, it, F);
			//====== Atualiza eta ======
			etaold = Parameters->SolverTolerance;
			Parameters->SolverTolerance = eta_newton(Fold, F, etaold, pricek, it, epsilon, eta0, Parameters);
			//========
			norm_F = sqrt(ddot(neq, F, F));
			
			delta = norm_F/norm_Fold;
			printf("\n |F| = %3.2E, |Fold| = %3.2E, |F|/|Fold| = %3.2E \n", norm_F, norm_Fold, delta);
			//========			
			if(norm_F < norm_Fold || w < c1*w_min){ 
				if(norm_F < norm_Fold && fd==1){
					w_max = (1.0 < c3*w_max)? 1.0 : c3*w_max; 
					w = (w_max < c2*w)? w_max : c2*w;
					qtcorr ++;
					printf("\n                 Soluccao aceita!");
					printf("\n Incremento no Fator de amortecimento para w = %f \n", w);
				}else{
					printf("\n                 Soluccao aceita!");
					printf("\n       Com Fator de amortecimento ** w = %f **\n", w);
				}
				cont = 0;
				//printf("\n AQUI 3\n");
			}else{
				//if(cont2==3){
					w = (w_min > c5*w)? w_min : c5*w;
					qtcorr ++;
					if(fd == 1){
						w_max = (w_min > c4*w_max)? w_min : c4*w_max;
						fd = 0;
						//printf("\n AQUI 5\n");
					}
					daxpy(neq, -1.*w_old, s, u);		//u = uold - w_old*s
					printf("\n Soluccao NAO aceita!");
					printf("\n Reduccao   no Fator de amortecimento para w = %f ", w);
					printf("\n ----------------------------------------------------------");
					cont2 = 1;				
				//}else		
				//	cont2 ++;	
			}
			//Fould = F
			//qtcorr ++;			
		}while(cont == 1);
		cont2 = 1;
		dcopy(neq, F, Fold);
		//
		//printf("\n AQUI 2\n");
		fprintf(OutFile2,"\t%d\t%f\t%d\n", it-1, w_old, qtcorr);
		//=======
		norms = sqrt(ddot(neq, s, s));
		normu = sqrt(ddot(neq, u, u));
		delta2 = norms/normu;
		
		
		printf("\n Eta: %3.2E, Iteracao Newton: %d, Iteracao GMRES: %d \n", Parameters->SolverTolerance, it-1, pricek-1);
		printf("\n     |s|/|u| = %3.2E  \n\n", delta2);

		//for (I=0;I<neq;I++)
		//	printf("s[%d]=%.15lf\n",I,s[I]);

					
	}

	printf("===================================================\n");
	printf("\n Iterações de Newton = %d  \n", it-1);
	Parameters-> NLiterations = it-1;
		
	Coef_Lift_Drag(Parameters, MatrixData, FemStructs, FemFunctions);

	//for(i=0;15;i++)
	//	printf("\n %f  %f  %f ", contGlobal[i][0], contGlobal[i][1], contGlobal[i][2]);  

	//SPARILU_clean(MatrixData->ILUp);	
	myfree(s);
	myfree(Fold);
	myfree(delta_old);
	myfree(normres_old);
	fclose(OutFile);
	fclose(OutFile2);

	return 0;
}


















