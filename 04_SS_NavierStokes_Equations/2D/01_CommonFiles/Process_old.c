#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Process(MatrixDataType *MatrixData, double *F,double *u, int **lm, NodeType *Node, ElementType *Element, ParametersType *Parameters)
{
	int i, it, itmax;
	int iter, iterold, pricek;	
	double neq, delta1, delta2, epsilon, tol, norms, normu, etaold, eta0;	
	double *s, *Fold;

	double (*ppresc)(double, double);
	double (*v1presc)(double, double);
	double (*v2presc)(double, double);
	void (*assembly)(ParametersType *, int , int,  double (*)[9], MatrixDataType *, int **);
	void (*end_assembly)(ParametersType *Parameters, int tag, MatrixDataType *MatrixData);
	int (*mv)(ParametersType *, MatrixDataType *, int, double *, double *, int **);
	int (*mvMK)(ParametersType *, MatrixDataType *, int, double *, double *, int **);
	int (*precond)(ParametersType *, MatrixDataType *, double *, double *);
	int (*precond_setup)(ParametersType *, MatrixDataType *, NodeType *, ElementType *, int, double *, int **);
	int (*solver)(ParametersType *, MatrixDataType *, double *, double *, int **, 
        int (*)(ParametersType *, MatrixDataType *, double *, double *),int (*)(ParametersType *, MatrixDataType *, int, double *, double *, int **));

	setProblem(Parameters, &v1presc, &v2presc, &ppresc);  
	setMatrixVectorProductType(Parameters, &assembly, &end_assembly, &mv, &mvMK);
	setSolver(Parameters,&solver);
	setPreconditioner(Parameters, &precond_setup,&precond);

//=========Iteração de Newton===========
	neq = Parameters->neq;
	s = (double*) mycalloc("s of 'Process'", neq+1, sizeof(double));
	Fold = (double*) mycalloc("Fold of 'Process'", neq+1, sizeof(double));

	for(i=0; i<neq+1; i++){
		u[i] = 0.0;
		s[i] = 0.0;
		Fold[i] = 0.0;
	}

	Build_K_F_SUPG_PSPG(Parameters, MatrixData, Node, Element, F, u, lm, assembly, end_assembly, ppresc, v1presc, v2presc);
	
	//====== Precondiona sistema ======
	precond_setup(Parameters, MatrixData, Node, Element, it, F, lm);

	delta1 = sqrt(ddot(neq, F, F));	
	delta2 = 1;	
	tol = 2e-1;	
	epsilon = tol*delta1;	
	printf("\n Norma de F_0. |F_0| = %3.2E === Saida de |F| = %3.2E  \n", delta1, epsilon);
	itmax = 1000;
	it = 0;
	iter = 0;
	//Parameters->SolverTolerance = eta_0 (0.1, 0.5, 0.9) 
	eta0 = Parameters->SolverTolerance;
	
	while((delta1 > epsilon || delta2 > tol) && it < itmax){
		it++;		
		iterold = iter;		
		dcopy(neq, F, Fold);   //Fould = F	
		printf("===================================================");	
		printf("\n\n Eta: %3.2E, Iteracao Newton: %d, Pricek: %d \n\n", Parameters->SolverTolerance, it, pricek);

		//====== Resolve sitema linear ======		
		//Parameters->SolverTolerance = 1e-6;
		solver(Parameters, MatrixData, F, s, lm, precond, mv);
		
		iter = Parameters->iterations;
		pricek = iter - iterold + 1; // tenho duvidas!!!!!
		daxpy(neq, 1, s, u);		//u = uold + s
				
		//====== Constroi matriz e vetor força ======		
		Build_K_F_SUPG_PSPG(Parameters, MatrixData, Node, Element, F, u, lm, assembly, end_assembly, ppresc, v1presc, v2presc);
		
		//====== Precondiona sistema ======
		precond_setup(Parameters, MatrixData, Node, Element, it, F, lm);

		//====== Atualiza eta ======
		etaold = Parameters->SolverTolerance;
		Parameters->SolverTolerance = eta_newton(Fold, F, etaold, pricek, it, epsilon, eta0, Parameters);
		
		//====== Condicoes de saida ======		
		delta1 = sqrt(ddot(neq, F, F));	
		norms = sqrt(ddot(neq, s, s));
		normu = sqrt(ddot(neq, u, u));
		delta2 = norms/normu;	
		printf("\n Norma de F. |F| = %3.2E  \n", delta1);		
		printf("\n         |s|/|u| = %3.2E  \n\n", delta2);
		
	}
	printf("===================================================\n");
	printf("\n Iterações de Newton = %d  \n", it);
		

	SPARILU_clean(MatrixData->ILUp);	
	free(s);
	free(Fold);
	
	return 0;
}


















