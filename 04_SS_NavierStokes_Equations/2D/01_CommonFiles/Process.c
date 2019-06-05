#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int i, neq, nel;
	int NonLinearMaxIter = Parameters->NonLinearMaxIter;
	double NonLinearTolerance = Parameters->NonLinearTolerance; 
	double delta2, norms, normu, normF;	
	double *s, *u, *F, *delta_old;
	
	setProblem(Parameters, FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters,FemOtherFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions);

	u = FemStructs->u;
	F = FemStructs->F;

	//=========Iteração de Newton===========
	neq = Parameters->neq;
	nel = Parameters->nel;
	s = (double*) mycalloc("s of 'Process'", neq+1, sizeof(double));
	delta_old = (double*) mycalloc("delta_old of 'Process'", nel, sizeof(double));
	FemStructs->delta_old = delta_old;

	i = 0;
	delta2 = NonLinearTolerance + 1.0;

	while(delta2 > NonLinearTolerance && i < NonLinearMaxIter){
		i++;		
		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, i, F);
		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, F, s);
		daxpy(neq, 1, s, u);		//u = uold + s
				
		norms = sqrt(ddot(neq, s, s));
		normu = sqrt(ddot(neq, u, u));
		normF = sqrt(ddot(neq, F, F));
		delta2 = norms/normu;	
		printf("===================================================");	
		printf("\n         |s|/|u| = %3.2E |s|=%lf |u|=%lf |F|=%lf (iteration %d)\n\n", delta2, norms, normu, normF, i);
		printf("===================================================");	
	}

	printf("===================================================\n");
	printf("\n Newton iterations = %d  \n", i);
		

	free(s);
	free(delta_old);
	
	return 0;
}


















