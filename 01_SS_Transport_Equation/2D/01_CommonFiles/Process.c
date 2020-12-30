#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int i, I, neq = Parameters->neq;
	int nel = Parameters->nel;
	double NonLinearTolerance = Parameters->NonLinearTolerance;
	double *u = FemStructs->u;
	double *uold, *diff, normu, normdiff;
	double error, *CbOld;

	setProblem(Parameters,FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters, FemOtherFunctions);
	setScaling(Parameters, FemFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions);
	
	CbOld = mycalloc("CbOld of 'Process'", nel, sizeof(double));
	for(I = 0; I < nel; I++)
		CbOld[I] = 0.0;
	
	FemStructs->CbOld = CbOld;
	
	dzero(neq+1,u);
	uold = mycalloc("uold",sizeof(double),neq+1);
	diff = mycalloc("diff",sizeof(double),neq+1);
	normdiff = NonLinearTolerance + 1;
	normu = 1/NonLinearTolerance;
	error = 100;
	i = 0;
	while (error > NonLinearTolerance && i <= 50){
		i++;
		
		printf("Iter Non Linear: %d\n", i);
		
		dcopy(neq, u, uold); // uold = u
		
		//for(I = 0; I < neq; I++){
		//	printf("uold[%d] = %lf\n", I, uold[I]);
		//}
		//getchar();
		
		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
//		printf("Before...\n");
//		for (i=0; i<neq; i++)
//			printf("F[%d]=%lf\n",i,FemStructs->F[i]);		
//		getchar();

	//	FemFunctions->scaling(Parameters, MatrixData, FemStructs);
	//	FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, i, FemStructs->F);
/*		printf("After...\n");
		for (i=0; i<neq; i++)
			printf("F[%d]=%lf\n",i,FemStructs->F[i]);
		exit(1);
*/
		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions);
	//	FemFunctions->unscaling(Parameters, MatrixData, FemStructs, FemStructs->u);
		//dcopy(neq, uold, diff); // diff = uold
		//daxpy(neq, -1, u, diff); // diff = diff - u
		for(I = 0; I < neq; I++){
			diff[I] = fabs(u[I] - uold[I]);
		}
		error = dmax(neq, diff);
		normdiff = sqrt(ddot(neq,diff,diff));
		normu = sqrt(ddot(neq,FemStructs->u,FemStructs->u));

		#ifdef debug
			printf("error = %.15lf normdiff/normu = %.15lf (normu=%lf)\n ",error, normdiff/normu, normu);
		#endif
	}

	free(uold);
	free(diff);
	free(CbOld);
		
	return 0;
}
