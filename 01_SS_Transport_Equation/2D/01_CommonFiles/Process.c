#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int i, neq = Parameters->neq;
	double NonLinearTolerance = Parameters->NonLinearTolerance;
	double *u = FemStructs->u;
	double *uold, normu, normdiff;

	setProblem(Parameters,FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters, FemOtherFunctions);
	setScaling(Parameters, FemFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemFunctions, FemOtherFunctions);

	uold = mycalloc("uold",sizeof(double),neq+1);
	normdiff = NonLinearTolerance + 1;
	normu = 1/NonLinearTolerance;
	i = 0;
	while (normdiff > normu*NonLinearTolerance){
		i++;
		dcopy(neq, u, uold);
		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
/*		printf("Before...\n");
		for (i=0; i<neq; i++)
			printf("F[%d]=%lf\n",i,FemStructs->F[i]);i=1;
*/
		FemFunctions->scaling(Parameters, MatrixData, FemStructs);
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, i, FemStructs->F);
/*		printf("After...\n");
		for (i=0; i<neq; i++)
			printf("F[%d]=%lf\n",i,FemStructs->F[i]);
		exit(1);
*/
		FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, FemStructs->F, FemStructs->u);
		FemFunctions->unscaling(Parameters, MatrixData, FemStructs, FemStructs->u);
		daxpy(neq, -1, u, uold);
		normdiff = sqrt(ddot(neq,uold,uold));
		normu = sqrt(ddot(neq,FemStructs->u,FemStructs->u));

		#ifdef debug
			printf("normdiff/normu = %.15lf (normu=%lf)\n ",normdiff/normu, normu);
		#endif
	}

	myfree(uold);

	return 0;
}
