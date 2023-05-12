#include "EulerEquations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Postprocess(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	FILE *OutFile;
	char FileName[2000];

	/*************************************************************/
	//		Paraview output to file
	/************************************************************/
	Paraview_Output(Parameters, FemStructs, FemFunctions);
	Paraview_Output_3D(Parameters, FemStructs, FemFunctions);
	/*************************************************************/


	/****************************************************************************************/
	// 			Printing final result
	/****************************************************************************************/
	printf("\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
	printf("Problem Title: %s\n", Parameters->ProblemTitle);
	printf("Number of nodes: %d\n", Parameters->nnodes);
	printf("Number of elements: %d\n", Parameters->nel);
	printf("Number of equations: %d\n", Parameters->neq);
	printf("Stabilization form used: %s\n", Parameters->StabilizationForm);
	printf("Discontinuities capture operator used: %s\n", Parameters->ShockCapture);
	printf("Stabilization coefficient tolerance: %E\n", Parameters->StabilizationTolerance);
	printf("Time integration method: %s\n", Parameters->TimeIntegration);
	printf("Time integration tolerance: %E\n", Parameters->TimeIntegrationTolerance);
	printf("Non linear tolerance: %E\n", Parameters->NonLinearTolerance);
	printf("Correction stopped: %s\n", Parameters->StopMulticorrection);
	printf("Maximum number of correction iteration: %d\n", Parameters->NonLinearMaxIter);
	printf("Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		printf("Reordering: %s (bandwidth before: %d) (bandwidth after: %d)\n", Parameters->reordering,
		Parameters->bandwidth_bef, Parameters->bandwidth_aft);
	}
	printf("Solver used: %s\n", Parameters->Solver);
	printf("Solver tolerance used: %E\n", Parameters->SolverTolerance);
	printf("Preconditioner used: %s\n", Parameters->Preconditioner);
	printf("Scaling used: %s\n", Parameters->Scaling);
	printf("Maximum number of solver iteration: %d\n", Parameters->LinearMaxIter);
	printf("Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
	printf("Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);
	printf("Alpha: %lf\t Step time: %E\t Final Time: %lf - Stopped at Steady State? %s (Current Time: %lf)\n", Parameters->Alpha, Parameters->DeltaT,
		Parameters->FinalTime, Parameters->StopAtSteadyState, Parameters->CurrentTime);
	printf("\n========================================================================\n\n");

	sprintf(FileName,"../../../../OUTPUT_DATA/%s_%s_%s_%s_%s_%s_%s_N%d_E%d.txt", Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->ShockCapture,
			Parameters->TimeIntegration, Parameters->MatrixVectorProductScheme, Parameters->Preconditioner, Parameters->nnodes, Parameters->nel);

	OutFile = myfopen(FileName,"w");
	fprintf(OutFile, "\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
	fprintf(OutFile, "Problem Title: %s\n", Parameters->ProblemTitle);
	fprintf(OutFile, "Number of nodes: %d\n", Parameters->nnodes);
	fprintf(OutFile, "Number of elements: %d\n", Parameters->nel);
	fprintf(OutFile, "Number of equations: %d\n", Parameters->neq);
	fprintf(OutFile, "Stabilization form used: %s\n", Parameters->StabilizationForm);
	fprintf(OutFile, "Discontinuities capture operator used: %s\n", Parameters->ShockCapture);
	fprintf(OutFile, "Stabilization coefficient tolerance: %E\n", Parameters->StabilizationTolerance);
	fprintf(OutFile, "Time integration method: %s\n", Parameters->TimeIntegration);
	fprintf(OutFile, "Time integration tolerance: %E\n", Parameters->TimeIntegrationTolerance);
	fprintf(OutFile, "Non linear tolerance: %E\n", Parameters->NonLinearTolerance);
	fprintf(OutFile, "Correction stopped: %s\n", Parameters->StopMulticorrection);
	fprintf(OutFile, "Maximum number of correction iteration: %d\n", Parameters->NonLinearMaxIter);
	fprintf(OutFile, "Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	fprintf(OutFile, "Reordering: %s\n", Parameters->reordering);
	fprintf(OutFile, "Solver used: %s\n", Parameters->Solver);
	fprintf(OutFile, "Solver tolerance used: %E\n", Parameters->SolverTolerance);
	fprintf(OutFile, "Maximum number of solver iteration: %d\n", Parameters->LinearMaxIter);
	fprintf(OutFile, "Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
	fprintf(OutFile, "Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);
	fprintf(OutFile, "Alpha: %lf\t Step time: %E\t Final Time: %lf - Stopped at Steady State? %s (Current Time: %lf)\n", Parameters->Alpha, Parameters->DeltaT,
		Parameters->FinalTime, Parameters->StopAtSteadyState, Parameters->CurrentTime);
	fprintf(OutFile, "\n========================================================================\n\n");
	fclose(OutFile);
	/****************************************************************************************/


	/***************************************************************************************/
	//				Memory deallocation
	/**************************************************************************************/
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3) == 0){
		myfree(MatrixData->A);
		myfree(MatrixData->Aaux);
		myfree(FemStructs->lmaux);
		myfree(FemStructs->lm);
	}
	else if (strncmp(Parameters->MatrixVectorProductScheme,"EDE",3) == 0){
		int I, nel;

		nel = Parameters->nel;
		myfree(MatrixData->A);
		myfree(MatrixData->Aaux);
		for (I = 0; I < nel; I++){
			myfree(MatrixData->Scheme_by_Element[I]);
			myfree(MatrixData->order[I]);
		}
		myfree(MatrixData->Scheme_by_Element);
		myfree(MatrixData->order);
		myfree(FemStructs->lmaux);
		myfree(FemStructs->lm);
		myfree(FemStructs->lm2aux);
		myfree(FemStructs->lm2);
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") == 0){
		int I, nel;

		nel = Parameters->nel;
		myfree(MatrixData->AA);
		myfree(MatrixData->JA);
		myfree(MatrixData->IA);
		myfree(MatrixData->Diag);
		myfree(MatrixData->invDiag);
		for (I = 0; I < nel; I++)
			myfree(MatrixData->Scheme_by_Element[I]);
		myfree(MatrixData->Scheme_by_Element);
		myfree(FemStructs->lmaux);
		myfree(FemStructs->lm);
		if ((strncmp(Parameters->Preconditioner,"ILU",3)==0)){
			SPARILU_clean(MatrixData->ILUp);
			SPARMAT_clean(MatrixData->mat);
			myfree(MatrixData->Ailu);
		}

	}

	if ((strcasecmp(Parameters->Preconditioner,"Jacobi")==0)||(strncmp(Parameters->Preconditioner,"SOR",3)==0)){
		myfree(MatrixData->invDe);
		myfree(MatrixData->invDeaux);
	}

	myfree(MatrixData);
	myfree(FemStructs->Node);
	myfree(FemStructs->Element);
	myfree(FemStructs->F);
	myfree(FemStructs->u);
	myfree(FemStructs->eqrho);
	myfree(FemStructs);
	myfree(FemFunctions);
	myfree(FemOtherFunctions);

	/***************************************************************************************/

	return 0;
}
