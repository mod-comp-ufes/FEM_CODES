#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "NavierStokesEquations.h"

int Postprocess(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	FILE *OutFile;
	char FileName[300];

	/*************************************************************/
	//		Paraview output to file
	/************************************************************/	
	Paraview_Output(Parameters, FemStructs, FemFunctions);
	//Paraview_Output_3D(Parameters, FemStructs, FemFunctions);
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
	//printf("Discontinuities capture operator used: %s\n", Parameters->ShockCapture);
	printf("Stabilization coefficient tolerance: %E\n", Parameters->CoefficientTolerance);
	printf("Time integration method: %s\n", Parameters->TimeIntegration);
	printf("Time integration tolerance: %E\n", Parameters->TimeIntegrationTolerance);
	printf("Correction tolerance: %E\n", Parameters->CorrectionTolerance);
	printf("Correction stopped: %s\n", Parameters->StopMulticorrection);
	printf("Maximum number of correction iteration: %d\n", Parameters->NumberCorrection);
	printf("Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		printf("Reordering: %s (bandwidth before: %d) (bandwidth after: %d)\n", Parameters->reordering, 
		Parameters->bandwidth_bef, Parameters->bandwidth_aft);
	}
	printf("Solver used: %s\n", Parameters->Solver);
	printf("Solver tolerance used: %E\n", Parameters->SolverTolerance);
	printf("Preconditioner used: %s\n", Parameters->Preconditioner);
	printf("Maximum number of solver iteration: %d\n", Parameters->LinearMaxIter);
	printf("Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
	printf("Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);
	printf("Alpha: %lf \t Step time: %E \t Final Time: %lf\n", Parameters->Alpha, Parameters->DeltaT, Parameters->FinalT);	
	printf("\n========================================================================\n\n");

	sprintf(FileName,"../03_output/%s_%s_%s_%s_%s_%s_N%d_E%d.txt", Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->ShockCapture, 
			Parameters->TimeIntegration, Parameters->MatrixVectorProductScheme, Parameters->nnodes, Parameters->nel);
	
	OutFile = myfopen(FileName,"w");
	fprintf(OutFile, "\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
	fprintf(OutFile, "Problem Title: %s\n", Parameters->ProblemTitle);
	fprintf(OutFile, "Number of nodes: %d\n", Parameters->nnodes);
	fprintf(OutFile, "Number of elements: %d\n", Parameters->nel);
	fprintf(OutFile, "Number of equations: %d\n", Parameters->neq);
	fprintf(OutFile, "Stabilization form used: %s\n", Parameters->StabilizationForm);
	//fprintf(OutFile, "Discontinuities capture operator used: %s\n", Parameters->ShockCapture);
	fprintf(OutFile, "Stabilization coefficient tolerance: %E\n", Parameters->CoefficientTolerance);
	fprintf(OutFile, "Time integration method: %s\n", Parameters->TimeIntegration);
	fprintf(OutFile, "Time integration tolerance: %E\n", Parameters->TimeIntegrationTolerance);
	fprintf(OutFile, "Correction tolerance: %E\n", Parameters->CorrectionTolerance);
	fprintf(OutFile, "Correction stopped: %s\n", Parameters->StopMulticorrection);
	fprintf(OutFile, "Maximum number of correction iteration: %d\n", Parameters->NumberCorrection);
	fprintf(OutFile, "Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	fprintf(OutFile, "Reordering: %s\n", Parameters->reordering);
	fprintf(OutFile, "Solver used: %s\n", Parameters->Solver);
	fprintf(OutFile, "Solver tolerance used: %E\n", Parameters->SolverTolerance);
	fprintf(OutFile, "Maximum number of solver iteration: %d\n", Parameters->LinearMaxIter);
	fprintf(OutFile, "Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
	fprintf(OutFile, "Number of %s iterations: %d\n", Parameters->Solver, Parameters->iterations);
	fprintf(OutFile, "Alpha: %lf \t Step time: %E \t Final Time: %lf", Parameters->Alpha, Parameters->DeltaT, Parameters->FinalT);	
	fprintf(OutFile, "\n========================================================================\n\n");
	fclose(OutFile);
	
	/****************************************************************************************/
	
	
	/***************************************************************************************/
	//				Memory deallocation	
	/**************************************************************************************/
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3) == 0){		
		free(MatrixData->A);
		free(MatrixData->Aaux);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
	}
	else if (strncmp(Parameters->MatrixVectorProductScheme,"EDE",3) == 0){		
		int I, nel;

		nel = Parameters->nel;
		free(MatrixData->A);
		free(MatrixData->Aaux);
		for (I = 0; I < nel; I++){
			free(MatrixData->Scheme_by_Element[I]);
			free(MatrixData->order[I]);
		}
		free(MatrixData->Scheme_by_Element);
		free(MatrixData->order);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
		free(FemStructs->lm2aux);
		free(FemStructs->lm2);
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") == 0){		
		int I, nel;

		nel = Parameters->nel;
		free(MatrixData->AA);
		free(MatrixData->JA);
		free(MatrixData->IA);
		free(MatrixData->Diag);
		free(MatrixData->invDiag);
		for (I = 0; I < nel; I++)
			free(MatrixData->Scheme_by_Element[I]);
		free(MatrixData->Scheme_by_Element);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
		SPARILU_clean(MatrixData->ILUp);
		
	}	


	
	free(MatrixData);

	free(FemStructs->Node);
	free(FemStructs->Element);
	free(FemStructs->F);
	free(FemStructs->u);
	free(FemStructs->eqpress);
	free(FemStructs);
	free(FemFunctions);
	free(FemOtherFunctions);

	


	/***************************************************************************************/	
	
	return 0;
}



