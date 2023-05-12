#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"

int Postprocess(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int I, nel, neq;
	char FileName[2000];
	FILE *OutFile;

	nel = Parameters->nel;
	neq = Parameters->neq;

	/*************************************************************/
	//		Paraview output to file
	/************************************************************/
	Paraview_Output(Parameters, FemStructs, FemFunctions);
	Data_Output(Parameters, FemStructs, FemFunctions);
	/*************************************************************/


	/****************************************************************************************/
		// 			Printing final result
	/****************************************************************************************/
	sprintf(FileName, "../03_output/%s_%s_%s_%s_%s_%s_%s_N%d_E%d.dat",Parameters->Experiments, Parameters->ProblemTitle,Parameters->StabilizationForm,Parameters->ShockCapture,
	        Parameters->h_Shock, Parameters->MatrixVectorProductScheme, Parameters->Preconditioner, Parameters->nnodes,Parameters->nel);
	OutFile = myfopen(FileName,"w");
	fprintf(OutFile, "\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
	fprintf(OutFile, "Problem Title: Steady State %s\n", Parameters->ProblemTitle);
	fprintf(OutFile, "Number of nodes: %d\n", Parameters->nnodes);
	fprintf(OutFile, "Number of elements: %d\n", nel);
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
		fprintf(OutFile, "Number of edges: %d\n", Parameters->nedge);
	}
	fprintf(OutFile, "Number of equations: %d\n", neq);
	fprintf(OutFile, "Stabilization form used: %s\n", Parameters->StabilizationForm);
	fprintf(OutFile, "Shock capture used: %s\n", Parameters->ShockCapture);
	fprintf(OutFile, "Parameter h for shock capture used: %s\n", Parameters->h_Shock);
	fprintf(OutFile, "Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		fprintf(OutFile,"Reordering: %s (bandwidth before: %d) (bandwidth after: %d)\n", Parameters->reordering,
		Parameters->bandwidth_bef, Parameters->bandwidth_aft);
	}
	fprintf(OutFile, "Solver used: %s\n", Parameters->Solver);
	fprintf(OutFile, "Preconditioner used: %s\n", Parameters->Preconditioner);
	fprintf(OutFile, "Scaling used: %s\n", Parameters->Scaling);
	fprintf(OutFile, "Solver tolerance used: %E\n", Parameters->SolverTolerance);
	fprintf(OutFile, "Non linear tolerance used: %E\n", Parameters->NonLinearTolerance);
	fprintf(OutFile, "Number of %s iterations: %d\n", Parameters->Solver, Parameters->SolverIterations);
	fprintf(OutFile, "\n========================================================================\n\n");
	fclose(OutFile);

	printf("\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
	printf("Problem Title: Steady State %s\n", Parameters->ProblemTitle);
	printf("Number of nodes: %d\n", Parameters->nnodes);
	printf("Number of elements: %d\n", nel);
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
		printf("Number of edges: %d\n", Parameters->nedge);
	}
	printf("Number of equations: %d\n", neq);
	printf("Stabilization form used: %s\n", Parameters->StabilizationForm);
	printf("Shock capture used: %s\n", Parameters->ShockCapture);
	printf("Parameter h for shock capture used: %s\n", Parameters->h_Shock);
	printf("Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		printf("Reordering: %s (bandwidth before: %d) (bandwidth after: %d)\n", Parameters->reordering,
		Parameters->bandwidth_bef, Parameters->bandwidth_aft);
	}
	printf("Solver used: %s\n", Parameters->Solver);
	printf("Preconditioner used: %s\n", Parameters->Preconditioner);
	printf("Scaling used: %s\n", Parameters->Scaling);
	printf("Solver tolerance used: %E\n", Parameters->SolverTolerance);
	printf("Non linear tolerance used: %E\n", Parameters->NonLinearTolerance);
	printf("Number of %s iterations: %d\n", Parameters->Solver, Parameters->SolverIterations);
	printf("\n========================================================================\n\n");

	/****************************************************************************************/



	/***************************************************************************************/
	//				Memory deallocation
	/**************************************************************************************/
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
		free(MatrixData->A);
		free(MatrixData->Aaux);
		free(FemStructs->lm);
		free(FemStructs->lmaux);
	}

	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
		free(MatrixData->A);
		free(MatrixData->Aaux);
		for (I = 0; I < nel; I++){
			free(MatrixData->Scheme_by_Element[I]);
			free(MatrixData->order[I]);
		}
		free(MatrixData->Scheme_by_Element);
		free(MatrixData->order);
		free(FemStructs->lm);
		free(FemStructs->lmaux);
		free(FemStructs->lm2);
		free(FemStructs->lm2aux);

	}

	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		free(MatrixData->AA);
		free(MatrixData->IA);
		free(MatrixData->JA);
		for (I = 0; I < nel; I++)
			free(MatrixData->Scheme_by_Element[I]);
		free(MatrixData->Scheme_by_Element);
		if (strncmp(Parameters->Preconditioner,"ILU",3)==0){
			SPARILU_clean(MatrixData->ILUp);
			SPARMAT_clean(MatrixData->mat);
			free(MatrixData->Ailu);
		}
	}

	if ((strcasecmp(Parameters->Preconditioner,"Jacobi")==0)||(strncmp(Parameters->Preconditioner,"SOR",3)==0)){
		free(MatrixData->invDe);
		free(MatrixData->invDeaux);
	}

	free(MatrixData->Diag);
	free(MatrixData->invDiag);
	free(MatrixData);
	free(FemStructs->Node);
	free(FemStructs->Element);
	free(FemStructs->F);
	free(FemStructs->u);
	free(FemStructs);
	free(FemFunctions);
	free(FemOtherFunctions);
	/***************************************************************************************/

	return 0;
}
