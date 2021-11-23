#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "SSTransportEquation3D.h"

int Postprocess(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	FILE *OutFile;
	char FileName[300];

	/*************************************************************/
	//		Paraview output to file
	/************************************************************/
	Paraview_Output(Parameters, FemStructs, FemFunctions);
	Gnuplot_Output(Parameters, FemStructs, FemFunctions);
	if (strcasecmp(Parameters->ExactSolution,"YES")==0){
		Paraview_Exact(Parameters, FemStructs, FemFunctions);
	}
	/*************************************************************/

	/****************************************************************************************/
	// 			Printing final result		
	/****************************************************************************************/
	printf("\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
	printf("Experiment: %s\n",Parameters->Experiments);
	printf("Problem Title: %s\n", Parameters->ProblemTitle);
	printf("Number Mesh: %d\n", Parameters->NumberMesh);
	printf("Mesh origin: %s\n", Parameters->OriginMesh);
	printf("Mesh type: %s\n", Parameters->TypeMesh);
	printf("Parameter Cst: %lf\n", Parameters->ConstApli);
	printf("Max Peclet Local: %lf\n", Parameters->MaxPecletLocal);
	printf("Min Peclet Local: %lf\n", Parameters->MinPecletLocal);
	printf("Media Peclet Local: %lf\n", Parameters->MedPecletLocal);
	printf("Max ReaDifRelation Local: %lf\n", Parameters->MaxReaDifRelationLocal);
	printf("Min ReaDifRelation Local: %lf\n", Parameters->MinReaDifRelationLocal);
	printf("Media ReaDifRelation Local: %lf\n", Parameters->MedReaDifRelationLocal);
	printf("Number of nodes: %d\n", Parameters->nnodes);
	printf("Number of elements: %d\n", Parameters->nel);
	printf("Number of equations: %d\n", Parameters->neq);
	printf("Problem has exact solution: %s\n", Parameters->ExactSolution);
	printf("What is the initial solution? %s\n", Parameters->InitialSolution);
	printf("Error type - Absolute Error (AE) or Relative Error (RE): %s\n", Parameters->ErrorType);
	printf("Stabilization form used: %s\n", Parameters->StabilizationForm);
	printf("Use Damping: %s\n", Parameters->UseDamping);
	printf("Compute Residual: %s\n", Parameters->ComputeResidual);
	printf("Use Peclet: %s\n", Parameters->UsePeclet);
	printf("Parameter macro scale: %s\n", Parameters->TauMacro);
	printf("Parameter micro scale: %s\n", Parameters->TauMicro);
	printf("Tolerance in GradU: %lf\n", Parameters->tolGradU);
	printf("Parameter wMacro: %lf\n", Parameters->wMacro);
	printf("Parameter wMicro: %lf\n", Parameters->wMicro);
	printf("Tetha value: %lf\n", Parameters->tetha);
	printf("Use output flow: %s\n", Parameters->OutputFlow);
	printf("Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	printf("Stop Multicorrection: %s\n", Parameters->StopMulticorrection);
	printf("NonLinear tolerance: %lf\n", Parameters->NonLinearTolerance);
	printf("NonLinearMaxIter: %d\n", Parameters->NonLinearMaxIter);
	printf("Tipo de h usado: %d\n", Parameters->TipoH);
	printf("Solver used: %s\n", Parameters->Solver);
	printf("Solver tolerance used: %E\n", Parameters->SolverTolerance);
	printf("Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
	printf("Maximum iterations for Solver: %d\n", Parameters->SolverMaxIter);
	printf("Number of %s iterations ex: %d\n", Parameters->Solver, Parameters->SolverIterations);
	printf("Preconditioner used: %s\n", Parameters->Preconditioner);	
	printf("Reordering used: %s\n", Parameters->Reordering);
	printf("xInf: %lf\n", Parameters->xInf);
	printf("xSup: %lf\n", Parameters->xSup);
	printf("yInf: %lf\n", Parameters->yInf);
	printf("ySup: %lf\n", Parameters->ySup);
	printf("zInf: %lf\n", Parameters->zInf);	
	printf("zSup: %lf\n", Parameters->zSup);
	printf("\n========================================================================\n\n");



	sprintf(FileName,"../03_Output/%s/%s_%s_%s_ExecutionData_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	
	OutFile = myfopen(FileName,"a");
	fprintf(OutFile,"\n\n======================= PROBLEM CHARACTERISTICS ========================\n\n");
	fprintf(OutFile, "Experiment: %s\n",Parameters->Experiments);
	fprintf(OutFile, "Problem Title: %s\n", Parameters->ProblemTitle);
	fprintf(OutFile, "Number of Mesh: %d\n", Parameters->NumberMesh);
	fprintf(OutFile, "Mesh Origin: %s\n", Parameters->OriginMesh);
	fprintf(OutFile, "Mesh Type: %s\n", Parameters->TypeMesh);
	fprintf(OutFile, "Parameter Cst: %lf\n", Parameters->ConstApli);
	fprintf(OutFile, "Max Peclet Local: %lf\n", Parameters->MaxPecletLocal);
	fprintf(OutFile, "Min Peclet Local: %lf\n", Parameters->MinPecletLocal);
	fprintf(OutFile, "Media Peclet Local: %lf\n", Parameters->MedPecletLocal);
	fprintf(OutFile, "Max ReaDifRelation Local: %lf\n", Parameters->MaxReaDifRelationLocal);
	fprintf(OutFile, "Min ReaDifRelation Local: %lf\n", Parameters->MinReaDifRelationLocal);
	fprintf(OutFile, "Media ReaDifRelation Local: %lf\n", Parameters->MedReaDifRelationLocal);
	fprintf(OutFile, "Number of nodes: %d\n", Parameters->nnodes);
	fprintf(OutFile, "Number of elements: %d\n", Parameters->nel);
	fprintf(OutFile, "Number of equations: %d\n", Parameters->neq);
	fprintf(OutFile, "Problem has exact solution: %s\n", Parameters->ExactSolution);
	fprintf(OutFile, "What is the initial solution? %s\n", Parameters->InitialSolution);
	fprintf(OutFile, "Error type - Absolute Error (AE) or Relative Error (RE): %s\n", Parameters->ErrorType);
	fprintf(OutFile, "Stabilization form used: %s\n", Parameters->StabilizationForm);
	fprintf(OutFile, "Use Damping: %s\n", Parameters->UseDamping);
	fprintf(OutFile, "Compute Residual: %s\n", Parameters->ComputeResidual);
	fprintf(OutFile, "Use Peclet: %s\n", Parameters->UsePeclet);
	fprintf(OutFile, "Parameter macro scale: %s\n", Parameters->TauMacro);
	fprintf(OutFile, "Parameter micro scale: %s\n", Parameters->TauMicro);
	fprintf(OutFile, "Tolerance in GradU: %lf\n", Parameters->tolGradU);
	fprintf(OutFile, "Parameter wMacro: %lf\n", Parameters->wMacro);
	fprintf(OutFile, "Parameter wMicro: %lf\n", Parameters->wMicro);
	fprintf(OutFile, "Tetha value: %lf\n", Parameters->tetha);
	fprintf(OutFile, "Use output flow: %s\n", Parameters->OutputFlow);
	fprintf(OutFile, "Matrix vector product scheme: %s\n", Parameters->MatrixVectorProductScheme);
	fprintf(OutFile, "Stop Multicorrection: %s\n", Parameters->StopMulticorrection);
	fprintf(OutFile, "NonLinear tolerance: %lf\n", Parameters->NonLinearTolerance);
	fprintf(OutFile, "NonLinearMaxIter: %d\n", Parameters->NonLinearMaxIter);
	fprintf(OutFile, "Tipo de h usado: %d\n", Parameters->TipoH);
	fprintf(OutFile, "Solver used: %s\n", Parameters->Solver);
	fprintf(OutFile, "Solver tolerance used: %E\n", Parameters->SolverTolerance);
	fprintf(OutFile, "Number of restart: %d\n", Parameters->KrylovBasisVectorsQuantity);
	fprintf(OutFile, "Maximum iterations for Solver: %d\n", Parameters->SolverMaxIter);
	fprintf(OutFile, "Number of %s iterations ex: %d\n", Parameters->Solver, Parameters->SolverIterations);
	fprintf(OutFile, "Preconditioner used: %s\n", Parameters->Preconditioner);	
	fprintf(OutFile, "Reordering used: %s\n", Parameters->Reordering);
	fprintf(OutFile, "xInf: %lf\n", Parameters->xInf);
	fprintf(OutFile, "xSup: %lf\n", Parameters->xSup);
	fprintf(OutFile, "yInf: %lf\n", Parameters->yInf);
	fprintf(OutFile, "ySup: %lf\n", Parameters->ySup);
	fprintf(OutFile, "zInf: %lf\n", Parameters->zInf);	
	fprintf(OutFile, "zSup: %lf\n", Parameters->zSup);
	fprintf(OutFile,"\n========================================================================\n\n");
	fclose(OutFile);
	
	/****************************************************************************************/
	
	/***************************************************************************************/
	//				Memory deallocation	
	/***************************************************************************************/
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3) == 0){		
		myfree(MatrixData->Aaux);
		myfree(MatrixData->A);
		myfree(FemStructs->lmaux);
		myfree(FemStructs->lm);
		//myfree(MatrixData->Gaux);
		//myfree(MatrixData->G);
		//free(MatrixData->Diag);
		//free(MatrixData->invDiag);
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") == 0){		
		int I, nel;

		nel = Parameters->nel;
		
		free(MatrixData->AA);
		free(MatrixData->IA);
		free(MatrixData->JA);
		free(MatrixData->Diag);
		free(MatrixData->invDiag);
		free(FemStructs->lmaux);
		free(FemStructs->lm);
		for (I = 0; I < nel; I++)
			free(MatrixData->Scheme_by_Element[I]);
		free(MatrixData->Scheme_by_Element);
		if (strncmp(Parameters->Preconditioner,"ILU",3)==0){
			SPARILU_clean(MatrixData->ILUp);
			SPARMAT_clean(MatrixData->mat);
			free(MatrixData->Ailu);
		}
	}
	
	myfree(MatrixData);
	myfree(FemStructs->Node);
	myfree(FemStructs->Element);
	myfree(FemStructs->CFF);
	myfree(FemStructs->F);
	myfree(FemStructs->u);
	myfree(FemStructs);
	myfree(FemFunctions);
	myfree(FemOtherFunctions);

	/***************************************************************************************/	
	
	return 0;
}


