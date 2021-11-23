#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

void Space_Algorithm_NonLinear(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int neq = Parameters->neq;
	int nel = Parameters->nel;
	int i;
	int NumStepSpace = 0;
	int tag = 1;
	double *uOld, *diff, normDiff, normU, *modU;
	double *CbOldMacro, *CbOldMicro;
	double *u, *F;
	double error = DBL_MAX;
	
	// Por padrão a solução inicial é zerada
	F = FemStructs->F;
	u = FemStructs->u;
	dzero(neq,u);
	
	uOld = mycalloc("uOld of 'Space_Algorithm_NonLinear'", neq+1, sizeof(double));
	diff = mycalloc("diff of 'Space_Algorithm_NonLinear'", neq, sizeof(double));
	modU = mycalloc("modU of 'Space_Algorithm_NonLinear'", neq, sizeof(double));
	
	if (strcasecmp(Parameters->InitialSolution,"Galerkin")==0){ // se usar Galerkin como solução inicial
		Build(Parameters, MatrixData, FemStructs, FemFunctions);
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag, F);
		tag++;
		FemOtherFunctions->Solver(Parameters, MatrixData, FemStructs, FemFunctions, F, u);
		//printf("NonLinearMaxIter: %d \t SolverMaxIter: %d \n", Parameters->NonLinearMaxIter, Parameters->SolverMaxIter);
		printf("IterGMRES: %d\n", Parameters->ContGMRES);
	}else if (strcasecmp(Parameters->InitialSolution,"VMS")==0){ // se usar Galerkin como solução inicial
		Build_VMS(Parameters, MatrixData, FemStructs, FemFunctions);
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag, F);
		tag++;
		FemOtherFunctions->Solver(Parameters, MatrixData, FemStructs, FemFunctions, F, u);
		//printf("NonLinearMaxIter: %d \t SolverMaxIter: %d \n", Parameters->NonLinearMaxIter, Parameters->SolverMaxIter);
		printf("IterGMRES: %d\n", Parameters->ContGMRES);
	}
	
	Parameters->SolverIterations = 0;
	
	CbOldMacro = mycalloc("CbOldMacro of 'Space_Algorithm_Nonlinear'", nel, sizeof(double));
	for(i = 0; i < nel; i++){
		CbOldMacro[i] = 0.0;
	}
	FemStructs->CbOldMacro = CbOldMacro;
	
	CbOldMicro = mycalloc("CbOldMicro of 'Space_Algorithm_Nonlinear'", nel, sizeof(double));
	for(i = 0; i < nel; i++){
		CbOldMicro[i] = 0.0;
	}
	FemStructs->CbOldMicro = CbOldMicro;
	
	FILE *OutFile;
	char FileName1[200];
	
	sprintf(FileName1,"../03_Output/%s/%s_%s_%s_Error_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	OutFile = myfopen(FileName1,"w");
	
	FILE *ResSisFile, *ResEqFile, *ResEqL2File;
	char FileName2[200], FileName3[200], FileName4[200];
	
	sprintf(FileName2,"../03_Output/%s/%s_%s_%s_ResidueSistem_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	sprintf(FileName3,"../03_Output/%s/%s_%s_%s_ResiduoEquation_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	sprintf(FileName4,"../03_Output/%s/%s_%s_%s_ResidueEquationL2_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	
	ResSisFile = myfopen(FileName2,"w");
	ResEqFile = myfopen(FileName3,"w");
	ResEqL2File = myfopen(FileName4,"w");
	
	do{
		NumStepSpace++;
		//if (NumStepSpace == 1 || NumStepSpace % 10 == 0)
			printf("\nNumStepSpace: %d\n", NumStepSpace);
		fprintf(OutFile, "\nNumStepSpace: %d\n", NumStepSpace);

		for(i = 0; i < neq; i++){
			uOld[i] = u[i];	
		}

		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);

		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag, F);
		tag++;

		FemOtherFunctions->Solver(Parameters, MatrixData, FemStructs, FemFunctions, F, u);
		printf("IterGMRES: %d\n\n", Parameters->ContGMRES);
		fprintf(OutFile, "IterGMRES: %d\n", Parameters->ContGMRES);

		double residue;
		strcpy(Parameters->ComputeResidual,"Sistem");
		residue =  Compute_Residual(Parameters, MatrixData, FemStructs, FemFunctions);
		//printf(" Residue Sistem step %d: \t %lf \t %lf\n", NumStepSpace, residue, log(residue));
		fprintf(ResSisFile, "%d: \t %lf \t %lf\n", NumStepSpace, residue, log(residue));
		strcpy(Parameters->ComputeResidual,"Equation");
		residue =  Compute_Residual(Parameters, MatrixData, FemStructs, FemFunctions);
		//printf(" Residue Equation step %d: \t %lf \t %lf\n", NumStepSpace, residue, log(residue));
		fprintf(ResEqFile, "%d: \t %lf \t %lf\n", NumStepSpace, residue, log(residue));
		strcpy(Parameters->ComputeResidual,"EquationL2");
		residue =  Compute_Residual(Parameters, MatrixData, FemStructs, FemFunctions);
		//printf(" Residue Equation L2 step %d: \t %lf \t %lf\n", NumStepSpace, residue, log(residue));
		fprintf(ResEqL2File, "%d: \t %lf \t %lf\n", NumStepSpace, residue, log(residue));
		
		for(i = 0; i < neq; i++){
			diff[i] = fabs(uOld[i] - u[i]);
		}
		for(i = 0; i < neq; i++){
			modU[i] = fabs(u[i]);
		}

		normDiff = sqrt(ddot(neq, diff, diff)); // dmax(neq, diff); // 
		normU = sqrt(ddot(neq,modU, modU)); // dmax(neq, modU); // 
		
		//printf("normDiff/normU: %e \t normDiff: %e \t normU: %e\n", normDiff/normU, normDiff, normU);
		fprintf(OutFile, "normDiff/normU: %e \t normDiff: %e \t normU: %e\n", normDiff/normU, normDiff, normU);
		
		if (strcasecmp(Parameters->ErrorType,"AE") == 0){
			error = normDiff;

		}else if (strcasecmp(Parameters->ErrorType,"RE") == 0){
			error = normDiff/normU;
	
		}else {
			printf("\nSpace_Algorithm_NonLinear Error: Error type is not defined correctly!\n");
			exit(1);
		}
		
	}while((NumStepSpace < Parameters->NonLinearMaxIter)&&(error >= Parameters->NonLinearTolerance)); 

	fclose(OutFile);

	fclose(ResSisFile);
	fclose(ResEqFile);
	fclose(ResEqL2File);

	FILE *OutFile2;
	char FileName[200];
	
	sprintf(FileName,"../03_Output/%s/%s_%s_%s_ExecutionData_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	
	OutFile2 = myfopen(FileName,"a");
	fprintf(OutFile2,"\n=======================================================================\n");
	fprintf(OutFile2, "Nonlinear Step: %d\n", NumStepSpace);
	fprintf(OutFile2, "normDiff/normU: %e \t normDiff: %e \t normU: %e\n", normDiff/normU, normDiff, normU);
	fprintf(OutFile2,"\n=======================================================================\n");
	fclose(OutFile2);

	printf("\n=======================================================================\n");
	printf("Nonlinear Step: %d\n", NumStepSpace);
	printf("normDiff/normU: %e \t normDiff: %e \t normU: %e\n", normDiff/normU, normDiff, normU);
	printf("\n=======================================================================\n");

	if (strcasecmp(Parameters->ExactSolution,"YES")==0){
		Calculating_Errors(Parameters, FemStructs, FemFunctions);
	}
	
	myfree(uOld);
	myfree(diff);
	myfree(modU);
	myfree(CbOldMacro);
	myfree(CbOldMicro);

}	
	
