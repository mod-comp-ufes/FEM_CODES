#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

void Space_Algorithm_Damping(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	int neq = Parameters->neq;
	int nel = Parameters->nel;
	int nnodes = Parameters->nnodes;
	int i, aux, first_damp;
	int NumStepSpace = 0;
	int tag = 1;
	int cont = 0;
	double *uOld, *uTil, *diff, normDiff, normU, *modU;
	double *CbOldMacro, *CbOldMicro;
	double *u, *F;
	double error = 0.0;
	double W, Wmin = 0.01, Wmax = 1.0, c1 = 1.001, c2 = 1.1, c3 = 1.001, c4 = 0.9;
	double Res = 0.0, ResOld = 0.0;
	
	F = FemStructs->F;
	u = FemStructs->u;
	dzero(neq+1,u);

	uOld = mycalloc("uOld of 'Space_Algorithm_Damping'", neq+1, sizeof(double));
	uTil = mycalloc("uOld of 'Space_Algorithm_Damping'", neq+1, sizeof(double));
	diff = mycalloc("diff of 'Space_Algorithm_Damping'", neq, sizeof(double));
	modU = mycalloc("modU of 'Space_Algorithm_Damping'", neq, sizeof(double));
	
	if (strcasecmp(Parameters->InitialSolution,"NULL") == 0){ // Solução nula como solução inicial
		
		double Fnode[nnodes];
		NodeType *Node = FemStructs->Node;

		dzero(nnodes,Fnode);

		for(i = 0; i < nnodes; i++){
			// F in each point in the mesh
			Fnode[i] = FemFunctions->f(Node[i].x, Node[i].y, Node[i].z, Parameters->ConstApli);
		}

		Res = sqrt(ddot(nnodes,Fnode,Fnode));
		printf("Res sol inicial: %lf\n", Res);

	}else if (strcasecmp(Parameters->InitialSolution,"Galekin") == 0){ // Galerkin como solução inicial
		
		Build(Parameters, MatrixData, FemStructs, FemFunctions);
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag, F);
		tag++;
		FemOtherFunctions->Solver(Parameters, MatrixData, FemStructs, FemFunctions, F, u);
		Res = Compute_Residual(Parameters, MatrixData, FemStructs, FemFunctions);
	
	}else if (strcasecmp(Parameters->InitialSolution,"VMS") == 0){ // VMS como solução inicial
		
		Build_VMS(Parameters, MatrixData, FemStructs, FemFunctions);
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag, F);
		tag++;
		FemOtherFunctions->Solver(Parameters, MatrixData, FemStructs, FemFunctions, F, u);
		Res = Compute_Residual(Parameters, MatrixData, FemStructs, FemFunctions);
	
	}

	Parameters->SolverIterations = 0;
	Parameters->ContGMRES = 0;
	
	CbOldMacro = mycalloc("CbOldMacro of 'Space_Algorithm_Damping'", nel, sizeof(double));
	for(i = 0; i < nel; i++){
		CbOldMacro[i] = 0.0;
	}
	FemStructs->CbOldMacro = CbOldMacro;
	
	CbOldMicro = mycalloc("CbOldMicro of 'Space_Algorithm_Damping'", nel, sizeof(double));
	for(i = 0; i < nel; i++){
		CbOldMicro[i] = 0.0;
	}
	FemStructs->CbOldMicro = CbOldMicro;
	
	FILE *OutFile, *ResGMRES, *ResFile;
	char FileName1[250], FileName2[250], FileName3[250];
	
	sprintf(FileName1,"../03_Output/%s/%s_%s_%s_Error_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	OutFile = myfopen(FileName1,"w");
	
	sprintf(FileName2,"../03_Output/%s/%s_%s_%s_ResGMRES_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	ResGMRES = myfopen(FileName2,"w");
	
	sprintf(FileName3,"../03_Output/%s/%s_%s_%s_ResFile_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	ResFile = myfopen(FileName3,"w");
	
	W = Wmax;

	do{
		NumStepSpace++;
		if (NumStepSpace == 1 || NumStepSpace % 10 == 0)
			printf("\nNumStepSpace: %d\n", NumStepSpace);
		fprintf(OutFile, "\nNumStepSpace: %d\n", NumStepSpace);
		
		for(i = 0; i < neq; i++){
			uOld[i] = u[i];
		}
		
		FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
		
		FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag, F);
		tag++;

		FemOtherFunctions->Solver(Parameters, MatrixData, FemStructs, FemFunctions, F, u);
		
		// Guarda a solução parcial
		for(i = 0; i < neq; i++){
			uTil[i] = u[i];
		}
		
		first_damp = 1;
		aux = 1;
		
		ResOld = Res;
		fprintf(ResFile, "%d: \t %e \t %lf\n", NumStepSpace, Res, log(Res));

		cont = 0;
		do{
			cont++;
			//printf("-> -> Passo: %d \t Ponderação: %d \t w = %lf \n", NumStepSpace, cont, W);
			
			//u = uOld + w*(uTil-uOld)
			for(i = 0; i < neq; i++){
				 u[i] = uOld[i] + W*(uTil[i] - uOld[i]);
			}
			
			Res = Compute_Residual(Parameters, MatrixData, FemStructs, FemFunctions);
			//printf("Res: %lf \t ResOld: %lf \t Res-ResOld: %lf\n", Res, ResOld, fabs(Res-ResOld));
			//getchar();

			if((Res < ResOld) || (W <= c1*Wmin)){ // A iteração de u foi aceita ou pq o novo resíduo é menor ou pq não é possível reduzir w
				if((Res < ResOld) && (first_damp == 1)){ // a iteração foi aceita sem que ouvesse rejeição, então o Wmax e w é aumentado
					Wmax = fmin(1.0, c3*Wmax);
					W = fmin(Wmax, c2*W);
				}
				aux = 0;
				//printf("IF -> Wmax = %lf \t W = %lf \t aux = %d\n", Wmax, W, aux);
			}else{ // rejeitou a iteração:w e Wmax são diminuidos
				W = fmax(Wmin, W/2.0);
				if(first_damp == 1){
					Wmax = fmax(Wmin, c4*Wmax);
					first_damp = 0;
				}
				//printf("ELSE -> W = %lf \t Wmax = %lf \t first_damp = %d\n", W, Wmax, first_damp);
			}
		}while(aux == 1);
		
		for(i = 0; i < neq; i++){
			diff[i] = fabs(uOld[i] - u[i]);
		}
		normDiff = dmax(neq, diff);

		for(i = 0; i < neq; i++){
			modU[i] = fabs(u[i]);
		}
		normU = dmax(neq, modU);

		//printf("IterGMRES: %d\n", Parameters->ContGMRES);
		fprintf(OutFile, "IterGMRES: %d\n", Parameters->ContGMRES);
		fprintf(ResGMRES, "%d: \t %e \t %lf\n", NumStepSpace, Parameters->ResGMRES, log(Parameters->ResGMRES));
		//printf("Res: %e \t Res-ResOld: %e \t normDiff/normU: %e \t normDiff: %e \t normU: %e\n", Res, fabs(Res-ResOld), normDiff/normU, normDiff, normU);
		fprintf(OutFile, "Res: %e \t Res-ResOld: %e \t normDiff/normU: %e \t normDiff: %e \t normU: %e\n", Res, fabs(Res-ResOld), normDiff/normU, normDiff, normU);
		
		if (strcasecmp(Parameters->ErrorType,"AE") == 0){
			error = normDiff;
		}else if (strcasecmp(Parameters->ErrorType,"RE") == 0) {
			error = normDiff/normU;
		}else{
			error = fabs(Res - ResOld);
		}
		
	}while( (error > Parameters->NonLinearTolerance) && (NumStepSpace < Parameters->NonLinearMaxIter) ); 
	
	fclose(OutFile);
	fclose(ResGMRES);
	fclose(ResFile);
	
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
	myfree(uTil);
	myfree(diff);
	myfree(modU);
	myfree(CbOldMacro);
	myfree(CbOldMicro);

}	
	
