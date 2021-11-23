#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

void Calculating_Errors(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{	
	int i;
	double *uh;
	double *uExact;
	double *error;
	double Cst = Parameters->ConstApli;
	//FILE *arq;
	//char ArqName[200];
	int nnodes = Parameters->nnodes;
	NodeType *Node = FemStructs->Node;
	
	//sprintf(ArqName,"../03_Output/%s/%s_%s_N%d_E%d.txt", Parameters->Experiments, Parameters->ProblemTitle, Parameters->Solver, Parameters->nnodes, Parameters->nel);
	//arq = myfopen(ArqName, "a");
		
	uh = mycalloc("uh of 'Calculating_Errors'", nnodes, sizeof(double));
	dzero(nnodes,uh);
	
	AproximationInEachNode(Parameters, FemStructs, FemFunctions, uh);
    
	uExact = mycalloc("uExact of 'Calculating_Errors'", nnodes, sizeof(double));
	dzero(nnodes,uExact);
		
	error = mycalloc("errorx of 'Calculating_Errors'", nnodes, sizeof(double));
	
	FemFunctions->ExactSolutionAllPoints(Node, nnodes, uExact, Cst);
	
	for(i = 0; i < nnodes; i++){
		error[i] = fabs(uh[i] - uExact[i]);
	}
	
	//fprintf(arq, "\nNode: ID U_aprox U_exata Erro\n");
	//for(i = 0; i < nnodes; i++){
	  // fprintf(arq, "%d:%d\t%lf\t%lf\t%lf\n", i, Node[i].id, uh[i], uExact[i], error[i]);
    //}
			
	double ErrorMax = dmax(nnodes,error);
	
	Norm_L2 (uh, Parameters, FemStructs, FemFunctions);
	Norm_H1 (uh, Parameters, FemStructs, FemFunctions);
	
	//double NormErroL2 = New_Norm_L2 (uh, Parameters->nnodes, Node, FemFunctions);
	
	//fprintf(arq, "\n NormMax = %lf\tNormL2 = %lf\tNewNormL2 = %lf\n", ErrorMax, Parameters->NormL2, NormErroL2);
	
	//fclose(arq);
	
	FILE *OutFile;
	char FileName[2000];
	
	sprintf(FileName,"../03_Output/%s/%s_%s_%s_ExecutionData_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	
	OutFile = myfopen(FileName,"a");
	fprintf(OutFile,"\n============================ Error Analise ============================\n\n");
	fprintf(OutFile,"Mesh Parameter h: %lf\n", Parameters->MaxVolume);
	fprintf(OutFile,"Mesh Parameter sqrt[3]{h}: %lf\n", pow(Parameters->MaxVolume,1.0/3.0));
	fprintf(OutFile,"Norm Max Error: %lf\n", ErrorMax);	
	fprintf(OutFile,"Norm L2 Error: %lf\n", Parameters->NormL2);
	fprintf(OutFile,"Norm H1 Error: %lf\n", Parameters->NormH1);
	fprintf(OutFile,"Norm Energy Error: %lf \t %e\n", Parameters->NormEnergy, Parameters->NormEnergy);
	//fprintf(OutFile,"New Norm L2 Error: %lf\n", NormErroL2);
	printf("Mesh Parameter h: %lf\n", Parameters->MaxVolume);	
	printf("Mesh Parameter sqrt[3]{h}: %lf\n", pow(Parameters->MaxVolume,1.0/3.0));
	printf("Norm Max Error: %lf\n", ErrorMax);	
	printf("Norm L2 Error: %lf\n", Parameters->NormL2);
	printf("Norm H1 Error: %lf\n", Parameters->NormH1);
	printf("Norm Energy Error: %lf \t %e\n", Parameters->NormEnergy, Parameters->NormEnergy);
	//printf("New Norm L2 Error: %lf\n", NormErroL2);	
	fprintf(OutFile,"\n=======================================================================\n");
	fclose(OutFile);
	
	myfree(uh);
	myfree(uExact);
	myfree(error);
	
}
