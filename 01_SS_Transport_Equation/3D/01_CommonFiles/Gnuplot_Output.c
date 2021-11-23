#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Gnuplot_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, nnodes;
	char FileName[2000];
	FILE *OutFile;
	double *v;
	double ySup = Parameters->ySup;
	NodeType *Node = FemStructs->Node;

	nnodes = Parameters->nnodes;

	v = (double*) mycalloc("v of 'Gnuplot_Output'", nnodes, sizeof(double));
	
	AproximationInEachNode(Parameters, FemStructs, FemFunctions, v);
	
	// Todo o grafico
	sprintf(FileName,"../03_Output/%s/%s_%s_%s_GNUPLOT_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	
	OutFile = myfopen(FileName,"w");
	
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"%.12lf\t%.12lf\t%.12lf\t%.12lf\n", Node[I].x, Node[I].y, Node[I].z, v[I]);

	fclose(OutFile);

	if (strcasecmp(Parameters->TypeMesh,"STRUCTURED") == 0){ // Se a malha é estruturada gera o arquivo de saída de corte do Gnuplot
	
		if(strcasecmp(Parameters->ProblemTitle,"CHANNEL") == 0){
			// Corte em y = Ly/2
			sprintf(FileName,"../03_Output/%s/%s_%s_%s_CORTE_GNUPLOT_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	
			OutFile = myfopen(FileName,"w");
	
			for (I = 0; I < nnodes; I++){
				if(fabs(Node[I].y - ySup/2.0) < 1e-6){
					fprintf(OutFile,"%.12lf\t%.12lf\t%.12lf\n", Node[I].z, Node[I].x, v[I]);
				}
			}
	
			fclose(OutFile);
		
		}else if(strcasecmp(Parameters->ProblemTitle,"PAREDE")==0){
			// Corte em z = 1.0
			sprintf(FileName,"../03_Output/%s/%s_%s_%s_CORTE_GNUPLOT_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	
			OutFile = myfopen(FileName,"w");
	
			for (I = 0; I < nnodes; I++){
				if(fabs(Node[I].z - 1.0) < 1e-6){
					fprintf(OutFile,"%.12lf\t%.12lf\t%.12lf\n", Node[I].x, Node[I].y, v[I]);
				}
			}
	
			fclose(OutFile);
		
		}else{
		// Corte em z = 0.5
		sprintf(FileName,"../03_Output/%s/%s_%s_%s_CORTE_GNUPLOT_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
		OutFile = myfopen(FileName,"w");
	
			for (I = 0; I < nnodes; I++){
				if(fabs(Node[I].z - 0.5) < 1e-6){
					fprintf(OutFile,"%.12lf\t%.12lf\t%.12lf\n", Node[I].x, Node[I].y, v[I]);
				}
			}
	
			fclose(OutFile);
			
		}
	
		if (strcasecmp(Parameters->ExactSolution,"YES") == 0){
			FemFunctions->ExactSolutionAllPoints(Node, nnodes, v, Parameters->ConstApli); // Calculate exact solution for each node
		
			// Corte em z = 0.5
			sprintf(FileName,"../03_Output/%s/%s_%s_%s_EXACT_CORTE_GNUPLOT_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
			OutFile = myfopen(FileName,"w");
	
			for (I = 0; I < nnodes; I++){
				if(fabs(Node[I].z - 0.5) < 1e-6){
					fprintf(OutFile,"%.12lf\t%.12lf\t%.12lf\n", Node[I].x, Node[I].y, v[I]);
				}
			}
	
			fclose(OutFile);
		}
	
	}
	
	free(v);

	return 0;
}
