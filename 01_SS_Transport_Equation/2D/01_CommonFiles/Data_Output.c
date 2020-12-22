#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Data_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, nnodes;
	char FileName[2000];
	FILE *OutFile;
	double *v;
	NodeType *Node = FemStructs->Node;
		
	nnodes = Parameters->nnodes;

	v = (double*) mycalloc("v of 'Data_Output'", nnodes, sizeof(double));
	
	AproximationInEachNode(Parameters, FemStructs, FemFunctions, v);
	
	sprintf(FileName,"../03_output/GNUPLOT_%s_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->nnodes, Parameters->nel);
	
	OutFile = myfopen(FileName,"w");
	
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"%.12lf\t%.12lf\t%.12lf\n", Node[I].x, Node[I].y, v[I]);

	fclose(OutFile);

	free(v);

	return 0;
}
