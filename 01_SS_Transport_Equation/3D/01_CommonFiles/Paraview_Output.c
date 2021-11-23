#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int Paraview_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, nnodes, nel;
	char FileName[2000];
	FILE *OutFile;
	double *v1;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	nnodes = Parameters->nnodes;
	nel = Parameters->nel;

	v1 = (double*) mycalloc("v1 of 'Paraview_Output'", nnodes, sizeof(double));
	
	AproximationInEachNode(Parameters, FemStructs, FemFunctions, v1);
	
	// 3D
	sprintf(FileName,"../03_Output/%s/%s_%s_%s_PARAVIEW_N%d_E%d.vtu", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	OutFile = myfopen(FileName,"w");

	fprintf(OutFile,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	fprintf(OutFile,"\t<UnstructuredGrid>\n");
	fprintf(OutFile,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nnodes, nel);
	fprintf(OutFile,"\t\t\t<PointData Scalars=\"Pressure\" Vectors=\"Velocity\">\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Velocity X\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", v1[I]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</PointData>\n");
	fprintf(OutFile,"\t\t\t<Points>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\t%.12lf\t%.12lf\n", Node[I].x, Node[I].y, Node[I].z);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Points>\n");
	fprintf(OutFile,"\t\t\t<Cells>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
		fprintf(OutFile,"\t\t\t\t   %d\t%d\t%d\t%d\n", Element[I].Vertex[0], Element[I].Vertex[1], Element[I].Vertex[2], Element[I].Vertex[3]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
		fprintf(OutFile,"\t\t\t\t   %d\n",(I+1)*4);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
			fprintf(OutFile,"\t\t\t\t   %d\n",10);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Cells>\n");
	fprintf(OutFile,"\t\t</Piece>\n");
	fprintf(OutFile,"\t</UnstructuredGrid>\n");
	fprintf(OutFile,"</VTKFile>\n");

	fclose(OutFile);
	
	free(v1);
	
	return 0;
}
