#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"


int Paraview_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, nel, nnodes;
	char FileName[2000];
	FILE *OutFile;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	double *U = FemStructs->u;

	nel = Parameters->nel;
	nnodes = Parameters->nnodes;

	sprintf(FileName, "../03_output/%s_%s_%s_%s_%s_%s_%s_N%d_E%d.vtu",Parameters->Experiments, Parameters->ProblemTitle,Parameters->StabilizationForm,Parameters->ShockCapture, Parameters->h_Shock,
		Parameters->MatrixVectorProductScheme, Parameters->Preconditioner, Parameters->nnodes,Parameters->nel); 	
	OutFile = myfopen(FileName,"w");

	fprintf(OutFile,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	fprintf(OutFile,"\t<UnstructuredGrid>\n");
	fprintf(OutFile,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nnodes, nel);
	fprintf(OutFile,"\t\t\t<PointData Scalars=\"scalars\">\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"solucao\" Format=\"ascii\">\n");


	for (I=0; I<nnodes; I++)
		if (Node[I].id>=0)
			fprintf(OutFile,"\t\t\t\t   %.12f\n",U[Node[I].id]);
		else
			fprintf(OutFile,"\t\t\t\t   %.12f\n", FemFunctions->upresc(Node[I].x,Node[I].y));

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</PointData>\n");
	fprintf(OutFile,"\t\t\t<Points>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");


        // Para visualizar em 3D//**************************************************************************************************
	for (I=0; I<nnodes; I++)
		if (Node[I].id>=0)
			fprintf(OutFile,"\t\t\t\t   %.12f\t%.12f\t%.12f\n", Node[I].x,Node[I].y, U[Node[I].id]);
		else

			fprintf(OutFile,"\t\t\t\t   %.12f\t%.12f\t%.12f\n", Node[I].x,Node[I].y, FemFunctions->upresc(Node[I].x,Node[I].y));

	//************************************************************************************************************************

	// Para visualizar em 2D************************************************************************************************
//	for (I=0; I<nnodes; I++)
//			fprintf(OutFile,"\t\t\t\t   %.12f\t%.12f\t%.1f\n", Node[I].x,Node[I].y, 0.0);
	//***********************************************************************************************************************

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Points>\n");
	fprintf(OutFile,"\t\t\t<Cells>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

	for (I=0;I<nel;I++)
		fprintf(OutFile,"\t\t\t\t   %d\t%d\t%d\n",Element[I].Vertex[0], Element[I].Vertex[1],Element[I].Vertex[2]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");

	for (I=0;I<nel;I++)
		fprintf(OutFile,"\t\t\t\t   %d\n",(I+1)*3);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");

	for (I=0;I<nel;I++)
			fprintf(OutFile,"\t\t\t\t   %d\n",5);


	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Cells>\n");
	fprintf(OutFile,"\t\t</Piece>\n");
	fprintf(OutFile,"\t</UnstructuredGrid>\n");
	fprintf(OutFile,"</VTKFile>\n");

	fclose(OutFile);

	return 0;
}
