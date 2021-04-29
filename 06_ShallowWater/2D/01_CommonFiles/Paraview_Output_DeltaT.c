#include "ShalowWater.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"


int Paraview_Output_DeltaT(double *U, NodeType *Node, ElementType *Element, ParametersType *Parameters, 
                           double (*hpresc)(double, double), double (*qxpresc)(double, double), double (*qypresc)(double, double), double t)
{
	int I, eq1, eq2, eq3, nnodes, nel;
	char FileName[2000];
	FILE *OutFile;
	double *h, *qx, *qy, X, Y;

	nnodes = Parameters->nnodes;
	nel = Parameters->nel;

	h = (double*) mycalloc("h of 'Paraview_Output'", nnodes, sizeof(double));
	qx = (double*) mycalloc("qx of 'Paraview_Output'", nnodes, sizeof(double));
	qy = (double*) mycalloc("qy of 'Paraview_Output'", nnodes, sizeof(double));

	for (I = 0; I < nnodes; I++){
		eq1 = Node[I].id[0];
		eq2 = Node[I].id[1];
		eq3 = Node[I].id[2];
		X = Node[I].x;
		Y = Node[I].y;

		if (eq1 >= 0)
			h[I] = U[eq1];
		else
			h[I] = hpresc(X, Y);

		if (eq2 >= 0)
			qx[I] = U[eq2];
		else
			qx[I] = qxpresc(X, Y);

	   	if (eq3 >= 0)
	   		qy[I] = U[eq3];
	   	else
			qy[I] = qypresc(X, Y);

	}
	sprintf(FileName,"%s%s_%s_%s_%s_%s_N%d_E%d_T%lf.vtu", Parameters->outPath, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->ShockCapture, 
			Parameters->MatrixVectorProductScheme, Parameters->nnodes, Parameters->nel, t);
	OutFile = myfopen(FileName,"w");

	fprintf(OutFile,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	fprintf(OutFile,"\t<UnstructuredGrid>\n");
	fprintf(OutFile,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nnodes, nel);

	fprintf(OutFile,"\t\t\t<PointData Scalars=\"scalars\">\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Height\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", h[I]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"XDischarges\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", qx[I]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"YDischarges\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", qy[I]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</PointData>\n");
	fprintf(OutFile,"\t\t\t<Points>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\t%.12lf\t%.12lf\n", Node[I].x,Node[I].y, 0.0);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Points>\n");
	fprintf(OutFile,"\t\t\t<Cells>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

	for (I = 0; I < nel; I++)
		fprintf(OutFile,"\t\t\t\t   %d\t%d\t%d\n", Element[I].Vertex[0], Element[I].Vertex[1], Element[I].Vertex[2]);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
	for (I = 0; I < nel; I++)
		fprintf(OutFile,"\t\t\t\t   %d\n",(I+1)*3);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
	for (I = 0; I < nel; I++)
			fprintf(OutFile,"\t\t\t\t   %d\n",5);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
	fprintf(OutFile,"\t\t\t</Cells>\n");
	fprintf(OutFile,"\t\t</Piece>\n");
	fprintf(OutFile,"\t</UnstructuredGrid>\n");
	fprintf(OutFile,"</VTKFile>\n");

	fclose(OutFile);

	myfree(h);
	myfree(qx);
	myfree(qy);

	return 0;
}
