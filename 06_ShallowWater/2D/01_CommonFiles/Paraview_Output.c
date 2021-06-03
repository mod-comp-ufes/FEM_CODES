#include "ShalowWater.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"


double analitical_h(double x, double y, ParametersType *Parameters);
double analitical_qx(double x, double y, ParametersType *Parameters);


int Paraview_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, eq1, eq2, eq3, nnodes, nel;
	char FileName[2000];
	FILE *OutFile;
	double *h, *qx, *qy, X, Y;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	double *U = FemStructs->u;
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
			h[I] = FemFunctions->hpresc(X, Y);

		if (eq2 >= 0)
			qx[I] = U[eq2];
		else
			qx[I] = FemFunctions->qxpresc(X, Y);

	   	if (eq3 >= 0)
	   		qy[I] = U[eq3];
	   	else
			qy[I] = FemFunctions->qypresc(X, Y);

	}
	sprintf(FileName,"%s%s_%s_%s_%s_%s_N%d_E%d.vtu", Parameters->outPath, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->ShockCapture, 
			Parameters->MatrixVectorProductScheme, Parameters->nnodes, Parameters->nel);
	OutFile = myfopen(FileName,"w");

	fprintf(OutFile,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	fprintf(OutFile,"\t<UnstructuredGrid>\n");
	fprintf(OutFile,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nnodes, nel);

	fprintf(OutFile,"\t\t\t<PointData Scalars=\"scalars\" Vectors = \"Velocity\">\n");	
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

	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"H Exact\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", analitical_h(Node[I].x, Node[I].y, Parameters));
	fprintf(OutFile,"\t\t\t\t</DataArray>\n");

	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Flux Exact\" Format=\"ascii\">\n");
	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", analitical_qx(Node[I].x, Node[I].y, Parameters));
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

double analitical_h(double x, double y, ParametersType *Parameters)
{
    double ksi = 4.183;
	double x13, x32, x20, 
           h0 = 1.0, h1 = 2.0, h2, h3,
           c0, c1, c2, c3,
           u0 = 0.0, u1 = 0.0, u2, u3;

    c0 = sqrt(9.81*h0);
    c1 = sqrt(9.81*h1);
	c2 = c0*sqrt(0.5*sqrt(1 + 8*(ksi/c0)*(ksi/c0)) - 0.5);
	c3 = 1.0/3.0*(2*c1 - x/Parameters->FinalTime);

	h2 = c2*c2/9.81;
	h3 = c3*c3/9.81;

    u2 = ksi - c0*c0/(4*ksi)*(1 + sqrt(1+8*(ksi/c0)*(ksi/c0)));
    u3 = 2.0/3.0*(2*c1 + x/Parameters->FinalTime);

	x13 = -c1*Parameters->FinalTime;
	x32 = (u2 - c2)*Parameters->FinalTime;
	x20 = ksi*Parameters->FinalTime;
	
	if(x < x13)
	    return h1;
	else if(x < x32)
		return h3;
	else if(x < x20)
		return h2;
	else
	    return h0;
}

double analitical_qx(double x, double y, ParametersType *Parameters)
{
    double ksi = 4.183;
	double x13, x32, x20, 
           h0 = 1.0, h1 = 2.0, h2, h3,
           c0, c1, c2, c3,
           u0 = 0.0, u1 = 0.0, u2, u3,
		   q0, q1, q2, q3;

    c0 = sqrt(9.81*h0);
    c1 = sqrt(9.81*h1);
	c2 = c0*sqrt(0.5*sqrt(1+8*(ksi/c0)*(ksi/c0))- 0.5);
	c3 = 1/3*(2*c1 - x/Parameters->FinalTime);

	h2 = c2*c2/9.81;
	h3 = c3*c3/9.81;

    u2 = ksi - c0*c0/(4*ksi)*(1 + sqrt(1+8*(ksi/c0)*(ksi/c0)));
    u3 = 2/3*(2*c1 + x/Parameters->FinalTime);

	x13 = -c1*Parameters->FinalTime;
	x32 = (u2 - c2)*Parameters->FinalTime;
	x20 = ksi*Parameters->FinalTime;

	q0 = h0*u0;
	q1 = h1*u1;
	q2 = h2*u2;
	q3 = h3*u3;
	
	if(x < x13)
	    return q1;
	else if(x < x32)
		return q3;
	else if(x < x20)
		return q2;
	else
	    return q0;
}