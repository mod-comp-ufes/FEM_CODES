#include "NavierStokesEquations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Paraview_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, eq1, eq2, eq3, nnodes, nel;
	char FileName[300];
	FILE *OutFile;
	double *v1, *v2, *pres, X, Y, t;// normu, normv, normp;
	double v1e, v2e, prese, *errov1, *errov2, *erropres, ev1, ev2, epres;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	double *U = FemStructs->u;
	nnodes = Parameters->nnodes;
	nel = Parameters->nel;
	t = Parameters->time;
	
	//rho = (double*) mycalloc("rho of 'Paraview_Output'", nnodes, sizeof(double));
	v1 = (double*) mycalloc("v1 of 'Paraview_Output'", nnodes, sizeof(double));
	v2 = (double*) mycalloc("v2 of 'Paraview_Output'", nnodes, sizeof(double));
	//e = (double*) mycalloc("e of 'Paraview_Output'", nnodes, sizeof(double));
	//temp = (double*) mycalloc("temp of 'Paraview_Output'", nnodes, sizeof(double));
	pres = (double*) mycalloc("pres of 'Paraview_Output'", nnodes, sizeof(double));
	
	
	errov1 = (double*) mycalloc("v1 of 'Paraview_Output'", nnodes, sizeof(double));
	errov2 = (double*) mycalloc("v2 of 'Paraview_Output'", nnodes, sizeof(double));
	erropres = (double*) mycalloc("pres of 'Paraview_Output'", nnodes, sizeof(double));
	

	for (I = 0; I < nnodes; I++){
		eq1 = Node[I].id[0];
		eq2 = Node[I].id[1];
		eq3 = Node[I].id[2];
		//eq4 = Node[I].id[3];
		X = Node[I].x;	
		Y = Node[I].y;	
		//cv = CV(X, Y);
		//gamma = Gamma(X, Y);

		if (eq1 >= 0)
			v1[I] = U[eq1];
		else
			v1[I] = FemFunctions->v1presc(X, Y);

		if (eq2 >= 0)
			v2[I] = U[eq2];
		else
			v2[I] = FemFunctions->v2presc(X, Y);
	   	
   		if (eq3 >= 0)
   			pres[I] = U[eq3];
   		else
			pres[I] = FemFunctions->ppresc(X, Y);


		//Sol Exata
		v1e = pow(1+t,2)*pow(X,2)*pow(1-X,2)*(2*Y-6*pow(Y,2)+4*pow(Y,3));
		v2e = pow(1+t,2)*pow(Y,2)*pow(1-Y,2)*(-2*X+6*pow(X,2)-4*pow(X,3));
		prese = pow(X,2) - pow(Y,2);	
		//Vetores erro		
		errov1[I] = fabs(v1[I] - v1e);
		errov2[I] = fabs(v2[I] - v2e);
		erropres[I] = fabs(pres[I] - prese);
	}
	
	// Calculo da norma de v e p
	//normu = sqrt(ddot(nnodes, v1, v1));
	//normv = sqrt(ddot(nnodes, v2, v2));
	//normp = sqrt(ddot(nnodes, pres, pres));
	//printf(" \n Normas: |v_x| = %.4lf, |v_y| = %.4lf, |p| = %.4lf (in Paraview_Output.c)\n", normu, normv, normp);
	
	// Calculo erro entre a solucao aprox e exata
	// Norma Euclidiana	
	//ev1 = sqrt(ddot(nnodes, errov1, errov1));
	//ev2 = sqrt(ddot(nnodes, errov2, errov2));
	//epres = sqrt(ddot(nnodes, erropres, erropres));
	// Norma do Maximo	
	ev1 = errov1[0];
	ev2 = errov2[0];
	epres = erropres[0];
	for (I = 1; I < nnodes; I++){
		if(ev1<errov1[I])
			ev1 = errov1[I];	
		if(ev2<errov2[I])
			ev2 = errov2[I];
		if(epres<erropres[I]){
			epres = erropres[I];
			//printf("\n No do erro para a pressao: %d, %f, %f", I, Node[I].x, Node[I].y);		
		}	
	}
	printf(" \n Normas do erro: |Ev_x| = %.4lf, |Ev_y| = %.4lf, |Ep| = %.4lf (in Paraview_Output.c)", ev1, ev2, epres); 


/*	sprintf(FileName,"../03_output/%s_%s_%s_%s_%s_%s_N%d_E%d.vtu", Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->ShockCapture, */
/*			Parameters->TimeIntegration,Parameters->MatrixVectorProductScheme,Parameters->nnodes, Parameters->nel);*/
	sprintf(FileName,"../03_output/%s_%s_%2.1lf_%s_N%d_E%d.vtu", Parameters->Experiments, Parameters->ProblemTitle, t, Parameters->MatrixVectorProductScheme,Parameters->nnodes, Parameters->nel);
	OutFile = myfopen(FileName,"w");

	fprintf(OutFile,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
	fprintf(OutFile,"\t<UnstructuredGrid>\n");
	fprintf(OutFile,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nnodes, nel);
	fprintf(OutFile,"\t\t\t<PointData Scalars=\"scalars\" Vectors = \"Velocity\">\n");
//	fprintf(OutFile,"\t\t\t<PointData Vectors = \"Velocity\">\n");
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf   %.12lf   %.12lf\n", v1[I], v2[I], 0.0);

	fprintf(OutFile,"\t\t\t\t</DataArray>\n");
/*	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n");*/

/*	for (I = 0; I < nnodes; I++)*/
/*		fprintf(OutFile,"\t\t\t\t   %.12lf\n", rho[I]);*/

/*	fprintf(OutFile,"\t\t\t\t</DataArray>\n");*/
/*	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Temperature\" Format=\"ascii\">\n");*/

/*	for (I = 0; I < nnodes; I++)*/
/*		fprintf(OutFile,"\t\t\t\t   %.12lf\n", temp[I]);*/

/*	fprintf(OutFile,"\t\t\t\t</DataArray>\n");*/
	fprintf(OutFile,"\t\t\t\t<DataArray type=\"Float32\" Name=\"Pressure\" Format=\"ascii\">\n");

	for (I = 0; I < nnodes; I++)
		fprintf(OutFile,"\t\t\t\t   %.12lf\n", pres[I]);

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

	//free(rho); 
	free(v1);
	free(v2);
	//free(e); 
	//free(temp);
	free(pres);
		
	return 0;
}


