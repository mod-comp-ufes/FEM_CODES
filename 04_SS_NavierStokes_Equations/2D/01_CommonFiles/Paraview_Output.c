#include "SSNavierStokesEquations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Paraview_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	int I, eq1, eq2, eq3, nnodes, nel;
	char FileName[2000];
	FILE *OutFile;
	double *v1, *v2, *pres, X, Y, normu, normv, normp;
	double twoArea;
	double y31, y12, x13, x21, X1, X2, X3, Y1, Y2, Y3, erv1, erv2, erv3, erv4, errovL2, erp1, erp2, erp3, erp4, erropL2, s, t;
	double u11, u12, u13, u21, u22, u23, p1, p2, p3, ve1, ve2, vh1, vh2, presh;
	double v1e, v2e, prese, *errov1, *errov2, *erropres;//, ev1, ev2, epres,
	int E, J1, J2, J3;
	
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	double *U = FemStructs->u;
	nnodes = Parameters->nnodes;
	nel = Parameters->nel;
	
	//rho = (double*) mycalloc("rho of 'Paraview_Output'", nnodes, sizeof(double));
	v1 = (double*) mycalloc("v1 of 'Paraview_Output'", nnodes, sizeof(double));
	v2 = (double*) mycalloc("v2 of 'Paraview_Output'", nnodes, sizeof(double));
	//e = (double*) mycalloc("e of 'Paraview_Output'", nnodes, sizeof(double));
	//temp = (double*) mycalloc("temp of 'Paraview_Output'", nnodes, sizeof(double));
	pres = (double*) mycalloc("pres of 'Paraview_Output'", nnodes, sizeof(double));
	

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

	}
	
	// Calculo da norma de v e p
	normu = sqrt(ddot(nnodes, v1, v1));
	normv = sqrt(ddot(nnodes, v2, v2));
	normp = sqrt(ddot(nnodes, pres, pres));
	printf(" \n Normas: |v_x| = %.4lf, |v_y| = %.4lf, |p| = %.4lf (in Paraview_Output.c)\n", normu, normv, normp);
	Parameters->normu = normu;
	Parameters->normv = normv;
	Parameters->normp = normp;
	
// Calculo erro entre a solucao aprox e exata
	if(strcasecmp(Parameters->ProblemTitle,"EXATA")==0){
	
		
		errov1 = (double*) mycalloc("v1 of 'Paraview_Output'", nnodes, sizeof(double));
		errov2 = (double*) mycalloc("v2 of 'Paraview_Output'", nnodes, sizeof(double));
		erropres = (double*) mycalloc("pres of 'Paraview_Output'", nnodes, sizeof(double));
		for (I = 0; I < nnodes; I++){
			X = Node[I].x;	
			Y = Node[I].y;
			//Sol Exata
			v1e = pow(X,2)*pow(1.0 - X,2)*(2*Y - 6*pow(Y,2) + 4*pow(Y,3));
			v2e = pow(Y,2)*pow(1.0 - Y,2)*(-2*X + 6*pow(X,2) - 4*pow(X,3));
			prese = pow(X,2) - pow(Y,2);	
			//Vetores erro		
			errov1[I] = fabs(v1[I] - v1e);
			errov2[I] = fabs(v2[I] - v2e);
			erropres[I] = fabs(pres[I] - prese);	
		}
		/*
		// Norma Euclidiana	
		ev1 = sqrt(ddot(nnodes, errov1, errov1));
		ev2 = sqrt(ddot(nnodes, errov2, errov2));
		epres = sqrt(ddot(nnodes, erropres, erropres));
		printf(" \n Normas do erro: |Ev_x| = %.6lf, |Ev_y| = %.6lf, |Ep| = %.6lf (in Paraview_Output.c)", ev1, ev2, epres); 			*/
		
		// Norma do Maximo	
/*		ev1 = errov1[0];
		ev2 = errov2[0];
		epres = erropres[0];
		for (I = 1; I < nnodes; I++){
			//if(fabs(Node[I].x-0.5)<1e-12)			
				if(ev1<errov1[I]){
					ev1 = errov1[I];
					//printf("\n No do erro V_x: %d, %f, %f", I, Node[I].x, Node[I].y);	
				}			
			//if(fabs(Node[I].y-0.5)<1e-12)
				if(ev2<errov2[I]){
					ev2 = errov2[I];
					//printf("\n No do erro V_y: %d, %f, %f", I, Node[I].x, Node[I].y);
				}
			if(epres<erropres[I]){
				epres = erropres[I];
				//printf("\n No do erro para a pressao: %d, %f, %f", I, Node[I].x, Node[I].y);		
			}	
		}
		printf(" \n Normas do erro: |Ev_x| = %E, |Ev_y| = %E, |Ep| = %E (in Paraview_Output.c)", ev1, ev2, epres); 
*/		
		// Norma do L2	
		errovL2 = 0.0;
		erropL2 = 0.0;
		prese = 0;		
		for (E=0; E<nel; E++){
			J1 = Element[E].Vertex[0];
			J2 = Element[E].Vertex[1];
			J3 = Element[E].Vertex[2];

			X1 = Node[J1].x;
			X2 = Node[J2].x;
			X3 = Node[J3].x;
			Y1 = Node[J1].y;
			Y2 = Node[J2].y;
			Y3 = Node[J3].y;

			//y23 = Y2 - Y3;
			y31 = Y3 - Y1;
			y12 = Y1 - Y2;

			//x32 = X3 - X2;	
			x13 = X1 - X3;
			x21 = X2 - X1;

			twoArea =  fabs(x21*y31 - x13*y12);
			//Area = 0.5*twoArea;
			//invArea = 1.0/Area;

			//==== velocidade-x aproximada nodal
			u11 = v1[J1]; 
			u12 = v1[J2];
			u13 = v1[J3];
			//==== velocidade-y aproximada nodal
			u21 = v2[J1];
			u22 = v2[J2];
			u23 = v2[J3];
			//==== pressao aproximada nodal
			p1 = pres[J1];
			p2 = pres[J2];
			p3 = pres[J3];

			//==== s e t: pontos de integracao de Gauss, trabalharemos com 4 pontos
			s = 1./3;
			t = 1./3;
			//==== transformada do plano (s,t) para o plano (x,y)
			X = X1 + (X2-X1)*s + (X3-X1)*t;
			Y = Y1 + (Y2-Y1)*s + (Y3-Y1)*t;	
			//==== solucao exata no ponto (s,t)
			ve1 = pow(X,2)*pow(1.0 - X,2)*(2*Y - 6*pow(Y,2) + 4*pow(Y,3));
			ve2 = pow(Y,2)*pow(1.0 - Y,2)*(-2*X + 6*pow(X,2) - 4*pow(X,3));
			prese = pow(X,2) - pow(Y,2);
			//==== solucao aproximada interpolada no ponto (s,t)
			vh1 = u11*(1-s-t) + u12*s + u13*t; 
			vh2 = u21*(1-s-t) + u22*s + u23*t;
			presh = p1*(1-s-t) + p2*s + p3*t;
			//==== integrando no primeiro ponto de Gauss			
			erv1 = ((ve1-vh1)*(ve1-vh1) + (ve2-vh2)*(ve2-vh2))*twoArea;
			erp1 = (prese - presh)*(prese - presh)*twoArea;	
		
			//==== 2 Ponto
			s = 1./5;
			t = 1./5;
			X = X1 + (X2-X1)*s + (X3-X1)*t;
			Y = Y1 + (Y2-Y1)*s + (Y3-Y1)*t;	
			ve1 = pow(X,2)*pow(1.0 - X,2)*(2*Y - 6*pow(Y,2) + 4*pow(Y,3));
			ve2 = pow(Y,2)*pow(1.0 - Y,2)*(-2*X + 6*pow(X,2) - 4*pow(X,3));
			prese = pow(X,2) - pow(Y,2);
			vh1 = u11*(1-s-t) + u12*s + u13*t; 
			vh2 = u21*(1-s-t) + u22*s + u23*t;
			presh = p1*(1-s-t) + p2*s + p3*t;
			erv2 = ((ve1-vh1)*(ve1-vh1) + (ve2-vh2)*(ve2-vh2))*twoArea;
			erp2 = (prese - presh)*(prese - presh)*twoArea;

			// ==== 3 Ponto
			s = 1./5;
			t = 3./5;
			X = X1 + (X2-X1)*s + (X3-X1)*t;
			Y = Y1 + (Y2-Y1)*s + (Y3-Y1)*t;	
			ve1 = pow(X,2)*pow(1.0 - X,2)*(2*Y - 6*pow(Y,2) + 4*pow(Y,3));
			ve2 = pow(Y,2)*pow(1.0 - Y,2)*(-2*X + 6*pow(X,2) - 4*pow(X,3));
			prese = pow(X,2) - pow(Y,2);
			vh1 = u11*(1-s-t) + u12*s + u13*t; 
			vh2 = u21*(1-s-t) + u22*s + u23*t;
			presh = p1*(1-s-t) + p2*s + p3*t;
			erv3 = ((ve1-vh1)*(ve1-vh1) + (ve2-vh2)*(ve2-vh2))*twoArea;
			erp3 = (prese - presh)*(prese - presh)*twoArea;

			// ==== 4 Ponto
			s = 3./5;
			t = 1./5;
			X = X1 + (X2-X1)*s + (X3-X1)*t;
			Y = Y1 + (Y2-Y1)*s + (Y3-Y1)*t;	
			ve1 = pow(X,2)*pow(1.0 - X,2)*(2*Y - 6*pow(Y,2) + 4*pow(Y,3));
			ve2 = pow(Y,2)*pow(1.0 - Y,2)*(-2*X + 6*pow(X,2) - 4*pow(X,3));
			prese = pow(X,2) - pow(Y,2);
			vh1 = u11*(1-s-t) + u12*s + u13*t; 
			vh2 = u21*(1-s-t) + u22*s + u23*t;
			presh = p1*(1-s-t) + p2*s + p3*t;
			erv4 = ((ve1-vh1)*(ve1-vh1) + (ve2-vh2)*(ve2-vh2))*twoArea;
			erp4 = (prese - presh)*(prese - presh)*twoArea;

			//
			errovL2 += -27./96*erv1 + 25./96*(erv2 + erv3 + erv4); 
			erropL2 += -27./96*erp1 + 25./96*(erp2 + erp3 + erp4); 
		} 
		printf("\n Norma L2 do erro: |u-u_h|_L2 = %E, |p-p_h|_L2 = %E (in Paraview_Output.c)", sqrt(errovL2), sqrt(erropL2)); 
	}


	sprintf(FileName,"../03_output/%s_%s_%s_%s_%s_N%d_E%d.vtu", Parameters->ProblemTitle, Parameters->Experiments, Parameters->StabilizationForm, Parameters->MatrixVectorProductScheme, Parameters->Preconditioner, Parameters->nnodes, Parameters->nel);
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
	myfree(v1);
	myfree(v2);
	//free(e); 
	//free(temp);
	myfree(pres);
		
	return 0;
}


