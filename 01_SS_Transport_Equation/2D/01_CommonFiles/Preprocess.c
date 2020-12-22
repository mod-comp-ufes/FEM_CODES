#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

int Preprocess(int narg, char **arguments, ParametersType **Parameters_out,  MatrixDataType **MatrixData_out,
			FemStructsType **FemStructs_out, FemFunctionsType **FemFunctions_out, FemOtherFunctionsType **FemOtherFunctions_out)
{
	int neq, nnodes, nedge, nel, I, J;
	int tag = 1; // Testing input error
	int **lm, *lmaux;
	int size = NDOF*NNOEL;
	int size2 = size*size;
	double *F, *u, *Diag, *invDiag;
	char FileName[1000], label[1000];
	FILE *InFile;
	NodeType *Node;
	ElementType *Element;
	ParametersType *Parameters;
	MatrixDataType *MatrixData;
	FemStructsType *FemStructs;
	FemFunctionsType *FemFunctions;
	FemOtherFunctionsType *FemOtherFunctions;

	/* **************************************************************************************************************************** */
	//												Testing initial parameters
	/* **************************************************************************************************************************** */
	if (narg!=2)
	{
		printf("Use ./SSTranspEquation2D <Parameters file according README>\n");
		exit(1);
	}
	/* **************************************************************************************************************************** */


	/***************************************************************************/
	//			Reading parameters from problem setting file
	/***************************************************************************/
	Parameters   = mycalloc("Parameters of 'Preprocess'",1,sizeof(ParametersType));
	MatrixData   = mycalloc("MatrixData of 'Preprocess'",1,sizeof(MatrixDataType));
	FemStructs   = mycalloc("FemStructs of 'Preprocess'",1,sizeof(FemStructsType));
	FemFunctions = mycalloc("FemFunctions of 'Preprocess'",1,sizeof(FemFunctionsType));
	FemOtherFunctions = mycalloc("FemOtherFunctions of 'Preprocess'",1,sizeof(FemOtherFunctionsType));

	InFile = myfopen(arguments[1], "r");
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->Experiments, label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->ProblemTitle, label);
	tag = fscanf(InFile, "%lf\t:%[^\n]\n", &(Parameters->SolverTolerance), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]\n", &(Parameters->NonLinearTolerance), label);
	tag = fscanf(InFile, "%d\t:%[^\n]\n", &(Parameters->SolverMaxIter), label);
	tag = fscanf(InFile, "%d\t:%[^\n]\n", &(Parameters->KrylovBasisVectorsQuantity), label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->Solver, label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->Preconditioner, label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->Scaling, label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->reordering, label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->MatrixVectorProductScheme, label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->StabilizationForm, label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->ShockCapture, label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->h_Shock, label);
	tag = fscanf(InFile, "%s\t:%[^\n]\n", Parameters->OutputFlow, label);
	tag = fscanf(InFile, "%d\t:%[^\n]\n",&(Parameters->nnodes), label);
	tag = fscanf(InFile, "%d\t:%[^\n]\n",&(Parameters->nel), label);
	fclose(InFile);
	/*****************************************************************************/


	/*****************************************************************************/
	//				Reading nodes
	/*****************************************************************************/
	sprintf(FileName,"../02_mesh/%s_%d_%d.dat", Parameters->ProblemTitle, Parameters->nnodes, Parameters->nel);
	InFile = myfopen(FileName, "r");
	tag =  fscanf(InFile, "%d", &nnodes);
	Node = (NodeType*) mycalloc("Node of 'Preprocess'",nnodes, sizeof(NodeType));
	for (I = 0, neq = 0; I < nnodes; I++)
	{
		tag = fscanf(InFile, "%lf%lf%d\n", &(Node[I].x), &(Node[I].y), &(Node[I].Type));

		if (Node[I].Type == 1)
			Node[I].id = neq++;
		else
			Node[I].id = -1;
	}
	/*****************************************************************************/
	
/*	for(I = 0; I < nnodes; I++){
		printf("%lf \t %lf\n", Node[I].x, Node[I].y);	
	}*/
	
/*	for(I = 0; I < nnodes; I++){
		printf("id = %d\n", Node[I].id);	
	}*/

	/************************************************************************************************************/
	//						Reading connection mesh
	/************************************************************************************************************/
	tag = fscanf(InFile, "%d", &nel);
	Element = (ElementType*) mycalloc("Element of 'Preprocess'",nel, sizeof(ElementType));
	for (I = 0; I < nel; I++)
		tag = fscanf(InFile, "%d%d%d%d\n", &(Element[I].Vertex[0]),  &(Element[I].Vertex[1]), &(Element[I].Vertex[2]), &(Element[I].Type));
	fclose(InFile);
	
/*	for(I = 0; I < nel; I++){
		printf("%d \t %d \t %d\n", Element[I].Vertex[0],  Element[I].Vertex[1], Element[I].Vertex[2]);	
	}
	getchar();*/
	
	int p, cont1, cont2, cont3, cont4;
	
	for(I = 0; I < nel; I++){
		cont1 = 0; cont2 = 0; cont3 = 0; cont4 = 0;
		for(J = 0; J < NNOEL; J++){
			p = Element[I].Vertex[J];
			if(Node[p].x == 0.0){ // reta 1: x = 0 
				cont1 = cont1 + 1;
			}
			if(Node[p].y == 0.0){ // reta 2: y = 0
				cont2 = cont2 + 1;
			}
			if(Node[p].x == 1.0){ // reta 3: x = 1
				cont3 = cont3 + 1;
			}
			if(Node[p].y == 1.0){ // reta 4: y = 1
				cont4 = cont4 + 1;
			}
		}
		
			
		if(cont1 == 2){
			if(cont1 == cont2){
				Element[I].Type = 5; // esquina entre x = 0 e y = 0
			}else if(cont1 == cont4){
				Element[I].Type = 8; // esquina entre y = 1 e x = 0
			}else{
				Element[I].Type = 1;
			}
		}else if(cont2 == 2){
			if(cont2 == cont3){
				Element[I].Type = 6; // esquina entre y = 0 e x = 1
			}else{
				Element[I].Type = 2;
			}
		}else if(cont3 == 2){
			if(cont3 == cont4){
				Element[I].Type = 7; // esquina entre x = 1 e y = 1
			}else{
				Element[I].Type = 3;
			}
		}else if(cont4 == 2){
				Element[I].Type = 4;
		}else{ // não está em nenhum contorno
			Element[I].Type = -1;
		}
	}
	/*************************************************************************************************************/


	/*****************************************************************************/
	//                      Memory allocations and Store strategies
        /*****************************************************************************/

	// Some variable inicializations

	MatrixData = (MatrixDataType *) mycalloc("MatrixData of 'Preprocess'",1, sizeof(MatrixDataType));
	F = (double*) mycalloc("F of 'Preprocess'",neq+1, sizeof(double));
	u = (double*) mycalloc("u of 'Preprocess'",neq+1, sizeof(double));
	Diag = (double*) mycalloc("Diag of 'Preprocess'",neq+1, sizeof(double));
	invDiag = (double*) mycalloc("invDiag of 'Preprocess'",neq+1, sizeof(double));
	lm = (int**) mycalloc("lm of 'Preprocess'",nel, sizeof(int*));
	lmaux = (int*) mycalloc("lmaux of 'Preprocess'",nel*size, sizeof(int));
	for (I=0;I<nel;I++)
		lm[I] = &lmaux[I*size];

	//Configuring equation according to variables and boundary conditions
	FemStructs->u = u;
	FemStructs->F = F;
	Fill_LM(neq, nel, lm, Node, Element);
	FemStructs->lm = lm;
	FemStructs->lmaux = lmaux;

	Parameters->neq = neq;

	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){

		double **K, *Kaux;

		K = (double**) mycalloc("K of 'Preprocess'", nel, sizeof(double*));
		Kaux = (double*) mycalloc("Kaux of 'Preprocess'", nel*size2,sizeof(double));

		for (I=0;I<nel;I++){
			K[I] = &Kaux[I*size2];
		}
		MatrixData->A = K;
		MatrixData->Aaux = Kaux;

	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){

		int **order, **EDGE_by_Element;
		double **K, *Kaux;

		order = mycalloc("order of 'Preprocess'", nel, sizeof(int*));
		for (I = 0; I < nel; I++)
			order[I] = (int*) mycalloc("order of 'Preprocess'", size, sizeof(int));

		ede_Initialization(Parameters, order, &lm, &lmaux, &EDGE_by_Element);
		nedge = Parameters->nedge;

		K = (double**) mycalloc("K of 'Preprocess'", nedge+1, sizeof(double*));
		Kaux = (double*) mycalloc("Kaux of 'Preprocess'", (nedge+1)*(NNOEL+1),sizeof(double));

		for (I = 0; I < nedge + 1; I++)
			K[I] = &Kaux[I*(NNOEL+1)];

		FemStructs->lm2 = lm;
		FemStructs->lm2aux = lmaux;
		MatrixData->A = K;
		MatrixData->Aaux = Kaux;
		MatrixData->Scheme_by_Element = EDGE_by_Element;
		MatrixData->order = order;

	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){

		int *IA, *JA, *perm, *invperm, **CSR_by_Element;
		double *K;

		csr_Initialization(Parameters, Node, &JA, &IA, &perm, &invperm, &lm, &lmaux, &CSR_by_Element);

		K = (double*) mycalloc("K of 'Preprocess'", Parameters->nnzero+1, sizeof(double));

		printf("nnzero=%d\n",Parameters->nnzero);

		MatrixData->AA = K;
		MatrixData->IA = IA;
		MatrixData->JA = JA;
		MatrixData->Scheme_by_Element = CSR_by_Element;
		MatrixData->Perm = perm;
	}
	else
		printf("Matrix vector produto scheme not defined!\n\n");

	/************************************************************************************/

	FemStructs->Node = Node;
	FemStructs->Element = Element;
	FemStructs->F = F;
	FemStructs->u = u;
	MatrixData->Diag = Diag;
	MatrixData->invDiag = invDiag;

	*Parameters_out = Parameters;
	*MatrixData_out = MatrixData;
	*FemStructs_out = FemStructs;
	*FemFunctions_out = FemFunctions;
	*FemOtherFunctions_out = FemOtherFunctions;

	if (tag<0){
		printf ("Error in some parameter\n");
		exit(1);
	}

	return 0;
}
