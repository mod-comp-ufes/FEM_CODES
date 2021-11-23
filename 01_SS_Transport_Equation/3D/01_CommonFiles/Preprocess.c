#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

int Preprocess(int narg, char **arguments, ParametersType **Parameters_out, MatrixDataType **MatrixData_out, FemStructsType **FemStructs_out, FemFunctionsType **FemFunctions_out, FemOtherFunctionsType **FemOtherFunctions_out)
{
	int neq, nnodes, nel, I, aux;
	int tag = 1; // Testing input error
	int **lm, *lmaux;
	int size = NDOF*NNOEL;
	int size2 = size*size;
	double *F, *u;
	char FileName[300], label[300];
	FILE *InFile;
	NodeType *Node;
	ElementType *Element;
	CoefFormFuncType *CFF;
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
		printf("Use ./PoissonEquations3D <Parameters file according README>\n");
		exit(1);
	}
	/* **************************************************************************************************************************** */

	/* ************************************************************************* */
	//			Reading parameters from problem setting file
	/* ************************************************************************* */
	Parameters   = mycalloc("Parameters of 'Preprocess'",1,sizeof(ParametersType));
	MatrixData   = mycalloc("MatrixData of 'Preprocess'",1,sizeof(MatrixDataType));
	FemStructs   = mycalloc("FemStructs of 'Preprocess'",1,sizeof(FemStructsType));
	FemFunctions = mycalloc("FemFunctions of 'Preprocess'",1,sizeof(FemFunctionsType));
	FemOtherFunctions = mycalloc("FemOtherFunctions of 'Preprocess'",1,sizeof(FemOtherFunctionsType));

	InFile = myfopen(arguments[1], "r");
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->Experiments, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->ProblemTitle, label);
	tag = fscanf(InFile, "%d\t:%[^\n]", &Parameters->NumberMesh, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->OriginMesh, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->TypeMesh, label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &Parameters->ConstApli, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->ExactSolution, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->InitialSolution, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->Solver, label);
	tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->KrylovBasisVectorsQuantity), label);
	tag = fscanf(InFile, "%d\t:%[^\n]", &Parameters->SolverMaxIter, label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->SolverTolerance), label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->Preconditioner, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->Reordering, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->MatrixVectorProductScheme, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->StabilizationForm, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->UseDamping, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->ComputeResidual, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->UsePeclet, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->TauMacro, label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->TauMicro, label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->tolGradU), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->wMacro), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->wMicro), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->tetha), label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->ErrorType, label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->NonLinearTolerance), label);
	tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->NonLinearMaxIter), label);
	tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->TipoH), label);
	tag = fscanf(InFile, "%s\t:%[^\n]", Parameters->OutputFlow, label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->Vmax), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->xInf), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->xSup), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->yInf), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->ySup), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->zInf), label);
	tag = fscanf(InFile, "%lf\t:%[^\n]", &(Parameters->zSup), label);
	tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->nnodes), label);
	tag = fscanf(InFile, "%d\t:%[^\n]", &(Parameters->nel), label);
	fclose(InFile);
	
	/* **************************************************************************************************************************** */
/*	printf("%s\n", Parameters->Experiments);
	printf("%s\n", Parameters->ProblemTitle);
	printf("%d\n", Parameters->NumberMesh);
	printf("%s\n", Parameters->OriginMesh);
	printf("%s\n", Parameters->TypeMesh);
	printf("%lf\n", Parameters->ConstApli);
	printf("%s\n", Parameters->ExactSolution);
	printf("%s\n", Parameters->InitialSolution);
	printf("%s\n", Parameters->Solver);
	printf("%d\n", Parameters->KrylovBasisVectorsQuantity);
	printf("%d\n", Parameters->SolverMaxIter);
	printf("%e\n", Parameters->SolverTolerance);
	printf("%s\n", Parameters->Preconditioner);
	printf("%s\n", Parameters->Reordering);
	printf("%s\n", Parameters->MatrixVectorProductScheme);
	printf("%s\n", Parameters->StabilizationForm);
	printf("%s\n", Parameters->UseDamping);
	printf("%s\n", Parameters->ComputeResidual);
	printf("%s\n", Parameters->UsePeclet);
	printf("%s\n", Parameters->TauMacro);
	printf("%s\n", Parameters->TauMicro);
	printf("%e\n", Parameters->tolGradU);
	printf("%lf\n", Parameters->wMacro);
	printf("%lf\n", Parameters->wMicro);
	printf("%lf\n", Parameters->tetha);
	printf("%s\n", Parameters->ErrorType);
	printf("%e\n", Parameters->NonLinearTolerance);
	printf("%d\n", Parameters->NonLinearMaxIter);
	printf("%d\n", Parameters->TipoH);
	printf("%s\n", Parameters->OutputFlow);
	printf("%lf\n", Parameters->Vmax);
	printf("%lf\n", Parameters->xInf);
	printf("%lf\n", Parameters->xSup);
	printf("%lf\n", Parameters->yInf);
	printf("%lf\n", Parameters->ySup);
	printf("%lf\n", Parameters->zInf);
	printf("%lf\n", Parameters->zSup);
	printf("%d\n",  Parameters->nnodes);
	printf("%d\n",  Parameters->nel);
	getchar();    */

	/* **************************************************************************************************************************** */
	//						Reading nodes
	/* **************************************************************************************************************************** */
	if(strcasecmp(Parameters->ProblemTitle,"CHANNEL")!=0){ // if problem is not channel
		sprintf(FileName,"../02_Mesh/%s/%s_%d_%d.dat", Parameters->ProblemTitle, Parameters->ProblemTitle, Parameters->nnodes, Parameters->nel);
	}else{ // if problem is channel
		sprintf(FileName,"../02_Mesh/%s/%s_M%d_%d_%d.dat", Parameters->ProblemTitle, Parameters->ProblemTitle, Parameters->NumberMesh, Parameters->nnodes, Parameters->nel);
	}
	
	InFile = myfopen(FileName, "r");
	tag = fscanf(InFile, "%d", &nnodes);
	Node = (NodeType*) mycalloc("Node of 'Preprocess'", nnodes, sizeof(NodeType));
	if (strcasecmp(Parameters->OriginMesh,"GMSH")==0){
		for (I = 0; I < nnodes; I++){
			// for mesh by GMSH
			tag = fscanf(InFile, "%lf%lf%lf%d\n", &(Node[I].x), &(Node[I].y), &(Node[I].z), &(Node[I].uType));
		}
	}else if(strcasecmp(Parameters->OriginMesh,"FREEFEM")==0){
		for (I = 0; I < nnodes; I++){
			// for mesh freefem
			tag = fscanf(InFile, "%lf%lf%lf%d\n", &(Node[I].x), &(Node[I].y), &(Node[I].z), &aux);
			if(aux == 0)//incognita
				Node[I].uType = 1;
			else
				Node[I].uType = 0;		
		} 
	}else{
		printf("Origin Mesh not defined!\n\n");
		exit(1);
	}

	Fill_ID(&neq, Node, nnodes);
	
	/* **************************************************************************************************************************** */

/*	printf("\n --- Coordenadas e Tipo dos 20 primeiros Nós --- \n");
	for(I = 0; I < 20; I++){
		printf("%d: %lf\t%lf\t%lf\t%d\n", I, Node[I].x, Node[I].y, Node[I].z, Node[I].uType);
	}
	getchar();   */
	
/*	printf("\n --- Matriz ID dos 20 primeiros nós ---\n");
	for(I = 0; I < 20; I++){
		printf("%d: \t %d\n", I, Node[I].id);
	}
	getchar();   */

	/* **************************************************************************************************************************** */
	//           				Reading connection mesh
	/* **************************************************************************************************************************** */
	int ver1, ver2, ver3, ver4;
	tag = fscanf(InFile, "%d", &nel);
	Element = (ElementType*) mycalloc("Element of 'Preprocess'", nel, sizeof(ElementType));
	if (strcasecmp(Parameters->OriginMesh,"GMSH")==0){
		for (I = 0; I < nel; I++){
			// for mesh by GMSH
			tag = fscanf(InFile, "%d%d%d%d%d", &(Element[I].Vertex[0]), &(Element[I].Vertex[1]), &(Element[I].Vertex[2]), &(Element[I].Vertex[3]), &(Element[I].Type));
		}
	}else if (strcasecmp(Parameters->OriginMesh,"FREEFEM")==0){
		for (I = 0; I < nel; I++){
			// for mesh freefem
			tag = fscanf(InFile, "%d%d%d%d%d\n", &(ver1), &(ver2), &(ver3), &(ver4), &aux);
			Element[I].Vertex[0] = ver1 - 1;
			Element[I].Vertex[1] = ver2 - 1;
			Element[I].Vertex[2] = ver3 - 1;
			Element[I].Vertex[3] = ver4 - 1;
		}
	}else{
		printf("Origin Mesh not defined!\n\n");
		exit(1);
	}

	fclose(InFile);
	
	ElementsType(Parameters, Node, Element);
	
	/* **************************************************************************************************************************** */
	
/*	printf("\n --- Conectividade dos 20 primeiros elementos --- \n");
	for(I = 0; I < 20; I++){
		printf("%d: %d\t%d\t%d\t%d\n", I, Element[I].Vertex[0], Element[I].Vertex[1], Element[I].Vertex[2], Element[I].Vertex[3]);
	}
	getchar();   */
	
/*	printf("\nConectividade\n");
	int j;
	for(I = 0; I < nel; I++){
		if (I == 35 || I == 147 || I == 191 || I == 222){
			printf("%d: \t %d \t %d \t %d \t %d \n", I, Element[I].Vertex[0], Element[I].Vertex[1], Element[I].Vertex[2], Element[I].Vertex[3]);
			printf("\nCoord dos nodes\n");
			for (j = 0; j < 4; j++){
				printf("%lf\t%lf\t%lf\t%d\n", Node[Element[I].Vertex[j]].x, Node[Element[I].Vertex[j]].y, Node[Element[I].Vertex[j]].z, Node[Element[I].Vertex[j]].id);
			}
			printf("\n");
		}
	}
	getchar(); */

	/* **************************************************************************************************************************** */
	//			          Memory allocations and Store strategies 
	/* **************************************************************************************************************************** */

	// Some variable initializations
	
	Parameters->neq = neq;
	Parameters->nel = nel;
	Parameters->nnodes = nnodes; 
	
	F = (double*) mycalloc("F of 'Preprocess'", neq+1, sizeof(double));
	u = (double*) mycalloc("u of 'Preprocess'", neq+1, sizeof(double));

	// Structure vector allocation that stores form function data in each element
	CFF = (CoefFormFuncType *) mycalloc("CFF of 'Preprocess'", nel, sizeof(CoefFormFuncType));
	FemStructs->CFF = CFF;

	// Configuring equation according to variables and boundary conditions
	lm = (int**) mycalloc("lm of 'Preprocess'", nel, sizeof(int*));
	lmaux = (int*) mycalloc("lmaux of 'Preprocess'", nel*size, sizeof(int));
	for (I = 0; I < nel; I++){
		lm[I] = &lmaux[I*size];
	}
	
	Fill_LM(neq, nel, lm, Node, Element);
	
/*	printf("\n --- Matriz LM dos 20 primeiros elementos --- \n");
	for(I = 0; I < 20; I++){
		printf("%d: \t %d \t %d \t %d \t %d \n", I, lm[I][0], lm[I][1], lm[I][2], lm[I][3]);
	}
	getchar(); 	*/
	
	FemStructs->lm = lm;
	FemStructs->lmaux = lmaux;
			
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3) == 0){
	
		double **M, *Maux;
		//double *Diag, *invDiag;

		M = (double**) mycalloc("M of 'Preprocess'", nel, sizeof(double*));
		Maux = (double*) mycalloc("Maux of 'Preprocess'", nel*size2,sizeof(double));
		//Diag = (double*) mycalloc("Diag of 'Preprocess'", neq+1, sizeof(double));
		//invDiag = (double*) mycalloc("invDiag of 'Preprocess'", neq+1, sizeof(double));
		
		for (I = 0; I < nel; I++){
			M[I] = &Maux[I*size2];
		}

		MatrixData->A = M;
		MatrixData->Aaux = Maux;
		//MatrixData->Diag = Diag;
		//MatrixData->invDiag = invDiag;
		
		//Para montar a matriz global 
	/*	double **N, *Naux;

		N = (double**) mycalloc("N of 'Preprocess'", neq, sizeof(double*));
		Naux = (double*) mycalloc("Naux of 'Preprocess'", neq*neq,sizeof(double));
		
		dzero(neq*neq,Naux);//tendo certeza que N está zerado
		
		for (I = 0; I < neq; I++){
			N[I] = &Naux[I*neq];
		}

		MatrixData->G = N;
		MatrixData->Gaux = Naux;*/
			
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR") == 0){
		int *IA, *JA, *perm, *invperm, **CSR_by_Element;
		double *AA, *Diag, *invDiag;
			
		csr_Initialization(Parameters, Node, &JA, &IA, &perm, &invperm, &lm, &lmaux, &CSR_by_Element);
		
		AA = (double*) mycalloc("AA of 'Preprocess'", Parameters->nnzero+1, sizeof(double));
		Diag = (double*) mycalloc("Diag of 'Preprocess'", neq+1, sizeof(double));
		invDiag = (double*) mycalloc("invDiag of 'Preprocess'", neq+1, sizeof(double));
		
		//printf("nnzero=%d\n",Parameters->nnzero);

		MatrixData->AA = AA;
		MatrixData->JA = JA;
		MatrixData->IA = IA;
		
/*		int i=0;
		printf("AA vector\n");
		for(i = 0; i < Parameters->nnzero; i++){
			printf("%E\t", AA[i]);
		}
		printf("\n\n");
		printf("JA vector\n");
		for(i = 0; i < Parameters->nnzero; i++){
			printf("%d\t", JA[i]);
		}
		printf("\n\n");
		printf("IA vector\n");
		for(i = 0; i < neq+1; i++){
			printf("%d\t", IA[i]);
		}*/
		
		MatrixData->Diag = Diag;
		MatrixData->invDiag = invDiag;
		MatrixData->Scheme_by_Element = CSR_by_Element;
		MatrixData->Perm = perm;
		MatrixData->invPerm = invperm;

	}
	else{
		printf("Matrix vector product scheme not defined!\n\n");
		exit(1);
	}
	
	/* **************************************************************************************************************************** */

	FemStructs->Node = Node;
	FemStructs->Element = Element;
	FemStructs->F = F;
	FemStructs->u = u;

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

