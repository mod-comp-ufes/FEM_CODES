#ifndef SSTransportEquation3D_h
#define SSTransportEquation3D_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"

#define NDIM 3            // number related to the dimension of the problem (2: two-dimensional, 3: three-dimensional)
#define NNOEL 4           // number of nodes per element
#define NDOF 1            // number of degrees of freedom
#define PI 3.141592654
//#define PI M_PI           // M_PI is the constant pi in the library math.h

typedef struct
{
	double x,y,z;	// coordinates of the node
	int uType; 		// mark of the node according to boundary conditions
	int id;         // identifies whether the node is prescribed or not for each property
}NodeType;             

typedef struct
{
	int Vertex[NNOEL];  // holds the global number of nodes that make up the element
	int Type;           // mark of the element
}ElementType;

typedef struct
{
	double volume;                 // Element volume 
	double a[4], b[4], c[4], d[4]; // Component of the form function
									  // 0: node 1
									  // 1: node 2
									  // 2: node 3
									  // 3: node 4
}CoefFormFuncType;

struct Node_List
{
	double value;
	int J;                  // vertice representando a cabeca do arco
	int count;              // counter to represent the number of a especific edge
	struct Node_List *next;
};
typedef struct Node_List NodeListType;

typedef struct
{
	char Experiments[20];				  // help organize experiments by name
	char ProblemTitle[20];                // problem name to be solved
	char TypeMesh[15];			// Type Mesh: STRUCTURED or UNSTRUCTURED
	char OriginMesh[10];		// Mesh origin: GMSH or FREEFEM
	char ExactSolution[10];				  // define if exist exact solution of the problem (YES or NOT)
	char Solver[10];                      // type of method used in the solution
	char MatrixVectorProductScheme[10];    // the global matrix storage form (EBE, CSR, EDE)
	char StabilizationForm[10];           // type of stabilization method
	char Preconditioner[10];              // preconditioners: yes - use or not - don't use
	char Reordering[10];			      // Reordering for CSR (NOT, Spectral, Weigthed Spectral (WSO) or RCM)
	char StopMulticorrection[10];
	char OutputFlow[10];			// Yes if use the flow outlet boundary to compute h
	char TauMacro[10];
	char TauMicro[10];
	char InitialSolution[10];		// NULL for zero and GA for Galerkin
	char UsePeclet[10];
	char ErrorType[15];			// Tipo de erro usado para sair do loop não linear: AE - Erro Absoluto, RE - Erro Relativo	
	char UseDamping[10];			// YES or NOT
	char ComputeResidual[15];		// Equation or Sistem	
	double SolverTolerance;               // tolerance for the solution method
	double NonLinearTolerance;
	double NormL2;
	double NormH1;
	double SemiNormH1;
	double NormEnergy;
	double MaxVolume;		// Volume Maximo
	double wMacro;    //Omega parâmetro de ponderação na estabilização
	double wMicro;
	double tetha;      // Parâmetro somado no denominador do cálculo de Cbtil no DD
	double ConstApli; // parametro A da aplicação 4, varia nos seguintes valores: 1, 5, 10, 25, 50, 100, 500, 1000
			  // parametro Re da aplicação 7, varia 1, 10, 100
			  // parametro Pe da aplicação do Romão
	double Vmax;
	double xInf, xSup, yInf, ySup, zInf, zSup;
	double ResGMRES;     // Resíduo do GMRES
	double tolGradU;	// tolerance for the gradiente u norm
	double MaxPecletLocal;
	double MinPecletLocal;
	double MedPecletLocal;
	double MaxReaDifRelationLocal;
	double MinReaDifRelationLocal;
	double MedReaDifRelationLocal;
	int NumberMesh;		// used to differentiate meshes in the channel problem (1,2,3,4,5,6,7) 
	int KrylovBasisVectorsQuantity;       // Krylov number of vectors in the basis for the restart
	int nnodes;                           // nnodes: number of nodes
	int nel;                              // nel: number of element
	int neq;                              // neq: number of equations to be solved
	int nedge;
	int nnzero;
	int SolverIterations;                 // total number of solver iteration
	int ContGMRES;				// Numero de iterações por loop não linear
	int SolverMaxIter;					  // Maximum number of solver iterations 
	int NonLinearMaxIter;
	int Step;				// guarda o número do passo não linear
	int TipoH;				// escolher 1: h = raiz cúbica de volume ou 2: h = norma de beta
	int size;
	int bandwidth_bef, bandwidth_aft;	// half bandwidth before and after reordering
}ParametersType;

typedef struct
{
	double **G;	// Armazena a matriz global montada
	double *Gaux;	
	double **A;   // all elementary matrix
	double *Aaux;
	double *AA;
	double *Diag;
	double *invDiag;
	double **invDe;
	double *invDeaux;
	double **LUe;
	double *LUeaux;
	int **Scheme_by_Element;
	int **order;
	int *JA;
	int *IA;
	SparILU *ILUp;
	SparMAT *mat;
	MAT *Ailu;
	int *Perm;
	int *invPerm;
}MatrixDataType;

typedef struct
{
	int **Id;	
	int **lm;
	int *lmaux;
	NodeType *Node;
	ElementType *Element;
	CoefFormFuncType *CFF;
	double *F;
	double *u;
	double *uOld;
	double *CbOldMacro;
	double *CbOldMicro;
	double *RES;
	double *RESold;
	double **Corte;
}FemStructsType;

typedef struct
{
	double (*upresc)(double, double, double, double, double);
	double (*f)(double, double, double, double);
	double (*ExactSolution)(double, double, double, double);
	void (*ExactSolutionAllPoints) (NodeType *, double, double *, double);
	void (*kappa) (double *, double *, double *, double, double, double, double);	
	void (*beta) (double *, double *, double *, double, double, double, double, double, double, double);
	double (*sigma) (double, double, double);
	double (*DuDx)(double, double, double, double);
	double (*DuDy)(double, double, double, double);
	double (*DuDz)(double, double, double, double);
	int (*ProductMatrixVector)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precondR)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond_setup)(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);	
	void (*MatrixAssembly)(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[4]);
}FemFunctionsType;

typedef struct
{
	int (*Solver) (ParametersType *, MatrixDataType *, FemStructsType*, FemFunctionsType*, double *, double *);
	int (*Build)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);
}FemOtherFunctionsType;

int Preprocess(int, char **, ParametersType **, MatrixDataType **, FemStructsType **, FemFunctionsType **, FemOtherFunctionsType **);

int Process(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Postprocess(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Fill_ID(int *, NodeType *, int);

int Fill_LM(int, int, int **, NodeType *, ElementType *);

int setProblem(ParametersType *, FemFunctionsType *);

int setMatrixVectorProductType(ParametersType *, FemFunctionsType *);

int setSolver(ParametersType *, FemOtherFunctionsType *);

int setPreconditioner(ParametersType *, FemFunctionsType *);

int setStabilizationForm(ParametersType *, FemOtherFunctionsType *);

int setZeros(ParametersType *, MatrixDataType *);

void Space_Algorithm_NonLinear(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

void Space_Algorithm_Damping(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

void Space_Algorithm_Galerkin(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

void Space_Algorithm_Frozen(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

void Space_Algorithm_VJ(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Build(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_SUPG(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_CAU(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_VMS(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_NSGS(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_DD(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

double Element_Configuration(int, ElementType *, NodeType *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

void ElementsType(ParametersType *, NodeType *, ElementType *);

double LengthMesh(int, ParametersType *, ElementType *, double, double, double, double, double *, double *, double *, double, double *);

void eval_U_Space(ParametersType *,FemStructsType *, FemFunctionsType *, double *);

void ebe_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[4]);

void csr_Initialization(ParametersType *, NodeType *, int **, int **, int **, int **, int ***, int **, int ***);

void csr_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int , double (*)[4]);

int AssemblyGlobalMatrix(MatrixDataType *, FemStructsType *, int, double (*)[4], int);

void font_assembly(int, double *, double (*)[4], FemFunctionsType *, FemStructsType *, int, double, double);

int AproximationInEachNode(ParametersType *, FemStructsType *, FemFunctionsType *, double *);

int Paraview_Output(ParametersType *, FemStructsType *, FemFunctionsType *);

int Paraview_Exact(ParametersType *, FemStructsType *, FemFunctionsType *);

int Gnuplot_Output(ParametersType *, FemStructsType *, FemFunctionsType *);

void Norm_L2 (double *, ParametersType *, FemStructsType *, FemFunctionsType *);

void Norm_H1 (double *, ParametersType *, FemStructsType *, FemFunctionsType *);

double New_Norm_L2 (double *, double, NodeType *, FemFunctionsType *);

void Calculating_Errors(ParametersType *, FemStructsType *, FemFunctionsType *);

double Compute_Residual(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

double TauMacro(int e, ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *Ue);

double TauMicro(int e, ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *Ue);

#endif

