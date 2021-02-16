#ifndef SSNavierStokesEquations_h
#define SSNavierStokesEquations_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/amg_precond.h"


#define NDIM 2            // number related to the dimension of the problem (2: two-dimensional, 3: three-dimensional)
#define NNOEL 3           // number of nodes per element
#define NDOF 3            // number of degrees of freedom
#define PI 3.14159265359

typedef struct
{
	double x,y;    // coordinates of the node
	int v1Type, v2Type, pType; // mark of the node according to boundary conditions
	int id[NDOF];  // vector that identifies whether the node is prescribed or not for each property
                         // 0: velocity in the x direction
                         // 1: velocity in the y direction
                         // 2: pressure
}NodeType;             

typedef struct
{
	int Vertex[NNOEL];  // holds the global number of nodes that make up the element
	int Type;           // mark of the element
}ElementType;

struct Node_List{
	double value;
	int J; // vertice representando a cabeca do arco
	int count; // counter to represent the number of a especific edge
	struct Node_List *next;
};

typedef struct Node_List NodeListType;

typedef struct
{
	char ProblemTitle[250];                // Problem name to be solved
	char Solver[200];                      // type of method used in the solution
	char TimeIntegration[200];             // Time integration method
	char MatrixVectorProductScheme[200];   // the global matrix storage form
	char StabilizationForm[200];           // type of stabilization method
	char ShockCapture[200];            // type discontinuities capture Operator
	char Preconditioner[200];           // precondicionators: yes - use or not - don't use
	char Experiments[200];
	char StopMulticorrection[200];         // Fala se o loop espacial para pela norma ou por um numero fixo de iteração - NORM: para pela norma; ITERATION: para pela iteração
	double ReynoldsNumber;		       // Reynolds Number
	double VelMax;			       //Velocidade maxima para o computo de Re local
	char reordering[200];		       // Reordering for CSR (NOT, Spectral, Weigthed Spectral (WSO) or RCM)
	double SolverTolerance;             // tolerance for the Linear solution method 
	int SolverToleranceCase;               // Case tolerance for the solution method 1(etapp), 2(etaewk), 3(etaewc), 4(etaglt)
	double SolverToleranceNonLin;          // tolerance for the Nonlinear solution method
	double CorrectionTolerance;            // tolerance for the loop of correction
	double TimeIntegrationTolerance;       // tolerance for the loop of time integration
	double CoefficientTolerance;           // tolerance for coefficient used in the stabilization form
	double Alpha;                          // parameter determining the stability control and accuracy time integration method
	double DeltaT;                         // time step
	double FinalT;                         // final time
	double invY[4];                        // vector that stores 1/U used in the capture operator YZBeta
	double normu, normv, normp;            // morm solutions vectors
	int KrylovBasisVectorsQuantity;        // Krylov number of vectors in the basis for the restart
	int nnodes;                            // nnodes: number of nodes
	int nel;                               // nel: number of element
	int neq;                               // neq: number of equations to be solved
	int nedge;
	int nnzero;
	int iterations;                        // iterations: total number of iteration 
	int gmres;
	int NLiterations;                        // NLiterations: total number of nonlinear iteration 
	int LinearMaxIter;                           // itermax: maximum number of iteration
	int NumberCorrection;                  // NumberCorrection: number of multicorrection
	int bandwidth_bef, bandwidth_aft;	// half bandwidth before and after reordering
}ParametersType;

typedef struct
{
	double **A;     
	double *Aaux;    
	double *AA;
	int *JA;
	int *IA;
	double *D;
	double *Diag;
	double *invDiag;
	double **BlockDiag;
	double *BlockDiagAux;
	double **invBlockDiag;
	double *invBlockDiagAux;
	int **Id;
	int *IdAux;
	int **Scheme_by_Element;
	int **order;
	SparILU *ILUp;
	SparMAT *mat;
	MAT *Ailu;
	int *Perm;
	int *invPerm;
	AMG_precond_data *amg_precond_data;
}MatrixDataType;

typedef struct{
	double **M2, **R1, **R2, *invN2;	
	double tolerance;
}AuxBuildStructuresType;

typedef struct
{
	int **lm;
	int *lmaux;
	int **lm2;	//Used only in EDE mode
	int *lm2aux;	//Used only in EDE mode
	int *eqrho;
	NodeType *Node;
	ElementType *Element;
	double *F;
	double *u;
	double *du;
	double *delta_old;
	double *normres_old;
	//double *deltam_old;
	double *uB;
	double *duB;
	AuxBuildStructuresType *AuxBuild;
}FemStructsType;

typedef struct
{
	double (*v1presc)(double, double);
	double (*v2presc)(double, double);
	double (*ppresc)(double, double);
	double (*f1ext)(double, double);
	double (*f2ext)(double, double);

	void (*assembly)(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[9]);
	//!!!!!!!!!!!!! VERIFICAR POIS ESTA DIFERENTE !!!!!!!!!
	int (*mv)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precondR)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond_setup)(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
}FemFunctionsType;

typedef struct
{
	int (*solver) (ParametersType *, MatrixDataType *, FemStructsType*, FemFunctionsType*, double *, double *);
	int (*Build)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);
}FemOtherFunctionsType;


int Preprocess(int, char **, ParametersType **, MatrixDataType **, FemStructsType **, FemFunctionsType **, FemOtherFunctionsType **); //ok

int Process(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *); //ok

int Postprocess(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *); //ok

int Fill_ID(int *, NodeType *, int); //ok

int Fill_LM(int, int, int **, NodeType *, ElementType *); //ok

int setProblem(ParametersType *, FemFunctionsType *); //ok

int setMatrixVectorProductType(ParametersType *, FemFunctionsType *);

int setSolver(ParametersType *, FemOtherFunctionsType *); //ok

int setPreconditioner(ParametersType *, FemFunctionsType *); //ok

int setStabilizationForm(ParametersType *,FemFunctionsType *, FemOtherFunctionsType *); //ok

void csr_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[9]); //ok

void csr_Initialization(ParametersType *, NodeType *, int **, int **, int **, int  **, int ***, int **, int ***); //ok

void csr_List_insertA(NodeListType **, int , int , int *);//ok

int csr_search(int, int, NodeListType *);//ok

void ebe_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[9]); //ok

int Build_K_F_SUPG(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *); //ok

int Build_K_F_MS(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *); //ok

int Build_K_F_MS_Time(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *); //ok

int Paraview_Output(ParametersType *, FemStructsType *, FemFunctionsType *); //ok

double eta_newton(double *, double *, double , int , int , double, double, ParametersType *); //ok

int Coef_Lift_Drag(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *); //ok

void eval_U(ParametersType *,FemStructsType *, FemFunctionsType *, double *);

void eval_U_dU(ParametersType *,FemStructsType *, FemFunctionsType *, double *,double *);

int setzeros(ParametersType *, MatrixDataType *);

#endif

