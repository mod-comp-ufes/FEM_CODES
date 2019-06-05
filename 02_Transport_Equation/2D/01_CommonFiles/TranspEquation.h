#ifndef TranspEquation_h
#define TranspEquation_h
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"

#define NDIM 2            // number related to the dimension of the problem (2: two-dimensional, 3: three-dimensional)
#define NNOEL 3           // number of nodes per element
#define NDOF 1            // number of degrees of freedom
#define PI 3.14159265359

typedef struct
{
	double x,y;
	int Type;
	int id;
}NodeType;

typedef struct
{
	int Vertex[NDIM+1];
	int Type;
}ElementType;

struct Node_List{
	int J; // vertice representando a cabeca do arco
	struct Node_List *next;
};

typedef struct Node_List NodeListType;


typedef struct
{
	char ProblemTitle[300];	
	char Solver[200];
	char Preconditioner[200];
	char Scaling[200];
	char reordering[200];
	char MatrixVectorProductScheme[200];
	char StabilizationForm[200];
	char ShockCapture[200];
	char TimeIntegration[200];
	char StopMulticorrection[200];
	char h_Shock[200];
	char StopAtSteadyState[10];
	double SolverTolerance, NonLinearTolerance, GradientTolerance, StabilizationTolerance, 
	       DeltaT, DeltaT_Build, Alpha, Alpha_Build, FinalTime, CurrentTime, TimeIntegrationTolerance;
	int KrylovBasisVectorsQuantity;
	int nnodes;
	int nel;
	int neq;
	int nedge;
	int nnzero;
	int iterations;
	int LinearMaxIter;
	int NonLinearMaxIter;
	int NumberCorrection;
	int bandwidth_bef, bandwidth_aft;	// half bandwidth before and after reordering
}ParametersType;

typedef struct
{
	double **A;   
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
}MatrixDataType;

typedef struct{
	double **M2, **R2, *invN2, *delta_old_NMV;	
	double tolerance;
}AuxBuildStructuresType;

typedef struct
{
	int **lm;
	int *lmaux;
	int **lm2;
	int *lm2aux;
	NodeType *Node;
	ElementType *Element;
	double *u;
	double *du;
	double *F;
	double *uB;
	double *duB;
	double *delta_old;
	AuxBuildStructuresType *AuxBuild;
}FemStructsType;

typedef struct
{	
	double (*Condutivity)(void);
	double (*Font)(double, double, double, double, double, double); 
	double (*Reaction)(void);
	void (*Velocity)(double, double, double []);
	double (*upresc)(double, double);
	int (*InitialSolution)(ParametersType *, NodeType*, double *);
	double (*h_shock)(double, double, double, double, double, double, double, double, double, double, double, double);
	double (*ShockCapture)(double,  double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);
	int (*StopCriteria)(ParametersType *, double, double, int);
	int (*StopTimeIntegration)(ParametersType *, double *, double *, double);

	void (*assembly)(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[3]);
	int (*mv)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precondR)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond_setup)(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
	int (*scaling)(ParametersType *, MatrixDataType *, FemStructsType *);
	int (*unscaling)(ParametersType *, MatrixDataType *, FemStructsType *, double *);
}FemFunctionsType;

typedef struct
{
	int (*solver) (ParametersType *, MatrixDataType *, FemStructsType*, FemFunctionsType *, double *, double *);
	int (*Build)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);
}FemOtherFunctionsType;

int Preprocess(int, char **, ParametersType **, MatrixDataType **, FemStructsType **, FemFunctionsType **, FemOtherFunctionsType **);

int Process(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Postprocess(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Paraview_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions);

double CAU_ShockCapture(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double,	double, double);

double CAU_DD_ShockCapture(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

double YZBeta_ShockCapture(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

double h_shock_2sqrtArea(double, double, double, double, double, double, double, double, double, double, double, double);

double h_shock_Option1(double, double, double, double, double, double, double, double, double, double, double, double);

double h_shock_Option2(double, double, double, double, double, double, double, double, double, double, double, double);

int Fill_LM(int, int, int **, NodeType *, ElementType *);

void ede_Initialization(ParametersType *, int **, int ***, int **, int ***);

void ede_List_insertA(NodeListType **, int, int, int *);

void csr_Initialization(ParametersType *, NodeType *, int **, int **, int **, int  **,	int ***, int **, int ***);

void csr_List_insertA(NodeListType **, int, int, int *);

int csr_search(int, int, NodeListType *);

int Build_M_K_R_SUPG(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_M_K_R_DD(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_M_F_DD(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int calculate_DaB(ParametersType *, FemStructsType *, FemFunctionsType *, double *, double *);

int uB_calculation(int, double *, double *, double (*)(double, double), NodeType *, ElementType *);

void csr_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[3]);

void ebe_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[3]);

void ede_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[3]);

int setzeros(ParametersType *, MatrixDataType *);

int setProblem(ParametersType *, FemFunctionsType *);

int setMatrixVectorProductType(ParametersType *, FemFunctionsType *);

int setSolver(ParametersType *, FemOtherFunctionsType *);

int setPreconditioner(ParametersType *, FemFunctionsType *);

int setScaling(ParametersType *, FemFunctionsType *);

int setStabilizationForm(ParametersType *, FemFunctionsType *, FemOtherFunctionsType *, 
			int (**)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *));

void eval_U_dU(ParametersType *,FemStructsType *, FemFunctionsType *FemFunctions, double *,double *);

int uB_InitialSolution(ParametersType *, FemStructsType *, FemFunctionsType *, double *, double *);


#endif

