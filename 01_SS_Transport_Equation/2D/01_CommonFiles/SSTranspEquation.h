#ifndef SSTranspEquation_h
#define SSTranspEquation_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/amg_precond.h"

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
	char h_Shock[200];
	double SolverTolerance, NonLinearTolerance;
	int KrylovBasisVectorsQuantity;
	int nnodes;
	int nel;
	int neq;
	int nedge;
	int nnzero;
	int iterations;
	int gmres;
	int LinearMaxIter;
	int bandwidth_bef, bandwidth_aft;	// half bandwidth before and after reordering
}ParametersType;

typedef struct
{
	double **A;   //A[0] refers to an M matrix and A[1] refers to an K matrix in transient problems
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
    AMG_precond_data *amg_precond_data;
}MatrixDataType;


typedef struct
{
	int **lm;
	int *lmaux;
	int **lm2;
	int *lm2aux;
	NodeType *Node;
	ElementType *Element;
	double *u;
	double *F;
}FemStructsType;

typedef struct
{
	double (*Condutivity)(void);
	double (*Font)(double, double, double, double, double, double);
	double (*Reaction)(void);
	void (*Velocity)(double, double, double []);
	double (*upresc)(double, double);
	double (*h_shock)(double, double, double, double, double, double, double, double, double, double, double, double);
	double (*ShockCapture)(double,  double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

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

double CAU_ShockCapture(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

double CAU_DD_ShockCapture(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

double YZBeta_ShockCapture(double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

double h_shock_2sqrtArea(double, double, double, double, double, double, double, double, double, double, double, double);

double h_shock_Option1(double, double, double, double, double, double, double, double, double, double, double, double);

double h_shock_Option2(double, double, double, double, double, double, double, double, double, double, double, double);

int Fill_LM(int, int, int **, NodeType *, ElementType *);

void ede_Initialization(ParametersType *, int **, int ***, int **, int ***);

void ede_List_insertA(NodeListType **, int, int, int *);

void csr_Initialization(ParametersType *, NodeType *, int **, int **, int **, int  **,	int ***, int **, int ***);

void csr_List_insertA(NodeListType **, int, int, int *);

int csr_search(int, int, NodeListType *);

int Build_K_F_SUPG(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_K_F_DD(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

void csr_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*ke)[3]);

void ebe_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*ke)[3]);

void ede_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*ke)[3]);

void F_assembly(NodeType *, int, int, int, double, double, double, double *, double [3][3],double [3], double [3], FemFunctionsType *);

int Paraview_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions);

int setzeros(ParametersType *, MatrixDataType *);

int setProblem(ParametersType *, FemFunctionsType *);

int setMatrixVectorProductType(ParametersType *, FemFunctionsType *);

int setSolver(ParametersType *, FemOtherFunctionsType *);

int setPreconditioner(ParametersType *, FemFunctionsType *);

int setScaling(ParametersType *Parameters, FemFunctionsType *FemFunctions);

int setStabilizationForm(ParametersType *, FemFunctionsType *, FemOtherFunctionsType *);

void eval_U(ParametersType *,FemStructsType *, FemFunctionsType *, double *);

//tools
int matrixData2octave(ParametersType *, MatrixDataType *, FemStructsType *);

#endif
