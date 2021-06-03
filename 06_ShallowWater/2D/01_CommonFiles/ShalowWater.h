#ifndef _ShalowWater_h_
#define _ShalowWater_h_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/amg_precond.h"


#define NDIM 2            // number related to the dimension of the problem (2: two-dimensional, 3: three-dimensional)
#define NNOEL 3           // number of nodes per element
#define NDOF 3            // number of degrees of freedom
#define PI 3.14159265359
#define TOL 1e-12         // tolerance

typedef struct
{
	double x,y;    // coordinates of the node
	int hType, qxType, qyType;           // mark of the node according to boundary conditions
	int id[NDOF];  // vector that identifies whether the node is prescribed or not for each property
                   // 0: h (height)
                   // 1: qx (discharge in x axis)
                   // 2: qy (discharge in y axis)
} NodeType;             

typedef struct
{
	int Vertex[NNOEL];  // holds the global number of nodes that make up the element
	int Type;           // mark of the element
}ElementType;

struct Node_List{
	double value;
	int J;  // vertice representando a cabeca do arco
	int count; // counter to represent the number of a especific edge
	struct Node_List *next;
};

typedef struct Node_List NodeListType;

typedef struct
{
	char outPath[200];
	char ProblemTitle[300];                // Problem name to be solved
	char Experiments[200];
	int Friction;						   // 1 friction, 0 frictionless
	double Viscosity;
	char Solver[200];                      // type of method used in the solution
	char TimeIntegration[200];             // Time integration method
	char MatrixVectorProductScheme[200];   // The global matrix storage form
	char StabilizationForm[200];           // Type of stabilization method
	char ShockCapture[200];                // Type discontinuities capture Operator
	char Preconditioner[200];              // Preconditioners: yes - use or not - don't use
	char StopMulticorrection[200];         // Fala se o loop espacial para pela norma ou por um numero fixo de iteracao - NORM: para pela norma; ITERATION: para pela iteracao
	char reordering[200];			       // Reordering for CSR (NOT, Spectral, Weigthed Spectral (WSO) or RCM)
	char StopAtSteadyState[200];		   // YES or NO if you want to stop in steady state or final time
	double SolverTolerance;                // Tolerance for the solution method
	double NonLinearTolerance;             // Tolerance for the loop of correction
	double TimeIntegrationTolerance;       // Tolerance for the loop of time integration
	double Alpha;                          // Parameter determining the stability control and accuracy time integration method
	double Alpha_Build;					   // Parameter Alpha to be used inside Build functions
	double StabilizationTolerance;         // Tolerance for coefficient used in the stabilization form
	double DeltaT;                         // Time step
	double DeltaT_Build;				   // Time step to be used inside Build functions
	double FinalTime;                      // Tinal time
	double CurrentTime;                    // Current time (to be used in steady state situation)
	int KrylovBasisVectorsQuantity;        // Krylov number of vectors in the basis for the restart
	int nnodes;                            // Number of nodes
	int nel;                               // Number of element
	int neq;                               // Number of equations to be solved
	int nedge;
	int nnzero;
	int iterations;                        // Total number of iteration 
	int LinearMaxIter;                     // Maximum number of iteration
	int NonLinearMaxIter;                  // Maximum nonlinear iterations: number of multicorrection
	int bandwidth_bef, bandwidth_aft;	   // Half bandwidth before and after reordering
	int SolverIterations;
	int SolverMaxIter;
	int ResGMRES;
	int ContGMRES;
} ParametersType;

typedef struct
{
	double **A;
	double *Aaux;
	double **BlockDiag;
	double *BlockDiagAux;
	double **invBlockDiag;
	double *invBlockDiagAux;
	double **invDe;
	double *invDeaux;
	double **LUe;
	double *LUeaux;
	int **Id;
	int *IdAux;
	int **Scheme_by_Element;
	int **order;
	double *AA;
	double *Diag;
	double *invDiag;
	int *JA;
	int *IA;
	SparILU *ILUp;
	SparMAT *mat;
	MAT *Ailu;
	int *Perm;
	int *invPerm;
	AMG_precond_data *amg_precond_data;
	double **G;
} MatrixDataType;

typedef struct{
	double **M2, **R2, *invN2, *delta_old_NMV;	
	double tolerance;
} AuxBuildStructuresType;

typedef struct
{
	int **lm;
	int *lmaux;
	int **lm2; // Used only in EDE mode
	int *lm2aux; // Used only in EDE mode
	int *eqrho;
	NodeType *Node;
	ElementType *Element;
	double *F;
	double *R;
	double *u;
	double *du;
	double *delta_old;
	double *uB;
	double *duB;
	AuxBuildStructuresType *AuxBuild;
} FemStructsType;

typedef struct
{
	double (*hpresc)(double, double);
	double (*qxpresc)(double, double);
	double (*qypresc)(double, double);
	double (*zb)(double, double);
	int (*InitialSolution)(ParametersType *, FemStructsType *);

	double (*ShockCapture)(double *, double *, double *, double *, double *,
                           double (*)[3], double (*)[3],
                           double, double, double, double, double, double, double,
                           double, double, double, double, double);
	int (*StopCriteria)(ParametersType *, double, double, int);
	int (*StopTimeIntegration)(ParametersType *, double *, double *, double);

	void (*assembly)(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[9]);
	int (*mv)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precondR)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond_setup)(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);

	int (*scaling)(ParametersType *, MatrixDataType *, FemStructsType *);
	int (*unscaling)(ParametersType *, MatrixDataType *, FemStructsType *, double *);

	int (*ProductMatrixVector)(ParametersType *, MatrixDataType *, FemStructsType *, double*, double*);
} FemFunctionsType;

typedef struct
{
	int (*solver) (ParametersType *, MatrixDataType *, FemStructsType*, FemFunctionsType*, double *, double *);
	int (*Build)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);
} FemOtherFunctionsType;


void A1_A2_calculations(double Ub[3], double (*A1)[3], double (*A2)[3], double g);

int Build_M_D_F_SUPG(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions);

void csr_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*Me)[9]);

void ebe_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*Me)[9]);

void csr_Initialization(ParametersType *Parameters, NodeType *Node, int **JA_out, int **IA_out, int **perm_out, int  **invperm_out,
			int ***lm_out, int **lmaux_out, int ***CSR_by_Element_out);

double deltaCAU(double *Ub, double *Sb, double *gradEta, double *gradUx, double *gradUy,
                double (*A1)[3], double (*A2)[3],
                double y23, double y31, double y12, double x32, double x13, double x21, double twoArea,
                double invhRef, double invqxRef, double invqyRef, double tau, double g);

double deltaYZB(double *Ub, double *Sb, double *gradEta, double *gradUx, double *gradUy,
                double (*A1)[3], double (*A2)[3],
                double y23, double y31, double y12, double x32, double x13, double x21, double twoArea,
                double invhRef, double invqxRef, double invqyRef, double tau, double g);

void eval_U_dU(ParametersType *Parameters,FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *U,double *dU);

int Fill_ID(int *neq_out, NodeType *Node, int nnodes);

int Fill_LM(int neq, int nel, int **lm, NodeType *Node, ElementType *Element);

int Paraview_Output_3D_DeltaT(double *U, NodeType *Node, ElementType *Element, ParametersType *Parameters, 
                              double (*hpresc)(double, double), double (*qxpresc)(double, double), double (*qypresc)(double, double), double t);

int Paraview_Output_3D(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions);

int Paraview_Output_DeltaT(double *U, NodeType *Node, ElementType *Element, ParametersType *Parameters, 
                           double (*hpresc)(double, double), double (*qxpresc)(double, double), double (*qypresc)(double, double), double t);

int Paraview_Output(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions);

int Postprocess(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions);

int Preprocess(int narg, char **arguments, ParametersType **Parameters_out, MatrixDataType **MatrixData_out, FemStructsType **FemStructs_out, FemFunctionsType **FemFunctions_out, FemOtherFunctionsType **FemOtherFunctions_out);

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions);

int setMatrixVectorProductType(ParametersType *Parameters, FemFunctionsType *FemFunctions);

int setPreconditioner(ParametersType *Parameters, FemFunctionsType *FemFunctions);

int setProblem(ParametersType *Parameters, FemFunctionsType *FemFunctions);

int setSolver(ParametersType *Parameters, FemOtherFunctionsType *FemOtherFunctions);

int setStabilizationForm(ParametersType *Parameters,FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions);

int setzeros(ParametersType *Parameters, MatrixDataType *MatrixData);

void F_assembly(int e, double *Fe, double (*De)[9], FemFunctionsType *FemFunctions, FemStructsType *FemStructs, int neq);

int PredictorMulticorrector(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
		FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions);

void setStopCriteria(ParametersType *, FemFunctionsType *);

int StopByIterations(ParametersType *, double, double, int);

int StopByNorm(ParametersType *, double, double, int);

int StopBySteadyState(ParametersType *, double *, double *, double);

int StopByTime(ParametersType *, double *, double *, double);


void printU(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions);

void checknull(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions);

int AssemblyGlobalMatrix(MatrixDataType *MatrixData, FemStructsType *FemStructs, int e, double (*Me)[9], int neq);

void print_CSR(MatrixDataType *MatrixData, ParametersType *Parameters);

void print_EBE(MatrixDataType *MatrixData, ParametersType *Parameters, int E);

double analitical(double x, double y, ParametersType *Parameters);

#endif
