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
	char ProblemTitle[300];                // Problem name to be solved
	char Solver[200];                      // type of method used in the solution
	char TimeIntegration[200];             // Time integration method
	char MatrixVectorProductScheme[200];   // the global matrix storage form
	char StabilizationForm[200];           // type of stabilization method
	char ShockCapture[200];            // type discontinuities capture Operator
	char Preconditioner[200];           // preconditioners: yes - use or not - don't use
	char Scaling[200];
	char Experiments[200];
	char StopMulticorrection[200];         // Fala se o loop espacial para pela norma ou por um numero fixo de iteracao - NORM: para pela norma; ITERATION: para pela iteracao
	char reordering[200];			// Reordering for CSR (NOT, Spectral, Weigthed Spectral (WSO) or RCM)
	char Dimensionless[200];		// Determines whether or not a problem is dimensionless
	char StopAtSteadyState[200];		// YES or NO if you want to stop in steady state or final time
	double SolverTolerance;                // tolerance for the solution method
	double NonLinearTolerance;            // tolerance for the loop of correction
	double TimeIntegrationTolerance;       // tolerance for the loop of time integration
	double StabilizationTolerance;           // tolerance for coefficient used in the stabilization form
	double Mach;				//Mach number of reference
	double Alpha;                          // parameter determining the stability control and accuracy time integration method
	double Alpha_Build;					// parameter Alpha to be used inside Build functions
	double DeltaT;                         // time step
	double DeltaT_Build;					// time step to be used inside Build functions
	double FinalTime;                         // final time
	double CurrentTime;                         // current time (to be used in steady state situation)
	double invY[4];                        // vector that stores 1/U used in the capture operator YZBeta
	int KrylovBasisVectorsQuantity;        // Krylov number of vectors in the basis for the restart
	int nnodes;                            // nnodes: number of nodes
	int nonpnodes;				//nonpnodes: number of nodes with no penetrability boundary conditions
	int nel;                               // nel: number of element
	int neq;                               // neq: number of equations to be solved
	int neqrho;                            // neqrho: number of equations relating the density
	int nedge;
	int nnzero;
	int iterations;                        // iterations: total number of iteration 
	int LinearMaxIter;                           // itermax: maximum number of iteration
	int NonLinearMaxIter;                  // Maximum nonlinear iterations: number of multicorrection
	int bandwidth_bef, bandwidth_aft;	// half bandwidth before and after reordering
		
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
	int (*InitialSolution)(ParametersType *, NodeType *, double *);

	double (*ShockCapture)(double *, double *, double *, double (*)[3], double (*)[3], double *, 
                           double, double, double, double, double, double, double, double, double);
	void (*A1_A2_calculations)(double *, double (*)[3], double (*)[3], double);
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


void A1_A2_calculations(double Ub[3], double A1[3][3], double A2[3][3], double g);

int Build_M_D_F_SUPG(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions);

void csr_assembly(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int E, double (*Me)[9]);

void csr_Initialization(ParametersType *Parameters, NodeType *Node, int **JA_out, int **IA_out, int **perm_out, int  **invperm_out,
			int ***lm_out, int **lmaux_out, int ***CSR_by_Element_out);

double delta91_MOD(double Ub[3], double gradUx[3], double gradUy[3], double A1[3][3], double A2[3][3], double Sb[3], 
                   double y23, double y31, double y12, double x32, double x13, double x21, double twoArea, 
				   double tau, double g);

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

int setStabilizationForm(ParametersType *Parameters,FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions,
						 int (**Predictor)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *));

int setzeros(ParametersType *Parameters, MatrixDataType *MatrixData);

void F_assembly(int e, double Fe[9], double De[9][9], FemFunctionsType *FemFunctions, FemStructsType *FemStructs, int neq);

#endif
