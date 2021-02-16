#ifndef EulerEquations_h
#define EulerEquations_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/ilup.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/amg_precond.h"

#define NDIM 2            // number related to the dimension of the problem (2: two-dimensional, 3: three-dimensional)
#define NNOEL 3           // number of nodes per element
#define NDOF 4            // number of degrees of freedom
#define PI 3.14159265359

typedef struct
{
	double x,y;    // coordinates of the node
	int rhoType, v1Type, v2Type, eType;           // mark of the node according to boundary conditions
	int id[NDOF];  // vector that identifies whether the node is prescribed or not for each property
                         // 0: density
                         // 1: velocity in the x direction
                         // 2: velocity in the y direction
                         // 3: energy
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
	char ProblemTitle[300];                // Problem name to be solved
	char Solver[200];                      // type of method used in the solution
	char TimeIntegration[200];             // Time integration method
	char MatrixVectorProductScheme[200];   // the global matrix storage form
	char StabilizationForm[200];           // type of stabilization method
	char ShockCapture[200];            // type discontinuities capture Operator
	char Preconditioner[1000];           // preconditioners: yes - use or not - don't use
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
}ParametersType;

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
}MatrixDataType;

typedef struct{
	double **M2, **R2, *invN2, *delta_old_NMV;	
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
	double *uB;
	double *duB;
	AuxBuildStructuresType *AuxBuild;
}FemStructsType;

typedef struct
{
	double (*gamma)(double, double);
	double (*cv)(double, double);
	double (*rhopresc)(double, double);
	double (*v1presc)(double, double);
	double (*v2presc)(double, double);
	double (*epresc)(double, double);
	int (*InitialSolution)(ParametersType *, NodeType *, double *);

	double (*ShockCapture)(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4],
					double *, double, double, double, double, double, double, double, int, double *, double*);
	void (*Ax_Ay_calculations)(double, double, double [4], double [4][4], double [4][4]);
	double (*BC_theta)(double, double);
	void (*BC_no_penetrability)(int , int, int, NodeType *, double, double, double, double, double, double,
					    double, double, double, double, double [4][4], double [4][4], double [12][12], double [12][4],
						double [4][12], double [12], double [4], double [12], double [12], double [4], double [4]);
	void (*BC_General_theta)(int, int, int, NodeType *, double [3], double (*)(double, double));
	int (*StopCriteria)(ParametersType *, double, double, int);
	int (*StopTimeIntegration)(ParametersType *, double *, double *, double);

	void (*assembly)(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[12]);
	int (*mv)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precondR)(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int (*precond_setup)(ParametersType *, MatrixDataType *, FemStructsType *, int, double *);
	int (*scaling)(ParametersType *, MatrixDataType *, FemStructsType *);
	int (*unscaling)(ParametersType *, MatrixDataType *, FemStructsType *, double *);
}FemFunctionsType;

typedef struct
{
	int (*solver) (ParametersType *, MatrixDataType *, FemStructsType*, FemFunctionsType*, double *, double *);
	int (*Build)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);
}FemOtherFunctionsType;


int Preprocess(int, char **, ParametersType **, MatrixDataType **, FemStructsType **, FemFunctionsType **, FemOtherFunctionsType **);

int Process(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Postprocess(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Fill_ID(int *, int *, NodeType *, int);

int Fill_LM(int, int, int **, NodeType *, ElementType *);

int setProblem(ParametersType *, FemFunctionsType *);

int setMatrixVectorProductType(ParametersType *, FemFunctionsType *);

int setSolver(ParametersType *, FemOtherFunctionsType *);

int setPreconditioner(ParametersType *, FemFunctionsType *);

int setScaling(ParametersType *, FemFunctionsType *);

int setStabilizationForm(ParametersType *,FemFunctionsType *, FemOtherFunctionsType *, int (**)(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *));

int setzeros(ParametersType *, MatrixDataType *);

double Delta_CAU(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4],
				double *, double, double, double, double, double, double, double, int, double *, double*);
				
double Delta_YZBeta(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4], 
				double *, double, double, double, double, double, double, double, int, double *, double *);

double Delta_YZBetaNMV(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4], 
				double *, double, double, double, double, double, double, double, int, double *, double *);

double Delta_DD(double, double *, double *, double *, double (*)[4], double (*)[4], double (*)[4],
				double *, double, double, double, double, double, double, double, int, double *, double *);

void NO_BC_no_penetrability(int, int, int, NodeType *, double, double, double, double, double, double, double, double, double, double, 
				    double [4][4], double [4][4], double [12][12], double [12][4], double [4][12],
		         		double [12], double [4], double [12], double [12], double [4], double [4]);

void no_penetrability(double, double, double, double, double, double, double, double, double, double, double, double, 
	                 double, double, double, double, double [4][4], double [4][4], double [12][12], double [12][4], double [4][12],
		         double [12], double [4], double [12], double [12], double [4], double [4]);

void set_BC_no_penetrability(ParametersType *, FemFunctionsType *);

void csr_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[12]);

void csr_Initialization(ParametersType *, NodeType *, int **, int **, int **, int  **, int ***, int **, int ***);

void csr_List_insertA(NodeListType **, int , int , int *);

int csr_search(int, int, NodeListType *);

void ebe_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[12]);

void ede_assembly(ParametersType *, MatrixDataType *, FemStructsType *, int, double (*)[12]);

void ede_Initialization(ParametersType *, NodeType *, ElementType *, int **, int ***, int **, int ***);

void edge_JI_assembly(int, int, int, double (*)[12], double **);

void ede_List_insertA(NodeListType **, int, int, int *);

int Build_M_K_F_SUPG(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_M_K_F_DD_Transiente(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

int Build_M_F_DD_Transiente(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *);

void eval_U_dU(ParametersType *,FemStructsType *, FemFunctionsType *, double *,double *);

void dimensionless_Ax_Ay_calculations(double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4]);

void dimensional_Ax_Ay_calculations(double gamma, double Mach, double Ub[4], double Ax[4][4], double Ay[4][4]);

void setDimensionlessness(ParametersType *, FemFunctionsType *);

int Predictor_Old(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_Old_BDF(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New_BDF(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New_TRBDF2(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_Old_TRBDF2(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Paraview_Output(ParametersType *, FemStructsType *, FemFunctionsType *);

int Paraview_Output_3D(ParametersType *, FemStructsType *, FemFunctionsType *);

int calculate_DaB(ParametersType *, FemStructsType *, FemFunctionsType *, double *, double *);

int uB_InitialSolution(ParametersType *, FemStructsType *, FemFunctionsType *, double *, double *);

void print_rho_Residue(FILE *outFile, int neqrho, double t,double *R, double *R_rho, int *eqrho);

void BC_theta_OK(int, int, int, NodeType *, double [3], double (*BC_theta)(double, double));

void BC_theta_NO(int, int, int, NodeType *, double [3], double (*BC_theta)(double, double));

void rotation(int tag, double theta, double M[12][12], double F[12]);

#endif

