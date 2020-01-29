#ifdef SSTranspEquation2D
	#include "../../01_SS_Transport_Equation/2D/01_CommonFiles/SSTranspEquation.h"
#endif
#ifdef TranspEquation2D
	#include "../../02_Transport_Equation/2D/01_CommonFiles/TranspEquation.h"
#endif
#ifdef EulerEquations2D
	#include "../../03_Euler_Equations/2D/01_CommonFiles/EulerEquations.h"
#endif
#ifdef SSNavierStokesEquations2D
	#include "../../04_SS_NavierStokes_Equations/2D/01_CommonFiles/SSNavierStokesEquations.h"
#endif
#ifdef NavierStokesEquations2D
	#include "../../05_NavierStokes_Equations/2D/01_CommonFiles/NavierStokesEquations.h"
#endif
#include "../BLAS_Operations/ourBLAS.h"
#include "../Allocation_Operations/allocations.h"
#include "amg_precond.h"

int pgmres (ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, double *, double *);

//int pcg (ParametersType *Parameters, MatrixDataType *MatrixData, double *B, double *X, int **lm, int (*precond)(ParametersType *, MatrixDataType *, double *, double *),
	// int (*mv)(ParametersType *, MatrixDataType *, int , double *, double *, int **));

