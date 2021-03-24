#ifndef _time_integration_h_
#define _time_integration_h_

#ifdef TranspEquation2D
	#include "../../02_Transport_Equation/2D/01_CommonFiles/TranspEquation.h"
#endif
#ifdef SSTransportEquation3D
	#include "../../01_SS_Transport_Equation/3D/01_CommonFiles/SSTransportEquation3D.h"
#endif
#ifdef EulerEquations2D
	#include "../../03_Euler_Equations/2D/01_CommonFiles/EulerEquations.h"
#endif
#ifdef SSNavierStokesEquations2D
	#include "../../04_SS_NavierStokes_Equations/2D/01_CommonFiles/SSNavierStokesEquations.h"
#endif
#ifdef SSNavierStokesEquations3D
	#include "../../04_SS_NavierStokes_Equations/3D/01_CommonFiles/SSNavierStokesEquations3D.h"
#endif
#ifdef NavierStokesEquations2D
	#include "../../05_NavierStokes_Equations/2D/01_CommonFiles/NavierStokesEquations.h"
#endif
#ifdef NavierStokesEquations3D
	#include "../../05_NavierStokes_Equations/3D/01_CommonFiles/NavierStokesEquations3D.h"
#endif
#ifdef ShalowWater
	#include "../../06_ShallowWater/2D/01_CommonFiles/ShalowWater.h"
#endif
#include "../BLAS_Operations/ourBLAS.h"
#include "../Allocation_Operations/allocations.h"
#include "../IO_Operations/io.h"
#include "../C_Operations/c.h"

int Predictor_Old(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_Old_BDF(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New_BDF(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_New_TRBDF2(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

int Predictor_Old_TRBDF2(ParametersType *, MatrixDataType *, FemStructsType *, FemFunctionsType *, FemOtherFunctionsType *);

void setStopCriteria(ParametersType *, FemFunctionsType *);

int StopByIterations(ParametersType *, double, double, int);

int StopByNorm(ParametersType *, double, double, int);

int StopBySteadyState(ParametersType *, double *, double *, double);

int StopByTime(ParametersType *, double *, double *, double);

#endif
