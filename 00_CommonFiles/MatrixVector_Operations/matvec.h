#ifdef SSTranspEquation2D
	#include "../../01_SS_Transport_Equation/2D/01_CommonFiles/SSTranspEquation.h"
	int ebemv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int edemv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
#endif
#ifdef SSTransportEquation3D
	#include "../../01_SS_Transport_Equation/3D/01_CommonFiles/SSTransportEquation3D.h"
#endif
#ifdef TranspEquation2D
	#include "../../02_Transport_Equation/2D/01_CommonFiles/TranspEquation.h"
	int ebemv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int edemv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
#endif
#ifdef EulerEquations2D
	#include "../../03_Euler_Equations/2D/01_CommonFiles/EulerEquations.h"
	int ebemvNDOF4(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int ebe2mvNDOF4(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int edemvNDOF4(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
	int ede2mvNDOF4(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
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
#ifdef SS_StokesEquations3D
	#include "../../06_SS_Stokes_Equations/3D/01_CommonFiles/SS_StokesEquations3D.h"
#endif
#ifdef SS7_StokesEquations3D
	#include "../../07_SS_Stokes_Equations/3D/01_CommonFiles/SS7_StokesEquations3D.h"
#endif
#ifdef PoissonEquation3D
	#include "../../08_Poisson_Equation/3D/01_CommonFiles/PoissonEquation3D.h"
#endif

int ebemv3Dndof1(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
int ebemvNDOF3(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
int ebemvNDOF4(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
int csrmv(ParametersType *, MatrixDataType *, FemStructsType *, double *, double *);
