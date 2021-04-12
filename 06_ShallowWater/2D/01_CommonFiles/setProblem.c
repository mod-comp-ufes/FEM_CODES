#include "ShalowWater.h"

extern double hpresc(double, double);
extern double qxpresc(double, double);
extern double qypresc(double, double);
extern double zb(double, double);
extern int InitialSolution(ParametersType *, FemStructsType *);


int setProblem(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	#ifndef onlyCompile

		FemFunctions->hpresc = hpresc;
		FemFunctions->qxpresc = qxpresc;
		FemFunctions->qypresc = qypresc;
		FemFunctions->zb = zb;
		FemFunctions->InitialSolution = InitialSolution;

	#endif

	return 0;
}
