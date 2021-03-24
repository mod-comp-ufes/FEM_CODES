#include "ShalowWater.h"


int setProblem(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->ProblemTitle,"OBLIQUO")==0)
	{
		//FemFunctions->rhopresc = OBLIQUO_rhopresc;
		//FemFunctions->v1presc = OBLIQUO_v1presc;
		//FemFunctions->v2presc = OBLIQUO_v2presc;
	}
	else{
		printf("Problem is not defined correctly!\n");
		exit(1);

	}
	return 0;
}
