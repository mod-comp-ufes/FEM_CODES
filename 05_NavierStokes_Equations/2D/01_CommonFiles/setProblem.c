#include "NavierStokesEquations.h"
#include "../CAVITY/cavity.h"
#include "../CHANNEL/channel.h"
#include "../EXATA/exata.h"

int setProblem(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->ProblemTitle,"CAVITY")==0){
		FemFunctions->v1presc = CAVITY_v1presc;
		FemFunctions->v2presc = CAVITY_v2presc;
		FemFunctions->ppresc = CAVITY_ppresc;
		FemFunctions->InitialSolution = CAVITY_InitialSolution;
		FemFunctions->f1ext = CAVITY_f1ext;
		FemFunctions->f2ext = CAVITY_f2ext;
	}else if (strcasecmp(Parameters->ProblemTitle,"CHANNEL")==0){
		FemFunctions->v1presc = CHANNEL_v1presc;
		FemFunctions->v2presc = CHANNEL_v2presc;
		FemFunctions->ppresc = CHANNEL_ppresc;
		FemFunctions->InitialSolution = CHANNEL_InitialSolution;	
		FemFunctions->f1ext = CHANNEL_f1ext;
		FemFunctions->f2ext = CHANNEL_f2ext;
	}else if (strcasecmp(Parameters->ProblemTitle,"EXATA")==0){
		FemFunctions->v1presc = EXATA_v1presc;
		FemFunctions->v2presc = EXATA_v2presc;
		FemFunctions->ppresc = EXATA_ppresc;
		FemFunctions->InitialSolution = EXATA_InitialSolution;
		FemFunctions->f1ext = EXATA_f1ext;
		FemFunctions->f2ext = EXATA_f2ext;
	}
	else{
		printf("Problem not defined!\n");
		exit(1);
	}
	return 0;
}
