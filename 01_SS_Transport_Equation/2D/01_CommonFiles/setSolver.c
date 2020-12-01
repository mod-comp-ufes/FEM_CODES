#include "SSTranspEquation.h"
#include "../../../00_CommonFiles/MatrixVector_Operations/matvec.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/solvers.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/CSR/AMG/amg_pgmres.h"

int setSolver(ParametersType *Parameters, FemOtherFunctionsType *FemOtherFunctions)
{
	if (strcasecmp(Parameters->Solver,"GMRES")==0){
		if (strncmp(Parameters->Preconditioner,"AMG",3)==0) {
			FemOtherFunctions->solver = pgmres;// AMG_GMRES;
		}
		else{
			FemOtherFunctions->solver = pgmres;
		}
	}
//	else if (strcasecmp(Parameters->Solver,"CG")==0){
//		FemOtherfunctions->solver = pcg;
//	}
	else{
		printf("Solver is not defined correctly!\n");
		exit(1);
	}

	return 0;
}



