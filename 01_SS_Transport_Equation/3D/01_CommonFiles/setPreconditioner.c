#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/preconditioners.h"

int setPreconditioner(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Preconditioner,"NOT")==0){
		FemFunctions->precond = NO_precond;
		FemFunctions->precondR = NO_precond;
		FemFunctions->precond_setup = NO_precond_setup;
	} else if (strcasecmp(Parameters->Preconditioner,"LU")==0){
		FemFunctions->precond_setup = LU_precond_EBE_setup_NNOEL4;
		FemFunctions->precondR = NO_precond;
		FemFunctions->precond = LU_precond_EBE_NNOEL4;
	} else if (strncmp(Parameters->Preconditioner,"ILU",3)==0){
		FemFunctions->precondR = NO_precond;
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
			FemFunctions->precond = ILUp_precond;
			FemFunctions->precond_setup = ILUp_precond_setup;
		} else {
			printf("Preconditioner definied only to CSR scheme\n");
			exit(1);
		}
	} else {
		printf("Preconditioner not defined!\n");
        exit(1);
	}

    return 0;
}
