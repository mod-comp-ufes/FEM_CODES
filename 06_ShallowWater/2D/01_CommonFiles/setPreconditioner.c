#include "ShalowWater.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/preconditioners.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/scaling.h"


int setPreconditioner(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Preconditioner,"NOT")==0)
	{
		FemFunctions->precond = NO_precond;
		FemFunctions->precondR = NO_precond;
		FemFunctions->precond_setup = NO_precond_setup;
	}
	/****************************************************************************/
	else if (strcasecmp(Parameters->Preconditioner,"Diag")==0)
	{
		FemFunctions->precond = Diag_precond;
		FemFunctions->precondR = NO_precond;
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0)
			FemFunctions->precond_setup = Diag_precond_CSR_setup;
		else
		{
			printf("Preconditioner definied only to CSR scheme\n");
			exit(1);
		}
	}
	/****************************************************************************/
	else if (strncmp(Parameters->Preconditioner,"ILU",3)==0)
	{
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0)
		{
			FemFunctions->precond = ILUp_precond;
			FemFunctions->precondR = NO_precond;
			FemFunctions->precond_setup = ILUp_precond_setup;
		}
		else
		{
			printf("Preconditioner definied only to CSR scheme\n");
			exit(1);
		}
	}
	/****************************************************************************/
	else if (strncmp(Parameters->Preconditioner,"AMG",3)==0)
	{
		FemFunctions->precondR = NO_precond;
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0)
		{
            FemFunctions->precond = AMG_precond;
			FemFunctions->precond_setup = AMG_precond_setup;
		}
		else
		{
			printf("Preconditioner defined only to CSR scheme\n");
			exit(1);
		}
	}
	
	/****************************************************************************/
	else
	{
		printf("Preconditioner is not defined correctly!\n");
		exit(1);
	}

	return 0;
}
