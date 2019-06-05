#include "EulerEquations.h"
#include "../../../00_CommonFiles/Solvers_and_Preconditioners/preconditioners.h"

int setPreconditioner(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Preconditioner,"NOT")==0){
		FemFunctions->precond = NO_precond;
		FemFunctions->precondR = NO_precond;
		FemFunctions->precond_setup = NO_precond_setup;
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE2")==0){
			FemFunctions->precond_setup = NOBlockDiag2_precond_EBE_setup;
		}
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE2")==0){
			FemFunctions->precond_setup = NOBlockDiag2_precond_EDE_setup;
		}
	}
	/****************************************************************************/
	else if (strcasecmp(Parameters->Preconditioner,"BlockDiag")==0){
		FemFunctions->precond = BlockDiag_precond;
		FemFunctions->precondR = NO_precond;
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
			FemFunctions->precond_setup = BlockDiag_precond_EBE_setup;
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
			FemFunctions->precond_setup = BlockDiag_precond_EDE_setup;
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE2")==0){
			FemFunctions->precond_setup = BlockDiag_precond_EBE_setup;
			FemFunctions->precond = NO_precond;
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE2")==0){
			FemFunctions->precond_setup = BlockDiag_precond_EDE_setup;
			FemFunctions->precond = NO_precond;
		}
		else{
			printf("Preconditioner definied only to EBE and EDE schemes\n");
			exit(1);
		}
	}
	/****************************************************************************/
	else if (strcasecmp(Parameters->Preconditioner,"SGSBlock")==0){
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
			FemFunctions->precond_setup = SGSBlock_precond_EBE_setup;
			FemFunctions->precond = SGSBlock_precond_EBE;
			FemFunctions->precondR = SGSBlock_precondR_EBE;
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
//			FemFunctions->precond_setup = JacobiBlockDOF4_precond_EDE_setup;
//			FemFunctions->precond = JacobiBlockDOF4_precond_EDE;
		}
		else{
			printf("Preconditioner definied only to EBE and EDE schemes\n");
			exit(1);
		}
	}
	/****************************************************************************/
	else if (strcasecmp(Parameters->Preconditioner,"LUBlock")==0){
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
			FemFunctions->precond_setup = LUBlock_precond_EBE_setup;
			FemFunctions->precond = LUBlock_precond_EBE;
			FemFunctions->precondR = NO_precond;
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
//			FemFunctions->precond_setup = JacobiBlockDOF4_precond_EDE_setup;
//			FemFunctions->precond = JacobiBlockDOF4_precond_EDE;
		}
		else{
			printf("Preconditioner definied only to EBE and EDE schemes\n");
			exit(1);
		}
	}

	/****************************************************************************/
	else if (strcasecmp(Parameters->Preconditioner,"JacobiBlock")==0){
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
			FemFunctions->precond_setup = JacobiBlockDOF4_precond_EBE_setup;
			FemFunctions->precond = JacobiBlockDOF4_precond_EBE;
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
//			FemFunctions->precond_setup = JacobiBlockDOF4_precond_EDE_setup;
//			FemFunctions->precond = JacobiBlockDOF4_precond_EDE;
		}
		else{
			printf("Preconditioner definied only to EBE and EDE schemes\n");
			exit(1);
		}
	}
	/****************************************************************************/
	else if (strcasecmp(Parameters->Preconditioner,"Jacobi")==0){
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
			FemFunctions->precond_setup = JacobiDOF4_precond_EBE_setup;
			FemFunctions->precond = JacobiDOF4_precond_EBE;
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
//			FemFunctions->precond_setup = JacobiDOF4_precond_EDE_setup;
//			FemFunctions->precond = JacobiDOF4_precond_EDE;
		}
		else{
			printf("Preconditioner definied only to EBE and EDE schemes\n");
			exit(1);
		}
	}
	/****************************************************************************/
	else if (strcasecmp(Parameters->Preconditioner,"Diag")==0){
		FemFunctions->precond = Diag_precond;
		FemFunctions->precondR = NO_precond;
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
			FemFunctions->precond_setup = Diag_precond_EBE_setup;
		}
		else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){
			FemFunctions->precond_setup = Diag_precond_EDE_setup;
		}
		else{
			FemFunctions->precond_setup = Diag_precond_CSR_setup;
		}
	}
	/****************************************************************************/
	else if (strncmp(Parameters->Preconditioner,"ILU",3)==0){
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
			FemFunctions->precond = ILUp_precond;
			FemFunctions->precondR = NO_precond;
			FemFunctions->precond_setup = ILUp_precond_setup;
		}
		else{
			printf("Preconditioner definied only to CSR scheme\n");
			exit(1);
		}
	}
	/****************************************************************************/
	else if (strncmp(Parameters->Preconditioner,"SORBlock",8)==0){
		if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3)==0){
    	FemFunctions->precond = SSORBlockDOF4_precond_EBE;
			FemFunctions->precond_setup = SSORBlockDOF4_precond_EBE_setup;
		}
		/*else if (strncmp(Parameters->MatrixVectorProductScheme,"EDE",3)==0){
      *precond = SSOR_precond_EDE;
			*precond_setup = SSOR_precond_EDE_setup;
		}*/
		else{
			printf("Preconditioner definied only to EBE or EDE schemes\n");
			exit(1);
		}
	}
  /****************************************************************************/
  else if (strncmp(Parameters->Preconditioner,"SOR",3)==0){
		if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3)==0){
    	FemFunctions->precond = SSORDOF4_precond_EBE;
			FemFunctions->precond_setup = SSORDOF4_precond_EBE_setup;
		}
		/*else if (strncmp(Parameters->MatrixVectorProductScheme,"EDE",3)==0){
      FemFunctions->precond = SSORDOF4_precond_EDE;
			FemFunctions->precond_setup = SSORDOF4_precond_EDE_setup;
		}*/
		else{
			printf("Preconditioner definied only to EBE or EDE schemes\n");
			exit(1);
		}
	}
	/****************************************************************************/
	else {
		printf("Preconditioner is not defined correctly!\n");
		exit(1);
	}

	return 0;
}
