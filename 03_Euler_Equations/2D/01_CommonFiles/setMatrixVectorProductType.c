#include "EulerEquations.h"
#include "../../../00_CommonFiles/MatrixVector_Operations/matvec.h"

int setMatrixVectorProductType(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->reordering,"NOT")!=0 && strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")!=0){
		printf("Reordering is applied only with CSR scheme!\n");
		exit(1);
	}
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3)==0){
		FemFunctions->assembly = ebe_assembly;
		FemFunctions->mv = ebemvNDOF4;
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE2")==0){
			FemFunctions->mv = ebe2mvNDOF4;
		}
	}
	else if (strncmp(Parameters->MatrixVectorProductScheme,"EDE",3)==0){
		FemFunctions->assembly = ede_assembly;
		FemFunctions->mv = edemvNDOF4;
		if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE2")==0){
			FemFunctions->mv = ede2mvNDOF4;
		}
	}
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		FemFunctions->assembly = csr_assembly;
		FemFunctions->mv = csrmv;
	}
	else{
		printf("\nMatrix vector product scheme is not defined correctly!\n\n");
		exit(1);
	}
	return 0;
}


