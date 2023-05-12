#include "NavierStokesEquations.h"
#include "../../../00_CommonFiles/MatrixVector_Operations/matvec.h"

int setMatrixVectorProductType(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->reordering,"NOT")!=0 && strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")!=0){
		printf("Reordering is applied only with CSR scheme!\n");
		exit(1);
	}
	if (strncmp(Parameters->MatrixVectorProductScheme,"EBE",3)==0){
		FemFunctions->assembly = ebe_assembly;
		FemFunctions->mv = ebemvNDOF3;
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


