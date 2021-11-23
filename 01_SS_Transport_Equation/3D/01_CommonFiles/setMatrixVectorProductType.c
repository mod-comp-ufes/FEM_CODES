#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/MatrixVector_Operations/matvec.h"

int setMatrixVectorProductType(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->Reordering,"NOT")!=0 && strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")!=0){
		printf("Reordering is applied only with CSR scheme!\n");
		exit(1);
	}
	if (strcasecmp(Parameters->MatrixVectorProductScheme,"EBE")==0){
		FemFunctions->MatrixAssembly = ebe_assembly;
		FemFunctions->ProductMatrixVector = ebemv3Dndof1;
		
	}
/*	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"EDE")==0){*/
/*		*assembly = ede_assembly;*/
/*		*end_assembly = ede_end_assembly;*/
/*		*mv = edemv;*/
/*		*mvMK = edemvMK;*/
/*	} */
	else if (strcasecmp(Parameters->MatrixVectorProductScheme,"CSR")==0){
		FemFunctions->MatrixAssembly = csr_assembly;
		FemFunctions->ProductMatrixVector = csrmv;
	}
	else{
		printf("\nMatrix vector product scheme not defined!\n\n");
		exit(1);
	}
	return 0;
}


