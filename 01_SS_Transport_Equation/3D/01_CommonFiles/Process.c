#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

int Process(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)
{
	
	setProblem(Parameters, FemFunctions);
	setMatrixVectorProductType(Parameters, FemFunctions);
	setSolver(Parameters,FemOtherFunctions);
	setPreconditioner(Parameters, FemFunctions);
	setStabilizationForm(Parameters, FemOtherFunctions);
	
	if ((strcasecmp(Parameters->StabilizationForm,"NSGS")==0)||(strcasecmp(Parameters->StabilizationForm,"DD")==0)||(strcasecmp(Parameters->StabilizationForm,"CAU")==0)){
		if(strcasecmp(Parameters->UseDamping,"YES")==0){
			Space_Algorithm_Damping(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
		}else{
			Space_Algorithm_NonLinear(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
		}
	}else if ((strcasecmp(Parameters->StabilizationForm,"Galerkin")==0)||(strcasecmp(Parameters->StabilizationForm,"VMS")==0)||(strcasecmp(Parameters->StabilizationForm,"SUPG")==0)){
		if(strcasecmp(Parameters->UseDamping,"YES")==0){
			Space_Algorithm_Damping(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
		}else{
			Space_Algorithm_Galerkin(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
		}
	}else{ 
		printf("\nProcess Error: Undefined combination between stabilization form and solution algorithm!\n\n");
		exit(1);
	}
	
	
	return 0;
	
}


