#include "../preconditioners.h"
#include "../../Allocation_Operations/allocations.h"

int SSOR_precond_EDE_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F)
{
	SSOR_precond_EDE (Parameters, MatrixData, FemStructs, F, F);

	return 0;
}
