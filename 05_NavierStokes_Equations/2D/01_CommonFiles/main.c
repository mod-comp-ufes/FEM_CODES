#include "NavierStokesEquations.h"
#include "../../../00_CommonFiles/CPU_Time_Operations/CPU_time.h"

int main(int argc, char **argv)
{
	ParametersType *Parameters;
	MatrixDataType *MatrixData;
	FemStructsType *FemStructs;
	FemFunctionsType *FemFunctions;
	FemOtherFunctionsType *FemOtherFunctions;
	struct timespec Start, End;
	double Preprocess_Time, Process_Time, Postprocess_Time;


	/* ******************************************************* Preprocess ****************************************************** */
	clock_gettime(CLOCK_MONOTONIC, &Start);
	Preprocess(argc, argv, &Parameters, &MatrixData, &FemStructs, &FemFunctions, &FemOtherFunctions);
	clock_gettime(CLOCK_MONOTONIC, &End);
	Preprocess_Time = End.tv_sec - Start.tv_sec + 1e-9*(End.tv_nsec - Start.tv_nsec);
	/* ************************************************************************************************************************* */

	/* ******************************************************** Process ******************************************************** */
	clock_gettime(CLOCK_MONOTONIC, &Start);
	Process(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
	clock_gettime(CLOCK_MONOTONIC, &End);
	Process_Time = End.tv_sec - Start.tv_sec + 1e-9*(End.tv_nsec - Start.tv_nsec);
	/* ************************************************************************************************************************* */
	
	/* ******************************************************* Postprocess ***************************************************** */

	clock_gettime(CLOCK_MONOTONIC, &Start);
	Postprocess(Parameters, MatrixData, FemStructs, FemFunctions, FemOtherFunctions);
	clock_gettime(CLOCK_MONOTONIC, &End);
	Postprocess_Time = End.tv_sec - Start.tv_sec + 1e-9*(End.tv_nsec - Start.tv_nsec);

	/* ************************************************************************************************************************* */

	/* ******************************************************* Total Time ****************************************************** */
	calculateTime(Preprocess_Time, Process_Time, Postprocess_Time, Parameters);
	/* ************************************************************************************************************************* */

	free(Parameters);
	
	return 0;
}


