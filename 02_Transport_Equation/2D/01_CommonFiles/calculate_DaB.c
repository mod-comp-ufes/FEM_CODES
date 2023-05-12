#include "TranspEquation.h"

int calculate_DaB(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *Da, double *DaB)
{
	int neq, nel, e, eq1, eq2, eq3;
	double M2Da;
	int **lm = FemStructs->lm;	
	double **M2, **R2, *invN2; 

	neq = Parameters->neq;
	nel = Parameters->nel;
	M2 = FemStructs->AuxBuild->M2;
	R2 = FemStructs->AuxBuild->R2;
	invN2 = FemStructs->AuxBuild->invN2;

	Da[neq] = 0;

	for (e=0; e<nel; e++)
	{
		eq1 = lm[e][0];
		eq2 = lm[e][1];
		eq3 = lm[e][2];

		M2Da = M2[e][0]*Da[eq1] + M2[e][1]*Da[eq2] + M2[e][2]*Da[eq3];

		DaB[e] = invN2[e]*(R2[e][0] - M2Da);
	}		
		
	return 0;
}


