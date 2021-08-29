#include "ShalowWater.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"

/*
 * Using to solve problem without bubble.
 * Call du/dt = a and Da = a_{t} - a_{t-dt}, the sistem solved is
 * Mtil*Da = R
 * where
 * Mtil = M + alpha*dt*K and R = F - Ma - Ku
*/

int PredictorMulticorrector(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
		FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)

{
	int I, i;
	int tag = 1;
	int neq, nnodes, steps;
	double t, dt, alpha, norm_a, norm_Da;
	double *a, *Da, *u, *u_old, *R; //Parametros do Preditor
	double kx, ky, kz;
	char FileName[2000];
	
	neq = Parameters->neq;
	nnodes = Parameters->nnodes;

	a = (double*) mycalloc("a of 'PreditorMulticorrector'", neq + 1, sizeof(double));
	Da = (double*) mycalloc("Da of 'PreditorMulticorrector'", neq + 1, sizeof(double));
	u_old = (double*) mycalloc("u_old of 'PreditorMulticorrector'", neq + 1, sizeof(double));

	dzero(neq + 1, a);
	dzero(neq + 1, Da);

	dt = Parameters->DeltaT;
	alpha = Parameters->Alpha;
	u = FemStructs->u;
	R = FemStructs->R;
	FemStructs->du = a;

	t = 0.0;
	Parameters->CurrentTime = 0.0;
	steps = 0;
	FemFunctions->InitialSolution(Parameters, FemStructs);

	do{
		printf("t = %lf\n",t);
		#ifdef printALL
			Paraview_Output_3D(Parameters, FemStructs, FemFunctions, t);
		#endif

		steps++;
		t += dt;
		Parameters->CurrentTime = t; 
		
		// PREDICAO
		i = 0;
		for(I=0; I<neq; I++)
		{
			u_old[I] = u[I];
			u[I] += (1.0 - alpha)*dt*a[I];
			a[I] = 0.0;
		}

		// MULTICORRECAO
		do{
			i++;
			FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);

			//FemFunctions->scaling(Parameters, MatrixData, FemStructs);

			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, R);

			FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, R, Da);

			//fprintf(OutFile, "MULTICORRECAO(%d) | iter GMRES=%d\n", i, Parameters->ContGMRES);
			printf("MULTICORRECAO(%d) | iter GMRES=%d\n", i, Parameters->ContGMRES);

			//FemFunctions->unscaling(Parameters, MatrixData, FemStructs, Da);

			daxpy(neq, 1, Da, a); // a = a + Da
			daxpy(neq, alpha*dt, Da, u); // u = u + alpha*dt*Da
			
			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = sqrt(ddot(neq, Da, Da));

		}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection

	}while(!FemFunctions->StopTimeIntegration(Parameters,u,u_old,t)); // end while time 

	myfree(a);
	myfree(Da);
	myfree(u_old);
	//fclose(OutFile);

	return 0;

}


int StopByIterations(ParametersType *Parameters, double norm_a, double norm_Da, int i)
{
	if (i>=Parameters->NonLinearMaxIter)
		return 1;
	else
		return 0;
}


int StopByNorm(ParametersType *Parameters, double norm_a, double norm_Da, int i)
{
	if (norm_Da < (Parameters->NonLinearTolerance)*norm_a || i>=Parameters->NonLinearMaxIter) {
		//printf("%d\n", i);
		return 1;
	}
	else
		return 0;
}


int StopBySteadyState(ParametersType *Parameters, double *u, double *u_old, double t)
{
	int i, neq;
	double norm_diff, norm_u;
	double *diff = mycalloc("diff of 'StopBySteadyState'",Parameters->neq,sizeof(double));
	
	neq = Parameters->neq;

	for (i=0; i<neq; i++)
		diff[i] = u[i] - u_old[i];

	norm_diff = sqrt(ddot(neq, diff, diff)); 
	norm_u = sqrt(ddot(neq,u,u));	

	free(diff);

	if((Parameters->TimeIntegrationTolerance)*norm_u > norm_diff || t >= Parameters->FinalTime){
		Parameters->CurrentTime = t;	
		return 1;
	}
	else
		return 0;
}


int StopByTime(ParametersType *Parameters, double *u, double *u_old, double t)
{
	if(t >= Parameters->FinalTime)
	{
		Parameters->CurrentTime = t;	
		return 1;
	}
	else
		return 0;
}
