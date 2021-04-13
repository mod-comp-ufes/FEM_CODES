#include "ShalowWater.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
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
	int neq, nnodes, steps;
	double t, dt, alpha, norm_a, norm_Da;
	double *a, *Da, *u, *u_old, *R; //Parametros do Preditor
	double kx, ky, kz;
	
	neq = Parameters->neq;
	nnodes = Parameters->nnodes;

	a = (double*) mycalloc("a of 'PreditorMulticorrector'", neq + 1, sizeof(double));
	Da = (double*) mycalloc("Da of 'PreditorMulticorrector'", neq + 1, sizeof(double));
	u_old = (double*) mycalloc("u_old of 'PreditorMulticorrector'", neq + 1, sizeof(double));

	dt = Parameters->DeltaT;
	alpha = Parameters->Alpha;
	u = FemStructs->u;
	R = FemStructs->F;
	FemStructs->du = a;
	
	t = 0.0;
	Parameters->CurrentTime = 0.0;
	steps = 0;
	int tag = 1;
	
	FemFunctions->InitialSolution(Parameters, FemStructs);
	
	do{
		steps++;
		//printf("\n*Predictor Steps: %d\n", steps);
		t += dt;
		Parameters->CurrentTime = t; 
		
		//PREDICAO
		i = 0;

		for(I = 0; I < neq; I++)
		{
			u_old[I] = u[I];
			u[I] += (1.0 - alpha)*dt*a[I];
			a[I] = 0.0;
		}

		/*for(I = 0; I < neq; I++){
			printf("uOld[%d] = %lf\n",I, u_old[I]);
			printf("u[%d] = %lf\n",I, u[I]);
			printf("a[%d] = %lf\n\n",I, a[I]);
		}
		getchar();*/

		// MULTICORRECAO
		do{
			i++;
			//printf("\n--Multicorrection: %d\n", i); 

			//FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
			Build_M_D_F_SUPG(Parameters, MatrixData, FemStructs, FemFunctions);

			FemFunctions->scaling(Parameters, MatrixData, FemStructs);

			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, R);

			FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, R, Da);

			FemFunctions->unscaling(Parameters, MatrixData, FemStructs, Da);

			daxpy(neq, 1, Da, a); // a = a + Da
			daxpy(neq, alpha*dt, Da, u); // u = u + alpha*dt*Da
			
			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = sqrt(ddot(neq, Da, Da));
			
			//#ifdef debug
				double normR, tol_correction;
				normR = sqrt(ddot(neq, R, R));
				tol_correction = Parameters->NonLinearTolerance;
				//printf("\n     Tol_correction = %e \t Tol*norm_a = %e \t Norm a = %e \t Norm Da = %e \t Norm_Res =%e \t i = %d \n", 
				//	tol_correction, tol_correction*norm_a, norm_a, norm_Da, normR, i);
				//getchar();
			//#endif
			
		}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection

		//#ifdef debug
			double diff[neq], normDiff, normU, tol_time;
			for(I = 0; I < neq; I++){
				diff[I] = u_old[i] - u[i];
			}
			normDiff = sqrt(ddot(neq, diff, diff));
			normU = sqrt(ddot(neq, u, u));
			tol_time = Parameters->TimeIntegrationTolerance;
			//printf("\n     Tol_time = %e \t Tol*norm_u = %e \t Norm U = %e \t NormDiff =%e \t t = %lf \n", 
			//		tol_time, tol_time*normU, normU, normDiff,t);
			//getchar();
			//printf("\n"); 
		//#endif
		
	}while(!FemFunctions->StopTimeIntegration(Parameters,u,u_old,t)); // end while time 

	myfree(a);
	myfree(Da);
	myfree(u_old);

	return 0;

}



