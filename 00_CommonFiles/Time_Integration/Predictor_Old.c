#include "time_integration.h"

int Predictor_Old(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs,
		FemFunctionsType *FemFunctions, FemOtherFunctionsType *FemOtherFunctions)

{
	int I, i;
	int nel, neq, passo;
	double t, dt, alpha, norm_a, norm_Da, tol_correction;
	double *a, *auxVec, *Da, *u, *u_old, *R; //Parametros do Preditor
	double *uB_old, *delta_old;
	AuxBuildStructuresType *AuxBuild;

	nel = Parameters->nel;
	neq = Parameters->neq;
	

	a = (double*) mycalloc("a of 'Preditor_Old'", neq + 1, sizeof(double));
	auxVec = (double*) mycalloc("auxVec of 'Preditor_Old'", neq + 1, sizeof(double));
	Da = (double*) mycalloc("Da of 'Preditor_Old'", neq + 1, sizeof(double));
	u_old = (double*) mycalloc("u_old of 'Preditor_Old'", neq + 1, sizeof(double));
	uB_old = (double*) mycalloc("uB_old of 'Preditor_Old'", nel*NDOF, sizeof(double));
	delta_old = (double*) mycalloc("delta_old of 'Preditor_Old'", nel, sizeof(double));

	dt = Parameters->DeltaT;
	alpha = Parameters->Alpha;
	Parameters->DeltaT_Build = dt;
	Parameters->Alpha_Build = alpha;
	tol_correction = Parameters->NonLinearTolerance;
	AuxBuild = (AuxBuildStructuresType*) mycalloc("AuxBuild of 'Predictor_Old'",1,sizeof(AuxBuildStructuresType));
	AuxBuild->tolerance = Parameters->StabilizationTolerance;
	FemStructs->AuxBuild = AuxBuild;
	u = FemStructs->u;
	R = FemStructs->F;
	FemStructs->du = a;
	FemStructs->uB = uB_old;
	FemStructs->duB = NULL;
	FemStructs->delta_old = delta_old;

	FemFunctions->InitialSolution(Parameters, FemStructs->Node, u);

	t = 0.0;
	passo = 0;
	int tag = 1;

	do{
		passo++;
		t += dt;

		#ifdef debug
			printf("\n\n Passo: %d\n", passo); 
		#endif

		//PREDICAO
		i = 0;

		for(I = 0; I < neq; I++)
		{
			u_old[I] = u[I];
			u[I] += (1.0 - alpha)*dt*a[I];
			a[I] = 0.0;
		}
		 
		norm_a = sqrt(ddot(neq, a, a));
		norm_Da = norm_a + 1.0;

		// MULTICORRECAO
		do{
			i++;

			FemOtherFunctions->Build(Parameters, MatrixData, FemStructs, FemFunctions);
			
			FemFunctions->scaling(Parameters, MatrixData, FemStructs);

			FemFunctions->precond_setup(Parameters, MatrixData, FemStructs, tag++, R);

			//---------------------------------------------------------------------------------	
		/*	int I, J;	
			FILE *Out;
			Out = fopen("Dados.m","w");
			fprintf(Out,"Ae=[\n");
			for (I=0; I<Parameters->nel; I++){
				for (J=0; J<144; J++)
					fprintf(Out,"%.15lf\t",MatrixData->A[I][J]);
				fprintf(Out,";\n");	
			}
			fprintf(Out,"];\n");
			fprintf(Out,"invDiag=[\n");
			for (I=0; I<Parameters->nel; I++){
				for (J=0; J<16; J++)
					fprintf(Out,"%.15lf\t",MatrixData->invBlockDiag[I][J]);
				fprintf(Out,";\n");
			}	
			fprintf(Out,"];\n");
			fprintf(Out,"lm=[\n");
			for (I=0; I<Parameters->nel; I++){
				for (J=0; J<12; J++)
					fprintf(Out,"%d\t",FemStructs->lm[I][J]+1);
				fprintf(Out,";\n");	
			}			
			fprintf(Out,"];\n");
			fprintf(Out,"F=[\n");
			for (I=0;I<Parameters->neq; I++)
				fprintf(Out,"%.15lf;\n",R[I]);
			fprintf(Out,"];\n");
			fclose(Out);	
			exit(1);
		*/
			//---------------------------------------------------------------------------------	

			FemOtherFunctions->solver(Parameters, MatrixData, FemStructs, FemFunctions, R, Da);

			FemFunctions->unscaling(Parameters, MatrixData, FemStructs, Da);

			daxpy(neq, 1, Da, a);
			daxpy(neq, alpha*dt, Da, u);
			
			norm_a = sqrt(ddot(neq, a, a));
			norm_Da = sqrt(ddot(neq, Da, Da));
			#ifdef debug
				double normR;
				normR = sqrt(ddot(neq, R, R));
				printf("Tol_correction = %lf \t  Norma_Res =%lf \t Norm a = %lf \t  Norma Da = %lf \t t = %lf \t i = %d \n", 
					tol_correction*norm_a, normR, norm_a, norm_Da, t, i);
			#endif
		}while(!FemFunctions->StopCriteria(Parameters,norm_a,norm_Da,i)); // end while multicorrection

		#ifdef debug
			printf("\n\n"); 
		#endif
		
	}while(!FemFunctions->StopTimeIntegration(Parameters,u,u_old,t)); // end while time

	
	free(a);
	free(Da);
	free(u_old);
	free(auxVec);
	free(uB_old);
	free(delta_old);
	free(AuxBuild);

	return 0;

}



