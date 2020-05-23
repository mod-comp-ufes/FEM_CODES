#include "CPU_time.h"

void calculateTime(double Preprocess_Time, double Process_Time, double Postprocess_Time, ParametersType *Parameters)
{
	int h, m;
	double s, Total_Time;
	FILE *OutFile,*CSV;
	char FileName[2000],CSVname[2000],CSVcontent[4000];

	Total_Time = Preprocess_Time + Process_Time + Postprocess_Time;
	h = (int) Total_Time;
	h = h/3600;
	m = (int) (Total_Time - 3600*h);
	m = m/60;
	s = Total_Time - (3600*h + 60*m);

	#ifdef SSTranspEquation2D
		sprintf(FileName, "../03_output/%s_%s_%s_%s_%s_%s_%s_N%d_E%d.dat",Parameters->ProblemTitle,Parameters->StabilizationForm,Parameters->ShockCapture, Parameters->h_Shock, 
		Parameters->MatrixVectorProductScheme,Parameters->Solver,Parameters->Preconditioner, Parameters->nnodes,Parameters->nel); 	
	#endif

	#ifdef TranspEquation2D
		sprintf(FileName, "../03_output/%s_%s_%s_%s_%s_%s_%s_%s_N%d_E%d.dat",Parameters->ProblemTitle,Parameters->TimeIntegration,Parameters->StabilizationForm,Parameters->ShockCapture,
	    Parameters->h_Shock, Parameters->MatrixVectorProductScheme,Parameters->Solver,Parameters->Preconditioner,Parameters->nnodes,Parameters->nel); 	
	#endif

	#ifdef EulerEquations2D
		sprintf(FileName,"../03_output/%s_%s_%s_%s_%s_%s_N%d_E%d.txt", Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->ShockCapture, 
		Parameters->TimeIntegration, Parameters->MatrixVectorProductScheme,Parameters->nnodes, Parameters->nel);
		sprintf(CSVname,"../../../../OUTPUT_DATA/EulerEquations2D.csv");
		sprintf(CSVcontent,"\n\"%lf\",\"%lf\",\"%lf\",\"%d\",\"%lf\",\"%lf\"", Preprocess_Time, Process_Time, Postprocess_Time, Parameters->iterations, /*media iteracoes*/Parameters->iterations*Parameters->DeltaT/(Parameters->NonLinearMaxIter*Parameters->FinalTime), /*media tempo*/Process_Time*Parameters->DeltaT/(Parameters->NonLinearMaxIter*Parameters->FinalTime);
		CSV = myfopen(CSVname,"a");
		fprintf(CSV,"%s", CSVcontent);
		fclose(CSV);
	#endif
	
	#ifdef SSNavierStokesEquations2D
		sprintf(FileName,"../03_output/%s_%s_%s_%s_%s_%s_N%d_E%d.txt", Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->ShockCapture,
		Parameters->TimeIntegration, Parameters->MatrixVectorProductScheme,Parameters->nnodes, Parameters->nel);
		sprintf(CSVname,"../03_output/SSNavierStokesEquations2D.csv");
		sprintf(CSVcontent,"\n\"%lf\",\"%lf\",\"%lf\",\"%d\",\"%d\",\"%lf\",\"%lf\"", Preprocess_Time, Process_Time, Postprocess_Time, Parameters->gmres, Parameters->iterations, /*media iteracoes*/1.0*Parameters->iterations/Parameters->gmres, /*media tempo*/Process_Time/Parameters->gmres);
		CSV = myfopen(CSVname,"a");
		fprintf(CSV,"%s", CSVcontent);
		fclose(CSV);
	#endif
	
	OutFile = myfopen(FileName,"a");
	fprintf(OutFile,"\n============================ Processing Time ============================\n\n");
	fprintf(OutFile,"\nPreprocess Time: %lf\n", Preprocess_Time);	
	fprintf(OutFile,"\nProcess Time: %lf\n", Process_Time);	
	fprintf(OutFile,"\nPostprocess Time: %lf\n", Postprocess_Time);	
	fprintf(OutFile,"\nTOTAL TIME: %lf (%dh %dm %lfs)\n", Total_Time, h, m, s);
	fprintf(OutFile,"\n=========================================================================\n\n");
	fclose(OutFile);
	printf("\n============================ Processing Time ============================\n\n");
	printf("\nPreprocess Time: %lf\n", Preprocess_Time);	
	printf("\nProcess Time: %lf\n", Process_Time);	
	printf("\nPostprocess Time: %lf\n", Postprocess_Time);	
	printf("\nTOTAL TIME: %lf (%dh %dm %lfs)\n", Total_Time, h, m, s);
	printf("\n=========================================================================\n\n");
	
}



