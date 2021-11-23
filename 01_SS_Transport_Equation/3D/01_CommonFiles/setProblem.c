#include "SSTransportEquation3D.h"
#include "../APLI4/apli4.h"
#include "../APLI6/apli6.h"
#include "../RAMPA2/Rampa2.h"
#include "../HIGH/high.h"
#include "../HEMKER/hemker.h"
#include "../DIFCONV/difconv.h"
#include "../DIFCONVREA/difconvrea.h"
#include "../DIFREA/difrea.h"
#include "../PAREDE/parede.h"
#include "../EXDIFREADOM/exdifreadom.h"

int setProblem(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{
	if (strcasecmp(Parameters->ProblemTitle,"HEMKER")==0){
		FemFunctions->upresc = HEMKER_upresc;
		FemFunctions->f = HEMKER_f;
		FemFunctions->kappa = HEMKER_kappa;
		FemFunctions->beta = HEMKER_beta;
		FemFunctions->sigma = HEMKER_sigma;
	}else if (strcasecmp(Parameters->ProblemTitle,"PAREDE")==0){
		FemFunctions->upresc = PAREDE_upresc;
		FemFunctions->f = PAREDE_f;
		FemFunctions->kappa = PAREDE_kappa;
		FemFunctions->beta = PAREDE_beta;
		FemFunctions->sigma = PAREDE_sigma;
	}else if (strcasecmp(Parameters->ProblemTitle,"RAMPA2")==0){
		FemFunctions->upresc = RAMPA2_upresc;
		FemFunctions->f = RAMPA2_f;
		FemFunctions->kappa = RAMPA2_kappa;
		FemFunctions->beta = RAMPA2_beta;
		FemFunctions->sigma = RAMPA2_sigma;
	}else if (strcasecmp(Parameters->ProblemTitle,"HIGH")==0){
		FemFunctions->upresc = HIGH_upresc;
		FemFunctions->f = HIGH_f;
		FemFunctions->kappa = HIGH_kappa;
		FemFunctions->beta = HIGH_beta;
		FemFunctions->sigma = HIGH_sigma;
	}else if (strcasecmp(Parameters->ProblemTitle,"EXDIFREADOM")==0){
		FemFunctions->upresc = EXDIFREADOM_upresc;
		FemFunctions->f = EXDIFREADOM_f;
		FemFunctions->kappa = EXDIFREADOM_kappa;
		FemFunctions->beta = EXDIFREADOM_beta;
		FemFunctions->sigma = EXDIFREADOM_sigma;
	}else if (strcasecmp(Parameters->ProblemTitle,"DIFCONV")==0){
		FemFunctions->upresc = DIFCONV_upresc;
		FemFunctions->f = DIFCONV_f;
		FemFunctions->ExactSolution = DIFCONV_ExactSolution;
		FemFunctions->ExactSolutionAllPoints = DIFCONV_ExactSolutionAllPoints;
		FemFunctions->kappa = DIFCONV_kappa;
		FemFunctions->beta = DIFCONV_beta;
		FemFunctions->sigma = DIFCONV_sigma;
		FemFunctions->DuDx = DIFCONV_DuDx;
		FemFunctions->DuDy = DIFCONV_DuDy;
		FemFunctions->DuDz = DIFCONV_DuDz;
	}else if (strcasecmp(Parameters->ProblemTitle,"DIFCONVREA")==0){
		FemFunctions->upresc = DIFCONVREA_upresc;
		FemFunctions->f = DIFCONVREA_f;
		FemFunctions->ExactSolution = DIFCONVREA_ExactSolution;
		FemFunctions->ExactSolutionAllPoints = DIFCONVREA_ExactSolutionAllPoints;
		FemFunctions->kappa = DIFCONVREA_kappa;
		FemFunctions->beta = DIFCONVREA_beta;
		FemFunctions->sigma = DIFCONVREA_sigma;
		FemFunctions->DuDx = DIFCONVREA_DuDx;
		FemFunctions->DuDy = DIFCONVREA_DuDy;
		FemFunctions->DuDz = DIFCONVREA_DuDz;
	}else if (strcasecmp(Parameters->ProblemTitle,"DIFREA")==0){
		FemFunctions->upresc = DIFREA_upresc;
		FemFunctions->f = DIFREA_f;
		FemFunctions->ExactSolution = DIFREA_ExactSolution;
		FemFunctions->ExactSolutionAllPoints = DIFREA_ExactSolutionAllPoints;
		FemFunctions->kappa = DIFREA_kappa;
		FemFunctions->beta = DIFREA_beta;
		FemFunctions->sigma = DIFREA_sigma;
		FemFunctions->DuDx = DIFREA_DuDx;
		FemFunctions->DuDy = DIFREA_DuDy;
		FemFunctions->DuDz = DIFREA_DuDz;
	}else if (strcasecmp(Parameters->ProblemTitle,"APLI4")==0){
		FemFunctions->upresc = APLI4_upresc;
		FemFunctions->f = APLI4_f;
		FemFunctions->ExactSolution = APLI4_ExactSolution;
		FemFunctions->ExactSolutionAllPoints = APLI4_ExactSolutionAllPoints;
		FemFunctions->kappa = APLI4_kappa;
		FemFunctions->beta = APLI4_beta;
		FemFunctions->sigma = APLI4_sigma;
		FemFunctions->DuDx = APLI4_DuDx;
		FemFunctions->DuDy = APLI4_DuDy;
		FemFunctions->DuDz = APLI4_DuDz;
	}else if (strcasecmp(Parameters->ProblemTitle,"APLI6")==0){
		FemFunctions->upresc = APLI6_upresc;
		FemFunctions->f = APLI6_f;
		FemFunctions->ExactSolution = APLI6_ExactSolution;
		FemFunctions->ExactSolutionAllPoints = APLI6_ExactSolutionAllPoints;
		FemFunctions->kappa = APLI6_kappa;
		FemFunctions->beta = APLI6_beta;
		FemFunctions->sigma = APLI6_sigma;
		FemFunctions->DuDx = APLI6_DuDx;
		FemFunctions->DuDy = APLI6_DuDy;
		FemFunctions->DuDz = APLI6_DuDz;
	}
	else{
		printf("Problem not defined!\n");
		exit(1);
	}
	return 0;
}



