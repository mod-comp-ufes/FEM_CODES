#include "SSTranspEquation.h"
#include "../PUDIM/pudim.h"
#include "../CARTOLA/cartola.h"
#include "../TESTE/teste.h"
#include "../HEMKER/hemker.h"
#include "../CONVECTION/convection.h"
#include "../REACTION/reaction.h"
#include "../RAMPA/rampa.h"
#include "../RAMPA2/rampa2.h"
#include "../VALIDACAO/validacao.h"

int setProblem(ParametersType *Parameters, FemFunctionsType *FemFunctions)
{		
	if (strcasecmp(Parameters->ProblemTitle,"RAMPA")==0){
		FemFunctions->Condutivity = RAMPA_Condutivity;	
		FemFunctions->Font = RAMPA_Font;
		FemFunctions->Reaction = RAMPA_Reaction;
		FemFunctions->Velocity = RAMPA_Velocity;
		FemFunctions->upresc = RAMPA_upresc;
	}else if (strcasecmp(Parameters->ProblemTitle,"RAMPA2")==0){
		FemFunctions->Condutivity = RAMPA2_Condutivity;	
		FemFunctions->Font = RAMPA2_Font;
		FemFunctions->Reaction = RAMPA2_Reaction;
		FemFunctions->Velocity = RAMPA2_Velocity;
		FemFunctions->upresc = RAMPA2_upresc;
	}else if (strcasecmp(Parameters->ProblemTitle,"PUDIM")==0){
		FemFunctions->Condutivity = PUDIM_Condutivity;	
		FemFunctions->Font = PUDIM_Font;
		FemFunctions->Reaction = PUDIM_Reaction;
		FemFunctions->Velocity = PUDIM_Velocity;
		FemFunctions->upresc = PUDIM_upresc;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"CARTOLA")==0){
		FemFunctions->Condutivity = CARTOLA_Condutivity;	
		FemFunctions->Font = CARTOLA_Font;
		FemFunctions->Reaction = CARTOLA_Reaction;
		FemFunctions->Velocity = CARTOLA_Velocity;
		FemFunctions->upresc = CARTOLA_upresc;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"TESTE")==0){
		FemFunctions->Condutivity = TESTE_Condutivity;	
		FemFunctions->Font = TESTE_Font;
		FemFunctions->Reaction = TESTE_Reaction;
		FemFunctions->Velocity = TESTE_Velocity;
		FemFunctions->upresc = TESTE_upresc;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"HEMKER")==0){
		FemFunctions->Condutivity = HEMKER_Condutivity;	
		FemFunctions->Font = HEMKER_Font;
		FemFunctions->Reaction = HEMKER_Reaction;
		FemFunctions->Velocity = HEMKER_Velocity;
		FemFunctions->upresc = HEMKER_upresc;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"CONVECTION")==0){
		FemFunctions->Condutivity = CONVECTION_Condutivity;	
		FemFunctions->Font = CONVECTION_Font;
		FemFunctions->Reaction = CONVECTION_Reaction;
		FemFunctions->Velocity = CONVECTION_Velocity;
		FemFunctions->upresc = CONVECTION_upresc;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"REACTION")==0){
		FemFunctions->Condutivity = REACTION_Condutivity;	
		FemFunctions->Font = REACTION_Font;
		FemFunctions->Reaction = REACTION_Reaction;
		FemFunctions->Velocity = REACTION_Velocity;
		FemFunctions->upresc = REACTION_upresc;
	}
	else if (strcasecmp(Parameters->ProblemTitle,"VALIDACAO")==0){
		FemFunctions->Condutivity = VALIDACAO_Condutivity;	
		FemFunctions->Font = VALIDACAO_Font;
		FemFunctions->Reaction = VALIDACAO_Reaction;
		FemFunctions->Velocity = VALIDACAO_Velocity;
		FemFunctions->upresc = VALIDACAO_upresc;
		FemFunctions->hflux = VALIDACAO_hflux;
	}
	else{
		printf("Problem not defined!\n");
		exit(1);
	}
	
	return 0;
}



