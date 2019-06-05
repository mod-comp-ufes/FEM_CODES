#include "EulerEquations.h"

int calculate_DaB(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *Da, double *DaB)
{

	int neq, nel, e, eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12;
	int eNDOF;
	double M2Da[4];
	int **lm = FemStructs->lm;	
	double **M2, **R2, *invN2; 
	
//	int J1, J2, J3;
//	double theta[3], proj1v1, proj1v2, proj2v1, proj2v2, proj3v1, proj3v2;
//	NodeType *Node = FemStructs->Node;
//	ElementType *Element = FemStructs->Element;

//	double DaJ1[4], DaJ2[4], DaJ3[4], DaBe[4];//variaveis debugger
	

	neq = Parameters->neq;
	nel = Parameters->nel;
	M2 = FemStructs->AuxBuild->M2;
	R2 = FemStructs->AuxBuild->R2;
	invN2 = FemStructs->AuxBuild->invN2;

	Da[neq] = 0.0;

	for(e = 0; e < nel; e++){
		
//		J1 = Element[e].Vertex[0];
//		J2 = Element[e].Vertex[1];
//		J3 = Element[e].Vertex[2];
//		FemFunctions->BC_General_theta(J1,J2,J3,Node,theta,FemFunctions->BC_theta);
		
		eq1 = lm[e][0];
		eq2 = lm[e][1];
		eq3 = lm[e][2];
		eq4 = lm[e][3];
		
		eq5 = lm[e][4];
		eq6 = lm[e][5];
		eq7 = lm[e][6];
		eq8 = lm[e][7];
		
		eq9 = lm[e][8];
		eq10 = lm[e][9];
		eq11 = lm[e][10];
		eq12 = lm[e][11];
		
/*		//Ajuste projecao J1
		if(Node[J1].v1Type < 0){
			proj1v1 = cos( theta[0] ) * Da[eq2];
			proj1v2 = sin( theta[0] ) * Da[eq2];
		}
		else{
			proj1v1 = Da[eq2];
			proj1v2 = Da[eq3];
		}
		//Ajuste projecao J2
		if(Node[J2].v1Type < 0){
			proj2v1 = cos( theta[1] ) * Da[eq6];
			proj2v2 = sin( theta[1] ) * Da[eq6];
		}
		else{
			proj2v1 = Da[eq6];
			proj2v2 = Da[eq7];
		}
		//Ajuste projecao J3
		if(Node[J3].v1Type < 0){
			proj3v1 = cos( theta[2] ) * Da[eq10];
			proj3v2 = sin( theta[2] ) * Da[eq10];
		}
		else{
			proj3v1 = Da[eq10];
			proj3v2 = Da[eq11];
		}
		
*/		
/*	M2Da[0] = M2[e][0]*Da[eq1]  + M2[e][1]*proj1v1		+ M2[e][2]*proj1v2  	+ M2[e][3]*Da[eq4]  
				+ M2[e][4]*Da[eq5]  + M2[e][5]*proj2v1  	+ M2[e][6]*proj2v2  	+ M2[e][7]*Da[eq8]
				+ M2[e][8]*Da[eq9]  + M2[e][9]*proj3v1  	+ M2[e][10]*proj3v2 	+ M2[e][11]*Da[eq12];

		M2Da[1] = M2[e][12]*Da[eq1] + M2[e][13]*proj1v1 	+ M2[e][14]*proj1v2 	+ M2[e][15]*Da[eq4]
				+ M2[e][16]*Da[eq5] + M2[e][17]*proj2v1 	+ M2[e][18]*proj2v2 	+ M2[e][19]*Da[eq8]
				+ M2[e][20]*Da[eq9] + M2[e][21]*proj3v1 	+ M2[e][22]*proj3v2 	+ M2[e][23]*Da[eq12];

		M2Da[2] = M2[e][24]*Da[eq1] + M2[e][25]*proj1v1 	+ M2[e][26]*proj1v2 	+ M2[e][27]*Da[eq4]
				+ M2[e][28]*Da[eq5] + M2[e][29]*proj2v1 	+ M2[e][30]*proj2v2 	+ M2[e][31]*Da[eq8]
				+ M2[e][32]*Da[eq9] + M2[e][33]*proj3v1 	+ M2[e][34]*proj3v2 	+ M2[e][35]*Da[eq12];

		M2Da[3] = M2[e][36]*Da[eq1] + M2[e][37]*proj1v1 	+ M2[e][38]*proj1v2 	+ M2[e][39]*Da[eq4]
				+ M2[e][40]*Da[eq5] + M2[e][41]*proj2v1 	+ M2[e][42]*proj2v2 	+ M2[e][43]*Da[eq8]
				+ M2[e][44]*Da[eq9] + M2[e][45]*proj3v1 	+ M2[e][46]*proj3v2 	+ M2[e][47]*Da[eq12];*/
		
		M2Da[0] = M2[e][0]*Da[eq1]  + M2[e][1]*Da[eq2]  + M2[e][2]*Da[eq3]  + M2[e][3]*Da[eq4]  
				+ M2[e][4]*Da[eq5]  + M2[e][5]*Da[eq6]  + M2[e][6]*Da[eq7]  + M2[e][7]*Da[eq8]
				+ M2[e][8]*Da[eq9]  + M2[e][9]*Da[eq10]  + M2[e][10]*Da[eq11] + M2[e][11]*Da[eq12];

		M2Da[1] = M2[e][12]*Da[eq1] + M2[e][13]*Da[eq2] + M2[e][14]*Da[eq3] + M2[e][15]*Da[eq4]
				+ M2[e][16]*Da[eq5] + M2[e][17]*Da[eq6] + M2[e][18]*Da[eq7] + M2[e][19]*Da[eq8]
				+ M2[e][20]*Da[eq9] + M2[e][21]*Da[eq10] + M2[e][22]*Da[eq11] + M2[e][23]*Da[eq12];

		M2Da[2] = M2[e][24]*Da[eq1] + M2[e][25]*Da[eq2] + M2[e][26]*Da[eq3] + M2[e][27]*Da[eq4]
				+ M2[e][28]*Da[eq5] + M2[e][29]*Da[eq6] + M2[e][30]*Da[eq7] + M2[e][31]*Da[eq8]
				+ M2[e][32]*Da[eq9] + M2[e][33]*Da[eq10] + M2[e][34]*Da[eq11] + M2[e][35]*Da[eq12];

		M2Da[3] = M2[e][36]*Da[eq1] + M2[e][37]*Da[eq2] + M2[e][38]*Da[eq3] + M2[e][39]*Da[eq4]
				+ M2[e][40]*Da[eq5] + M2[e][41]*Da[eq6] + M2[e][42]*Da[eq7] + M2[e][43]*Da[eq8]
				+ M2[e][44]*Da[eq9] + M2[e][45]*Da[eq10] + M2[e][46]*Da[eq11] + M2[e][47]*Da[eq12];
		
		
		eNDOF = e*NDOF;

		DaB[eNDOF]   = invN2[e]*(R2[e][0] - M2Da[0]); 
		DaB[eNDOF+1] = invN2[e]*(R2[e][1] - M2Da[1]);
		DaB[eNDOF+2] = invN2[e]*(R2[e][2] - M2Da[2]);
		DaB[eNDOF+3] = invN2[e]*(R2[e][3] - M2Da[3]);
		
		//para debug
/*		DaBe[0] = DaB[eNDOF]; 
		DaBe[1] = DaB[eNDOF+1];
		DaBe[2] = DaB[eNDOF+2];
		DaBe[3] = DaB[eNDOF+3];
		
		DaJ1[0] = Da[eq1];
		DaJ1[1] = Da[eq2];
		DaJ1[2] = Da[eq3];
		DaJ1[3] = Da[eq4];

		DaJ2[0] = Da[eq5];
		DaJ2[1] = Da[eq6];
		DaJ2[2] = Da[eq7];
		DaJ2[3] = Da[eq8];
	
		DaJ3[0] = Da[eq9];
		DaJ3[1] = Da[eq10];
		DaJ3[2] = Da[eq11];
		DaJ3[3] = Da[eq12];*/
		
	}

	return 0;
}


