#include "EulerEquations.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"

int uB_InitialSolution(ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions, double *u, double *uB)
{
	int i, e, nel, eNDOF, J1, J2, J3;
	double x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21, twoArea, Area, gamma, deltaNMV, tolerance;
	double third = 1.0/3.0, ninefortieth = 9./40.;
	double Ue[12], dUe[12], Ub[4], dUb[4], gradUx[4], gradUy[4];
	double FB[4]={0.,0.,0.,0.};
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	tolerance = Parameters->StabilizationTolerance;
	nel = Parameters->nel;
	
	double *U = (double*) mycalloc("U of 'Build_M_F_DD_Transiente'", 4*Parameters->nnodes, sizeof(double));
	double *dU = (double*) mycalloc("dU of 'Build_M_F_DD_Transiente'", 4*Parameters->nnodes, sizeof(double));
	eval_U_dU(Parameters,FemStructs,FemFunctions,U,dU);

	for (e = 0; e < nel; e++){
		// *** Global node that composes the element
		J1 = Element[e].Vertex[0];
		J2 = Element[e].Vertex[1];
		J3 = Element[e].Vertex[2];

		// *** Nodal coordinates and differential operator
		x1 = Node[J1].x;
		x2 = Node[J2].x;
		x3 = Node[J3].x;
		
		y1 = Node[J1].y;
		y2 = Node[J2].y;
		y3 = Node[J3].y;
		
		gamma = FemFunctions->gamma(third*(x1+x2+x3), third*(y1+y2+y3));

		y23 = y2 - y3;
		y31 = y3 - y1;
		y12 = y1 - y2;

		x32 = x3 - x2;
		x13 = x1 - x3;
		x21 = x2 - x1;

		// *** Element Area
		twoArea = fabs(x21 * y31 - x13 * y12);
		Area = 0.5*twoArea; 
		
		//Depois desse ponto vc utilizar os valores de theta[0],  theta[1] ou theta[2] convenientimente

		//Density
		Ue[0] = U[4*J1];
		Ue[4] = U[4*J2];
		Ue[8] = U[4*J3];

		dUe[0] = dU[4*J1];
		dUe[4] = dU[4*J2];
		dUe[8] = dU[4*J3];

		//x Velocity
		Ue[1] = U[4*J1+1];
		Ue[5] = U[4*J2+1];
		Ue[9] = U[4*J3+1];

		dUe[1] = dU[4*J1+1];
		dUe[5] = dU[4*J2+1];
		dUe[9] = dU[4*J3+1];

		//y Velocity
		Ue[2] = U[4*J1+2];
		Ue[6] = U[4*J2+2];
		Ue[10] = U[4*J3+2];

		dUe[2] = dU[4*J1+2];
		dUe[6] = dU[4*J2+2];
		dUe[10] = dU[4*J3+2];

		//Energy
		Ue[3] = U[4*J1+3];
		Ue[7] = U[4*J2+3];
		Ue[11] = U[4*J3+3];

		dUe[3] = dU[4*J1+3];
		dUe[7] = dU[4*J2+3];
		dUe[11] = dU[4*J3+3];


		// *** The triangle centroid
		for (i = 0; i < 4; i++){
			Ub[i] = (Ue[i] + Ue[i+4] + Ue[i+8]) * third;
			dUb[i] = (dUe[i] + dUe[i+4] + dUe[i+8]) * third;
		}

		// *** Calculation of Gradient (gradu = Bu)
		for (i = 0; i < 4; i++)
		{
			gradUx[i] = (Ue[i]*y23 + Ue[i+4]*y31 + Ue[i+8]*y12) / twoArea;
			gradUy[i] = (Ue[i]*x32 + Ue[i+4]*x13 + Ue[i+8]*x21) / twoArea;
		}

		
		// *** Jacobian matrices Ax and Ay calculations
		double Ax[4][4];
		double Ay[4][4];
	
		FemFunctions->Ax_Ay_calculations(gamma,Parameters->Mach,Ub, Ax, Ay);

		/*---------------------------------------------ALTERACAO Y^-1)-------------------------------------------------------*/
		double A0[4][4];		
		//------------------------------------------------------------------------------
		//  MONTAGENS DAS MATRIZES
		//------------------------------------------------------------------------------

		//Coeficientes utilizados na matriz Khh e KBh

		// Matriz Ka = y23Ax + x32Ay
		double Ka11, Ka12, Ka13, Ka14, Ka21, Ka22, Ka23, Ka24, Ka31, Ka32, Ka33, Ka34, Ka41, Ka42, Ka43, Ka44;
		Ka11 = y23 * Ax[0][0] + x32 * Ay[0][0];
		Ka12 = y23 * Ax[0][1] + x32 * Ay[0][1];
		Ka13 = y23 * Ax[0][2] + x32 * Ay[0][2];
		Ka14 = y23 * Ax[0][3] + x32 * Ay[0][3];
		Ka21 = y23 * Ax[1][0] + x32 * Ay[1][0];
		Ka22 = y23 * Ax[1][1] + x32 * Ay[1][1];
		Ka23 = y23 * Ax[1][2] + x32 * Ay[1][2];
		Ka24 = y23 * Ax[1][3] + x32 * Ay[1][3];
		Ka31 = y23 * Ax[2][0] + x32 * Ay[2][0];
		Ka32 = y23 * Ax[2][1] + x32 * Ay[2][1];
		Ka33 = y23 * Ax[2][2] + x32 * Ay[2][2];
		Ka34 = y23 * Ax[2][3] + x32 * Ay[2][3];
		Ka41 = y23 * Ax[3][0] + x32 * Ay[3][0];
		Ka42 = y23 * Ax[3][1] + x32 * Ay[3][1];
		Ka43 = y23 * Ax[3][2] + x32 * Ay[3][2];
		Ka44 = y23 * Ax[3][3] + x32 * Ay[3][3];

		// Matriz Kb = y31Ax + x13Ay
		double Kb11, Kb12, Kb13, Kb14, Kb21, Kb22, Kb23, Kb24, Kb31, Kb32, Kb33, Kb34, Kb41, Kb42, Kb43, Kb44;
		Kb11 = y31 * Ax[0][0] + x13 * Ay[0][0];
		Kb12 = y31 * Ax[0][1] + x13 * Ay[0][1];
		Kb13 = y31 * Ax[0][2] + x13 * Ay[0][2];
		Kb14 = y31 * Ax[0][3] + x13 * Ay[0][3];
		Kb21 = y31 * Ax[1][0] + x13 * Ay[1][0];
		Kb22 = y31 * Ax[1][1] + x13 * Ay[1][1];
		Kb23 = y31 * Ax[1][2] + x13 * Ay[1][2];
		Kb24 = y31 * Ax[1][3] + x13 * Ay[1][3];
		Kb31 = y31 * Ax[2][0] + x13 * Ay[2][0];
		Kb32 = y31 * Ax[2][1] + x13 * Ay[2][1];
		Kb33 = y31 * Ax[2][2] + x13 * Ay[2][2];
		Kb34 = y31 * Ax[2][3] + x13 * Ay[2][3];
		Kb41 = y31 * Ax[3][0] + x13 * Ay[3][0];
		Kb42 = y31 * Ax[3][1] + x13 * Ay[3][1];
		Kb43 = y31 * Ax[3][2] + x13 * Ay[3][2];
		Kb44 = y31 * Ax[3][3] + x13 * Ay[3][3];

		// Matriz Kc = y12Ax + x21Ay
		double Kc11, Kc12, Kc13, Kc14, Kc21, Kc22, Kc23, Kc24, Kc31, Kc32, Kc33, Kc34, Kc41, Kc42, Kc43, Kc44;
		Kc11 = y12 * Ax[0][0] + x21 * Ay[0][0];
		Kc12 = y12 * Ax[0][1] + x21 * Ay[0][1];
		Kc13 = y12 * Ax[0][2] + x21 * Ay[0][2];
		Kc14 = y12 * Ax[0][3] + x21 * Ay[0][3];
		Kc21 = y12 * Ax[1][0] + x21 * Ay[1][0];
		Kc22 = y12 * Ax[1][1] + x21 * Ay[1][1];
		Kc23 = y12 * Ax[1][2] + x21 * Ay[1][2];
		Kc24 = y12 * Ax[1][3] + x21 * Ay[1][3];
		Kc31 = y12 * Ax[2][0] + x21 * Ay[2][0];
		Kc32 = y12 * Ax[2][1] + x21 * Ay[2][1];
		Kc33 = y12 * Ax[2][2] + x21 * Ay[2][2];
		Kc34 = y12 * Ax[2][3] + x21 * Ay[2][3];
		Kc41 = y12 * Ax[3][0] + x21 * Ay[3][0];
		Kc42 = y12 * Ax[3][1] + x21 * Ay[3][1];
		Kc43 = y12 * Ax[3][2] + x21 * Ay[3][2];
		Kc44 = y12 * Ax[3][3] + x21 * Ay[3][3];

		// *** Matriz de Rigidez KBh 4x12
		double KBh11, KBh12, KBh13, KBh14, KBh21, KBh22, KBh23, KBh24, KBh31, KBh32, KBh33, KBh34, KBh41, KBh42, KBh43, KBh44;
		KBh11 = ninefortieth * Ka11;
		KBh12 = ninefortieth * Ka12;
		KBh13 = ninefortieth * Ka13;
		KBh14 = ninefortieth * Ka14;
		KBh21 = ninefortieth * Ka21;
		KBh22 = ninefortieth * Ka22;
		KBh23 = ninefortieth * Ka23;
		KBh24 = ninefortieth * Ka24;
		KBh31 = ninefortieth * Ka31;
		KBh32 = ninefortieth * Ka32;
		KBh33 = ninefortieth * Ka33;
		KBh34 = ninefortieth * Ka34;
		KBh41 = ninefortieth * Ka41;
		KBh42 = ninefortieth * Ka42;
		KBh43 = ninefortieth * Ka43;
		KBh44 = ninefortieth * Ka44;
		
		double KBh15, KBh16, KBh17, KBh18, KBh25, KBh26, KBh27, KBh28, KBh35, KBh36, KBh37, KBh38, KBh45, KBh46, KBh47, KBh48;
		KBh15 = ninefortieth * Kb11;
		KBh16 = ninefortieth * Kb12;
		KBh17 = ninefortieth * Kb13;
		KBh18 = ninefortieth * Kb14;
		KBh25 = ninefortieth * Kb21;
		KBh26 = ninefortieth * Kb22;
		KBh27 = ninefortieth * Kb23;
		KBh28 = ninefortieth * Kb24;
		KBh35 = ninefortieth * Kb31;
		KBh36 = ninefortieth * Kb32;
		KBh37 = ninefortieth * Kb33;
		KBh38 = ninefortieth * Kb34;
		KBh45 = ninefortieth * Kb41;
		KBh46 = ninefortieth * Kb42;
		KBh47 = ninefortieth * Kb43;
		KBh48 = ninefortieth * Kb44;
		
		double KBh19, KBh110, KBh111, KBh112, KBh29, KBh210, KBh211, KBh212, KBh39, KBh310, KBh311, KBh312, KBh49, KBh410, KBh411, KBh412;
		KBh19 = ninefortieth * Kc11;
		KBh110 = ninefortieth * Kc12;
		KBh111 = ninefortieth * Kc13;
		KBh112 = ninefortieth * Kc14;
		KBh29 = ninefortieth * Kc21;
		KBh210 = ninefortieth * Kc22;
		KBh211 = ninefortieth * Kc23;
		KBh212 = ninefortieth * Kc24;
		KBh39 = ninefortieth * Kc31;
		KBh310 = ninefortieth * Kc32;
		KBh311 = ninefortieth * Kc33;
		KBh312 = ninefortieth * Kc34;
		KBh49 = ninefortieth * Kc41;
		KBh410 = ninefortieth * Kc42;
		KBh411 = ninefortieth * Kc43;
		KBh412 = ninefortieth * Kc44;
		

		deltaNMV = Delta_YZBetaNMV(tolerance, FemStructs->delta_old, gradUx, gradUy, Ax, Ay, A0, dUb, y23, y31, y12, x32, x13, x21, twoArea, e, Parameters->invY, Ub);
		// *** Matriz de Rigidez KBB 4x4
		double KBB, inv_KBB;
		
		//delta = delta/Cstab;
		KBB = (81.0*deltaNMV*(y23*y23 + y31*y31 + y12*y12 + y23*y31 + y23*y12 + y31*y12 + x32*x32 + x13*x13 + x21*x21 + x32*x13 + x32*x21 + x13*x21))/(40.0*Area);
					
				
		// Agora a parte resolvida: R2 = R2 - MBh*du - KBh*u;
		double KBhUe[4];
		KBhUe[0]  = KBh11*Ue[0] + KBh12*Ue[1] + KBh13*Ue[2] + KBh14*Ue[3] + KBh15*Ue[4] + KBh16*Ue[5] 
			  + KBh17*Ue[6] + KBh18*Ue[7] + KBh19*Ue[8] + KBh110*Ue[9] + KBh111*Ue[10] + KBh112*Ue[11];

		KBhUe[1]  = KBh21*Ue[0] + KBh22*Ue[1] + KBh23*Ue[2] + KBh24*Ue[3] + KBh25*Ue[4] + KBh26*Ue[5] 
			  + KBh27*Ue[6] + KBh28*Ue[7] + KBh29*Ue[8] + KBh210*Ue[9] + KBh211*Ue[10] + KBh212*Ue[11];

		KBhUe[2]  = KBh31*Ue[0] + KBh32*Ue[1] + KBh33*Ue[2] + KBh34*Ue[3] + KBh35*Ue[4] + KBh36*Ue[5] 
                          + KBh37*Ue[6] + KBh38*Ue[7] + KBh39*Ue[8] + KBh310*Ue[9] + KBh311*Ue[10] + KBh312*Ue[11];

		KBhUe[3]  = KBh41*Ue[0] + KBh42*Ue[1] + KBh43*Ue[2] + KBh44*Ue[3] + KBh45*Ue[4] + KBh46*Ue[5] 
                          + KBh47*Ue[6] + KBh48*Ue[7] + KBh49*Ue[8] + KBh410*Ue[9] + KBh411*Ue[10] + KBh412*Ue[11];
	
		// Valor de uB e duB no elemento
		eNDOF = e*NDOF;
		if (fabs(KBB)<1e-10){
			uB[eNDOF]     =	0;
			uB[eNDOF + 1] = 0;
			uB[eNDOF + 1] = 0;
			uB[eNDOF + 1] = 0;
		}
		else{
			inv_KBB = 1.0/KBB;
			uB[eNDOF]     =	inv_KBB*(FB[0] - KBhUe[0]);
			uB[eNDOF + 1] = inv_KBB*(FB[1] - KBhUe[1]);
			uB[eNDOF + 1] = inv_KBB*(FB[1] - KBhUe[1]);
			uB[eNDOF + 1] = inv_KBB*(FB[1] - KBhUe[1]);
		}
	}//for elemento
	
		
	// Liberacao dos espacos alocados
	free(U);
	free(dU);
	

	return 0;
}// end build



