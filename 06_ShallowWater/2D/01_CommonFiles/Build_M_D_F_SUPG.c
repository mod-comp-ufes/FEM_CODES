#include "ShalowWater.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"


int Build_M_D_F_SUPG(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	double mu = 0.01;
	double rho = 1000;
	double g = 9.8;
	double n = 0.018;
	double Cf, gamma;
	int e, i, j, r;
	int J1, J2, J3, nel, neq;
	double x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21, zb1, zb2, zb3;
	double Area, twoArea, invArea, second = 0.5, third = 1.0/3.0, forth = 0.25, sixth = 1.0/6.0, twelve = 1.0/12.0;
	double u[3], v[3];
	double tau, delta;
	double *U, *dU;
	double Me[9][9], De[9][9], Ue[9], dUe[9], Ub[3], dUb[3], Fb[3], gradUx[3], gradUy[3];
	double S1[3], S2[3], S3[3], Sb[3], Fe[9], gradzbx, gradzby;
	double alpha = Parameters->Alpha_Build;
	double delta_t = Parameters->DeltaT_Build;
	double *R = FemStructs->F;
	int **lm = FemStructs->lm;

	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	
	nel = Parameters->nel;
	neq = Parameters->neq;
	
	dzero(neq+1, R);
	setzeros(Parameters, MatrixData);

	U = (double*) mycalloc("U of 'Build_M_F_SUPG_Transiente'", 3*Parameters->nnodes, sizeof(double));
	dU = (double*) mycalloc("dU of 'Build_M_F_SUPG_Transiente'", 3*Parameters->nnodes, sizeof(double));

	eval_U_dU(Parameters,FemStructs,FemFunctions,U,dU);

	for (e = 0; e < nel; e++){
		// *** Nos do elemento
		J1 = Element[e].Vertex[0];
		J2 = Element[e].Vertex[1];
		J3 = Element[e].Vertex[2];

		x1 = Node[J1].x;
		x2 = Node[J2].x;
		x3 = Node[J3].x;
		
		y1 = Node[J1].y;
		y2 = Node[J2].y;
		y3 = Node[J3].y;

		y23 = y2 - y3;
		y31 = y3 - y1;
		y12 = y1 - y2;

		x32 = x3 - x2;
		x13 = x1 - x3;
		x21 = x2 - x1;

		// *** Area do elemento
		twoArea = fabs(x21 * y31 - x13 * y12);
		Area = 0.5*twoArea;
		invArea = 1.0/Area;

		// Fill 'Ue' with the values ​​of the previous solution or workaround for each degree of liberty of the element node
		// height
		Ue[0] = U[3*J1];
		Ue[3] = U[3*J2];
		Ue[6] = U[3*J3];

		dUe[0] = dU[3*J1];
		dUe[3] = dU[3*J2];
		dUe[6] = dU[3*J3];

		// discharge in X
		Ue[1] = U[3*J1+1];
		Ue[4] = U[3*J2+1];
		Ue[7] = U[3*J3+1];

		dUe[1] = dU[3*J1+1];
		dUe[4] = dU[3*J2+1];
		dUe[7] = dU[3*J3+1];

		// discharge in Y
		Ue[2] = U[3*J1+2];
		Ue[5] = U[3*J2+2];
		Ue[8] = U[3*J3+2];

		dUe[2] = dU[3*J1+2];
		dUe[5] = dU[3*J2+2];
		dUe[8] = dU[3*J3+2];

		// velocities
		u[0] = Ue[1]/Ue[0];
		u[1] = Ue[4]/Ue[3];
		u[2] = Ue[7]/Ue[6];

		v[0] = Ue[2]/Ue[0];
		v[1] = Ue[5]/Ue[3];
		v[2] = Ue[8]/Ue[6];

		// Baricentro do triangulo
		for (i = 0; i < 3; i++){
			Ub[i] = (Ue[i] + Ue[i+3] + Ue[i+6]) * third;
			dUb[i] = (dUe[i] + dUe[i+3] + dUe[i+6]) * third;
		}

		// Calculo do Gradiente
		for (i = 0; i < 3; i++)
		{
			gradUx[i] = (Ue[i]*y23 + Ue[i+3]*y31 + Ue[i+6]*y12) / twoArea;
			gradUy[i] = (Ue[i]*x32 + Ue[i+3]*x13 + Ue[i+6]*x21) / twoArea;
		}

		// *** A = [A1 A2]
		double A1[3][3];
		double A2[3][3];
		A1_A2_calculations(Ue, A1, A2, g);
		
		double a11[3][3], a12[3][3], a21[3][3], a22[3][3];

		a11[0][0] = A1[1][0];
		a11[0][1] = A1[1][1];
		a11[0][2] = 0.0;
		a11[1][0] = A1[1][0]*A1[1][0] + A1[1][1]*A1[2][0];
		a11[1][1] = A1[1][0] + A1[1][1]*A1[1][1];
		a11[1][2] = 0.0;
		a11[2][0] = A1[2][1]*A1[1][0] + A1[2][2]*A1[2][0];
		a11[2][1] = A1[2][0] + A1[2][1]*A1[1][1] + A1[2][2]*A1[2][1];
		a11[2][2] = A1[2][2]*A1[2][2];

		a12[0][0] = A2[1][0];
		a12[0][1] = A2[1][1];
		a12[0][2] = A2[1][2];
		a12[1][0] = A1[1][1]*A2[1][0];
		a12[1][1] = A1[1][1]*A2[1][1];
		a12[1][2] = A1[1][0] + A1[1][1]*A2[1][2];
		a12[2][0] = A1[2][1]*A2[1][0] + A1[2][2]*A2[2][0];
		a12[2][1] = A1[2][1]*A2[1][1];
		a12[2][2] = A1[2][0] + A1[2][1]*A2[1][2] + A1[2][2]*A2[2][2];

		a21[0][0] = A1[2][0];
		a21[0][1] = A1[2][1];
		a21[0][2] = A1[2][2];
		a21[1][0] = A2[1][1]*A1[1][0] + A2[1][2]*A1[2][0];
		a21[1][1] = A2[1][0] + A2[1][1]*A1[1][1] + A2[2][1]*A1[2][1];
		a21[1][2] = A2[1][2]*A2[2][2];
		a21[2][0] = A2[2][2]*A1[2][0];
		a21[2][1] = A2[2][0] + A2[2][2]*A1[2][1];
		a21[2][2] = A2[2][2]*A2[2][2];

		a22[0][0] = A2[2][0];
		a22[0][1] = 0.0;
		a22[0][2] = A2[2][2];
		a22[1][0] = A2[1][1]*A2[1][0] + A2[1][2]*A2[2][0];
		a22[1][1] = A2[1][1]*A2[1][1];
		a22[1][2] = A2[1][0] + A2[1][1]*A2[2][0] + A2[2][1]*A2[2][2];
		a22[2][0] = A2[2][0];
		a22[2][1] = 0.0;
		a22[2][2] = A2[2][0] + A2[2][2]*A2[2][2];

		// Source

		zb1 = FemFunctions->zb(x1, y1);
		zb2 = FemFunctions->zb(x2, y2);
		zb3 = FemFunctions->zb(x3, y3);

		gradzbx = (zb1*y23 + zb2*y31 + zb3*y12)/(2*Area);
		gradzby = (zb1*x32 + zb2*x13 + zb3*x21)/(2*Area);

		Cf = (Ub[0]>=0) ? g*n*n*pow(Ub[0], -7/3) : 0.0;
		r = Ub[1]*Ub[1] + Ub[2]*Ub[2];
		gamma = (r>=0) ? Cf*sqrt(r) : 0.0;
		
		S1[0] = 0.0;
		S2[0] = 0.0;
		S3[0] = 0.0;

		S1[1] = -g*gradzbx*Ue[0] - gamma*Ue[1];
		S2[1] = -g*gradzbx*Ue[3] - gamma*Ue[4];
		S3[1] = -g*gradzbx*Ue[6] - gamma*Ue[7];

		S1[2] = -g*gradzby*Ue[0] - gamma*Ue[2];
		S2[2] = -g*gradzby*Ue[3] - gamma*Ue[5];
		S3[2] = -g*gradzby*Ue[6] - gamma*Ue[8];

		Sb[0] = third*(S1[0] + S2[0] + S3[0]);
		Sb[1] = third*(S1[1] + S2[2] + S3[1]);
		Sb[2] = third*(S1[2] + S2[2] + S3[2]);

		// *** Calculo do parametro de estabilizacao (tau) do SUPG            
		double tauSUGN1 = 0.0, tauSUGN2 = 0.0;
		double jbold[2], normjbold = 0.0;

		jbold[0] = second*invArea*(y23*(Ue[0] + zb1) + y31*(Ue[3] + zb2) + y12*(Ue[6] + zb3));
		jbold[1] = second*invArea*(x32*(Ue[0] + zb1) + x13*(Ue[3] + zb2) + x21*(Ue[6] + zb3));
		r = jbold[0]*jbold[0] + jbold[1]*jbold[1];
		normjbold = (r>=0) ? sqrt(r) : 0.0;
		jbold[0] /= normjbold;
		jbold[1] /= normjbold;

		r = g*Ub[0];
		
		tauSUGN1 = (r>=0) ? 1.0/(second*invArea*sqrt(r)*(abs(jbold[0]*y23 + jbold[1]*x32) + abs(Ub[1]*y23 + Ub[2]*x32) +
		                                                 abs(jbold[0]*y31 + jbold[1]*x13) + abs(Ub[1]*y31 + Ub[2]*x13) +
												         abs(jbold[0]*y12 + jbold[1]*x21) + abs(Ub[2]*y12 + Ub[2]*x21))) : 0.0;
		tauSUGN2 = delta_t/2.0;

		tau = (tauSUGN1>=0 && tauSUGN2>=0) ? 1.0/sqrt(1.0/(tauSUGN1*tauSUGN1) + 1.0/(tauSUGN2*tauSUGN2)) : 0.0;

		delta = FemFunctions->ShockCapture(Ub, gradUx, gradUy, A1, A2, Sb, y23, y31, y12, x32, x13, x21, twoArea, tau, g);

		//------------------------------------------------------------------------------
		//  MONTAGENS DAS MATRIZES
		//------------------------------------------------------------------------------

		// *** Matriz de Massa do Galerkin
		double Mg1, Mg2;
		Mg1 = Area * sixth;
		Mg2 = Area * twelve;

		// *** Matriz de Massa do SUPG 
		double ms1[3][3], ms2[3][3], ms3[3][3];

		ms1[0][0] = 0.0;
		ms1[0][1] = y23;
		ms1[0][2] =                x32;
		ms1[1][0] = y23*A1[1][0] + x32*A2[1][0];
		ms1[1][1] = y23*A1[1][1] + x32*A2[1][1];
		ms1[1][2] =                x32*A2[1][2];
		ms1[2][0] = y23*A1[2][0] + x32*A2[2][0];
		ms1[2][1] = y23*A1[2][1];
		ms1[2][2] = y23*A1[2][2] + x32*A2[2][2];

		ms2[0][0] = 0.0;
		ms2[0][1] = y31;
		ms2[0][2] =                x13;
		ms2[1][0] = y31*A1[1][0] + x13*A2[1][0];
		ms2[1][1] = y31*A1[1][1] + x13*A2[1][1];
		ms2[1][2] =                x13*A2[1][2];
		ms2[2][0] = y31*A1[2][0] + x13*A2[2][0];
		ms2[2][1] = y31*A1[2][1];
		ms2[2][2] = y31*A1[2][2] + x13*A2[2][2];

		ms3[0][0] = 0.0;
		ms3[0][1] = y12;
		ms3[0][2] =                x21;
		ms3[1][0] = y12*A1[1][0] + x21*A2[1][0];
		ms3[1][1] = y12*A1[1][1] + x21*A2[1][1];
		ms3[1][2] =                x21*A2[1][2];
		ms3[2][0] = y12*A1[2][0] + x21*A2[2][0];
		ms3[2][1] = y12*A1[2][1];
		ms3[2][2] = y12*A1[2][2] + x21*A2[2][2];


		// *** Coeficientes da Matriz de Massa [Me]

		Me[0][0] = Mg1 + tau*sixth*ms1[0][0];
		Me[0][1] = Mg1 + tau*sixth*ms1[0][1];
		Me[0][2] = Mg1 + tau*sixth*ms1[0][2];
		Me[1][0] = Mg1 + tau*sixth*ms1[1][0];
		Me[1][1] = Mg1 + tau*sixth*ms1[1][1];
		Me[1][2] = Mg1 + tau*sixth*ms1[1][2];
		Me[2][0] = Mg1 + tau*sixth*ms1[2][0];
		Me[2][1] = Mg1 + tau*sixth*ms1[2][1];
		Me[2][2] = Mg1 + tau*sixth*ms1[2][2];

		Me[3][3] = Mg1 + tau*sixth*ms2[0][0];
		Me[3][4] = Mg1 + tau*sixth*ms2[0][1];
		Me[3][5] = Mg1 + tau*sixth*ms2[0][2];
		Me[4][3] = Mg1 + tau*sixth*ms2[1][0];
		Me[4][4] = Mg1 + tau*sixth*ms2[1][1];
		Me[4][5] = Mg1 + tau*sixth*ms2[1][2];
		Me[5][3] = Mg1 + tau*sixth*ms2[2][0];
		Me[5][4] = Mg1 + tau*sixth*ms2[2][1];
		Me[5][5] = Mg1 + tau*sixth*ms2[2][2];

		Me[6][6] = Mg1 + tau*sixth*ms3[0][0];
		Me[6][7] = Mg1 + tau*sixth*ms3[0][1];
		Me[6][8] = Mg1 + tau*sixth*ms3[0][2];
		Me[7][6] = Mg1 + tau*sixth*ms3[1][0];
		Me[7][7] = Mg1 + tau*sixth*ms3[1][1];
		Me[7][8] = Mg1 + tau*sixth*ms3[1][2];
		Me[8][6] = Mg1 + tau*sixth*ms3[2][0];
		Me[8][7] = Mg1 + tau*sixth*ms2[2][1];
		Me[8][8] = Mg1 + tau*sixth*ms3[2][2];

		Me[0][6] = Me[0][3] = Mg2 + tau*sixth*ms1[0][0];
		Me[0][7] = Me[0][4] = Mg2 + tau*sixth*ms1[0][1];
		Me[0][8] = Me[0][5] = Mg2 + tau*sixth*ms1[0][2];
		Me[1][6] = Me[1][3] = Mg2 + tau*sixth*ms1[1][0];
		Me[1][7] = Me[1][4] = Mg2 + tau*sixth*ms1[1][1];
		Me[1][8] = Me[1][5] = Mg2 + tau*sixth*ms1[1][2];
		Me[2][6] = Me[2][3] = Mg2 + tau*sixth*ms1[2][0];
		Me[2][7] = Me[2][4] = Mg2 + tau*sixth*ms1[2][1];
		Me[2][8] = Me[2][5] = Mg2 + tau*sixth*ms1[2][2];

		Me[3][6] = Me[3][0] = Mg2 + tau*sixth*ms2[0][0];
		Me[3][7] = Me[3][1] = Mg2 + tau*sixth*ms2[0][1];
		Me[3][8] = Me[3][2] = Mg2 + tau*sixth*ms2[0][2];
		Me[4][6] = Me[4][0] = Mg2 + tau*sixth*ms2[1][0];
		Me[4][7] = Me[4][1] = Mg2 + tau*sixth*ms2[1][1];
		Me[4][8] = Me[4][2] = Mg2 + tau*sixth*ms2[1][2];
		Me[5][6] = Me[5][0] = Mg2 + tau*sixth*ms2[2][0];
		Me[5][7] = Me[5][1] = Mg2 + tau*sixth*ms2[2][1];
		Me[5][8] = Me[5][2] = Mg2 + tau*sixth*ms2[2][2];

		Me[6][3] = Me[6][0] = Mg2 + tau*sixth*ms3[0][0];
		Me[6][4] = Me[6][1] = Mg2 + tau*sixth*ms3[0][1];
		Me[6][5] = Me[6][2] = Mg2 + tau*sixth*ms3[0][2];
		Me[7][3] = Me[7][0] = Mg2 + tau*sixth*ms3[1][0];
		Me[7][4] = Me[7][1] = Mg2 + tau*sixth*ms3[1][1];
		Me[7][5] = Me[7][2] = Mg2 + tau*sixth*ms3[1][2];
		Me[8][3] = Me[8][0] = Mg2 + tau*sixth*ms3[2][0];
		Me[8][4] = Me[8][1] = Mg2 + tau*sixth*ms3[2][1];
		Me[8][5] = Me[8][2] = Mg2 + tau*sixth*ms3[2][2];

		// *** Matriz de Conveccao de Galerkin
		double d12, d13, d23;

		d12 = forth*invArea*(y23*y31*mu/rho + x32*x13*mu/rho);
		d13 = forth*invArea*(y23*y12*mu/rho + x32*x21*mu/rho);
		d23 = forth*invArea*(y31*y12*mu/rho + x13*x21*mu/rho);

		// *** Matriz de Conveccao do SUPG
		double ds11[3][3], ds12[3][3], ds13[3][3], ds21[3][3], ds22[3][3], ds23[3][3], ds31[3][3], ds32[3][3], ds33[3][3];

		ds11[0][0] = y23*(y23*a11[0][0] + x32*a21[0][0]) + x32*(y23*a12[0][0] + x32*a22[0][0]);
		ds11[0][1] = y23*(y23*a11[0][1] + x32*a21[0][1]) + x32*(y23*a12[0][1] + x32*a22[0][1]);
		ds11[0][2] = y23*(y23*a11[0][2] + x32*a21[0][2]) + x32*(y23*a12[0][2] + x32*a22[0][2]);
		ds11[1][0] = y23*(y23*a11[1][0] + x32*a21[1][0]) + x32*(y23*a12[1][0] + x32*a22[1][0]);
		ds11[1][1] = y23*(y23*a11[1][1] + x32*a21[1][1]) + x32*(y23*a12[1][1] + x32*a22[1][1]);
		ds11[1][2] = y23*(y23*a11[1][2] + x32*a21[1][2]) + x32*(y23*a12[1][2] + x32*a22[1][2]);
		ds11[2][0] = y23*(y23*a11[2][0] + x32*a21[2][0]) + x32*(y23*a12[2][0] + x32*a22[2][0]);
		ds11[2][1] = y23*(y23*a11[2][1] + x32*a21[2][1]) + x32*(y23*a12[2][1] + x32*a22[2][1]);
		ds11[2][2] = y23*(y23*a11[2][2] + x32*a21[2][2]) + x32*(y23*a12[2][2] + x32*a22[2][2]);

		ds12[0][0] = y31*(y23*a11[0][0] + x32*a21[0][0]) + x13*(y23*a12[0][0] + x32*a22[0][0]);
		ds12[0][1] = y31*(y23*a11[0][1] + x32*a21[0][1]) + x13*(y23*a12[0][1] + x32*a22[0][1]);
		ds12[0][2] = y31*(y23*a11[0][2] + x32*a21[0][2]) + x13*(y23*a12[0][2] + x32*a22[0][2]);
		ds12[1][0] = y31*(y23*a11[1][0] + x32*a21[1][0]) + x13*(y23*a12[1][0] + x32*a22[1][0]);
		ds12[1][1] = y31*(y23*a11[1][1] + x32*a21[1][1]) + x13*(y23*a12[1][1] + x32*a22[1][1]);
		ds12[1][2] = y31*(y23*a11[1][2] + x32*a21[1][2]) + x13*(y23*a12[1][2] + x32*a22[1][2]);
		ds12[2][0] = y31*(y23*a11[2][0] + x32*a21[2][0]) + x13*(y23*a12[2][0] + x32*a22[2][0]);
		ds12[2][1] = y31*(y23*a11[2][1] + x32*a21[2][1]) + x13*(y23*a12[2][1] + x32*a22[2][1]);
		ds12[2][2] = y31*(y23*a11[2][2] + x32*a21[2][2]) + x13*(y23*a12[2][2] + x32*a22[2][2]);

		ds13[0][0] = y12*(y23*a11[0][0] + x32*a21[0][0]) + x21*(y23*a12[0][0] + x32*a22[0][0]);
		ds13[0][1] = y12*(y23*a11[0][1] + x32*a21[0][1]) + x21*(y23*a12[0][1] + x32*a22[0][1]);
		ds13[0][2] = y12*(y23*a11[0][2] + x32*a21[0][2]) + x21*(y23*a12[0][2] + x32*a22[0][2]);
		ds13[1][0] = y12*(y23*a11[1][0] + x32*a21[1][0]) + x21*(y23*a12[1][0] + x32*a22[1][0]);
		ds13[1][1] = y12*(y23*a11[1][1] + x32*a21[1][1]) + x21*(y23*a12[1][1] + x32*a22[1][1]);
		ds13[1][2] = y12*(y23*a11[1][2] + x32*a21[1][2]) + x21*(y23*a12[1][2] + x32*a22[1][2]);
		ds13[2][0] = y12*(y23*a11[2][0] + x32*a21[2][0]) + x21*(y23*a12[2][0] + x32*a22[2][0]);
		ds13[2][1] = y12*(y23*a11[2][1] + x32*a21[2][1]) + x21*(y23*a12[2][1] + x32*a22[2][1]);
		ds13[2][2] = y12*(y23*a11[2][2] + x32*a21[2][2]) + x21*(y23*a12[2][2] + x32*a22[2][2]);

		ds21[0][0] = y23*(y31*a11[0][0] + x13*a21[0][0]) + x32*(y31*a12[0][0] + x13*a22[0][0]);
		ds21[0][1] = y23*(y31*a11[0][1] + x13*a21[0][1]) + x32*(y31*a12[0][1] + x13*a22[0][1]);
		ds21[0][2] = y23*(y31*a11[0][2] + x13*a21[0][2]) + x32*(y31*a12[0][2] + x13*a22[0][2]);
		ds21[1][0] = y23*(y31*a11[1][0] + x13*a21[1][0]) + x32*(y31*a12[1][0] + x13*a22[1][0]);
		ds21[1][1] = y23*(y31*a11[1][1] + x13*a21[1][1]) + x32*(y31*a12[1][1] + x13*a22[1][1]);
		ds21[1][2] = y23*(y31*a11[1][2] + x13*a21[1][2]) + x32*(y31*a12[1][2] + x13*a22[1][2]);
		ds21[2][0] = y23*(y31*a11[2][0] + x13*a21[2][0]) + x32*(y31*a12[2][0] + x13*a22[2][0]);
		ds21[2][1] = y23*(y31*a11[2][1] + x13*a21[2][1]) + x32*(y31*a12[2][1] + x13*a22[2][1]);
		ds21[2][2] = y23*(y31*a11[2][2] + x13*a21[2][2]) + x32*(y31*a12[2][2] + x13*a22[2][2]);

		ds22[0][0] = y31*(y31*a11[0][0] + x13*a21[0][0]) + x13*(y31*a12[0][0] + x13*a22[0][0]);
		ds22[0][1] = y31*(y31*a11[0][1] + x13*a21[0][1]) + x13*(y31*a12[0][1] + x13*a22[0][1]);
		ds22[0][2] = y31*(y31*a11[0][2] + x13*a21[0][2]) + x13*(y31*a12[0][2] + x13*a22[0][2]);
		ds22[1][0] = y31*(y31*a11[1][0] + x13*a21[1][0]) + x13*(y31*a12[1][0] + x13*a22[1][0]);
		ds22[1][1] = y31*(y31*a11[1][1] + x13*a21[1][1]) + x13*(y31*a12[1][1] + x13*a22[1][1]);
		ds22[1][2] = y31*(y31*a11[1][2] + x13*a21[1][2]) + x13*(y31*a12[1][2] + x13*a22[1][2]);
		ds22[2][0] = y31*(y31*a11[2][0] + x13*a21[2][0]) + x13*(y31*a12[2][0] + x13*a22[2][0]);
		ds22[2][1] = y31*(y31*a11[2][1] + x13*a21[2][1]) + x13*(y31*a12[2][1] + x13*a22[2][1]);
		ds22[2][2] = y31*(y31*a11[2][2] + x13*a21[2][2]) + x13*(y31*a12[2][2] + x13*a22[2][2]);

		ds23[0][0] = y12*(y31*a11[0][0] + x13*a21[0][0]) + x21*(y31*a12[0][0] + x13*a22[0][0]);
		ds23[0][1] = y12*(y31*a11[0][1] + x13*a21[0][1]) + x21*(y31*a12[0][1] + x13*a22[0][1]);
		ds23[0][2] = y12*(y31*a11[0][2] + x13*a21[0][2]) + x21*(y31*a12[0][2] + x13*a22[0][2]);
		ds23[1][0] = y12*(y31*a11[1][0] + x13*a21[1][0]) + x21*(y31*a12[1][0] + x13*a22[1][0]);
		ds23[1][1] = y12*(y31*a11[1][1] + x13*a21[1][1]) + x21*(y31*a12[1][1] + x13*a22[1][1]);
		ds23[1][2] = y12*(y31*a11[1][2] + x13*a21[1][2]) + x21*(y31*a12[1][2] + x13*a22[1][2]);
		ds23[2][0] = y12*(y31*a11[2][0] + x13*a21[2][0]) + x21*(y31*a12[2][0] + x13*a22[2][0]);
		ds23[2][1] = y12*(y31*a11[2][1] + x13*a21[2][1]) + x21*(y31*a12[2][1] + x13*a22[2][1]);
		ds23[2][2] = y12*(y31*a11[2][2] + x13*a21[2][2]) + x21*(y31*a12[2][2] + x13*a22[2][2]);

		ds31[0][0] = y23*(y12*a11[0][0] + x21*a21[0][0]) + x32*(y12*a12[0][0] + x21*a22[0][0]);
		ds31[0][1] = y23*(y12*a11[0][1] + x21*a21[0][1]) + x32*(y12*a12[0][1] + x21*a22[0][1]);
		ds31[0][2] = y23*(y12*a11[0][2] + x21*a21[0][2]) + x32*(y12*a12[0][2] + x21*a22[0][2]);
		ds31[1][0] = y23*(y12*a11[1][0] + x21*a21[1][0]) + x32*(y12*a12[1][0] + x21*a22[1][0]);
		ds31[1][1] = y23*(y12*a11[1][1] + x21*a21[1][1]) + x32*(y12*a12[1][1] + x21*a22[1][1]);
		ds31[1][2] = y23*(y12*a11[1][2] + x21*a21[1][2]) + x32*(y12*a12[1][2] + x21*a22[1][2]);
		ds31[2][0] = y23*(y12*a11[2][0] + x21*a21[2][0]) + x32*(y12*a12[2][0] + x21*a22[2][0]);
		ds31[2][1] = y23*(y12*a11[2][1] + x21*a21[2][1]) + x32*(y12*a12[2][1] + x21*a22[2][1]);
		ds31[2][2] = y23*(y12*a11[2][2] + x21*a21[2][2]) + x32*(y12*a12[2][2] + x21*a22[2][2]);

		ds32[0][0] = y31*(y12*a11[0][0] + x21*a21[0][0]) + x13*(y12*a12[0][0] + x21*a22[0][0]);
		ds32[0][1] = y31*(y12*a11[0][1] + x21*a21[0][1]) + x13*(y12*a12[0][1] + x21*a22[0][1]);
		ds32[0][2] = y31*(y12*a11[0][2] + x21*a21[0][2]) + x13*(y12*a12[0][2] + x21*a22[0][2]);
		ds32[1][0] = y31*(y12*a11[1][0] + x21*a21[1][0]) + x13*(y12*a12[1][0] + x21*a22[1][0]);
		ds32[1][1] = y31*(y12*a11[1][1] + x21*a21[1][1]) + x13*(y12*a12[1][1] + x21*a22[1][1]);
		ds32[1][2] = y31*(y12*a11[1][2] + x21*a21[1][2]) + x13*(y12*a12[1][2] + x21*a22[1][2]);
		ds32[2][0] = y31*(y12*a11[2][0] + x21*a21[2][0]) + x13*(y12*a12[2][0] + x21*a22[2][0]);
		ds32[2][1] = y31*(y12*a11[2][1] + x21*a21[2][1]) + x13*(y12*a12[2][1] + x21*a22[2][1]);
		ds32[2][2] = y31*(y12*a11[2][2] + x21*a21[2][2]) + x13*(y12*a12[2][2] + x21*a22[2][2]);

		ds33[0][0] = y12*(y12*a11[0][0] + x21*a21[0][0]) + x21*(y12*a12[0][0] + x21*a22[0][0]);
		ds33[0][1] = y12*(y12*a11[0][1] + x21*a21[0][1]) + x21*(y12*a12[0][1] + x21*a22[0][1]);
		ds33[0][2] = y12*(y12*a11[0][2] + x21*a21[0][2]) + x21*(y12*a12[0][2] + x21*a22[0][2]);
		ds33[1][0] = y12*(y12*a11[1][0] + x21*a21[1][0]) + x21*(y12*a12[1][0] + x21*a22[1][0]);
		ds33[1][1] = y12*(y12*a11[1][1] + x21*a21[1][1]) + x21*(y12*a12[1][1] + x21*a22[1][1]);
		ds33[1][2] = y12*(y12*a11[1][2] + x21*a21[1][2]) + x21*(y12*a12[1][2] + x21*a22[1][2]);
		ds33[2][0] = y12*(y12*a11[2][0] + x21*a21[2][0]) + x21*(y12*a12[2][0] + x21*a22[2][0]);
		ds33[2][1] = y12*(y12*a11[2][1] + x21*a21[2][1]) + x21*(y12*a12[2][1] + x21*a22[2][1]);
		ds33[2][2] = y12*(y12*a11[2][2] + x21*a21[2][2]) + x21*(y12*a12[2][2] + x21*a22[2][2]);

		// *** Matriz de Conveccao da captura de choque
		double sh12, sh13, sh23;

		sh12 = y23*y31 + x32*x13;
		sh13 = y23*y12 + x32*x21;
		sh23 = y31*y12 + x13*x21;

		// *** Coeficientes da Matriz de Rigidez [De]

		De[0][0] = sixth*ms1[0][0]                 + tau*forth*invArea*ds11[0][0] + delta*forth*invArea*(-sh12 -sh13);
		De[0][1] = sixth*ms1[0][1]                 + tau*forth*invArea*ds11[0][1] + delta*forth*invArea*(-sh12 -sh13);
		De[0][2] = sixth*ms1[0][2]                 + tau*forth*invArea*ds11[0][2] + delta*forth*invArea*(-sh12 -sh13);
		De[1][0] = sixth*ms1[1][0]                 + tau*forth*invArea*ds11[1][0] + delta*forth*invArea*(-sh12 -sh13);
		De[1][1] = sixth*ms1[1][1] + (- d12 - d13) + tau*forth*invArea*ds11[1][1] + delta*forth*invArea*(-sh12 -sh13);
		De[1][2] = sixth*ms1[1][2]                 + tau*forth*invArea*ds11[1][2] + delta*forth*invArea*(-sh12 -sh13);
		De[2][0] = sixth*ms1[2][0]                 + tau*forth*invArea*ds11[2][0] + delta*forth*invArea*(-sh12 -sh13);
		De[2][1] = sixth*ms1[2][1]                 + tau*forth*invArea*ds11[2][1] + delta*forth*invArea*(-sh12 -sh13);
		De[2][2] = sixth*ms1[2][2] + (- d12 - d13) + tau*forth*invArea*ds11[2][2] + delta*forth*invArea*(-sh12 -sh13);

		De[0][3] = sixth*ms2[0][0]       + tau*forth*invArea*ds12[0][0] + delta*forth*invArea*sh12;
		De[0][4] = sixth*ms2[0][1]       + tau*forth*invArea*ds12[0][1] + delta*forth*invArea*sh12;
		De[0][5] = sixth*ms2[0][2]       + tau*forth*invArea*ds12[0][2] + delta*forth*invArea*sh12;
		De[1][3] = sixth*ms2[1][0]       + tau*forth*invArea*ds12[1][0] + delta*forth*invArea*sh12;
		De[1][4] = sixth*ms2[1][1] + d12 + tau*forth*invArea*ds12[1][1] + delta*forth*invArea*sh12;
		De[1][5] = sixth*ms2[1][2]       + tau*forth*invArea*ds12[1][2] + delta*forth*invArea*sh12;
		De[2][3] = sixth*ms2[2][0]       + tau*forth*invArea*ds12[2][0] + delta*forth*invArea*sh12;
		De[2][4] = sixth*ms2[2][1]       + tau*forth*invArea*ds12[2][1] + delta*forth*invArea*sh12;
		De[2][5] = sixth*ms2[2][2] + d12 + tau*forth*invArea*ds12[2][2] + delta*forth*invArea*sh12;

		De[0][6] = sixth*ms3[0][0]       + tau*forth*invArea*ds13[0][0] + delta*forth*invArea*sh13;
		De[0][7] = sixth*ms3[0][1]       + tau*forth*invArea*ds13[0][1] + delta*forth*invArea*sh13;
		De[0][8] = sixth*ms3[0][2]       + tau*forth*invArea*ds13[0][2] + delta*forth*invArea*sh13;
		De[1][6] = sixth*ms3[1][0]       + tau*forth*invArea*ds13[1][0] + delta*forth*invArea*sh13;
		De[1][7] = sixth*ms3[1][1] + d13 + tau*forth*invArea*ds13[1][1] + delta*forth*invArea*sh13;
		De[1][8] = sixth*ms3[1][2]       + tau*forth*invArea*ds13[1][2] + delta*forth*invArea*sh13;
		De[2][6] = sixth*ms3[2][0]       + tau*forth*invArea*ds13[2][0] + delta*forth*invArea*sh13;
		De[2][7] = sixth*ms3[2][1]       + tau*forth*invArea*ds13[2][1] + delta*forth*invArea*sh13;
		De[2][8] = sixth*ms3[2][2] + d13 + tau*forth*invArea*ds13[2][2] + delta*forth*invArea*sh13;

		De[3][0] = sixth*ms1[0][0]       + tau*forth*invArea*ds12[0][0] + delta*forth*invArea*sh12;
		De[3][1] = sixth*ms1[0][1]       + tau*forth*invArea*ds12[0][1] + delta*forth*invArea*sh12;
		De[3][2] = sixth*ms1[0][2]       + tau*forth*invArea*ds12[0][2] + delta*forth*invArea*sh12;
		De[4][0] = sixth*ms1[1][0]       + tau*forth*invArea*ds12[1][0] + delta*forth*invArea*sh12;
		De[4][1] = sixth*ms1[1][1] + d13 + tau*forth*invArea*ds12[1][1] + delta*forth*invArea*sh12;
		De[4][2] = sixth*ms1[1][2]       + tau*forth*invArea*ds12[1][2] + delta*forth*invArea*sh12;
		De[5][0] = sixth*ms1[2][0]       + tau*forth*invArea*ds12[2][0] + delta*forth*invArea*sh12;
		De[5][1] = sixth*ms1[2][1]       + tau*forth*invArea*ds12[2][1] + delta*forth*invArea*sh12;
		De[5][2] = sixth*ms1[2][2] + d13 + tau*forth*invArea*ds12[2][2] + delta*forth*invArea*sh12;

		De[3][3] = sixth*ms2[0][0]                 + tau*forth*invArea*ds22[0][0] + delta*forth*invArea*(-sh12 -sh23);
		De[3][4] = sixth*ms2[0][1]                 + tau*forth*invArea*ds22[0][1] + delta*forth*invArea*(-sh12 -sh23);
		De[3][5] = sixth*ms2[0][2]                 + tau*forth*invArea*ds22[0][2] + delta*forth*invArea*(-sh12 -sh23);
		De[4][3] = sixth*ms2[1][0]                 + tau*forth*invArea*ds22[1][0] + delta*forth*invArea*(-sh12 -sh23);
		De[4][4] = sixth*ms2[1][1] + (- d12 - d13) + tau*forth*invArea*ds22[1][1] + delta*forth*invArea*(-sh12 -sh23);
		De[4][5] = sixth*ms2[1][2]                 + tau*forth*invArea*ds22[1][2] + delta*forth*invArea*(-sh12 -sh23);
		De[5][3] = sixth*ms2[2][0]                 + tau*forth*invArea*ds22[2][0] + delta*forth*invArea*(-sh12 -sh23);
		De[5][4] = sixth*ms2[2][1]                 + tau*forth*invArea*ds22[2][1] + delta*forth*invArea*(-sh12 -sh23);
		De[5][5] = sixth*ms2[2][2] + (- d12 - d13) + tau*forth*invArea*ds22[2][2] + delta*forth*invArea*(-sh12 -sh23);

		De[3][6] = sixth*ms3[0][0]       + tau*forth*invArea*ds23[0][0] + delta*forth*invArea*sh23;
		De[3][7] = sixth*ms3[0][1]       + tau*forth*invArea*ds23[0][1] + delta*forth*invArea*sh23;
		De[3][8] = sixth*ms3[0][2]       + tau*forth*invArea*ds23[0][2] + delta*forth*invArea*sh23;
		De[4][6] = sixth*ms3[1][0]       + tau*forth*invArea*ds23[1][0] + delta*forth*invArea*sh23;
		De[4][7] = sixth*ms3[1][1] + d23 + tau*forth*invArea*ds23[1][1] + delta*forth*invArea*sh23;
		De[4][8] = sixth*ms3[1][2]       + tau*forth*invArea*ds23[1][2] + delta*forth*invArea*sh23;
		De[5][6] = sixth*ms3[2][0]       + tau*forth*invArea*ds23[2][0] + delta*forth*invArea*sh23;
		De[5][7] = sixth*ms3[2][1]       + tau*forth*invArea*ds23[2][1] + delta*forth*invArea*sh23;
		De[5][8] = sixth*ms3[2][2] + d23 + tau*forth*invArea*ds23[2][2] + delta*forth*invArea*sh23;

		De[6][0] = sixth*ms1[0][0]       + tau*forth*invArea*ds31[0][0] + delta*forth*invArea*sh13;
		De[6][1] = sixth*ms1[0][1]       + tau*forth*invArea*ds31[0][1] + delta*forth*invArea*sh13;
		De[6][2] = sixth*ms1[0][2]       + tau*forth*invArea*ds31[0][2] + delta*forth*invArea*sh13;
		De[7][0] = sixth*ms1[1][0]       + tau*forth*invArea*ds31[1][0] + delta*forth*invArea*sh13;
		De[7][1] = sixth*ms1[1][1] + d13 + tau*forth*invArea*ds31[1][1] + delta*forth*invArea*sh13;
		De[7][2] = sixth*ms1[1][2]       + tau*forth*invArea*ds31[1][2] + delta*forth*invArea*sh13;
		De[8][0] = sixth*ms1[2][0]       + tau*forth*invArea*ds31[2][0] + delta*forth*invArea*sh13;
		De[8][1] = sixth*ms1[2][1]       + tau*forth*invArea*ds31[2][1] + delta*forth*invArea*sh13;
		De[8][2] = sixth*ms1[2][2] + d13 + tau*forth*invArea*ds31[2][2] + delta*forth*invArea*sh13;

		De[6][3] = sixth*ms2[0][0]       + tau*forth*invArea*ds32[0][0] + delta*forth*invArea*sh23;
		De[6][4] = sixth*ms2[0][1]       + tau*forth*invArea*ds32[0][1] + delta*forth*invArea*sh23;
		De[6][5] = sixth*ms2[0][2]       + tau*forth*invArea*ds32[0][2] + delta*forth*invArea*sh23;
		De[7][3] = sixth*ms2[1][0]       + tau*forth*invArea*ds32[1][0] + delta*forth*invArea*sh23;
		De[7][4] = sixth*ms2[1][1] + d23 + tau*forth*invArea*ds32[1][1] + delta*forth*invArea*sh23;
		De[7][5] = sixth*ms2[1][2]       + tau*forth*invArea*ds32[1][2] + delta*forth*invArea*sh23;
		De[8][3] = sixth*ms2[2][0]       + tau*forth*invArea*ds32[2][0] + delta*forth*invArea*sh23;
		De[8][4] = sixth*ms2[2][1]       + tau*forth*invArea*ds32[2][1] + delta*forth*invArea*sh23;
		De[8][5] = sixth*ms2[2][2] + d23 + tau*forth*invArea*ds32[2][2] + delta*forth*invArea*sh23;

		De[6][6] = sixth*ms3[0][0]                 + tau*forth*invArea*ds33[0][0] + delta*forth*invArea*(-sh13 -sh23);
		De[6][7] = sixth*ms3[0][1]                 + tau*forth*invArea*ds33[0][1] + delta*forth*invArea*(-sh13 -sh23);
		De[6][8] = sixth*ms3[0][2]                 + tau*forth*invArea*ds33[0][2] + delta*forth*invArea*(-sh13 -sh23);
		De[7][6] = sixth*ms3[1][0]                 + tau*forth*invArea*ds33[1][0] + delta*forth*invArea*(-sh13 -sh23);
		De[7][7] = sixth*ms3[1][1] + (- d12 - d23) + tau*forth*invArea*ds33[1][1] + delta*forth*invArea*(-sh13 -sh23);
		De[7][8] = sixth*ms3[1][2]                 + tau*forth*invArea*ds33[1][2] + delta*forth*invArea*(-sh13 -sh23);
		De[8][6] = sixth*ms3[2][0]                 + tau*forth*invArea*ds33[2][0] + delta*forth*invArea*(-sh13 -sh23);
		De[8][7] = sixth*ms3[2][1]                 + tau*forth*invArea*ds33[2][1] + delta*forth*invArea*(-sh13 -sh23);
		De[8][8] = sixth*ms3[2][2] + (- d12 - d23) + tau*forth*invArea*ds33[2][2] + delta*forth*invArea*(-sh13 -sh23);

		// *** Matriz do termo fonte

		Fe[0] = Area*twelve*(2*S1[0] + S2[0] + S3[0]) + tau*second*(ms1[0][0]*Sb[0] + ms1[0][1]*Sb[1] + ms1[0][2]*Sb[2]);
		Fe[1] = Area*twelve*(2*S1[1] + S2[1] + S3[1]) + tau*second*(ms1[1][0]*Sb[0] + ms1[1][1]*Sb[1] + ms1[1][2]*Sb[2]);
		Fe[2] = Area*twelve*(2*S1[2] + S2[2] + S3[2]) + tau*second*(ms1[2][0]*Sb[0] + ms1[2][1]*Sb[1] + ms1[2][2]*Sb[2]);
		Fe[3] = Area*twelve*(S1[0] + 2*S2[0] + S3[0]) + tau*second*(ms2[0][0]*Sb[0] + ms2[0][1]*Sb[1] + ms2[0][2]*Sb[2]);
		Fe[4] = Area*twelve*(S1[1] + 2*S2[1] + S3[1]) + tau*second*(ms2[1][0]*Sb[0] + ms2[1][1]*Sb[1] + ms2[1][2]*Sb[2]);
		Fe[5] = Area*twelve*(S1[2] + 2*S2[2] + S3[2]) + tau*second*(ms2[2][0]*Sb[0] + ms2[2][1]*Sb[1] + ms2[2][2]*Sb[2]);
		Fe[6] = Area*twelve*(S1[0] + S2[0] + 2*S3[0]) + tau*second*(ms3[0][0]*Sb[0] + ms3[0][1]*Sb[1] + ms3[0][2]*Sb[2]);
		Fe[7] = Area*twelve*(S1[1] + S2[1] + 2*S3[1]) + tau*second*(ms3[1][0]*Sb[0] + ms3[1][1]*Sb[1] + ms3[1][2]*Sb[2]);
		Fe[8] = Area*twelve*(S1[2] + S2[2] + 2*S3[2]) + tau*second*(ms3[2][0]*Sb[0] + ms3[2][1]*Sb[1] + ms3[2][2]*Sb[2]);

		F_assembly(e, Fe, De, FemFunctions, FemStructs, neq);

		//*****************************************

		//Fill local Re and Global R
		double Ae[9][9], Re[9], MedUe[9], DeUe[9];
		
		for(i=0; i<9; i++){
			MedUe[i] = 0;
			DeUe[i] = 0;
			for (j=0; j<9; j++){
				MedUe[i] += Me[i][j]*dUe[j];
				DeUe[i] += De[i][j]*Ue[j];
				Ae[i][j] = Me[i][j] + alpha*delta_t*De[i][j];
			}
			Re[i] = - MedUe[i] - DeUe[i];  // F já está sendo adicionado ao R via F_assembly() com tratamento de contorno
		}

		for (i=0; i<9; i++)
			R[lm[e][i]] += Re[i];
		R[neq] = 0;
		
		FemFunctions->assembly(Parameters, MatrixData, FemStructs, e, Ae);

	}

	myfree(U);
	myfree(dU);

	return 0;
}
