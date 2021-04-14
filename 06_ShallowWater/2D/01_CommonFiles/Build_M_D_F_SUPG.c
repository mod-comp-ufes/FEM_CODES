#include "ShalowWater.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"


int Build_M_D_F_SUPG(ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	// Counts
	int e, i, j;

	// Constants
	double g = 9.81;
	double n = 0.018;
	double mu = 0.001;
	double rho = 1000;
	double Cf, gamma;
	double tau, delta;
	double second = 0.5, third = 1.0/3.0, forth = 0.25, sixth = 1.0/6.0, twelve = 1.0/12.0;

	// Geometry
	int J1, J2, J3, nel, neq;
	double x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21, zb1, zb2, zb3;
	double Area, twoArea, invArea, invtwoArea;

	// Variables
	double u[3], v[3];	
	double *U, *dU;
	double Me[9][9], De[9][9], Ue[9], dUe[9], Ub[3], dUb[3], Fb[3], gradUx[3], gradUy[3];
	double S1[3], S2[3], S3[3], Sb[3], Fe[9], gradzbx, gradzby;
	double alpha = Parameters->Alpha_Build;
	double delta_t = Parameters->DeltaT_Build;
	double *R = FemStructs->F;
	int **lm = FemStructs->lm;
	double A1[3][3], A2[3][3];

	// Auxiliars
	double r;
	double Mg1, Mg2;
	double sh12, sh13, sh23;
	double tauSUGN1, tauSUGN2;
	double jbold[2], normjbold;
	double d12[2], d13[2], d23[2];
	double ms1[3][3], ms2[3][3], ms3[3][3];
	double a11[3][3], a12[3][3], a21[3][3], a22[3][3];
	double ds11[3][3], ds12[3][3], ds13[3][3], ds21[3][3], ds22[3][3], ds23[3][3], ds31[3][3], ds32[3][3], ds33[3][3];

	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;

	nel = Parameters->nel;
	neq = Parameters->neq;

	dzero(neq+1, R);
	setzeros(Parameters, MatrixData);

	U = (double*) mycalloc("U of 'Build_M_F_SUPG_Transiente'", 3*Parameters->nnodes, sizeof(double));
	dU = (double*) mycalloc("dU of 'Build_M_F_SUPG_Transiente'", 3*Parameters->nnodes, sizeof(double));

	eval_U_dU(Parameters,FemStructs,FemFunctions,U,dU);

	for(e=0; e<nel; e++)
	{
		// *Nos do elemento
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

		// Area do elemento
		twoArea = fabs(x21 * y31 - x13 * y12);
		Area = 0.5*twoArea;
		invArea = 1.0/Area;
		invtwoArea = 1.0/(2*Area);

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

		for(i=0; i<3; i++)
		{
			// Baricentro do triangulo
			Ub[i] = (Ue[i] + Ue[i+3] + Ue[i+6])*third;
			dUb[i] = (dUe[i] + dUe[i+3] + dUe[i+6])*third;

			// Calculo do Gradiente
			gradUx[i] = (Ue[i]*y23 + Ue[i+3]*y31 + Ue[i+6]*y12)*invtwoArea;
			gradUy[i] = (Ue[i]*x32 + Ue[i+3]*x13 + Ue[i+6]*x21)*invtwoArea;
		}

		// A = [A1 A2]
		A1_A2_calculations(Ue, A1, A2, g);

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
		a21[1][1] = A2[1][0] + A2[1][1]*A1[1][1] + A2[1][2]*A1[2][1];
		a21[1][2] = A2[1][2]*A2[2][2];
		a21[2][0] = A2[2][2]*A1[2][0];
		a21[2][1] = A2[2][0] + A2[2][2]*A1[2][1];
		a21[2][2] = A2[2][2]*A2[2][2];

		a22[0][0] = A2[2][0];
		a22[0][1] = 0.0;
		a22[0][2] = A2[2][2];
		a22[1][0] = A2[1][1]*A2[1][0] + A2[1][2]*A2[2][0];
		a22[1][1] = A2[1][1]*A2[1][1];
		a22[1][2] = A2[1][0] + A2[1][1]*A2[2][0] + A2[1][2]*A2[2][2];
		a22[2][0] = A2[2][0];
		a22[2][1] = 0.0;
		a22[2][2] = A2[2][0] + A2[2][2]*A2[2][2];

		// Source
		zb1 = FemFunctions->zb(x1, y1);
		zb2 = FemFunctions->zb(x2, y2);
		zb3 = FemFunctions->zb(x3, y3);

		gradzbx = (zb1*y23 + zb2*y31 + zb3*y12)*invtwoArea;
		gradzby = (zb1*x32 + zb2*x13 + zb3*x21)*invtwoArea;

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

		// Parametro de estabilizacao (tau) do SUPG
		jbold[0] = (y23*(Ue[0] + zb1) + y31*(Ue[3] + zb2) + y12*(Ue[6] + zb3))*invtwoArea;
		jbold[1] = (x32*(Ue[0] + zb1) + x13*(Ue[3] + zb2) + x21*(Ue[6] + zb3))*invtwoArea;
		r = jbold[0]*jbold[0] + jbold[1]*jbold[1];
		normjbold = (r>=0) ? sqrt(r) : 0.0;
		jbold[0] /= normjbold;
		jbold[1] /= normjbold;

		r = g*Ub[0];
		tauSUGN1 = (r>0) ? 1.0/(sqrt(r)*(abs(jbold[0]*y23 + jbold[1]*x32)*invtwoArea + abs(Ub[1]*y23 + Ub[2]*x32)*invtwoArea +
		                                 abs(jbold[0]*y31 + jbold[1]*x13)*invtwoArea + abs(Ub[1]*y31 + Ub[2]*x13)*invtwoArea +
								         abs(jbold[0]*y12 + jbold[1]*x21)*invtwoArea + abs(Ub[2]*y12 + Ub[2]*x21)*invtwoArea)) : 0.0;
		tauSUGN2 = delta_t/2.0;
		tau = (tauSUGN1>0 && tauSUGN2>0) ? 1.0/sqrt(1.0/(tauSUGN1*tauSUGN1) + 1.0/(tauSUGN2*tauSUGN2)) : 0.0;

		delta = FemFunctions->ShockCapture(Ub, gradUx, gradUy, A1, A2, Sb, y23, y31, y12, x32, x13, x21, twoArea, tau, g);

		// Matriz de Massa do Galerkin
		Mg1 = Area*sixth;
		Mg2 = Area*twelve;

		// Matriz de Massa do SUPG 
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

		// Coeficientes da Matriz de Massa [Me]
		for(i=0; i<3; i++)
		{
			for(j=0; j<3; j++)
			{
				Me[  i][  j] = Mg1 + tau*sixth*ms1[i][j];
				Me[3+i][3+j] = Mg1 + tau*sixth*ms2[i][j];
				Me[6+i][6+j] = Mg1 + tau*sixth*ms3[i][j];
			}
			for(j=0; j<3; j++)
			{
				Me[  i][6+j] = Me[  i][3+j] = Mg2 + tau*sixth*ms1[i][j];
				Me[3+i][6+j] = Me[3+i][  j] = Mg2 + tau*sixth*ms2[i][j];
				Me[6+i][3+j] = Me[6+i][  j] = Mg2 + tau*sixth*ms3[i][j];
			}
		}

		// Matriz de Conveccao de Galerkin
		d12[0] = d13[0] = d23[0] = 0;
		d12[1] = y23*y31*mu/rho + x32*x13*mu/rho;
		d13[1] = y23*y12*mu/rho + x32*x21*mu/rho;
		d23[1] = y31*y12*mu/rho + x13*x21*mu/rho;

		// Matriz de Conveccao do SUPG
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
			{
				ds11[i][j] = y23*(y23*a11[i][j] + x32*a21[i][j]) + x32*(y23*a12[i][j] + x32*a22[i][j]);
				ds12[i][j] = y31*(y23*a11[i][j] + x32*a21[i][j]) + x13*(y23*a12[i][j] + x32*a22[i][j]);
				ds13[i][j] = y12*(y23*a11[i][j] + x32*a21[i][j]) + x21*(y23*a12[i][j] + x32*a22[i][j]);
				ds21[i][j] = y23*(y31*a11[i][j] + x13*a21[i][j]) + x32*(y31*a12[i][j] + x13*a22[i][j]);
				ds22[i][j] = y31*(y31*a11[i][j] + x13*a21[i][j]) + x13*(y31*a12[i][j] + x13*a22[i][j]);
				ds23[i][j] = y12*(y31*a11[i][j] + x13*a21[i][j]) + x21*(y31*a12[i][j] + x13*a22[i][j]);
				ds31[i][j] = y23*(y12*a11[i][j] + x21*a21[i][j]) + x32*(y12*a12[i][j] + x21*a22[i][j]);
				ds32[i][j] = y31*(y12*a11[i][j] + x21*a21[i][j]) + x13*(y12*a12[i][j] + x21*a22[i][j]);
				ds33[i][j] = y12*(y12*a11[i][j] + x21*a21[i][j]) + x21*(y12*a12[i][j] + x21*a22[i][j]);
			}

		// Matriz de Conveccao da captura de choque
		sh12 = y23*y31 + x32*x13;
		sh13 = y23*y12 + x32*x21;
		sh23 = y31*y12 + x13*x21;

		// Coeficientes da Matriz de Rigidez [De]
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
			{
				De[  i][  j] = sixth*ms1[i][j] + forth*invArea*(- d12[i&(i==j)] - d13[i&(i==j)]) + tau*forth*invArea*ds11[i][j] + delta*forth*invArea*(-sh12 -sh13);
				De[  i][3+j] = sixth*ms2[i][j] + forth*invArea*d12[i&(i==j)]                     + tau*forth*invArea*ds12[i][j] + delta*forth*invArea*sh12;
				De[  i][6+j] = sixth*ms3[i][j] + forth*invArea*d13[i&(i==j)]                     + tau*forth*invArea*ds13[i][j] + delta*forth*invArea*sh13;
				De[3+i][  j] = sixth*ms1[i][j] + forth*invArea*d12[i&(i==j)]                     + tau*forth*invArea*ds12[i][j] + delta*forth*invArea*sh12;
				De[3+i][3+j] = sixth*ms2[i][j] + forth*invArea*(- d12[i&(i==j)] - d13[i&(i==j)]) + tau*forth*invArea*ds22[i][j] + delta*forth*invArea*(-sh12 -sh23);
				De[3+i][6+j] = sixth*ms3[i][j] + forth*invArea*d23[i&(i==j)]                     + tau*forth*invArea*ds23[i][j] + delta*forth*invArea*sh23;
				De[6+i][  j] = sixth*ms1[i][j] + forth*invArea*d13[i&(i==j)]                     + tau*forth*invArea*ds31[i][j] + delta*forth*invArea*sh13;
				De[6+i][3+j] = sixth*ms2[i][j] + forth*invArea*d23[i&(i==j)]                     + tau*forth*invArea*ds32[i][j] + delta*forth*invArea*sh23;
				De[6+i][6+j] = sixth*ms3[i][j] + forth*invArea*(- d12[i&(i==j)] - d23[i&(i==j)]) + tau*forth*invArea*ds33[i][j] + delta*forth*invArea*(-sh13 -sh23);
			}

		// Matriz do termo fonte
		for(i=0; i<3; i++)
		{
			Fe[  i] = Area*twelve*(2*S1[i] + S2[i] + S3[i]) + tau*second*(ms1[i][0]*Sb[0] + ms1[i][1]*Sb[1] + ms1[i][2]*Sb[2]);
			Fe[3+i] = Area*twelve*(S1[i] + 2*S2[i] + S3[i]) + tau*second*(ms2[i][0]*Sb[0] + ms2[i][1]*Sb[1] + ms2[i][2]*Sb[2]);
			Fe[6+i] = Area*twelve*(S1[i] + S2[i] + 2*S3[i]) + tau*second*(ms3[i][0]*Sb[0] + ms3[i][1]*Sb[1] + ms3[i][2]*Sb[2]);
		}
		F_assembly(e, Fe, De, FemFunctions, FemStructs, neq);

		//Fill local Re and Global R
		double Ae[9][9], Re[9], MedUe[9], DeUe[9];
		
		for(i=0; i<9; i++)
		{
			MedUe[i] = 0;
			DeUe[i] = 0;
			for(j=0; j<9; j++)
			{
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
