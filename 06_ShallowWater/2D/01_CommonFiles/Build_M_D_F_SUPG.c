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
	double tau, delta, Cf;
	double second = 0.5, third = 1.0/3.0, forth = 0.25, sixth = 1.0/6.0, twelve = 1.0/12.0;

	// Geometry
	int J1, J2, J3, nel, neq;
	double x1, x2, x3, y1, y2, y3, y23, y31, y12, x32, x13, x21, zb1, zb2, zb3;
	double Area, twoArea, invArea, invtwoArea;

	// Variables
	double h, u, v;	// baricentro
	double *U, *dU;
	double Me[9][9], De[9][9], Ue[9], dUe[9], Ub[3], dUb[3], Fb[3], gradUx[3], gradUy[3], Rbold[3];
	double S1[3], S2[3], S3[3], Sb[3], Fe[9], gradzbx, gradzby;
	double alpha = Parameters->Alpha;
	double delta_t = Parameters->DeltaT;
	double *R = FemStructs->F;
	int **lm = FemStructs->lm;
	double A1[3][3], A2[3][3];

	// Auxiliars
	double r;
	double Mg1, Mg2;
	double sh12[3][3], sh13[3][3], sh23[3][3];
	double tauSUGN1, tauSUGN2;
	double gradEta[2], normgradEta;
	double d12[3][3], d13[3][3], d23[3][3];
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
		// Nos do elemento
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
		twoArea = fabs(x21*y31 - x13*y12);
		Area = 0.5*twoArea;
		invArea = 1.0/Area;
		invtwoArea = 1.0/twoArea;

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

		for(i=0; i<3; i++)
		{
			// Baricentro do triangulo
			Ub[i] = (Ue[i] + Ue[i+3] + Ue[i+6])*third;
			dUb[i] = (dUe[i] + dUe[i+3] + dUe[i+6])*third;

			// Calculo do Gradiente
			gradUx[i] = (Ue[i]*y23 + Ue[i+3]*y31 + Ue[i+6]*y12)*invtwoArea;
			gradUy[i] = (Ue[i]*x32 + Ue[i+3]*x13 + Ue[i+6]*x21)*invtwoArea;
		}

		// Baricentro do triangulo
		h = Ub[0];
		u = Ub[1]/h;
		v = Ub[2]/h;

		// A = [A1 A2]
		A1[0][0] = 0.0;
		A1[0][1] = 1.0;
		A1[0][2] = 0.0;
		A1[1][0] = g*h - u*u;
		A1[1][1] = 2.0*u;
		A1[1][2] = 0.0;
		A1[2][0] = -u*v;
		A1[2][1] = v;
		A1[2][2] = u;

		// A2 coefficients
		A2[0][0] = 0.0;
		A2[0][1] = 0.0;
		A2[0][2] = 1.0;
		A2[1][0] = -u*v;
		A2[1][1] = v;
		A2[1][2] = u;
		A2[2][0] = g*h - v*v;
		A2[2][1] = 0.0;
		A2[2][2] = 2*v;

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
		a21[2][2] = A2[2][2]*A1[2][2];

		a22[0][0] = A2[2][0];
		a22[0][1] = 0.0;
		a22[0][2] = A2[2][2];
		a22[1][0] = A2[1][1]*A2[1][0] + A2[1][2]*A2[2][0];
		a22[1][1] = A2[1][1]*A2[1][1];
		a22[1][2] = A2[1][0] + A2[1][1]*A2[1][2] + A2[1][2]*A2[2][2];
		a22[2][0] = A2[2][0];
		a22[2][1] = 0.0;
		a22[2][2] = A2[2][0] + A2[2][2]*A2[2][2];

		// Source
		zb1 = FemFunctions->zb(x1, y1);
		zb2 = FemFunctions->zb(x2, y2);
		zb3 = FemFunctions->zb(x3, y3);

		gradzbx = (zb1*y23 + zb2*y31 + zb3*y12)*invtwoArea;
		gradzby = (zb1*x32 + zb2*x13 + zb3*x21)*invtwoArea;

		Cf = g*n*n*pow(Ub[0], -7/3);

		S1[0] = 0.0;
		S2[0] = 0.0;
		S3[0] = 0.0;

		r = FemFunctions->gammaBed(Cf, x2, y2);
		S1[1] = -g*gradzbx*Ue[0] - r*Ue[1];
		S2[1] = -g*gradzbx*Ue[3] - r*Ue[4];
		S3[1] = -g*gradzbx*Ue[6] - r*Ue[7];

		r = FemFunctions->gammaBed(Cf, x3, y3);
		S1[2] = -g*gradzby*Ue[0] - r*Ue[2];
		S2[2] = -g*gradzby*Ue[3] - r*Ue[5];
		S3[2] = -g*gradzby*Ue[6] - r*Ue[8];

		Sb[0] = third*(S1[0] + S2[0] + S3[0]);
		Sb[1] = third*(S1[1] + S2[2] + S3[1]);
		Sb[2] = third*(S1[2] + S2[2] + S3[2]);
		printf("(%.2lf, %.2lf)\t(%.2lf, %.2lf)\t(%.2lf, %.2lf)\n", x1, y1, x2, y2, x3, y3);
		// Parametro de estabilizacao (tau) do SUPG
		gradEta[0] = (y23*(Ue[0] + zb1) + y31*(Ue[3] + zb2) + y12*(Ue[6] + zb3))*invtwoArea;
		gradEta[1] = (x32*(Ue[0] + zb1) + x13*(Ue[3] + zb2) + x21*(Ue[6] + zb3))*invtwoArea;
		normgradEta = sqrt(gradEta[0]*gradEta[0] + gradEta[1]*gradEta[1]);

		if(fabs(normgradEta) >= TOL) {
			gradEta[0] /= normgradEta;
			gradEta[1] /= normgradEta;
		}
		else {
			gradEta[0] = 0.0;
			gradEta[1] = 0.0;
		}

		r = sqrt(g*Ue[0])*fabs(gradEta[0]*y23 + gradEta[1]*x32) + fabs(u*y23 + v*x32) +
		    sqrt(g*Ue[3])*fabs(gradEta[0]*y31 + gradEta[1]*x13) + fabs(u*y31 + v*x13) +
			sqrt(g*Ue[6])*fabs(gradEta[0]*y12 + gradEta[1]*x21) + fabs(u*y12 + v*x21);
		tauSUGN1 = (fabs(r) >= 1e-12) ? twoArea/r : 0.0;

		tauSUGN2 = delta_t/2.0;

		if(fabs(tauSUGN1) >= TOL)
			tau = 1.0/sqrt(1.0/(tauSUGN1*tauSUGN1) + 1.0/(tauSUGN2*tauSUGN2));
		else
			tau = 1.0/sqrt(1.0/(tauSUGN2*tauSUGN2));

		for(i=0; i<3; i++)
		{
			// Calculo de R = AxgradUx + AygradUy - S
			Rbold[i] = A1[i][0]*gradUx[0] + A1[i][1]*gradUx[1] + A1[i][2]*gradUx[2] + 
			           A2[i][0]*gradUy[0] + A2[i][1]*gradUy[1] + A2[i][2]*gradUy[2] - Sb[i];
		}
		delta = FemFunctions->ShockCapture(Ub, Rbold, gradUx, gradUy, y23, y31, y12, x32, x13, x21, invtwoArea, tau, g);

		// Matriz de Massa do Galerkin
		Mg1 = Area*sixth;
		Mg2 = Area*twelve;

		// Matriz de Massa do SUPG
		for(i = 0; i < 3; i++)
			for(j = 0; j < 3; j++) {
				ms1[i][j] = y23*A1[i][j] + x32*A2[i][j];
				ms2[i][j] = y31*A1[i][j] + x13*A2[i][j];
				ms3[i][j] = y12*A1[i][j] + x21*A2[i][j];
			}

		// Coeficientes da Matriz de Massa [Me]
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
			{
				Me[  i][  j] = Mg1 + tau*sixth*ms1[i][j];
				Me[3+i][3+j] = Mg1 + tau*sixth*ms2[i][j];
				Me[6+i][6+j] = Mg1 + tau*sixth*ms3[i][j];

				Me[  i][3+j] = Mg2 + tau*sixth*ms1[i][j];
				Me[3+i][  j] = Mg2 + tau*sixth*ms2[i][j];
				Me[6+i][  j] = Mg2 + tau*sixth*ms3[i][j];

				Me[  i][6+j] = Mg2 + tau*sixth*ms1[i][j];
				Me[3+i][6+j] = Mg2 + tau*sixth*ms2[i][j];
				Me[6+i][3+j] = Mg2 + tau*sixth*ms3[i][j];
			}

		// Matriz de Conveccao de Galerkin
		for(i=0; i<3; i++) {
			for(j=0; j<3; j++) {
				d12[i][j] = 0.0;
				d13[i][j] = 0.0;
				d23[i][j] = 0.0;
			}
			if(i>0) {
				d12[i][i] = y23*y31*mu/rho + x32*x13*mu/rho;
				d13[i][i] = y23*y12*mu/rho + x32*x21*mu/rho;
				d23[i][i] = y31*y12*mu/rho + x13*x21*mu/rho;
			}
		}

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
		for(i=0; i<3; i++) {
			for(j=0; j<3; j++) {
				sh12[i][j] = 0.0;
				sh13[i][j] = 0.0;
				sh23[i][j] = 0.0;
			}
			sh12[i][i] = y23*y31 + x32*x13;
			sh13[i][i] = y23*y12 + x32*x21;
			sh23[i][i] = y31*y12 + x13*x21;
		}

		// Coeficientes da Matriz de Rigidez [De]	
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
			{
				De[  i][  j] = sixth*ms1[i][j] + forth*invArea*(-d12[i][j] -d13[i][j]) + tau*forth*invArea*ds11[i][j] + delta*forth*invArea*(-sh12[i][j] -sh13[i][j]);
				De[  i][3+j] = sixth*ms2[i][j] + forth*invArea*  d12[i][j]             + tau*forth*invArea*ds12[i][j] + delta*forth*invArea*sh12[i][j];
				De[  i][6+j] = sixth*ms3[i][j] + forth*invArea*  d13[i][j]             + tau*forth*invArea*ds13[i][j] + delta*forth*invArea*sh13[i][j];
				De[3+i][  j] = sixth*ms1[i][j] + forth*invArea*  d12[i][j]             + tau*forth*invArea*ds21[i][j] + delta*forth*invArea*sh12[i][j];
				De[3+i][3+j] = sixth*ms2[i][j] + forth*invArea*(-d12[i][j] -d23[i][j]) + tau*forth*invArea*ds22[i][j] + delta*forth*invArea*(-sh12[i][j] -sh23[i][j]);
				De[3+i][6+j] = sixth*ms3[i][j] + forth*invArea*  d23[i][j]             + tau*forth*invArea*ds23[i][j] + delta*forth*invArea*sh23[i][j];
				De[6+i][  j] = sixth*ms1[i][j] + forth*invArea*  d13[i][j]             + tau*forth*invArea*ds31[i][j] + delta*forth*invArea*sh13[i][j];
				De[6+i][3+j] = sixth*ms2[i][j] + forth*invArea*  d23[i][j]             + tau*forth*invArea*ds32[i][j] + delta*forth*invArea*sh23[i][j];
				De[6+i][6+j] = sixth*ms3[i][j] + forth*invArea*(-d13[i][j] -d23[i][j]) + tau*forth*invArea*ds33[i][j] + delta*forth*invArea*(-sh13[i][j] -sh23[i][j]);
			}

		// Matriz do termo fonte
		for(i=0; i<3; i++)
		{
			Fe[  i] = Area*twelve*(2*S1[i] + S2[i] + S3[i]) + tau*sixth*(ms1[i][0]*Sb[0] + ms1[i][1]*Sb[1] + ms1[i][2]*Sb[2]);
			Fe[3+i] = Area*twelve*(S1[i] + 2*S2[i] + S3[i]) + tau*sixth*(ms2[i][0]*Sb[0] + ms2[i][1]*Sb[1] + ms2[i][2]*Sb[2]);
			Fe[6+i] = Area*twelve*(S1[i] + S2[i] + 2*S3[i]) + tau*sixth*(ms3[i][0]*Sb[0] + ms3[i][1]*Sb[1] + ms3[i][2]*Sb[2]);
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
