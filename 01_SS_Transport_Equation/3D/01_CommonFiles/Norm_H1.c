#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/BLAS_Operations/ourBLAS.h"
#include "../../../00_CommonFiles/Allocation_Operations/allocations.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

void Norm_H1 (double *U, ParametersType *Parameters, FemStructsType *FemStructs, FemFunctionsType *FemFunctions){
	
	int e, i, nel, J1, J2, J3, J4;
	double detJ, Volume, DUeDx = 0.0, DUeDy= 0.0, DUeDz=0.0, DUhDx, DUhDy, DUhDz, Xi[5][3], w[5];
	double x[4], y[4], z[4], X, Y, Z, erroX, erroY, erroZ, SomaInt, normL2;
	double b1, b2, b3, b4, c1, c2, c3 ,c4, d1, d2, d3, d4;
	double kappa, sigma, aux=0.0;
	double Cst = Parameters->ConstApli;
	NodeType *Node = FemStructs->Node;
	ElementType *Element = FemStructs->Element;
	CoefFormFuncType *CFF = FemStructs->CFF;
	
	FemFunctions->kappa(&kappa, &aux, &aux, aux, aux, aux, Cst);
	sigma = FemFunctions->sigma(aux, aux, aux);
	
	nel = Parameters->nel;
	
	normL2 = Parameters->NormL2;
	
	// Fills in with the gaussian quadrature weights and points
	Xi[0][0] = Xi[0][1] = Xi[0][2] = 1.0/6.0;  //ponto 1
	Xi[1][0] = Xi[1][1] = 1.0/6.0; Xi[1][2] = 1.0/2.0;  //ponto 2
	Xi[2][0] = 1.0/6.0;  Xi[2][1] = 1.0/2.0; Xi[2][2] = 1.0/6.0; //ponto 3
	Xi[3][0] = 1.0/2.0;  Xi[3][1] = Xi[3][2] = 1.0/6.0;  //ponto 4
	Xi[4][0] = Xi[4][1] = Xi[4][2] = 1.0/4.0;  //ponto 5
	
	w[0] = w[1] = w[2] = w[3] = 9.0/20.0; // 3.0/40.0;
	w[4] = -16.0/20.0; // -4.0/30.0; 
	
	SomaInt = 0.0;
	
	for(e = 0; e < nel; e++){
		Volume = CFF[e].volume;
		
		b1 = CFF[e].b[0];
		b2 = CFF[e].b[1];
		b3 = CFF[e].b[2];
		b4 = CFF[e].b[3];
		c1 = CFF[e].c[0];
		c2 = CFF[e].c[1];
		c3 = CFF[e].c[2];
		c4 = CFF[e].c[3];
		d1 = CFF[e].d[0];
		d2 = CFF[e].d[1];
		d3 = CFF[e].d[2];
		d4 = CFF[e].d[3];
		
		J1 = Element[e].Vertex[0];
		J2 = Element[e].Vertex[1];
		J3 = Element[e].Vertex[2];
		J4 = Element[e].Vertex[3];
		
		// Nodal coordinates
		x[0] = Node[J1].x;
		x[1] = Node[J2].x;
		x[2] = Node[J3].x;
		x[3] = Node[J4].x;
		
		y[0] = Node[J1].y;
		y[1] = Node[J2].y;
		y[2] = Node[J3].y;
		y[3] = Node[J4].y;
		
		z[0] = Node[J1].z;
		z[1] = Node[J2].z;
		z[2] = Node[J3].z;
		z[3] = Node[J4].z;
		
/*		x21 = x[1] - x[0];
		x31 = x[2] - x[0];
		x41 = x[3] - x[0];
		y21 = y[1] - y[0];
		y31 = y[2] - y[0];
		y41 = y[3] - y[0];
		z21 = z[1] - z[0];
		z31 = z[2] - z[0];
		z41 = z[3] - z[0];   */
		
		detJ = 6.0*Volume;
		
		//--------------------------------------------------------------
/*		double **C, **invC, **M, **gradX;
		
		C = (double**) mycalloc("C line of 'Norm_H1'", 4, sizeof(double*));
		invC = (double**) mycalloc("invC line of 'Norm_H1'", 4, sizeof(double*));
		M = (double**) mycalloc("M line of 'Norm_H1'", 4, sizeof(double*));
		gradX = (double**) mycalloc("gradX line of 'Norm_H1'", 4, sizeof(double*));
		for (i = 0; i < 4; i++){
			C[i] = (double*) mycalloc("C colum of 'Norm_H1'", 4, sizeof(double));
			invC[i] = (double*) mycalloc("invC colum of 'Norm_H1'", 4, sizeof(double));
			M[i] = (double*) mycalloc("M colum of 'Norm_H1'", 3, sizeof(double));
			gradX[i] = (double*) mycalloc("gradX colum of 'Norm_H1'", 3, sizeof(double));
		}
		
		C[0][0] = 1.0;  C[0][1] = 1.0;  C[0][2] = 1.0;  C[0][3] = 1.0;
		C[1][0] = x[0]; C[1][1] = x[1]; C[1][2] = x[2]; C[1][3] = x[3];
		C[2][0] = y[0]; C[2][1] = y[1]; C[2][2] = y[2]; C[2][3] = y[3];
		C[3][0] = z[0]; C[3][1] = z[1]; C[3][2] = z[2]; C[3][3] = z[3];
		
		inverse_matrix(C, invC, 4);
		
		M[0][0] = 0.0; M[0][1] = 0.0; M[0][2] = 0.0;
		M[1][0] = 1.0; M[1][1] = 0.0; M[1][2] = 0.0;
		M[2][0] = 0.0; M[2][1] = 1.0; M[2][2] = 0.0;
		M[3][0] = 0.0; M[3][1] = 0.0; M[3][2] = 1.0;
		
		matrix_multiplication(invC, 4, 4, M, 4, 3, gradX); 
		
		printf("Elemento %d\n", e+1);
		printf("gradX[0][0] = %lf  DN1/dx = %lf\n", gradX[0][0], CFF[e].b[0]/detJ);
		printf("gradX[0][1] = %lf  DN1/dy = %lf\n", gradX[0][1], CFF[e].c[0]/detJ);
		printf("gradX[0][2] = %lf  DN1/dz = %lf\n\n", gradX[0][2], CFF[e].d[0]/detJ);
		
		printf("gradX[1][0] = %lf  DN2/dx = %lf\n", gradX[1][0], CFF[e].b[1]/detJ);
		printf("gradX[1][1] = %lf  DN2/dy = %lf\n", gradX[1][1], CFF[e].c[1]/detJ);
		printf("gradX[1][2] = %lf  DN2/dz = %lf\n\n", gradX[1][2], CFF[e].d[1]/detJ);
		
		printf("gradX[2][0] = %lf  DN3/dx = %lf\n", gradX[2][0], CFF[e].b[2]/detJ);
		printf("gradX[2][1] = %lf  DN3/dy = %lf\n", gradX[2][1], CFF[e].c[2]/detJ);
		printf("gradX[2][2] = %lf  DN3/dz = %lf\n\n", gradX[2][2], CFF[e].d[2]/detJ);
		
		printf("gradX[3][0] = %lf  DN4/dx = %lf\n", gradX[3][0], CFF[e].b[3]/detJ);
		printf("gradX[3][1] = %lf  DN4/dy = %lf\n", gradX[3][1], CFF[e].c[3]/detJ);
		printf("gradX[3][2] = %lf  DN4/dz = %lf\n\n", gradX[3][2], CFF[e].d[3]/detJ);
		getchar(); */
		
		// gradient form function of default element considering N1 = 1 - xi - eta - zeta, N2 = xi, N3 = eta, N4 = zeta
/*		gradXi[0][0] = -1.0; gradXi[0][1] = 1.0; gradXi[0][2] = 0.0; gradXi[0][3] = 0.0;
		gradXi[1][0] = -1.0; gradXi[1][1] = 0.0; gradXi[1][2] = 1.0; gradXi[1][3] = 0.0;
		gradXi[2][0] = -1.0; gradXi[2][1] = 0.0; gradXi[2][2] = 0.0; gradXi[2][3] = 1.0;
		
		// Jacobian Matriz
		J[0][0] = x21; J[0][1] = x31; J[0][2] = x41;
		J[1][0] = y21; J[1][1] = y31; J[1][2] = y41;
		J[2][0] = z21; J[2][1] = z31; J[2][2] = z41;*/
		
/*		double invJt[3][3];
		invJt[0][0] = ((y31*z41) - (z31*y41))/detJ;
		invJt[1][0] = (- (x31*z41) + (z31*x41))/detJ;
		invJt[2][0] = ((x31*y41) - (y31*x41))/detJ; 
		invJt[0][1] = (- (y21*z41) + (z21*y41))/detJ; 
		invJt[1][1] = ((x21*z41) - (z21*x41))/detJ;
		invJt[2][1] = (- (x21*y41) + (y21*x41))/detJ;
		invJt[0][2] = ((y21*z31) - (z21*y31))/detJ;
		invJt[1][2] = (- (x21*z31) + (z21*x31))/detJ; 
		invJt[2][2] = ((x21*y31) - (x31*y21))/detJ;*/
		
	/*	for(i=0;i<3;i++){
			printf("%lf\t%lf\t%lf\n",invJt[i][0], invJt[i][1], invJt[i][2]);
		}
		getchar(); */
		
		// gradient form function of real element
/*		gradX[0][0] = invJt[0][0]*gradXi[0][0] + invJt[0][1]*gradXi[1][0] + invJt[0][2]*gradXi[2][0]; 
		gradX[1][0] = invJt[1][0]*gradXi[0][0] + invJt[1][1]*gradXi[1][0] + invJt[1][2]*gradXi[2][0];
		gradX[2][0] = invJt[2][0]*gradXi[0][0] + invJt[2][1]*gradXi[1][0] + invJt[2][2]*gradXi[2][0];
                                                                                        
		gradX[0][1] = invJt[0][0]*gradXi[0][1] + invJt[0][1]*gradXi[1][1] + invJt[0][2]*gradXi[2][1]; 
		gradX[1][1] = invJt[1][0]*gradXi[0][1] + invJt[1][1]*gradXi[1][1] + invJt[1][2]*gradXi[2][1];
		gradX[2][1] = invJt[2][0]*gradXi[0][1] + invJt[2][1]*gradXi[1][1] + invJt[2][2]*gradXi[2][1]; 
                                                                                         
		gradX[0][2] = invJt[0][0]*gradXi[0][2] + invJt[0][1]*gradXi[1][2] + invJt[0][2]*gradXi[2][2];
		gradX[1][2] = invJt[1][0]*gradXi[0][2] + invJt[1][1]*gradXi[1][2] + invJt[1][2]*gradXi[2][2];
		gradX[2][2] = invJt[2][0]*gradXi[0][2] + invJt[2][1]*gradXi[1][2] + invJt[2][2]*gradXi[2][2];
                                                                                       
		gradX[0][3] = invJt[0][0]*gradXi[0][3] + invJt[0][1]*gradXi[1][3] + invJt[0][2]*gradXi[2][3];
		gradX[1][3] = invJt[1][0]*gradXi[0][3] + invJt[1][1]*gradXi[1][3] + invJt[1][2]*gradXi[2][3];
		gradX[2][3] = invJt[2][0]*gradXi[0][3] + invJt[2][1]*gradXi[1][3] + invJt[2][2]*gradXi[2][3]; */
		
	/*	printf("gradX[0][0] = %lf \t DN1/dx = %lf\n", gradX[0][0], CFF[e].b[0]/detJ);
		printf("gradX[1][0] = %lf \t DN1/dx = %lf\n", gradX[1][0], CFF[e].c[0]/detJ);
		printf("gradX[2][0] = %lf \t DN1/dx = %lf\n", gradX[2][0], CFF[e].d[0]/detJ);
		getchar();*/
		
		//--------------------------------------------------------------
		
		for(i = 0; i < 5; i++){
			// Transformed coordinates
			X = (1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2])*x[0] + Xi[i][0]*x[1] + Xi[i][1]*x[2] + Xi[i][2]*x[3];
			Y = (1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2])*y[0] + Xi[i][0]*y[1] + Xi[i][1]*y[2] + Xi[i][2]*y[3];
			Z = (1.0 - Xi[i][0] - Xi[i][1] - Xi[i][2])*z[0] + Xi[i][0]*z[1] + Xi[i][1]*z[2] + Xi[i][2]*z[3];
			
			// Derived from the exact solucion
			DUeDx = FemFunctions->DuDx(X, Y, Z, Cst);
			DUeDy = FemFunctions->DuDy(X, Y, Z, Cst);
			DUeDz = FemFunctions->DuDz(X, Y, Z, Cst);
			
			// Derived from aproximate solution
			DUhDx = (b1*U[J1] + b2*U[J2] + b3*U[J3] + b4*U[J4])/detJ;
			DUhDy = (c1*U[J1] + c2*U[J2] + c3*U[J3] + c4*U[J4])/detJ;
			DUhDz = (d1*U[J1] + d2*U[J2] + d3*U[J3] + d4*U[J4])/detJ;
			
			// Local Error
			erroX = DUeDx - DUhDx;
			erroY = DUeDy - DUhDy;
			erroZ = DUeDz - DUhDz;
			
			SomaInt = SomaInt + Volume*w[i]*(pow(erroX,2.0) + pow(erroY,2.0) + pow(erroZ,2.0)); // pow norm grad error
		
		}//end for integração
		
		// 1ª tentativa
	/*	for(i = 0; i < 4; i++){
			myfree(C[i]);
			myfree(invC[i]);
			myfree(M[i]);
			myfree(gradX[i]);
		}
		myfree(C);
		myfree(invC);
		myfree(M);
		myfree(gradX); */
	
	}//end for elemento
	
	Parameters->SemiNormH1 = sqrt(SomaInt);
	Parameters->NormH1 = sqrt(pow(normL2,2) + SomaInt);
	Parameters->NormEnergy = sqrt(kappa*SomaInt + sigma*pow(normL2,2));
	
	FILE *OutFile;
	char FileName[2000];
	
	sprintf(FileName,"../03_Output/%s/%s_%s_%s_ExecutionData_N%d_E%d.txt", Parameters->ProblemTitle, Parameters->Experiments, Parameters->ProblemTitle, Parameters->StabilizationForm, Parameters->nnodes, Parameters->nel);
	OutFile = myfopen(FileName,"a");
	fprintf(OutFile, "pow||E||L2 = %e \t pow||gradE||L2 = %e\n", pow(normL2,2), SomaInt);	
	printf("\n\npow||E||L2 = %e \t pow||gradE||L2 = %e\n\n", pow(normL2,2), SomaInt);
	fclose(OutFile);
	
}//end Norm_H1
