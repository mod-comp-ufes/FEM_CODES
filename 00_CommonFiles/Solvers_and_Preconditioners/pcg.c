#include "solvers.h"
#include "preconditioners.h"

int pcg (ParametersType *Parameters,	MatrixDataType *MatrixData, FemStructsType *FemStructs,
			FemFunctionsType *FemFunctions)
{
	int iter, itermax, neq;
	double *B, *X, *R, *P, *Q, gammaOld, gammaNew, tau, alpha, beta, r0, Tol;

	B = FemStructs->F;
	X = FemStructs->u; // x no alg
	neq = Parameters->neq;

	// Some CG inicializations
	itermax = Parameters->SolverMaxIter;
	Tol = Parameters->SolverTolerance;

	/*************************************************************/
	//		Memory allocations for CG 
	/************************************************************/
	P = mycalloc("P of 'cg'", (neq+1),sizeof(double)); // d no alg
	Q = mycalloc("Q of 'cg'", (neq+1),sizeof(double)); // v no alg
	R = mycalloc("R of 'cg'", neq,sizeof(double)); // r no alg
	//Z = mycalloc("Z of 'cg'", neq,sizeof(double));
	/************************************************************/

	dzero(neq,X); // x0 = 0

	//memcpy(R, B, neq*sizeof(double)); //R <-- B
	dcopy(neq,B,R); // r0 = b
	dcopy(neq,B,P); // d0 = b
	
	//FemFunctions->precond(Parameters, MatrixData, FemStructs, R, Z); //precond(Parameters, MatrixData, R, Z); // Preconditioning Z = invM * R

	//memcpy(P, Z, neq*sizeof(double)); //P <-- Z
	//dcopy(neq,R,P);
	
	gammaNew = ddot(neq, R, R); // gammaNew = R dot R
	
	r0 = gammaNew; // gamma0 no alg

//	if (gammaNew < Tol*Tol*r0)
//		return 0;
	
	while((gammaNew > Tol*Tol*r0)||(iter <= itermax)){
	//for (iter=0; iter<itermax; iter++){
	
		FemFunctions->ProductMatrixVector(Parameters, MatrixData, FemStructs, P, Q); // Q <-- A*P     vi = Adi
	
		tau = ddot(neq, P, Q); // di^t*vi
	
		alpha = gammaNew/tau; // lambda_i = gammaNew/di^t*vi

		daxpy(neq, alpha, P, X); // X <-- X + alpha*P

		daxpy(neq, -alpha, Q, R); // R <-- R - alpha*Q

		//FemFunctions->precond(Parameters, MatrixData, FemStructs, R, Z); //precond(Parameters, MatrixData, R, Z); // Preconditioning Z = invM * R

		gammaOld = gammaNew;

		gammaNew = ddot(neq, R, R); // ri dot ri 

//		#ifdef debug			
//			printf("Gamma = %.15lf (iter = %d)\n", gammaNew,iter);
//		#endif

//		if (gammaNew < Tol*Tol*r0)	
//			break;
		
		beta = gammaNew/gammaOld;
		
		// di+1 = ri+1 + betai+1*di 
		dscal(neq, beta, P); // d = beta*d
		daxpy(neq, 1.0, R, P); // d = d + r
		iter++;
	}

	#ifdef debug
		printf("cg iterations = %d\n",iter);
	#endif

	Parameters->SolverIterations = iter;

	free(P);
	free(Q);
	free(R);
	//free(Z);
	return 0;
}


