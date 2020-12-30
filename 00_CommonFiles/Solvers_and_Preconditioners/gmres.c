#include "solvers.h"

int gmres (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, FemFunctionsType *FemFunctions)
{
	
	//----------------------------------------------
	// NÃO FOI IMPLEMENTADO O PRECONDICIONADOR
	//----------------------------------------------
	 
	int i, j, k, l, cont, kmax, lmax, neq;
	double normb, normb2, ui2, uip12, r2, eps, r, rho, soma, h1, h2;
	double tol;
	double **u, *u_aux, **h, *h_aux;
	double *c, *s, *y, *e, *B, *X;
	
	kmax = Parameters->KrylovBasisVectorsQuantity;
	lmax = Parameters->SolverMaxIter;
	neq = Parameters->neq;
	tol = Parameters->SolverTolerance;
	B = FemStructs->F;
	X = FemStructs->u;
	
/*	printf("\nVetor B=F Global dentro do solver\n");
		for(i = 0; i <= neq; i++){
			printf("%d: %lf\n", i, B[i]);
		}
	getchar(); */
	
	u = (double**) mycalloc("u of 'gmres'",kmax+1,sizeof(double));
	u_aux = (double*) mycalloc("u_aux of 'gmres'",(kmax+1)*(neq+1),sizeof(double));
	for (i = 0; i < kmax+1; i++)
		u[i] = &u_aux[i*(neq+1)];
	
	h = (double**) mycalloc("h of 'gmres'",kmax+1,sizeof(double*));
	h_aux = (double*) mycalloc("h_aux of 'gmres'", (kmax+1)*(kmax+1) ,sizeof(double*));
	for (i = 0; i < (kmax+1); i++)
		h[i] = &h_aux[i*(kmax+1)];
		
	e = (double*) mycalloc("e of 'gmres'",kmax+1,sizeof(double));
	c = (double*) mycalloc("c of 'gmres'",kmax,sizeof(double));
	s = (double*) mycalloc("s of 'gmres'",kmax,sizeof(double)); 
	y = (double*) mycalloc("y of 'gmres'",kmax,sizeof(double));

	// Calcula ||b||_2
	normb2 = ddot(neq,B,B);
	normb = sqrt(normb2);
//printf("\t norma B: %lf \n",normb);

	// Calcula eps = tol*||b||_2 
	eps = tol*normb;
//printf("\t eps: %e \n",eps);getchar();

	// Inicializa matrizes e vetores com zero
	dzero(neq+1, X);
   
	l = 0;
	cont = 0;

	do
	{
		//printf("Entrou no primeiro DO!\n");
			
		i = 0;
		for(j = 0; j < kmax+1; j++)
			dzero(neq, u[j]);

		// ui = AX
		FemFunctions->ProductMatrixVector(Parameters, MatrixData, FemStructs, X, u[i]);

	/*	printf("Eq    Vetor    Solucao\n");
		int w;
		for(w = 0; w <= neq; w++){
			printf("%d: %lf\t%lf\n", w, X[w], u[0][w]);
		}
		getchar(); */

		daxpy(neq, -1.0, B, u[i]);     // ui = ui - b
		dscal(neq, -1.0, u[i]);        // ui = - ui = (b - ui) = (b - Ax)
      
		// ei = ||ui||_2
		ui2 = ddot(neq, u[i], u[i]);
		e[i] = sqrt(ui2);

		// ui = ui / ei
		dscal(neq,1.0/e[i],u[i]);

		rho = e[i];

		do
		{
			//printf("Entrou no segundo DO!\n");
			
			cont++;
			
			// uj = ui+1 = Aui
			FemFunctions->ProductMatrixVector(Parameters, MatrixData, FemStructs, u[i], u[i+1]);
	
			// Ortogonalizacao de Gram-Schmidt
			for (j = 0; j <= i; j++)
			{
				// hji = ui+1*uj (Produto Interno)
				h[j][i] = ddot(neq, u[j], u[i+1]);

				// ui+1 = ui+1 - hji*uj
				daxpy(neq, -h[j][i], u[j], u[i+1]);
			}

			// hi+1,i = ||ui+1||
			uip12 = ddot(neq, u[i+1], u[i+1]);
			h[i+1][i] = sqrt(uip12);

			// ui+1 = ui+1 / hi+1,i
			dscal(neq, 1.0/h[i+1][i], u[i+1]);

			// Algoritmo QR                                            
			for (j = 0; j <= i-1; j++)
			{                                     
				// hji = cj*hji + sj*hj+1,i                                          
				h1 =  c[j]*h[j][i] + s[j]*h[j+1][i];                  

				// hj+1,i = -sj*hji + cj*hj+1,i                                          
				h2 = -s[j]*h[j][i] + c[j]*h[j+1][i];

				h[j][i] = h1;
				h[j+1][i] = h2;
			}

			// r = sqrt((hii)^2 + (hi+1,i)^2)
			r2 = pow(h[i][i],2) + pow(h[i+1][i],2);            
			r = sqrt(r2);

			// ci = hii / r
			c[i] = h[i][i] / r;

			// si = hi+1,i / r
			s[i] = h[i+1][i] / r;

			// hii = r
			h[i][i] = r;                                                   

			// hi+1,i = 0.0
			h[i+1][i] = 0.0;

			// ei+1 = -si*ei
			e[i+1] = -s[i]*e[i];

			// ei = ci*ei
			e[i]   =  c[i]*e[i];

			// rho = |ei+1|
			rho = fabs(e[i+1]);  

			i++;
		
		}while ((rho > eps)&&(i < kmax));
		
		i--;

		y[i] = e[i] / h[i][i];

		for (j = i-1; j >= 0; j--)
		{
			soma = 0.0;
			for (k = j+1; k <= i; k++)
				soma = soma + h[j][k]*y[k];
			y[j] = (e[j] - soma)/h[j][j];
		}

		for (j = 0; j <= i; j++)
			for (k = 0; k < neq; k++)
				X[k] = X[k] + u[j][k]*y[j];

		l++;

	}while((rho > eps)&&(l<lmax));

	Parameters->ResGMRES = rho;
	Parameters->ContGMRES = cont;
	Parameters->SolverIterations += cont;
	
	free(u_aux);
	free(u);
	free(h_aux);
	free(h);
	free(e); 
	free(c); 
	free(s);
	free(y);
	
	return 0;
}



