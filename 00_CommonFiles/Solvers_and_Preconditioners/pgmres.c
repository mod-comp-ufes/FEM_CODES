#include "solvers.h"

int pgmres (ParametersType *Parameters,	MatrixDataType *MatrixData, FemStructsType *FemStructs,
			FemFunctionsType *FemFunctions,	double *B, double *X)
{
	
	int i, j, k, l, cont, kmax, lmax, neq;
	double normb, eps, r, rho, soma, h1, h2;
	double tol;
	double **u, *u_aux, **h, *h_aux;
	double *c, *s, *y, *e, *z, *v;

	kmax = Parameters->KrylovBasisVectorsQuantity;
	lmax = Parameters->LinearMaxIter;
	neq = Parameters->neq;
	tol = Parameters->SolverTolerance;

	u = (double**) mycalloc("u of 'gmres'",kmax,sizeof(double));
	u_aux = (double*) mycalloc("u_aux of 'gmres'",kmax*(neq+1),sizeof(double));
	for (i = 0; i < kmax; i++)
		u[i] = &u_aux[i*(neq+1)];

	h = (double**) mycalloc("h of 'gmres'",kmax,sizeof(double*));
	h_aux = (double*) mycalloc("h_aux of 'gmres'", kmax*kmax ,sizeof(double*));
	for (i = 0; i < kmax; i++)
		h[i] = &h_aux[i*kmax];

	e = (double*) mycalloc("e of 'gmres'",kmax,sizeof(double));
	c = (double*) mycalloc("c of 'gmres'",kmax,sizeof(double));
	s = (double*) mycalloc("s of 'gmres'",kmax,sizeof(double)); 
	y = (double*) mycalloc("y of 'gmres'",kmax,sizeof(double));
	v = (double*) mycalloc("z of 'gmres'",neq + 1,sizeof(double));
	z = (double*) mycalloc("z of 'gmres'",neq + 1,sizeof(double));

	// Inicializa matrizes e vetores com zero
	dzero(neq, X);

	// Calcula ||b||_2
	normb = ddot(neq,B,B);
	normb = sqrt(normb);

	// Calcula eps = tol*||b||_2 
	eps = tol*normb;

	l = 0;
	cont = 0;

	do
	{
		i = 0;
		for(j = 0; j < kmax; j++)
			dzero(neq, u[j]);
		
		dcopy(neq, X, v);
				
		FemFunctions->precondR(Parameters,MatrixData,FemStructs,v,v);

		// ui = AX
		FemFunctions->mv(Parameters, MatrixData, FemStructs, v, u[i]);

		// Preconditioning M u = z
		FemFunctions->precond(Parameters, MatrixData, FemStructs, u[i], u[i]);

		daxpy(neq, -1.0, B, u[i]);     // ui = ui - b
		dscal(neq, -1.0, u[i]);          // ui = - ui = (b - ui) = (b - Ax)

		// ei = ||ui||_2
		e[i] = ddot(neq, u[i], u[i]);
		e[i] = sqrt(e[i]);

		// ui = ui / ei
		dscal(neq,1.0/e[i],u[i]);

		rho = e[i];


		do
		{
			cont++;
			
			dcopy(neq, u[i], v);

			FemFunctions->precondR(Parameters,MatrixData,FemStructs,v,v);

			// uj = ui+1 = A z
			FemFunctions->mv(Parameters, MatrixData, FemStructs, v, u[i+1]);

			// Preconditioning M ui+1 = z 
			FemFunctions->precond(Parameters, MatrixData, FemStructs, u[i+1], u[i+1]);

			// Ortogonalizacao de Gram-Schmidt
			for (j = 0; j <= i; j++)
			{
				// hji = ui+1*uj (Produto Interno)
				h[j][i] = ddot(neq, u[j], u[i+1]);

				// ui+1 = ui+1 - hji*uj
				daxpy(neq, -h[j][i], u[j], u[i+1]);
			}

			// hi+1,i = ||ui+1||
			h[i+1][i] = ddot(neq, u[i+1], u[i+1]);
			h[i+1][i] = sqrt(h[i+1][i]);

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
			r = h[i][i]*h[i][i] + h[i+1][i]*h[i+1][i];            
			r = sqrt(r);

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
		}
		while ((rho > eps)&&(i < kmax-1));

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

	#ifdef debug 	
		printf(" Iteracoes GMRES: %d \n", cont);
	#endif
	Parameters->iterations += cont;

	FemFunctions->precondR (Parameters, MatrixData, FemStructs, X, X);

	free(u_aux);
	free(u);
	free(h_aux);
	free(h);
	free(e); 
	free(c); 
	free(s);
	free(y);
	free(v);
	free(z);

	return 0;
}



