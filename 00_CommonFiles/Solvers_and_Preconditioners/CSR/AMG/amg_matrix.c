# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include "amg_matrix.h"
# include "amg_list.h"
# include "../../../Allocation_Operations/allocations.h"

matrix * create_m (int m, int n, int nnz) {
	matrix *A = (matrix *) mycalloc("A of create_m",1,sizeof(matrix));
	A->m = m; A->n = n; A->nnz = nnz;
	A->val = (double *) mycalloc("A->val of create_m",nnz,sizeof(double));
	A->diag = (double *) mycalloc("A->diag of create_m", n, sizeof(double));
	A->col_ind = (int *) mycalloc("A->col_ind of create_m",nnz,sizeof(int));
	A->row_ptr = (int *) mycalloc("A->row_ptr of create_m",(m+1),sizeof(int));
	return A;
};

matrix * copy_m(matrix * A) {
	matrix * B = create_m(A->n, A->m, A->nnz);
	int i;

	B->m = A->m;
	B->n = A->n;
	B->nnz = A->nnz;

	for(i=0; i<A->nnz; i++) {
		B->val[i] = A->val[i];
		B->col_ind[i] = A->col_ind[i];
	}

	for(i=0; i<A->n; i++)
		B->diag[i] = A->diag[i];

	for(i=0; i<(A->m+1); i++)
		B->row_ptr[i] = A->row_ptr[i];

	return B;
}

void destroy_m (matrix *A) {
	myfree(A->val); myfree(A->diag);
	myfree(A->col_ind); myfree(A->row_ptr);
	myfree(A);
};

void print_m (matrix *A) {
	int i, j, n=A->m;
	for (i=0; i<n; i++) {
		printf("%d:", i+1);
		for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) printf(" %.6f (%d);", A->val[j], A->col_ind[j]+1);
		printf("\n");
	}
};

double get_Aij (matrix *A, int i, int j) {
	int k;
	for (k=A->row_ptr[i]; k<A->row_ptr[i+1]; k++) {
		if (A->col_ind[k] == j) return A->val[k];
	}
	return 0.0;
};

matrix * csr_transpose (matrix *A) {
	int i, j, m=A->m, n=A->n;
	matrix *A_t;
	Node *aux;
	list **row_elem = (list **) mycalloc("row_elem of csr_transpose",n,sizeof(list *));
	for (i=0; i<n; i++) row_elem[i] = create_l();
	for (i=0; i<m; i++) {
		for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) insert_l_tail(row_elem[A->col_ind[j]], i, A->val[j]);
	}
	A_t = create_m(n, m, A->nnz);
	for (A_t->row_ptr[0]=0, j=0, i=0; i<A_t->m; i++) {
		for (aux=row_elem[i]->head; aux; aux=aux->next, j++) {
			A_t->val[j] = aux->val;
			A_t->col_ind[j] = aux->elem;
		}
		A_t->diag[i] = A->diag[i];
		A_t->row_ptr[i+1] = j;
	}
	for (i=0; i<n; i++) destroy_l(row_elem[i]);
	myfree(row_elem);
	return A_t;
};

matrix * sum(matrix *A, matrix *A_t) {
	matrix *A_s;
	Node *aux;
	int i, j, j_t, nnz=0;
	list **row_elem = (list **) mycalloc("row_elem of sum",A->n,sizeof(list *));
	for (i=0; i<A->n; i++) row_elem[i] = create_l();
	for (i=0; i<A->n; i++)
	{		
		j 	= A->	row_ptr[i+1]-1;
		j_t 	= A_t->	row_ptr[i+1]-1;

		while((j >= A->row_ptr[i]) || (j_t >= A_t->row_ptr[i])) { /*se esta na linha*/
			
			if((j >= A->row_ptr[i]) && (j_t >= A_t->row_ptr[i]) && 
                           (A->col_ind[j] == A_t->col_ind[j_t])) {
				insert_l_head(row_elem[i], A->col_ind[j], A->val[j] + A_t->val[j_t]); //fabs(A->val[j]) + fabs(A_t->val[j_t]));
				j--;
				j_t--;
			} else
			if(((j >= A->row_ptr[i]) &&
                            (A->col_ind[j] > A_t->col_ind[j_t])) ||
                           (j_t < A_t->row_ptr[i])) {
				insert_l_head(row_elem[i], A->col_ind[j], A->val[j]); //fabs(A->val[j]));
				j--;
			} else
			{
				insert_l_head(row_elem[i], A_t->col_ind[j_t], A_t->val[j_t]); //fabs(A_t->val[j_t]));
				j_t--;
			}
			
			nnz++;
		}
	}

	A_s = create_m(A->n, A->m, nnz);
	for (A_s->row_ptr[0]=0, j=0, i=0; i<A_s->m; i++) {
		for (aux=row_elem[i]->head; aux; aux=aux->next, j++) {
			A_s->val[j] = aux->val;
			A_s->col_ind[j] = aux->elem;
		}
		//A_s->diag[i] = 0; //nao precisa devido ao calloc
		A_s->row_ptr[i+1] = j;
	}
	for (i=0; i<A->n; i++) destroy_l(row_elem[i]);
	myfree(row_elem);
	return A_s;
}

matrix * sum_abs(matrix *A, matrix *A_t) {
	matrix *A_s;
	Node *aux;
	int i, j, j_t, nnz=0;
	list **row_elem = (list **) mycalloc("row_elem of sum_abs",A->n,sizeof(list *));
	for (i=0; i<A->n; i++) row_elem[i] = create_l();
	for (i=0; i<A->n; i++)
	{		
		j 	= A->	row_ptr[i+1]-1;
		j_t 	= A_t->	row_ptr[i+1]-1;

		while((j >= A->row_ptr[i]) || (j_t >= A_t->row_ptr[i])) { /*se esta na linha*/
			
			if((j >= A->row_ptr[i]) && (j_t >= A_t->row_ptr[i]) && 
                           (A->col_ind[j] == A_t->col_ind[j_t])) {
				insert_l_head(row_elem[i], A->col_ind[j], fabs(A->val[j]) + fabs(A_t->val[j_t]));
				j--;
				j_t--;
			} else
			if(((j >= A->row_ptr[i]) &&
                            (A->col_ind[j] > A_t->col_ind[j_t])) ||
                           (j_t < A_t->row_ptr[i])) {
				insert_l_head(row_elem[i], A->col_ind[j], fabs(A->val[j]));
				j--;
			} else
			{
				insert_l_head(row_elem[i], A_t->col_ind[j_t], fabs(A_t->val[j_t]));
				j_t--;
			}
			
			nnz++;
		}
	}

	A_s = create_m(A->n, A->m, nnz);
	for (A_s->row_ptr[0]=0, j=0, i=0; i<A_s->m; i++) {
		for (aux=row_elem[i]->head; aux; aux=aux->next, j++) {
			A_s->val[j] = aux->val;
			A_s->col_ind[j] = aux->elem;
		}
		//A_s->diag[i] = 0; //nao precisa devido ao calloc
		A_s->row_ptr[i+1] = j;
	}
	for (i=0; i<A->n; i++) destroy_l(row_elem[i]);
	myfree(row_elem);
	return A_s;
}

matrix * sum_upper(matrix *A, matrix *A_t) {
	matrix *A_s;
	Node *aux;
	int i, j, j_t, nnz=0;
	list **row_elem = (list **) mycalloc("row_elem of sum_upper",A->n,sizeof(list *));
	for (i=0; i<A->n; i++) row_elem[i] = create_l();
	for (i=0; i<A->n; i++)
	{		
		j 	= A->	row_ptr[i+1]-1;
		j_t 	= A_t->	row_ptr[i+1]-1;

		while(((j >= A->row_ptr[i]) || (j_t >= A_t->row_ptr[i])) && /*se esta na linha*/
		((A->col_ind[j] > i) || (A_t->col_ind[j_t] > i))) { /*se esta na triangular superior*/
			
			if(A->col_ind[j] == A_t->col_ind[j_t]) {
				insert_l_head(row_elem[i], A->col_ind[j], fabs(A->val[j]) + fabs(A_t->val[j_t]));
				j--;
				j_t--;
			} else
			if(A->col_ind[j] > A_t->col_ind[j_t]) {
				insert_l_head(row_elem[i], A->col_ind[j], fabs(A->val[j]));
				j--;
			} else
			{
				insert_l_head(row_elem[i], A_t->col_ind[j_t], fabs(A_t->val[j_t]));
				j_t--;
			}
			
			nnz++;
		}
	}

	A_s = create_m(A->n, A->m, nnz);
	for (A_s->row_ptr[0]=0, j=0, i=0; i<A_s->m; i++) {
		for (aux=row_elem[i]->head; aux; aux=aux->next, j++) {
			A_s->val[j] = aux->val;
			A_s->col_ind[j] = aux->elem;
		}
		//A_s->diag[i] = 0; //nao precisa devido ao calloc
		A_s->row_ptr[i+1] = j;
	}
	for (i=0; i<A->n; i++) destroy_l(row_elem[i]);
	myfree(row_elem);
	return A_s;
}

double maxCoeff_abs (matrix *A) {
	int i, j, n=A->m;
	double max = 0.0;
	for (i=0; i<n; i++) {
		for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) { if (fabs(A->val[j])>max) max = fabs(A->val[j]); }
	}
	return max;
};

matrix * matmat_sparseproduct (matrix *A, matrix *B) {
	int *xb, p=A->m, r=B->n, ip, i, jp, j, kp, k;
	double *x;
	matrix *C;
	Node *aux;
	list **row_elem = (list **) mycalloc("row_elem of matmat_sparseproduct",p,sizeof(list *));
	xb = (int *) mycalloc("xb of matmat_sparseproduct",r, sizeof(int));
	x = (double *) mycalloc("x of matmat_sparseproduct",r,sizeof(double));
	ip = 0;
	for (i=0; i<p; i++) { // each iteration fills a row in matrix C
		row_elem[i] = create_l();
		for (jp=A->row_ptr[i]; jp<A->row_ptr[i+1]; jp++) {
			j = A->col_ind[jp];
			for (kp=B->row_ptr[j]; kp<B->row_ptr[j+1]; kp++) {
				k = B->col_ind[kp];
				if (xb[k]!=i+1) {
					insert_l_tail(row_elem[i], k, 0.0);
					ip++;
					xb[k] = i+1;
					x[k] = A->val[jp]*B->val[kp];
				} else x[k] += A->val[jp]*B->val[kp];
			}
		}
		for (aux=row_elem[i]->head; aux; aux=aux->next) aux->val = x[aux->elem];
	}
	myfree(xb); myfree(x);
	C = create_m(p, r, ip);
	for (C->row_ptr[0]=0, j=0, i=0; i<p; i++) {
		for (aux=row_elem[i]->head; aux; aux=aux->next, j++) {
			C->val[j] = aux->val;
			C->col_ind[j] = aux->elem; if (aux->elem==i) C->diag[i] = aux->val;
		}
		C->row_ptr[i+1] = j;
	}
	for (i=0; i<p; i++) destroy_l(row_elem[i]);
	myfree(row_elem);
	return C;
};

void mat_vec (double *p, matrix *A, double *v) {
	int i, j, n=A->m;
	for (i=0; i<n; i++) {
		p[i] = 0.0;
		for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) p[i] += A->val[j]*v[A->col_ind[j]];
	}
};

void mat_vec_plus (double *p, matrix *A, double *v) {
	int i, j, n=A->m;
	for (i=0; i<n; i++) {
		for (j=A->row_ptr[i]; j<A->row_ptr[i+1]; j++) p[i] += A->val[j]*v[A->col_ind[j]];
	}
};

