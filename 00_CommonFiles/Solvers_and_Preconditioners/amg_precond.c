# include "CSR/AMG/amg.h"
# include "CSR/AMG/amg_dpa.h"
# include "CSR/AMG/amg_util.h"
# include "amg_precond.h"
# include "preconditioners.h"
# include "../Allocation_Operations/allocations.h"

void writeFile(matrix* A, double *F) {
    FILE *fp;
    fp = fopen("matrix.mtx", "w+");
    fprintf(fp,"%d %d %d\n",A->n,A->n,A->nnz);
    for(int i=0; i<A->n; i++)
        for(int j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
	   fprintf(fp,"%d %d %la\n",A->col_ind[j]+1,i+1,A->val[j]);
    fclose(fp);

    fp = fopen("matrix.F", "w+");
    fprintf(fp,"%d\n",A->n);
    for(int i=0; i<A->n; i++)
	fprintf(fp,"%la\n",F[i]);
    fclose(fp);

    exit (0);
}

int AMG_precond_setup (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, int tag, double *F) {
    
    AMG_precond_data *data = MatrixData->amg_precond_data = (AMG_precond_data*)mycalloc("AMG_precond_data in AMG_precond_setup",1,sizeof(AMG_precond_data));
    int     n = Parameters->neq;

    int     precond;
    int     NCL;
    double  str_thr;
    double  omega;
    int     nr;
    int     aggress;
    int     refinePass;
    int     truncq;
    double  trunc_fact;
    precondAMG *AMG_data = NULL;
    precondDPA *DPA_data = NULL;

    sscanf(Parameters->Preconditioner,"AMG %d %d %lf %lf %d",&precond,&NCL,&str_thr,&omega,&nr);
    if (precond == 1) {
    	sscanf(Parameters->Preconditioner,"AMG %*d %*d %*f %*f %*d %d %d %d %lf",&aggress, &refinePass, &truncq, &trunc_fact);
        data->AMG_data = AMG_data  = (precondAMG*)mycalloc("AMG_data in AMG_precond_setup",1,sizeof(precondAMG));
		data->AMG_data->aggress    = aggress;
		data->AMG_data->refinePass = refinePass;
		data->AMG_data->truncq     = truncq;
		data->AMG_data->trunc_fact = trunc_fact;
    }
    else if (precond == 2) {
        data->DPA_data = DPA_data = (precondDPA*)mycalloc("DPA_data in AMG_precond_setup",1,sizeof(precondDPA));
    }

    data->precond    = precond;
    data->NCL        = NCL;
    data->str_thr    = str_thr;
    data->omega      = omega;
    data->nr         = nr;

    matrix* A   = (matrix*)mycalloc("Matrix A in AMG_precond_setup",1,sizeof(matrix)); //se usar o create, criará vetores vazios desnecessários
    A->val      = MatrixData->AA;
    A->col_ind  = MatrixData->JA;
    A->row_ptr  = MatrixData->IA;
    A->n        = Parameters->neq;
    A->m        = Parameters->neq;
    A->nnz      = Parameters->nnzero;
    A->diag     = (double*)mycalloc("A->diag in AMG_precond_setup",n,sizeof(double));
    for(int i = 0; i < n; i++)
        for(int j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++)
            if(A->col_ind[j] == i)
                A->diag[i] = A->val[j];

    //writeFile(A,F); exit(0);

    double  *f2;
    //double  *u      = (double *)mycalloc("u in AMG_precond_setup",A->n, sizeof(double));

    if (precond == 1) {
        AMG_setup(A, F, NULL, NCL, str_thr, AMG_data->aggress, AMG_data->refinePass, AMG_data->truncq, AMG_data->trunc_fact, &(data->A), &(data->f), &(data->u), &(data->r), &(AMG_data->I_cf), &(AMG_data->I_cf_t));
    } else
    if (precond == 2) {
        dpa_setup(A, F, NULL, NCL, &(data->A), &(data->f), &(data->u), &(data->r), &(DPA_data->match), str_thr);
    }
    if (precond) { // modify right-hand side
        f2 = (double *) mycalloc("f2 in AMG_precond_setup",n, sizeof (double));
        precond_sol(F, f2, precond, data);
        for (int i = 0; i < n; i++) F[i] = f2[i];
        myfree(f2);
    }

    return 0;
}

int AMG_precond (ParametersType *Parameters, MatrixDataType *MatrixData, FemStructsType *FemStructs, double *p, double *z) {
    int n = MatrixData->amg_precond_data->A[0]->n;
    double* aux = (double*)mycalloc("aux in AMG_precond",n,sizeof(double));
    for (int i = 0; i < n; i++)
        aux[i] = p[i];
    precond_sol(aux, z, MatrixData->amg_precond_data->precond, MatrixData->amg_precond_data);
    myfree(aux);
    return 0;
}

void AMG_precond_data_destroy(AMG_precond_data *data) {
	if (data==NULL)
		return;

	if(data->precond == 1) {
		AMG_destroy(data->NCL, data->A, data->f, data->u, data->r, data->AMG_data->I_cf, data->AMG_data->I_cf_t);
		myfree(data->AMG_data);
		
	}
	else if(data->precond == 2) {
		dpa_destroy(data->NCL, data->A, data->f, data->u, data->r, data->DPA_data->match);
		myfree(data->DPA_data);
	}

	myfree(data);
}
