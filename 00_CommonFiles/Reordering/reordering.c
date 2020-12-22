#include "reordering.h"

void reordering(ParametersType *Parameters, int *JA, int *IA, int *perm, int *PermCSR)
{
	Parameters->bandwidth_bef = MATRIX_bandwidth(Parameters,JA,IA);
	Parameters->bandwidth_aft = Parameters->bandwidth_bef;

	if (strcasecmp(Parameters->reordering,"NOT")==0){
		return;
	}
/*	else if (strcasecmp(Parameters->reordering,"Spectral")==0)
		REORDERING_SPECTRAL (Parameters, JA,  IA, perm, PermCSR);
*/
	else if (strcasecmp(Parameters->reordering,"SYMRCM")==0)
		REORDERING_SYMRCM (Parameters, JA,  IA, perm, PermCSR);

/*	else if (strncmp(Parameters->reordering,"Sloan",5)==0 || strncmp(Parameters->reordering,"RCM",3)==0){
		REORDERING_RCM_or_SLOAN (Parameters, JA,  IA, perm, PermCSR);
	}*/
	else{
		printf("Reordering scheme not defined!\n");
		exit(1);
	}
	Parameters->bandwidth_aft = MATRIX_bandwidth(Parameters,JA,IA);
}

/*----------------------------------------------------------------------------
 * SPECTRAL reordering
 *--------------------------------------------------------------------------*/
/*void REORDERING_SPECTRAL (ParametersType *Parameters, int *ja, int *ia, int *p, int *pT)
{
	int i;
	int n = Parameters->neq;
	
	double *fvector = calloc (n ,sizeof(double));
	int    *list    = calloc (n ,sizeof(int));
	int    *info    = calloc (10,sizeof(int));
	int   lirn = Parameters->nnzero;
	int    nnz = Parameters->nnzero;
*/	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
/*	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] += 1;
	}

	mc73_fiedler(&n,&lirn,ja,ia,list,fvector,info,NULL);
*/	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
/*	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] -= 1;
	}
	
	ARRAY2* R = malloc(n*sizeof(ARRAY2));
	for (i = 0; i < n; ++i)
	{
		R[i].arr1 = fvector[i];
		R[i].arr2 = i;		
	}	
	
	qsort (R,n,sizeof(ARRAY2),COMPARE_eig); 
	
	for (i = 0; i < n; ++i){ 
		p[i] = R[i].arr2; 
	}

	free(R);
	free(info);
	free(list);
	free(fvector);
	
	MATRIX_ROW_permutation (Parameters, ja, ia, p, pT);
	MATRIX_COL_permutation (Parameters, ja, ia, p, pT);
	
	return;
}
*/
/*----------------------------------------------------------------------------
 * RCM reordering
 *--------------------------------------------------------------------------*/
/*void REORDERING_RCM_or_SLOAN (ParametersType *Parameters, int *ja, int *ia, int *p, int *pT)
{
	int 	i;
	int 	n = Parameters->neq;
	int 	nnz = Parameters->nnzero;
	int 	lirn = nnz;
	int 	nsup = n;
	int    	*svar    = calloc (n,sizeof(int));
 	int    	*vars    = calloc (n,sizeof(int));
	int    	*info    = calloc (4,sizeof(int));
	int   	*icptr	= calloc (n+1,sizeof(int));	
	int   	*irn	= calloc (lirn,sizeof(int));	
	int   	jcntl[2];		#include "../../../00_CommonFiles/MatrixVector_Operations/matvec.h"
	int    	*iw    = calloc (2*n+2,sizeof(int));
	int 	control;
	double  weight[2];		
	double 	*rinfo    = calloc (4,sizeof(double));
*/
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
/*	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] += 1;
	}

	if (strncmp(Parameters->reordering,"Sloan",5)==0){
		jcntl[0]=0; //Sloan reordering
		control = atoi(&(Parameters->reordering[5])); 
		if (control==0)
			jcntl[1]=0; //Automatic choice of pseudoperipheral pairs
		else if (control==1)
			jcntl[1]=1; //Pseudoperipheral pairs specified in PAIR
		else if (control==2)
			jcntl[1]=2; //Global priority vector given in PERMSV (Sloan's algorithm only).
		else{
			printf("Sloan control algorithm not defined!\n");
			exit(1);
		}
		
	}
	else if (strncmp(Parameters->reordering,"RCM",3)==0){
		jcntl[0]=1; //RCM reordering 
		control = atoi(&(Parameters->reordering[3])); 
		if (control==0)
			jcntl[1]=0; //Automatic choice of pseudoperipheral pairs
		else if (control==1) 
			jcntl[1]=1; //Pseudoperipheral pairs specified in PAIR
		else{
			printf("RCM control algorithm not defined!\n");
			exit(1);
		}

	}

	weight[0] = 2.0;
	weight[1] = 1.0;

//	mc60bd_(&n, &lirn,irn,icptr, &nsup, svar, vars, iw);

	free(iw);
	nsup = n;
	
	int    *possv    = calloc (nsup,sizeof(int));
	int    *permsv    = calloc (nsup,sizeof(int));
	int    **pair;
	double  *w = calloc(nsup,sizeof(double));		
	iw    = calloc (3*nsup+1,sizeof(int));

	pair = calloc(2,sizeof(int*));
	pair[0] = calloc(nsup/2,sizeof(int));
	pair[1] = calloc(nsup/2,sizeof(int));
	
	for (i=0; i<n; i++)
		vars[i] = 1;


	mc60cd_(&n,&nsup,&lirn,irn,icptr,vars,jcntl,p,weight,pair,info,iw,w);

	mc60fd_(&n,&nsup,&lirn,irn,icptr,vars,p,iw,rinfo); 


//	mc60dd_(&n,&nsup,svar,vars,permsv,p,possv);
*/
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
/*	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] -= 1;
	}

	for (i = 0; i < n; i++){ 
		p[i] -= 1;
		printf("p[%d]=%d\n",i,p[i]);
	}
	printf("rinfo[0]=%lf\n",rinfo[0]);
	printf("rinfo[1]=%lf\n",rinfo[1]);
	printf("rinfo[2]=%lf\n",rinfo[2]);
	printf("rinfo[3]=%lf\n",rinfo[3]);


	free(info);
	free(rinfo);
	free(vars);
	free(svar);
	free(permsv);
	free(possv);
	free(pair[0]);
	free(pair[1]);
	free(pair);
	free(iw);
	free(w);
	free(icptr);
	free(irn);	

	MATRIX_ROW_permutation (Parameters, ja, ia, p, pT);
	MATRIX_COL_permutation (Parameters, ja, ia, p, pT);

	return;
}
*/
int COMPARE_eig (const void * a, const void * b)
{
	if (((ARRAY2*)a)->arr1 > ((ARRAY2*)b)->arr1) return  1;
	if (((ARRAY2*)a)->arr1 < ((ARRAY2*)b)->arr1) return -1;
	return 0;
}

int COMPARE_array (const void * a, const void * b)
{
	if (((ARRAY*)a)->arr3 <  ((ARRAY*)b)->arr3) return -1;
	if (((ARRAY*)a)->arr3 >  ((ARRAY*)b)->arr3) return  1;
	if (((ARRAY*)a)->arr3 == ((ARRAY*)b)->arr3)
	{
		if (((ARRAY*)a)->arr2 < ((ARRAY*)b)->arr2) return -1;
		if (((ARRAY*)a)->arr2 > ((ARRAY*)b)->arr2) return  1;
	}
	return 0;
}

/*----------------------------------------------------------------------------
 * Perform the colunm permutation
 *--------------------------------------------------------------------------*/
void MATRIX_COL_permutation (ParametersType * Parameters, int *JA, int *IA, int *p, int *pT)
{
	int i, j, k;
	int n   = Parameters->neq;
	int nz  = Parameters->nnzero;  

	ARRAY* a = calloc (nz,sizeof(ARRAY));
	int*   q = calloc (n ,sizeof(int));
	
	for (i = 0; i < n; ++i) 
		q[p[i]] = i; 

	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = IA[i]; j <= IA[i+1] - 1; ++j)
		{
			a[k].arr1 = pT[j];
			a[k].arr2 = q[JA[j]];
			a[k].arr3 = i;
				k = k + 1;
		}
		IA[i+1] = k;    
	}

	qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	
	for (i = 0; i < nz; ++i)
	{
		pT[i] = a[i].arr1;
		JA[i] = a[i].arr2;
	}

	free(a);
	free(q);
}

/*----------------------------------------------------------------------------
 * Perform the colunm permutation
 *--------------------------------------------------------------------------*/
void MATRIX_ROW_permutation (ParametersType * Parameters, int *JA, int *IA, int *p, int *pT)
{
	int i, j, k;
	int n   = Parameters->neq;
	int nz  = Parameters->nnzero;  

	int* 	auxpT = malloc( nz  *sizeof(int));
	int*    auxJA = malloc( nz  *sizeof(int));
	int*    auxIA = malloc((n+1)*sizeof(int));
  
	auxIA[0] = 0;
	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = IA[p[i]]; j <= IA[p[i]+1] - 1; ++j)
		{
			auxpT[k] = pT[j];
			auxJA[k] = JA[j];
			      k  = k + 1;
		}
		auxIA[i+1] = k;    
	}

	memcpy(&pT[0],&auxpT[0],nz*sizeof(int));
	memcpy(&JA[0],&auxJA[0],nz*sizeof(int));
	memcpy(&IA[0],&auxIA[0],(n+1)*sizeof(int));

	free(auxpT);
	free(auxJA);
	free(auxIA);
}

int MATRIX_bandwidth (ParametersType *Parameters, int *JA, int *IA)
{
	int i;
	int bandl, bandr;
	int n = Parameters->neq;
	int bandwidth=0;
	
	bandl = 0;
	bandr = 0;
	for (i = 0; i < n ; i++)
	{
		if (fabs(i - JA[IA[i]]) > bandl)
			bandl= i - JA[IA[i]];
		if (fabs(JA[IA[i+1]-1]-i) > bandr)	
			bandr = JA[IA[i+1]-1]-i; 
	}
	if (bandl>bandr)
		bandwidth = bandl;
	else
		bandwidth = bandr;

	return bandwidth;
}

