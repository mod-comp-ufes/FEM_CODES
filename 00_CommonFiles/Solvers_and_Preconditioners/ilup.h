/*----------------------------------------------------------------------------
 * MATRIX HEADER FILE
 *--------------------------------------------------------------------------*/
#ifndef PRECONDITIONERS_H
#define PRECONDITIONERS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define min(a,b) (((a)>(b))?(b):(a))
#define max(a,b) (((a)>(b))?(a):(b))

typedef struct
{
	int n, nz;
	double *D;
	double *AA;
	int *JA;
	int *IA;	
} MAT;

typedef struct
{
	int           n;
	int*    nzcount;  /* length of each row                          */
	int**        ja;  /* pointer-to-pointer to store column indices  */
	double**     ma;  /* pointer-to-pointer to store nonzero entries */
} SparMAT;

typedef struct
{
	int         n;
	SparMAT*    L;   /* L part elements   */
	double*     D;   /* diagonal elements */
	SparMAT*    U;   /* U part elements   */
	int*     work;   /* working buffer    */
} SparILU;

void SPARMAT_setup  (SparMAT* mat, int n);
void SPARILU_setup  (SparILU* lu, int n);
void SPARIUL_setup  (SparILU* ul, int n);
void SPARILU_row    (SparILU* lu, int nrow);
void SPARIUL_row    (SparILU* ul, int nrow);
void CSRto_SPARMAT (MAT* A, SparMAT* mat);
void CSRto_SPARMAT_setup (MAT* A, SparMAT* mat);
void SPARILU_toCSR (SparILU* lu, MAT* L, MAT* U);
int  LEVEL_OF_FILL_lu  (SparMAT* csmat, SparILU* lu, int p);
int  LEVEL_OF_FILL_ul  (SparMAT* csmat, SparILU* ul, int p);
void ILUP           (SparMAT* csmat, SparILU* lu, int p);
void IULP           (SparMAT* csmat, SparILU* lu, int p);
void SPARMAT_clean  (SparMAT* mat);
void SPARILU_clean  (SparILU* lu);
void SPARILU_print  (SparILU* lu);
void QSPLIT         (double *a, int *ind, int n, int Ncut);
void ILUT           (SparMAT* csmat, SparILU* lu, int lfil, double tol);

#endif /* PRECONDITIONERS_H */
