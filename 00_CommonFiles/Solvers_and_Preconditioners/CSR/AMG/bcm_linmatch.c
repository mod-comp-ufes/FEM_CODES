/*
                BootCMatch
     Bootstrap AMG based on Compatible weighted Matching, version 0.9
    (C) Copyright 2017
                       Pasqua D'Ambra         IAC-CNR, IT
                       Panayot S. Vassilevski Portland State University, OR USA

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions, and the following disclaimer in the
       documentation and/or other materials provided with the distribution.
    3. The name of the BootCMatch group or the names of its contributors may
       not be used to endorse or promote products derived from this
       software without specific written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BootCMatch GROUP OR ITS CONTRIBUTORS
  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.

 */

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "amg_matrix.h"
# include "../../../Allocation_Operations/allocations.h"

/*----------------------------------------------------------------------------
 * bcm_CSRMatrixHMatch compute an half approximate maximum product matching 
 * in the graph of a sparse matrix by applying the algorithm in:
 * Preis, Linear Time 1/2-approximation algorithm for maximum weighted matching
 * in General Graph, in STACS'99, LNCS, vol.1563 (1999).
 *
 * Note: we expect a square matrix with symmetric sparsity pattern and null diagonal;
 *       Indeed we work with a strictly upper triangular matrix
 *--------------------------------------------------------------------------*/

int pairs; // contabiliza quantas arestas tem o matching


int bcm_trymatch(int rowindex, int colindex, matrix *W, int *jrowindex,
        int ljrowindex, int *jcolindex, int ljcolindex, int *rmatch);

int * bcm_CSRMatrixHMatch(matrix *A);

/*----------------------------------------------------------------*/

int bcm_trymatch(int rowindex, int colindex, matrix *W, int *jrowindex,
        int ljrowindex, int *jcolindex, int ljcolindex, int *rmatch) {

    int tryrowmatch, trycolmatch;
    int i, j, k, nzrow_W, startj, kindex;
    double cweight, nweight;

    int *W_i = W->row_ptr;
    int *W_j = W->col_ind;
    int nrows_W = W->n;
    double *W_data = W->val;

    k = -1;
    i = 0;
    while (i < ljrowindex && k == -1) {
        if (jrowindex[i] == colindex) k = i;
        i++;
    }
    if (k >= 0) {
        ljrowindex = ljrowindex - 1;
        for (i = k; i < ljrowindex; ++i) jrowindex[i] = jrowindex[i + 1];
    }

    k = -1;
    i = 0;
    while (i < ljcolindex && k == -1) {
        if (jcolindex[i] == rowindex) k = i;
        i++;
    }
    if (k >= 0) {
        ljcolindex = ljcolindex - 1;
        for (i = k; i < ljcolindex; ++i) jcolindex[i] = jcolindex[i + 1];
    }

    nzrow_W = W_i[rowindex + 1] - W_i[rowindex];
    startj = W_i[rowindex];
    kindex = -1;
    k = 0;
    while (k < nzrow_W && kindex == -1) {
        if (W_j[startj + k] == colindex) kindex = startj + k;
        k++;
    }
    cweight = W_data[kindex];

    while ((rmatch[rowindex] == -1 && rmatch[colindex] == -1)
            && (ljrowindex != 0 || ljcolindex != 0)) {

        if (rmatch[rowindex] == -1 && ljrowindex != 0) {
            tryrowmatch = jrowindex[0];
            nzrow_W = W_i[rowindex + 1] - W_i[rowindex];
            startj = W_i[rowindex];
            k = 0;
            kindex = -1;
            while (k < nzrow_W && kindex == -1) {
                if (W_j[startj + k] == tryrowmatch) kindex = startj + k;
                k++;
            }
            nweight = W_data[kindex];

            ljrowindex = ljrowindex - 1;
            for (i = 0; i < ljrowindex; ++i) jrowindex[i] = jrowindex[i + 1];

            if (nweight > cweight && rmatch[tryrowmatch] == -1) {

                nzrow_W = W_i[tryrowmatch + 1] - W_i[tryrowmatch];
                int *trymatchindexrow;

                trymatchindexrow = (int *) mycalloc("trymatchindexrow of bcm_trymatch",nzrow_W, sizeof (int));

                startj = W_i[tryrowmatch];
                for (k = 0; k < nzrow_W; ++k) trymatchindexrow[k] = W_j[startj + k];

                bcm_trymatch(rowindex, tryrowmatch, W, jrowindex, ljrowindex,
                        trymatchindexrow, nzrow_W, rmatch);
                myfree(trymatchindexrow);
            }
        }

        if (rmatch[colindex] == -1 && ljcolindex != 0) {
            trycolmatch = jcolindex[0];
            nzrow_W = W_i[colindex + 1] - W_i[colindex];
            startj = W_i[colindex];
            k = 0;
            kindex = -1;
            while (k < nzrow_W && kindex == -1) {
                if (W_j[startj + k] == trycolmatch) kindex = startj + k;
                k++;
            }
            nweight = W_data[kindex];

            ljcolindex = ljcolindex - 1;
            for (i = 0; i < ljcolindex; ++i) jcolindex[i] = jcolindex[i + 1];

            if (nweight > cweight && rmatch[trycolmatch] == -1) {
                nzrow_W = W_i[trycolmatch + 1] - W_i[trycolmatch];
                int *trymatchindexcol;
                trymatchindexcol = (int *) mycalloc("trymatchindexcol of bcm_trymatch",nzrow_W, sizeof (int));

                startj = W_i[trycolmatch];
                for (k = 0; k < nzrow_W; ++k) trymatchindexcol[k] = W_j[startj + k];

                bcm_trymatch(colindex, trycolmatch, W, jcolindex, ljcolindex,
                        trymatchindexcol, nzrow_W, rmatch);
                myfree(trymatchindexcol);
            }
        }
    }

    if (rmatch[rowindex] == -1 & rmatch[colindex] == -1) {
        rmatch[rowindex] = colindex;
        rmatch[colindex] = rowindex;
        pairs++;
    }

    return 0;
}

int * bcm_CSRMatrixHMatch(matrix *A) {
    matrix *B = A;

    int i, j, k, *rmatch;
    double *c, alpha;
    int jbp, nzrows_B;
    double tmp = 0.0;
    int rno, cno, nzrows_cno, startj, ljrowindex, ljjrowindex;
    int *jrowindex, *jjrowindex;

    int *B_i = B->row_ptr;
    int *B_j = B->col_ind;
    double *B_data = B->val;
    int nrows_B = B->n;
    int ncols_B = B->m;
    int nnz_B = B->nnz;

    nzrows_B = B_i[1] - B_i[0];

    matrix * W = copy_m(B);

    int *W_i = W->row_ptr;
    int *W_j = W->col_ind;
    int nrows_W = W->n;
    double *W_data = W->val;

    /* needed for computing weights as in Bora code */

    double min_B = fabs(B_data[0]);
    for (i = 1; i < nnz_B; i++) {
        if (min_B > fabs(B_data[i])) min_B = fabs(B_data[i]);
    }

    /* end of weights as in Bora code */

    for (i = 0; i < nrows_W; ++i) {
        for (j = W_i[i]; j < W_i[i + 1]; ++j) {
            /* As in Duff paper */
            //if(fabs(B_data[j])>0.) W_data[j]=-log(c[B_j[j]])+log(fabs(B_data[j]));
            // else W_data[j]=FLT_MAX; 
            /* As in Scott paper for maximum cardinality */
            // if(fabs(B_data[j])>0.) W_data[j]=alpha+log(fabs(B_data[j]))+(alpha-c[B_j[j]]);
            //   else W_data[j]=FLT_MAX; 
            /* As Bora code */
            W_data[j] = log(fabs(B_data[j]) / (0.999 * min_B)); /* This seems to be generally better */
        }
    }

    rmatch = (int *) mycalloc("rmatch of bcm_CSRMatrixHMatch",(nrows_B + 1), sizeof (int));

    for (i = 0; i < nrows_B; ++i) rmatch[i] = -1;

    jrowindex = (int *) mycalloc("jrowindex of bcm_CSRMatrixHMatch",nrows_B, sizeof (int));
    jjrowindex = (int *) mycalloc("jjrowindex of bcm_CSRMatrixHMatch",nrows_B, sizeof (int));

    pairs = 0;

    jbp = 0;
    for (i = 0; i < nrows_B; ++i) {
        nzrows_B = B_i[i + 1] - B_i[i];

        for (j = 0; j < nzrows_B; ++j) jrowindex[j] = B_j[jbp + j];
        for (j = 0; j < nzrows_B; ++j) {
            rno = i;
            cno = B_j[jbp + j];

            startj = B_i[cno];
            nzrows_cno = B_i[cno + 1] - startj;

            for (k = 0; k < nzrows_cno; ++k) jjrowindex[k] = B_j[startj + k];

            if (rmatch[rno] == -1 && rmatch[cno] == -1)
                bcm_trymatch(rno, cno, W, jrowindex, nzrows_B, jjrowindex, nzrows_cno, rmatch);

        }

        jbp = jbp + nzrows_B;
    }

    rmatch[nrows_B] = pairs;

    destroy_m(W);
    myfree(jrowindex);
    myfree(jjrowindex);
    /* myfree(c); */

    return rmatch;
}
