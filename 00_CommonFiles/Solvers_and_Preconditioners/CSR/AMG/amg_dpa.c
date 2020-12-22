#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "amg_list.h"
#include "amg_matrix.h"
#include "amg_solvers.h"
#include "amg_util.h"
#include "amg_dpa_counting.h"
#include "bcm_linmatch.h"
# include "../../../Allocation_Operations/allocations.h"

int pairs; // contabiliza quantas arestas tem o matching

/*----------------------------------------------------------------*/

int * matching_c(matrix * A, double beta);
matrix * dpa_galerkin(matrix * A, int * rmatch);
void dpa_fine2coarse(int * match, double * fine, int n, double * coarse);
void dpa_coarse2fine(int * match, double * coarse, int n, double * fine);
void dpa_coarse2fine_pure(int * match, double * coarse, int n, double * fine);
void dpa_setup(matrix *A_o, double *f_o, double *u_o, int NCL, matrix ***A, double ***f, double ***u, double **r, int ***match, double beta);
void dpa_Vcycle(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match);
void dpa_Vcycle_precond(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match);
void dpa_destroy(int NCL, matrix **A, double **f, double **u, double *r, int **match);
int dpa_AMG(matrix *A_o, double *f_o, double *u_o, int NCL, double omega, int nr, double tol, int lmax, double beta);

/*----------------------------------------------------------------*/

int * matching_c(matrix * A, double beta/*, int checkDD*/) {
    matrix * A_t = csr_transpose(A);
    A = sum_abs(A, A_t);
    destroy_m(A_t);

    int neg = 0;

    int i, j, max_j, max_j_exist, nc = 0, n = A->n;
    int * row_ptr = A->row_ptr, * col_ind = A->col_ind;
    double * val = A->val;
    Node* aux;

    int * rmatch = (int*) mycalloc("rmatch of matching_c" , (n + 1) , sizeof (int));

    counting * c = create_c(n);

    /*double * max;
    if (checkDD) {
            max = (double*)calloc(n,sizeof(double));
            for(i=0, line=0; i < nnz; i++) {
                    if (i == row_ptr[line+1]) line++;
                    if (val[i] > max[line]) max[line] = val[i];
            }
    }*/

    list ** s = (list**) mycalloc("s of matching_c",n , sizeof (list*));
    for (i = 0; i < n; i++) s[i] = create_l();


    /*TODO: checkDD*/

    /*	Encontra maior coeficiente de cada variavel [max]	*/
    double * max;
    if (neg) {
        max = (double*) mycalloc("max of matching_c",n, sizeof (double));
        for (i = 0; i < n; i++) {
            for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
                if (val[j] < 0) {
                    if (fabs(val[j]) > max[i]) {
                        /*printf("\n %f > %f",fabs(val[j]),max[i]);*/
                        max[i] = fabs(val[j]);
                    }
                    /*else
                        printf("\n %f <= %f",fabs(val[j]),max[i]);*/
                }
                /*else
                    printf("\n %f >= 0",val[j]);*/
        }
    } else {
        max = (double*) mycalloc("max of matching_c",n, sizeof (double));
        for (i = 0; i < n; i++) {
            for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
                if (fabs(val[j]) > max[i]) {
                    max[i] = fabs(val[j]);
                }
        }
    }

    /*	Insere vizinhos fortemente conectados em s	*/
    if (neg) {
        for (i = 0; i < n; i++) {
            for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
                if (i != col_ind[j])
                    if (val[j] < -beta * max[i]) {
                        //if (s[i]->n == 0) printf("\n [%d]",i);
                        insert_l_tail(s[i], 0, col_ind[j]);
                    }
        }
    } else {
        for (i = 0; i < n; i++) {
            for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
                if (i != col_ind[j])
                    if (val[j] > beta * max[i])
                        insert_l_tail(s[i], 0, col_ind[j]);
        }
    }

    /*	insere no na estrutura do counting sort	*/
    for (i = 0; i < n; i++) {
        insert_c(c, i, (s[i]->n)+(c->n));
    }
    int k = 0;
    while (c->min_m < (2 * n)) {
        k++;
        /* Passo 1 */
        i = c->count[c->min_m]->head->i;
        nc++;

        /* Passo 2 */
        max_j_exist = 0;
        max_j = row_ptr[i];
        for (j = row_ptr[i]; j < row_ptr[i + 1]; j++)
            if (i != col_ind[j])
                if (c->nodes[col_ind[j]]) {
                    if ((val[j] >= val[max_j]) || (max_j_exist==0)) {
                        max_j = j;
                    }
                    max_j_exist = 1;
                }
        j = col_ind[max_j];

        /* passo 3 */
        if (check_in_l(s[i], j) && max_j_exist) { // se j E Si
            /* Passo 4 */
            rmatch[i] = j;
            rmatch[j] = i;
            remove_c(c, i);
            remove_c(c, j);

            /* Passo 5 */
            for (aux = s[i]->head; aux; aux = aux->next) {
                if (c->nodes[(int) (aux->val)]) {
                    move_c(c, (int) (aux->val), c->m[(int) (aux->val)] - 1);
                }
            }
            for (aux = s[j]->head; aux; aux = aux->next) {
                if (c->nodes[(int) (aux->val)]) {
                    move_c(c, (int) (aux->val), c->m[(int) (aux->val)] - 1);
                }
            }
        } else {
            /* Passo 4 */
            rmatch[(int) (i)] = -1;
            remove_c(c, i);

            /* Passo 5 */
            for (aux = s[i]->head; aux; aux = aux->next) {
                if (c->nodes[(int) (aux->val)]) {
                    move_c(c, (int) (aux->val), c->m[(int) (aux->val)] - 1);
                }
            }
        }
    }
    
    myfree(max);
    for (i = 0; i < n; i++) destroy_l(s[i]);
    myfree(s);
    destroy_c(c);
    destroy_m(A);

    rmatch[n] = n - nc;

    return rmatch;
}

matrix * dpa_galerkin(matrix * A, int * rmatch) {
    /*
            O operador de galerkin consiste no produto matriz*matriz*matriz P_t*A*P
		
            O prolongador P consiste numa matriz com altura igual a da matriz A e largura igual ao numero de agregados
            onde cada linha i tem exatamente um numero 1 na coluna j indicando que a variavel i pertence ao agregado j
		
            Uma vez que a matriz P tem exatamente um numero 1 por linha, e no maximo dois numeros 1 por coluna, nao ha
            necessidade de efetuar o produto matriz*matriz completo, que eh de ordem cubica.
		
            Esta funcao usa a seguinte estrategia: somam-se as linhas que fazem parte de um mesmo agregado, e em cada
            linha, somam-se os elementos que fazem parte de um mesmo agregado. O algoritmo soma elementos das linhas
            ao mesmo tempo em que soma as linhas, logo, a matriz eh lida apenas uma vez e sequencialmente, e portanto
            a funcao eh de ordem nnz (numero de nao nulos).
		
            Como nao e possivel prever o numero de elementos nao nulos da matriz resultante, primeiro a matriz eh
            criada numa lista de adjacencias, e depois transferida para CSR.
     */

    matrix * resulting_matrix;
    Node *aux;
    /* actual_line armazena qual linha da matriz resultante esta sendo preenchida atualmente */
    int actual_line = 0, i, nnz = 0, n = A->n, resulting_n = n - rmatch[n];
    /*
            line_1_ini e line_2_ini indicam os indices onde comecam as duas linhas a serem somadas, nos vetores val e
            col_ind do CSR da matriz original.
            Do mesmo modo, line_1_end e line_2_end indicam o termino e line_1_actual e line_2_actual indicam os
            elementos sendo lidos atualmente.		
     */
    int line_1_end, line_1_actual, line_2_end, line_2_actual;
    list **row_elem = (list **) mycalloc("row_elem of dpa_galerkin",resulting_n , sizeof (list *));
    for (i = 0; i < resulting_n; i++) row_elem[i] = create_l();
    /*
            Quando se insere um elemento tal que seu par nao foi inserido ainda, o ponteiro para o novo 'node'
            inserido eh guardado no vetor pair no indice do par. Quando o outro elemento do par for lido, em vez de
            inserir um novo 'node' o seu valor sera somado ao 'node' apontado no vetor.
     */
    Node **pair = (Node **) mycalloc("pair of dpa_galerkin", n , sizeof (Node *));
    /*
            Para saber se o par ja foi inserido ou nao, eh so verificar se pair_aux na posicao do 'node' contem o numero
            da linha. Apos inserir um elemento, armazena-se o numero da linha na posicao do 'node' no pair_aux para
            ajudar mais tarde.
     */
    int *pair_aux = (int*) mycalloc("pair_aux of dpa_galerkin",n , sizeof (int));
    for (i = 0; i < n; i++) pair_aux[i] = -1;
    /*
            Pair_ind indica a qual agregado cada variavel pertence, o que eh necessario para saber o col_ind ao gerar
            o CSR da matriz resultante. 
     */
    int *pair_ind = (int*) mycalloc("pair_ind of dpa_galerkin",n, sizeof (int));
    /*	atalhos		*/
    double * val = A->val;
    double * diag = A->diag;
    int * col_ind = A->col_ind;
    int * row_ptr = A->row_ptr;

    for (i = 0; i < n; i++) {
        /*
                Percorre as linhas da matriz original. Se a linha nao tiver par, eh transferida diretamente para
                a matriz resultante, pois nao ha nada a se somar a ela. Se ela tiver par, somam-se as duas linhas

                Como estamos lendo linha por linha, iremos passar pelas duas linhas de cada par, ou seja, cada par
                sera visitado duas vezes. Para processar cada par apenas uma vez, eh checado se o indice da linha
                eh menor que o indice do par. Se nao for, ignora. Assim cada par eh processado exatamente uma vez.
         */
        if (rmatch[i] == -1) {
            /*	rmatch[i] == -1 significa que i nao tem par	*/

            /*	calcula-se o inicio e o fim da linha a ser transferida		*/
            line_1_actual = row_ptr[i];
            line_1_end = row_ptr[i + 1];

            while (line_1_actual < line_1_end) {
                /*	line_1_actual < line_1_end significa que a linha ainda nao acabou	*/
                if (pair_aux[col_ind[line_1_actual]] == i) {
                    /*	se o par ja foi inserido, soma o valor e limpa o vetor	*/
                    pair[col_ind[line_1_actual]]->val += val[line_1_actual];
                } else {
                    /*	se o par nao foi inserido ainda, insere novo 'node' e contabiliza	*/
                    if (rmatch[col_ind[line_1_actual]] == -1)
                        insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                    else if (col_ind[line_1_actual] > rmatch[col_ind[line_1_actual]])
                        insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_1_actual]], val[line_1_actual]);
                    else {
                        insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                        /*	se houver par, grava o ponteiro para novo 'node' para posterior soma dos valores	*/
                        pair[rmatch[col_ind[line_1_actual]]] = row_elem[actual_line]->tail;
                        pair_aux[rmatch[col_ind[line_1_actual]]] = i;
                    }
                    nnz++;
                }
                /*	avanca na leitura da linha	*/
                line_1_actual++;
            }
            /*	o indice do elemento na nova matriz eh armazenada em pair_ind, para obter o col_ind mais tarde	*/
            pair_ind[i] = actual_line;
            actual_line++;
        } else
            if (i < rmatch[i]) {
            /*	i < rmatch[i] significa que o indice da linha i eh menor que o indice do par de i	*/

            /*	calcula-se inicio e fim das duas linhas a serem somadas		*/
            line_1_actual = row_ptr[i];
            line_1_end = row_ptr[i + 1];
            line_2_actual = row_ptr[rmatch[i]];
            line_2_end = row_ptr[rmatch[i] + 1];

            /*	enquanto houver coisa para transferir	*/
            while ((line_1_actual < line_1_end) || (line_2_actual < line_2_end)) {
                /*	Se a linha 2 acabou, certamente a linha 1 nao acabou, pois senao o while teria dado false.
                        Logo,  transfere a linha 1 sem somar nada pois nao ha nada a somar	*/
                if (line_2_actual == line_2_end) {
                    if (pair_aux[col_ind[line_1_actual]] == i) {
                        pair[col_ind[line_1_actual]]->val += val[line_1_actual];
                    } else {
                        if (rmatch[col_ind[line_1_actual]] == -1)
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                        else if (col_ind[line_1_actual] > rmatch[col_ind[line_1_actual]])
                            insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_1_actual]], val[line_1_actual]);
                        else {
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                            pair[rmatch[col_ind[line_1_actual]]] = row_elem[actual_line]->tail;
                            pair_aux[rmatch[col_ind[line_1_actual]]] = i;
                        }
                        nnz++;
                    }
                    line_1_actual++;
                } else
                    /*	Se a linha 2 nao acabou e a linha 1 sim, transfere a linha 2 sem somar nada	*/
                    /*	Se as duas linhas nao tiverem acabado, mas o indice dos elementos lidos forem diferentes,
                            transfere o de indice menor, a fim de manter a ordem	*/
                    if ((line_1_actual == line_1_end) || (col_ind[line_1_actual] > col_ind[line_2_actual])) {
                    if (pair_aux[col_ind[line_2_actual]] == i) {
                        pair[col_ind[line_2_actual]]->val += val[line_2_actual];
                    } else {
                        if (rmatch[col_ind[line_2_actual]] == -1)
                            insert_l_tail(row_elem[actual_line], col_ind[line_2_actual], val[line_2_actual]);
                        else if (col_ind[line_2_actual] > rmatch[col_ind[line_2_actual]])
                            insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_2_actual]], val[line_2_actual]);
                        else {
                            insert_l_tail(row_elem[actual_line], col_ind[line_2_actual], val[line_2_actual]);
                            pair[rmatch[col_ind[line_2_actual]]] = row_elem[actual_line]->tail;
                            pair_aux[rmatch[col_ind[line_2_actual]]] = i;
                        }
                        nnz++;
                    }
                    line_2_actual++;
                } else
                    if (col_ind[line_1_actual] < col_ind[line_2_actual]) {
                    if (pair_aux[col_ind[line_1_actual]] == i) {
                        pair[col_ind[line_1_actual]]->val += val[line_1_actual];
                    } else {
                        if (rmatch[col_ind[line_1_actual]] == -1)
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                        else if (col_ind[line_1_actual] > rmatch[col_ind[line_1_actual]])
                            insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_1_actual]], val[line_1_actual]);
                        else {
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual]);
                            pair[rmatch[col_ind[line_1_actual]]] = row_elem[actual_line]->tail;
                            pair_aux[rmatch[col_ind[line_1_actual]]] = i;
                        }
                        nnz++;
                    }
                    line_1_actual++;
                } else
                    /*	O unico caso que resta eh o que as duas linhas nao acabaram e os dois elementos lidos tem
                            o mesmo indice. Neste caso, transfere a soma. Assim, teremos somado as duas linhas	*/
                {
                    if (pair_aux[col_ind[line_1_actual]] == i) {
                        pair[col_ind[line_1_actual]]->val += val[line_1_actual] + val[line_2_actual];
                    } else {
                        if (rmatch[col_ind[line_1_actual]] == -1)
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual] + val[line_2_actual]);
                        else if (col_ind[line_1_actual] > rmatch[col_ind[line_1_actual]])
                            insert_l_sort(row_elem[actual_line], rmatch[col_ind[line_1_actual]], val[line_1_actual] + val[line_2_actual]);
                        else {
                            insert_l_tail(row_elem[actual_line], col_ind[line_1_actual], val[line_1_actual] + val[line_2_actual]);
                            pair[rmatch[col_ind[line_1_actual]]] = row_elem[actual_line]->tail;
                            pair_aux[rmatch[col_ind[line_1_actual]]] = i;
                        }
                        nnz++;
                    }
                    line_1_actual++;
                    line_2_actual++;
                }
            }
            /*	o indice do elemento na nova matriz eh armazenada em pair_ind, para obter o col_ind mais tarde	*/
            pair_ind[i] = pair_ind[rmatch[i]] = actual_line;
            actual_line++;
        }
    }

    resulting_matrix = create_m(resulting_n, resulting_n, nnz);
    val = resulting_matrix->val;
    diag = resulting_matrix->diag;
    col_ind = resulting_matrix->col_ind;
    row_ptr = resulting_matrix->row_ptr;
    int count = 0;
    row_ptr[0] = 0;
    for (i = 0; i < resulting_n; i++) {
        for (aux = row_elem[i]->head; aux; aux = aux->next, count++) {
            val[count] = aux->val;
            col_ind[count] = pair_ind[aux->elem];
            if (pair_ind[aux->elem] == i) diag[i] = aux->val;
            //if(i==1838)printf("\n[%d][%d/%d->%d]=%f",i,aux->elem,rmatch[aux->elem],pair_ind[aux->elem],aux->val);
            //printf("\n[%d][%d]=%f",i,pair_ind[aux->elem],aux->val); // descomente esta linha para printar a matriz
        }
        row_ptr[i + 1] = count;
    }

    for (i = 0; i < resulting_n; i++) destroy_l(row_elem[i]);
    myfree(row_elem);
    myfree(pair);
    myfree(pair_aux);
    myfree(pair_ind);

    return resulting_matrix;
}

void dpa_fine2coarse(int * match, double * fine, int n, double * coarse) {

    for (int i = 0, j = 0; i < n; i++)
        if (match[i] == -1)
            coarse[j++] = fine[i];
        else if (i < match[i])
            coarse[j++] = fine[i] + fine[match[i]];
}

void dpa_coarse2fine(int * match, double * coarse, int n, double * fine) {
    for (int i = 0, j = 0; i < n; i++)
        if (match[i] == -1)
            fine[i] += coarse[j++];
        else if (i < match[i]) {
            fine[i] += coarse[j];
            fine[match[i]] += coarse[j++];
        }
}

// doesn't sum arrays
void dpa_coarse2fine_pure(int * match, double * coarse, int n, double * fine) {
    for (int i = 0, j = 0; i < n; i++)
        if (match[i] == -1)
            fine[i] = coarse[j++];
        else if (i < match[i]) {
            fine[i] = coarse[j];
            fine[match[i]] = coarse[j++];
        }
}

void dpa_setup(matrix *A_o, double *f_o, double *u_o, int NCL, matrix ***A, double ***f, double ***u, double **r, int ***match, double beta) {
    int i;
    NCL *= 2; // DPA
    *A = (matrix **) mycalloc("A of dpa setup",(NCL + 1) , sizeof (matrix *));
    *f = (double **) mycalloc("f of dpa setup",(NCL + 1) , sizeof (double *));
    *u = (double **) mycalloc("u of dpa setup",(NCL + 1) , sizeof (double *));
    *match = (int **) mycalloc("match of dpa setup",NCL , sizeof (int *));
    (*A)[0] = A_o;
    (*f)[0] = f_o;
    (*u)[0] = u_o;
    matrix *A_s, *A_t;
    
    /*struct mc64_control mc_cntrl;
    struct mc64_info    mc_info;*/  
    
    // SETUP PHASE
    
    // gera A[i+1] com base em A[i]
    if (beta == -1) 
        for (i = 0; i < NCL; i++) {
            A_t = csr_transpose((*A)[i]);
            A_s = sum_upper((*A)[i], A_t);
            (*match)[i] = bcm_CSRMatrixHMatch(A_s);
            (*A)[i + 1] = dpa_galerkin((*A)[i], (*match)[i]);
            destroy_m(A_t);
            destroy_m(A_s);
        }
    else
        for (i = 0; i < NCL; i++) {
            (*match)[i] = matching_c((*A)[i], beta);
            (*A)[i + 1] = dpa_galerkin((*A)[i], (*match)[i]);
        }

    for (i = 1; i < (NCL + 1); i++) {
        (*f)[i] = (double *) mycalloc("f[i] of dpa setup",(*A)[i]->m , sizeof (double));
        (*u)[i] = (double *) mycalloc("u[i] of dpa setup",(*A)[i]->n , sizeof (double));
    }
    *r = (double *) mycalloc("r of dpa setup",A_o->m , sizeof (double));
};

// na descida, primeiro faz correcao por erro, depois NAO faz
void dpa_Vcycle(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match) {

    int i, j;
    NCL*=2;
    for (i = 0; i < NCL; i+=2) {
        SOR_relax(A[i], f[i], u[i], omega, nr);
        residual(r, f[i], A[i], u[i], 0, 0);
        dpa_fine2coarse(match[i],   r,      A[i]->n  , f[i+1]); // f[i+1] <- r
        dpa_fine2coarse(match[i+1], f[i+1], A[i+1]->n, f[i+2]); // f[i+2] <- f[i+1]
        for (j = 0; j < A[i+1]->n; j++) u[i+1][j] = 0.0;
        for (j = 0; j < A[i+2]->n; j++) u[i+2][j] = 0.0;
    }
    SOR_relax(A[i], f[i], u[i], omega, 2 * nr);
    for (; i > 0; i-=2) {
        dpa_coarse2fine_pure(match[i-1], u[i],   A[i-1]->n, u[i-1]); // u[i-1] <- u[i]
        dpa_coarse2fine     (match[i-2], u[i-1], A[i-2]->n, u[i-2]); // u[i-2] <- u[i-1]
        SOR_relax(A[i-2], f[i-2], u[i-2], omega, nr);
    }
};

void dpa_Vcycle_precond(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match) {

    int i, j;
    NCL*=2;
    
    for (i=0; i<A[0]->m; i++) u[0][i] = 0.0;
    for (i = 0; i < NCL; i+=2) {
        SOR_relax(A[i], f[i], u[i], omega, nr);
        residual(r, f[i], A[i], u[i], 0, 0);
        dpa_fine2coarse(match[i],   r,      A[i]->n  , f[i+1]); // f[i+1] <- r
        dpa_fine2coarse(match[i+1], f[i+1], A[i+1]->n, f[i+2]); // f[i+2] <- f[i+1]
        for (j = 0; j < A[i+1]->n; j++) u[i+1][j] = 0.0;
        for (j = 0; j < A[i+2]->n; j++) u[i+2][j] = 0.0;
    }
    SOR_relax(A[i], f[i], u[i], omega, nr);
    SOR_relax_rev(A[i], f[i], u[i], omega, nr);
    for (; i > 0; i-=2) {
        dpa_coarse2fine_pure(match[i-1], u[i],   A[i-1]->n, u[i-1]); // u[i-1] <- u[i]
        dpa_coarse2fine     (match[i-2], u[i-1], A[i-2]->n, u[i-2]); // u[i-2] <- u[i-1]
        SOR_relax_rev(A[i-2], f[i-2], u[i-2], omega, nr);
    }
};

void xdpa_precond(int NCL, matrix **A, double **f, double **u, double *r, double omega, int nr, int **match) {

    int i, j;
    NCL*=2;
    
    for (i=0; i<A[0]->m; i++) u[0][i] = 0.0;
    for (i = 0; i < NCL; i+=2) {
        SOR_relax(A[i], f[i], u[i], omega, nr);
        residual(r, f[i], A[i], u[i], 0, 0);
        dpa_fine2coarse(match[i],   r,      A[i]->n  , f[i+1]); // f[i+1] <- r
        dpa_fine2coarse(match[i+1], f[i+1], A[i+1]->n, f[i+2]); // f[i+2] <- f[i+1]
        for (j = 0; j < A[i+1]->n; j++) u[i+1][j] = 0.0;
        for (j = 0; j < A[i+2]->n; j++) u[i+2][j] = 0.0;
    }
    GMRES(A[i], f[i], u[i], 100000, 1e-8, 1);
    for (; i > 0; i-=2) {
        dpa_coarse2fine_pure(match[i-1], u[i],   A[i-1]->n, u[i-1]); // u[i-1] <- u[i]
        dpa_coarse2fine     (match[i-2], u[i-1], A[i-2]->n, u[i-2]); // u[i-2] <- u[i-1]
	GMRES(A[i-2], f[i-2], u[i-2], 10, 1e-8, 1);
    }
};

void dpa_destroy(int NCL, matrix **A, double **f, double **u, double *r, int **match) {
    int i;
    NCL*=2;
    for (i = 1; i < (NCL + 1); i++) {
        myfree(f[i]);
        myfree(u[i]);
    }
    myfree(f);
    myfree(u);
    myfree(r);
    for (i = 1; i < (NCL + 1); i++) destroy_m(A[i]);
    // matriz mais fina tem dados de outras structs
    myfree(A[0]->diag);
    myfree(A[0]);
    for (i = 0; i < NCL; i++) myfree(match[i]);
    myfree(A);
    myfree(match);
};

int dpa_AMG(matrix *A_o, double *f_o, double *u_o, int NCL, double omega, int nr, double tol, int lmax, double beta) {
    int k;
    double **f, **u, *r, delta, norm_f;
    matrix **A/*, **I_cf, **I_cf_t*/;
    int **match;

    // AMG_setup generates A matrices, u and r arrays, and prolongator matrices
    dpa_setup(A_o, f_o, u_o, NCL, &A, &f, &u, &r, &match, beta);
    residual(r, f_o, A_o, u_o, 0, 0);
    norm_f = norm_inf(f_o, A_o->m);
    delta = norm_inf(r, A_o->m) / norm_f;
    k = 0;
    while (k < lmax && delta > tol) {
        // multiply matrices by u for intergrid transfer
        dpa_Vcycle(NCL, A, f, u, r, omega, nr, match);
        residual(r, f_o, A_o, u_o, 0, 0);
        delta = norm_inf(r, A_o->m) / norm_f;
        k++;
    }
    dpa_destroy(NCL, A, f, u, r, match);
    return (delta <= tol) ? 1 : 0;
};
