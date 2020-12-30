#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include "gmsh2Mesh.h"

extern int setmark_according_boundary(NodeType *node, int *mark, int nNode);

int print_bound[17][4] = {{0,0,0,0}/*mark 0 */, {0,0,0,1}/*mark 1 */, {0,0,1,0}/*mark 2 */, {0,0,1,1}/*mark 3 */,
				         {0,1,0,0}/*mark 4 */, {0,1,0,1}/*mark 5 */, {0,1,1,0}/*mark 6 */, {0,1,1,1}, /*mark 7 */
						 {1,0,0,0}/*mark 8 */, {1,0,0,1}/*mark 9 */, {1,0,1,0}/*mark 10 */, {1,0,1,1}/*mark 11 */,
				         {1,1,0,0}/*mark 12 */, {1,1,0,1}/*mark 13 */, {1,1,1,0}/*mark 14 */, {1,1,1,1} /*mark 15 */,
						 {2,1,1,1}/*mark 16 */ };

char* read_line(FILE *arq) {
    size_t n = 0, r;
    char *buf = NULL;

    r = getline(&buf, &n, arq);
    if(r < 0)
        return(NULL);

    if(buf && r >= UINT_MAX) {
        free(buf);
        return(NULL);
    }
    return(buf);
}


int main(int argc, char *argv[]) {
	size_t i, j; // counter
	char *line, // arbtrary line to read
		 filename_out[255];
	int tmp[4], m, n, // readers from file in
		nNode, nElement, 
		NDOF,
		*mark;
	NodeType *node;
	ElementType *element;
	FILE *in, *out;

	if (argc<3){
		printf("Use ./gmsh2Mesh <msh file> <problem>\n");
		//Exemplo ./gmsh2Mesh cavity.msh CAVITY
		exit(1);
	}

	if((in = fopen(argv[1], "r")) == NULL) {
		printf("File %s not found!\n", argv[1]);
		exit(1);
	}
	
    while(!feof(in)) {
        line = read_line(in);
        if(!line) continue; // inválida?

		// Nós
        if(strncmp(line, "$Nodes", 6) == 0) {
			free(line);
			// Leitura de algo + quantidade de nós
			fscanf(in, "%d %d %d %d\n", &tmp[0], &nNode, &tmp[1], &tmp[2]);
			node = (NodeType*) malloc(nNode*sizeof(NodeType));
			for(i = 0; i < nNode; ) {
				// Leitura das informações dos blocos
				fscanf(in, "%d %d %d %d\n", &tmp[0], &tmp[1], &tmp[2], &m);
				// Leitura do id e das coordenadas dos nós
				for(j = 0; j < m; j++)
					fscanf(in, "%d\n", &node[i+j].id);
				// Leitura das coordenadas dos nós
				for(j = 0; j < m; j++, i++)
					fscanf(in, "%lf %lf %lf\n", &node[i].x, &node[i].y, &node[i].z);
			}
			// Chegou ao final do bloco
			line = read_line(in);
			if(strncmp(line, "$EndNodes", 9) != 0) { // Alguma coisa deu errado
				perror("Problema durante leitura dos Nós!");
				exit(1);
			}
			free(line);
			continue;
        }

		// Elementos
        else if(strncmp(line, "$Elements", 9) == 0) {
			free(line);
			// Leitura de algo + quantidade de nós
			fscanf(in, "%d %d %d %d\n", &tmp[0], &n, &tmp[1], &tmp[2]);
			for(i = 0; i < n; ) {
				// Leitura das informações dos blocos
				fscanf(in, "%d %d %d %d\n", &tmp[0], &tmp[1], &tmp[2], &m);
				// Ignora os pontos
				if(tmp[2] == 15)
					for(j = 0; j < m; j++, i++)
						fscanf(in, "%d %d\n", &tmp[0], &tmp[1]);
				
				// Ignora as retas
				else if(tmp[2] == 1)
					for(j = 0; j < m; j++, i++)
						fscanf(in, "%d %d %d\n", &tmp[0], &tmp[1], &tmp[2]);
				
				// Elementos triangulares (p.199, doc)
				else if(tmp[2] == 2) {
					nElement = m;
					element = (ElementType*) malloc(nElement*sizeof(ElementType));
					for(j = 0; j < nElement; j++, i++) {
						fscanf(in, "%d %d %d %d\n", &element[j].id, &element[j].Vertex[0], &element[j].Vertex[1], &element[j].Vertex[2]);
						element[j].Type = -1; // por hora, nenhum caso possui diferenciação
					}
				}
				// Alguma coisa deu errado
				else
					break;
			}
			// Chegou ao final do bloco
			line = read_line(in);
			if(strncmp(line, "$EndElements", 12) != 0) { // Alguma coisa deu errado
				perror("Problema durante leitura dos Elementos!");
				exit(1);
			}
			free(line);
			continue;
        }
		else {
			free(line);
			continue;
		}
    }
	fclose(in);

	// Escrita do arquivo .dat
	sprintf(filename_out,"%s_%d_%d.dat", argv[2], nNode, nElement);
	if((out = fopen(filename_out, "w")) == NULL) {
		printf("File %s not found!\n", filename_out);
		exit(1);
	}

	// Nós
	mark = (int*) malloc(nNode*sizeof(int));
	NDOF = setmark_according_boundary(node, mark, nNode);
	fprintf(out, "%d\n", nNode);
	for(i = 0; i < nNode; i++) {
		fprintf(out, "%.14lf\t%.14lf", node[i].x, node[i].y);
		for(j = 0; j < NDOF; j++)
			fprintf(out, "\t%d", print_bound[mark[i]][j]);
		fprintf(out,"\n");
	}

	// Elementos
	fprintf(out, "%d\n", nElement);
	for(i = 0; i < nElement; i++)
		fprintf(out, "%d\t%d\t%d\t%d\n", element[i].Vertex[0]-1, element[i].Vertex[1]-1, element[i].Vertex[2]-1, element[i].Type);

	fclose(out);

	// Desalocações
	free(node);
	free(element);
	free(mark);
	
	return 0;
}
