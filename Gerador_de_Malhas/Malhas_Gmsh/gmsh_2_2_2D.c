#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gmsh2Mesh.h"

extern int setmark_according_boundary(NodeType *node, int *mark, int nNode);

int print_bound[18][4] = {{0,0,0,0}/*mark 0 */, {0,0,0,1}/*mark 1 */, {0,0,1,0}/*mark 2 */, {0,0,1,1}/*mark 3 */,
				         {0,1,0,0}/*mark 4 */, {0,1,0,1}/*mark 5 */, {0,1,1,0}/*mark 6 */, {0,1,1,1}, /*mark 7 */
						 {1,0,0,0}/*mark 8 */, {1,0,0,1}/*mark 9 */, {1,0,1,0}/*mark 10 */, {1,0,1,1}/*mark 11 */,
				         {1,1,0,0}/*mark 12 */, {1,1,0,1}/*mark 13 */, {1,1,1,0}/*mark 14 */, {1,1,1,1}, /*mark 15 */
						 // special cases
						 // cylinder
						 {-1,-1,1,0}/*mark 16 */, {-2,-2,1,0}/*mark 17 */};


int main(int argc, char **argv)
{
	int nel, nnodes, I, index, aux1, aux2, aux3, aux4, aux5;
	double z;	
	NodeType *Node;
	ElementType *Element;
	char label[300], filename[300];
	char *aux;
	const char delimiter[3]=" \n";
	FILE *In, *Out;

	if (argc<3){
		printf("Use ./gmsh2Mesh <msh file> <problem>\n");
		//Exemplo ./gmsh2Mesh cavity.msh CAVITY
		exit(1);
	}

	if ((In=fopen(argv[1],"r"))==NULL){
		printf("File %s not found!\n", argv[1]);
		exit(1);
	}		
	
	//Reading header information
	do{
		fscanf(In,"%s",label);
	}while(strcmp(label, "$Nodes")!=0);

	fscanf(In,"%d",&nnodes);
	
	if ((Node = calloc(nnodes,sizeof(NodeType)))==NULL){
		printf("Memmory allocation error in Node!\n");
		exit(1);
	}

	//Reading nodes
	for (I=0; I<nnodes; I++)
		fscanf(In,"%d%lf%lf%lf", &index, &Node[I].x, &Node[I].y, &z);
	
	//Reading header information
	do{
		fscanf(In,"%s",label);
	}while(strcmp(label, "$Elements")!=0);

	fscanf(In,"%d\n",&nel);

	do{
		fgets(label, 300, In);
		aux = strtok(label, delimiter);
		aux = strtok(NULL, delimiter);
		nel--;
	}while (strcmp(aux, "2")!=0);

	nel++;
	
	if ((Element = calloc(nel,sizeof(ElementType)))==NULL){
		printf("Memmory allocation error in Element!\n");
		exit(1);
	}

	aux = strtok(NULL, delimiter);
	aux = strtok(NULL, delimiter);
	aux = strtok(NULL, delimiter);
	aux = strtok(NULL, delimiter);
	Element[0].Vertex[0] = atoi(aux);
	aux = strtok(NULL, delimiter);
	Element[0].Vertex[1] = atoi(aux);
	aux = strtok(NULL, delimiter);
	Element[0].Vertex[2] = atoi(aux);
	
	for(I=1; I<nel;I++)
		fscanf(In, "%d%d%d%d%d%d%d%d", &aux1, &aux2, &aux3, &aux4, &aux5, &Element[I].Vertex[0], &Element[I].Vertex[1], &Element[I].Vertex[2]);  	
	
	fclose(In);

    sprintf(filename,"%s_%d_%d.dat",argv[2],nnodes,nel);
    
	if((Out = fopen(filename, "w")) == NULL) {
		printf("File %s not found!\n", filename);
		exit(1);
	}

	// Nodes
	int *mark = (int*) malloc(nnodes*sizeof(int));
	int NDOF = setmark_according_boundary(Node, mark, nnodes);
	fprintf(Out, "%d\n", nnodes);

	for(I=0;I<nnodes;I++) {
		fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
		for(int j=0;j<NDOF;j++)
			fprintf(Out, "\t%d", print_bound[mark[I]][j]);
		fprintf(Out,"\n");
	}

	// Elements
	fprintf(Out,"%d\n",nel);
		for(I=0;I<nel;I++)
			fprintf(Out, "%d\t%d\t%d\t-1\n", Element[I].Vertex[0]-1, Element[I].Vertex[1]-1, Element[I].Vertex[2]-1);

	fclose(Out);

	return 0;
}