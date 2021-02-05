#include "../../gmsh2Mesh.h"
#include <math.h>
	

int setmark_according_boundary(NodeType *Node, int *mark, int nNode) {
	int I;
	int NDOF;

	NDOF = 3;
	for(I=0;I<nNode;I++){
		if(fabs(Node[I].x-0.) <= 1e-15 || fabs(Node[I].y-0.) <= 1e-15 || fabs(Node[I].y-0.41) <= 1e-15 || fabs(Node[I].y-2.2) <= 1e-15)
			mark[I] = 16; //fprintf(Out, "-1\t-1\t1\n"); // (= 1) incognita, diferente c.c. -1 borda externa. -2 borda do cilindro 
		else if(fabs(((Node[I].x-0.2)*(Node[I].x-0.2) + (Node[I].y-0.2)*(Node[I].y-0.2)) - 0.0025) <= 1e-6 )
			mark[I] = 17; //fprintf(Out, "-2\t-2\t1\n"); // (= 1) incognita, diferente c.c. -1 borda externa. -2 borda do cilindro 
		else 
			mark[I] = 15; //fprintf(Out, "1\t1\t1\n");   // (= 1) incognita, diferente c.c. -1 borda externa. -2 borda do cilindro 
	}
	return NDOF;
}