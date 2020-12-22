#include "../../gmsh2Mesh.h"
#include <math.h>


int setmark_according_boundary(NodeType *Node, int *mark, int nNode) {
	int I;
	int NDOF;

	NDOF = 4;
	for(I=0;I<nNode;I++){
		if (fabs(Node[I].x)<=1e-15 || fabs(Node[I].y)<=1e-15 || fabs(Node[I].x-2)<=1e-15 || fabs(Node[I].y-2)<=1e-15)
			mark[I] = 0; //fprintf(Out, "0\t0\t0\t0\n");
		else
			mark[I] = 15; //fprintf(Out, "1\t1\t1\t1\n");
	}
	
	return NDOF;
}