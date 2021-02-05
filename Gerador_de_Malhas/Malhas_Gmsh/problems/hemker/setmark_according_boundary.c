#include "../../gmsh2Mesh.h"
#include <math.h>

int setmark_according_boundary(NodeType *Node, int *mark, int nNode) {
	int I;
	int NDOF;

	NDOF = 1;
	for(I=0;I<nNode;I++){
		if (fabs(Node[I].x*Node[I].x + Node[I].y*Node[I].y - 1.0)<=1e-15 || fabs(Node[I].x + 3.0)<=1e-15)
			mark[I] = 0; //fprintf(Out, "0\n");
		else
			mark[I] = 15; //fprintf(Out, "1\n");
	}
	
	return NDOF;
}