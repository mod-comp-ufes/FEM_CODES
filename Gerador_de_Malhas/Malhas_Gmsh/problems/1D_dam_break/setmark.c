#include "../../gmsh2Mesh.h"
#include <math.h>
	

int setmark_according_boundary(NodeType *Node, int *mark, int nNode)
{
	int NDOF = 3;
	// 0 prescrito, 1 incognita
	for(int I=0;I<nNode;I++) {
		if((fabs(50-fabs(Node[I].x)) <= 1e-8) && (fabs(1-Node[I].y) <= 1e-8 || fabs(Node[I].y) <= 1e-8))
			mark[I] = 0; // 0 0 0
		else if(fabs(50-fabs(Node[I].x)) <= 1e-8)
			mark[I] = 2; // 0 0 1
		else if(fabs(1-Node[I].y) <= 1e-8 || fabs(Node[I].y) <= 1e-8)
			mark[I] = 12; // 1 1 0
		else
			mark[I] = 15; // 1 1 1
	}
	return NDOF;
}
