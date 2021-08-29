#include "../../gmsh2Mesh.h"
#include <math.h>
	

int setmark_according_boundary(NodeType *Node, int *mark, int nNode)
{
	int NDOF = 3;
	// 0 prescrito, 1 incognita
	for(int I=0; I < nNode; I++) {
		if(fabs(1000000.0 - sqrt(Node[I].x*Node[I].x + Node[I].y*Node[I].y)) <= 1e-8)
			mark[I] = 8; // 1 0 0
		else
			mark[I] = 15; // 1 1 1
	}
	return NDOF;
}
