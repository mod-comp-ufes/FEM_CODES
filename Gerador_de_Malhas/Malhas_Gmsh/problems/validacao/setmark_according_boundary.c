#include <math.h>
#include "../../gmsh2Mesh.h"

int setmark_according_boundary(NodeType *node, int *mark, int nNode) {
	int i;
	int NDOF;

	NDOF = 1;
	for(i = 0; i < nNode; i++) {
		if(fabs(node[i].y) == 0 || fabs(node[i].x) == 1 || fabs(node[i].y) == 1)
			mark[i] = 0;
		else if(fabs(node[i].x) == 0)
			mark[i] = 16;
		else
			mark[i] = 15;
	}
	
	return NDOF;
}