#include "../../gmsh2Mesh.h"
#include <math.h>


int setmark_according_boundary(NodeType *Node, int *mark, int nNode) {
	int I;
	int NDOF;

	NDOF = 3;
	for(I=0;I<nNode;I++){
		if(fabs(Node[I].x-0.)<=1e-15 || fabs(Node[I].y-0.)<=1e-15 || fabs(Node[I].y-1.5)<=1e-15)
		    mark[I] = 2; //fprintf(Out, "0\t0\t1\n");
		else if (fabs(Node[I].y-0.75)<=1e-15 && Node[I].x <= 15.0 )
		    mark[I] = 2; //fprintf(Out, "0\t0\t1\n");
		else if (fabs(Node[I].x-15)<=1e-15 && Node[I].y <= 0.75 )
		    mark[I] = 2; //fprintf(Out, "0\t0\t1\n");
		else if (fabs(Node[I].x-240.0)<=1e-15 )
		    mark[I] = 10; //fprintf(Out, "1\t0\t1\n");
		else 
		    mark[I] = 15; //fprintf(Out, "1\t1\t1\n");
	}
	return NDOF;
}