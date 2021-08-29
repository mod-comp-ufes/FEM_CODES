#include "../../gmsh2Mesh.h"
#include <math.h>

//{{0,0,0,0}/*mark 0 */, {0,0,0,1}/*mark 1 */, {0,0,1,0}/*mark 2 */, {0,0,1,1}/*mark 3 */,
// {0,1,0,0}/*mark 4 */, {0,1,0,1}/*mark 5 */, {0,1,1,0}/*mark 6 */, {0,1,1,1}, /*mark 7 */
// {1,0,0,0}/*mark 8 */, {1,0,0,1}/*mark 9 */, {1,0,1,0}/*mark 10 */, {1,0,1,1}/*mark 11 */,
// {1,1,0,0}/*mark 12 */, {1,1,0,1}/*mark 13 */, {1,1,1,0}/*mark 14 */, {1,1,1,1} /*mark 15 */,
// {2,1,1,1}/*mark 16 */ };

int setmark_according_boundary(NodeType *Node, int *mark, int nNode)
{
	int NDOF = 3;
	// 0 prescrito, 1 incognita
	for(int I=0;I<nNode;I++) {
		if((fabs(40-fabs(Node[I].x)) <= 1e-8) && (fabs(40-fabs(Node[I].y)) <= 1e-8))
			mark[I] = 8; // 1 0 0
		else if(fabs(40-fabs(Node[I].x)) <= 1e-8)
			mark[I] = 10; // 1 0 1
		else if(fabs(40-fabs(Node[I].y)) <= 1e-8)
			mark[I] = 12; // 1 1 0
		else
			mark[I] = 15; // 1 1 1
	}
	return NDOF;
}
