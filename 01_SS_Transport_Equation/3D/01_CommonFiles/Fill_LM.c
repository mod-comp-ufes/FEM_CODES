#include "SSTransportEquation3D.h"

int Fill_LM(int neq, int nel, int **lm, NodeType *Node, ElementType *Element)
{
	int e, i, j, l;

	for (e = 0; e < nel; e++)   //'for' what through all the elements
	{
		l = 0;   //'l' assists in the assembly lm
		for (i = 0; i < 4; i++)   //'for' what through all the nodes of elements
		{
			j = Element[e].Vertex[i];   //'j' stores the global value of the node
			if(Node[j].id == -1){
				lm[e][l] = neq;
				l++;
			}else{
				lm[e][l] = Node[j].id;
				l++;
			}
		}
	}

	return 0;
}

