#include "NavierStokesEquations.h"

int Fill_LM(int neq, int nel, int **lm, NodeType *Node, ElementType *Element)
{
	int e, i, j, k, l;

	for (e = 0; e < nel; e++)   //'for' what through all the elements
	{
		l = 0;   //'l' assists in the assembly lm
		for (i = 0; i < 3; i++)   //'for' what through all the nodes of elements
		{
			j = Element[e].Vertex[i];   //'j' stores the global value of the node
			for(k = 0; k < 3; k++)   //'for' what through all the degrees of freedom
			{
				if(Node[j].id[k] == -1){
					lm[e][l] = neq;
					l++;
				}else{
					lm[e][l] = Node[j].id[k];
					l++;
				}
			}
		}
	}

	return 0;
}


