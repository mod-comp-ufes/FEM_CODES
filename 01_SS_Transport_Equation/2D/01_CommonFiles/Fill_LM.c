#include "SSTranspEquation.h"

int Fill_LM(int neq, int nel, int **lm, NodeType *Node, ElementType *Element)
{
	int J,K,JJ;

	for (K=0; K<nel; K++)
	{
		lm[K][0]=neq;
		lm[K][1]=neq;
		lm[K][2]=neq;
		for (J=0; J<3; J++)
		{
			JJ = Element[K].Vertex[J];

			if (Node[JJ].Type>=1)
				lm[K][J] = Node[JJ].id;

		}
	}

	return 0;

}


