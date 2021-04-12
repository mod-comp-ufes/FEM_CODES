#include "../01_CommonFiles/ShalowWater.h"

int InitialSolution(ParametersType *Parameters, FemStructsType *FemStructs)
{
    int i, j, id, nnodes = Parameters->nnodes;
	double *u = FemStructs->u;
    NodeType *Node = FemStructs->Node;

    for(i = 0; i<nnodes; i++)
        for(j=0; j<3; j++)
        {
            id = Node[i].id[j];
            if(id != -1)
                u[id] = 0;
        }
}