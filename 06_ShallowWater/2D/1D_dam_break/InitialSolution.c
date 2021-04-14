#include "../01_CommonFiles/ShalowWater.h"


int InitialSolution(ParametersType *Parameters, FemStructsType *FemStructs)
{
    int i, id, nnodes = Parameters->nnodes;
	double *u = FemStructs->u;
    NodeType *Node = FemStructs->Node;

    for(i=0; i<nnodes; i++)
    {
        id = Node[i].id[0];
        if(id != -1)
            u[id] = 0.0;

        id = Node[i].id[1];
        if(id != -1)
            u[id] = 0.0;

        id = Node[i].id[2];
        if(id != -1)
            u[id] = 10e10;
    }
}