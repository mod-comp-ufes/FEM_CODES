#include "../01_CommonFiles/ShalowWater.h"


int InitialSolution(ParametersType *Parameters, FemStructsType *FemStructs)
{
    int i, id, nnodes = Parameters->nnodes;
    double x, y;
	double *u = FemStructs->u;
    NodeType *Node = FemStructs->Node;

    for(i=0; i<nnodes; i++)
    {
        x = Node[i].x;
        y = Node[i].y;

        id = Node[i].id[0];
        if(id!=-1) {
            u[id] = ((x-20.0)*(x-20.0) + (y-20.0)*(y-20.0)) <= 6.25 ? 12.5 : 0.5;
        }

        id = Node[i].id[1];
        if(id!=-1) {
            u[id] = 0.0;
        }

        id = Node[i].id[2];
        if(id!=-1) {
            u[id] = 0.0;
        }
    }
}