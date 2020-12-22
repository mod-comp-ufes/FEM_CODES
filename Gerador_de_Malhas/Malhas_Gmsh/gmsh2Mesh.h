#ifndef _GMSH2MESH_H_
#define _GMSH2MESH_H_

typedef struct
{
	int id;
	double x,y,z;
} NodeType;

typedef struct
{
	int id;
	int Type;
	int Vertex[3];
} ElementType;


#endif