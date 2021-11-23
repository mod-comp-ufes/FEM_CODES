#include "SSTransportEquation3D.h"
#include "../../../00_CommonFiles/IO_Operations/io.h"

double Element_Configuration(int e, ElementType *Element, NodeType *Node, double *X, double *Y, double *Z, double *a, double *b, double *c, double *d, double *Ue, double *uSpace, double *uSpaceB){
	
	int J1, J2, J3, J4;
	double Volume;
	
	// Global node that composes the element
	J1 = Element[e].Vertex[0];
	J2 = Element[e].Vertex[1];
	J3 = Element[e].Vertex[2];
	J4 = Element[e].Vertex[3];
    
	// Nodal coordinates
	X[0] = Node[J1].x;
	X[1] = Node[J2].x;
	X[2] = Node[J3].x;
	X[3] = Node[J4].x;
	
	Y[0] = Node[J1].y;
	Y[1] = Node[J2].y;
	Y[2] = Node[J3].y;
	Y[3] = Node[J4].y;
	
	Z[0] = Node[J1].z;
	Z[1] = Node[J2].z;
	Z[2] = Node[J3].z;
	Z[3] = Node[J4].z;
	
	// Component of the form function
	a[0] =  (X[1]*Y[2]*Z[3] + X[3]*Y[1]*Z[2] + X[2]*Y[3]*Z[1]) - (X[3]*Y[2]*Z[1] + X[1]*Y[3]*Z[2] + X[2]*Y[1]*Z[3]);
	a[1] = -(X[0]*Y[2]*Z[3] + X[3]*Y[0]*Z[2] + X[2]*Y[3]*Z[0]) + (X[3]*Y[2]*Z[0] + X[0]*Y[3]*Z[2] + X[2]*Y[0]*Z[3]);
	a[2] =  (X[0]*Y[1]*Z[3] + X[3]*Y[0]*Z[1] + X[1]*Y[3]*Z[0]) - (X[3]*Y[1]*Z[0] + X[0]*Y[3]*Z[1] + X[1]*Y[0]*Z[3]);
	a[3] = -(X[0]*Y[1]*Z[2] + X[2]*Y[0]*Z[1] + X[1]*Y[2]*Z[0]) + (X[2]*Y[1]*Z[0] + X[0]*Y[2]*Z[1] + X[1]*Y[0]*Z[2]);
	
	b[0] = -(Y[1]-Y[3])*(Z[2]-Z[3]) + (Y[2]-Y[3])*(Z[1]-Z[3]);
	b[1] =  (Y[0]-Y[3])*(Z[2]-Z[3]) - (Y[2]-Y[3])*(Z[0]-Z[3]);
	b[2] = -(Y[0]-Y[3])*(Z[1]-Z[3]) + (Y[1]-Y[3])*(Z[0]-Z[3]);
	b[3] =  (Y[0]-Y[2])*(Z[1]-Z[2]) - (Y[1]-Y[2])*(Z[0]-Z[2]);
	
	c[0] =  (X[1]-X[3])*(Z[2]-Z[3]) - (X[2]-X[3])*(Z[1]-Z[3]);
	c[1] = -(X[0]-X[3])*(Z[2]-Z[3]) + (X[2]-X[3])*(Z[0]-Z[3]);
	c[2] =  (X[0]-X[3])*(Z[1]-Z[3]) - (X[1]-X[3])*(Z[0]-Z[3]);
	c[3] = -(X[0]-X[2])*(Z[1]-Z[2]) + (X[1]-X[2])*(Z[0]-Z[2]);
	
	d[0] = -(X[1]-X[3])*(Y[2]-Y[3]) + (X[2]-X[3])*(Y[1]-Y[3]);
	d[1] =  (X[0]-X[3])*(Y[2]-Y[3]) - (X[2]-X[3])*(Y[0]-Y[3]);
	d[2] = -(X[0]-X[3])*(Y[1]-Y[3]) + (X[1]-X[3])*(Y[0]-Y[3]);
	d[3] =  (X[0]-X[2])*(Y[1]-Y[2]) - (X[1]-X[2])*(Y[0]-Y[2]); 

	//printf("Dentro do Element\n");
	//printf("b1+b2+b3+b4 = %lf\n", b[1]+b[2]+b[3]+b[0]);
	//printf("c1+c2+c3+c4 = %lf\n", c[1]+c[2]+c[3]+c[0]);
	//printf("d1+d2+d3+d4 = %lf\n", d[1]+d[2]+d[3]+d[0]);

	// Volume
	Volume = (a[0] + X[0]*b[0] + Y[0]*c[0] + Z[0]*d[0])/6.0;
	
/*	FILE *OutFile;
	char FileName[2000];
	
	sprintf(FileName,"../03_Output/DIFCONV/US5_ElementVolume_N76889_E431380.txt");
	
	OutFile = myfopen(FileName,"a");
	
	fprintf(OutFile, "Elemento %d \t Volume %lf\n", e, Volume);
	printf("Elemento %d \t Volume %lf\n", e, Volume);
	
	fclose(OutFile);
	
	if(e%10000 == 0){
		getchar();
	}*/
		
	// calculates the delayed u in space
	Ue[0] = uSpace[J1];
	Ue[1] = uSpace[J2];
	Ue[2] = uSpace[J3];
	Ue[3] = uSpace[J4];
	
	// tetrahedro centroid
	*uSpaceB = (Ue[0] + Ue[1] + Ue[2] + Ue[3])/4.0;
	
/*	double x21 = X[1] - X[0];
	double x31 = X[2] - X[0];
	double x41 = X[3] - X[0];
	double y21 = Y[1] - Y[0];
	double y31 = Y[2] - Y[0];
	double y41 = Y[3] - Y[0];
	double z21 = Z[1] - Z[0];
	double z31 = Z[2] - Z[0];
	double z41 = Z[3] - Z[0];
		
		// Jacobian Matrix (dXi/dx)		
		double invJ[3][3];
		double sixV = 6.0*Volume;
		
		invJ[0][0] = ((y31*z41) - (z31*y41))/sixV;
		invJ[1][0] = (- (x31*z41) + (z31*x41))/sixV;
		invJ[2][0] = ((x31*y41) - (y31*x41))/sixV; 
		invJ[0][1] = (- (y21*z41) + (z21*y41))/sixV; 
		invJ[1][1] = ((x21*z41) - (z21*x41))/sixV;
		invJ[2][1] = (- (x21*y41) + (y21*x41))/sixV;
		invJ[0][2] = ((y21*z31) - (z21*y31))/sixV;
		invJ[1][2] = (- (x21*z31) + (z21*x31))/sixV; 
		invJ[2][2] = ((x21*y31) - (x31*y21))/sixV;


	printf("[0][0]= %lf \t b2/6V = %lf\n", invJ[0][0], b[1]/sixV);
	printf("[1][0]= %lf \t b3/6V = %lf\n", invJ[1][0], b[2]/sixV);
	printf("[2][0]= %lf \t b4/6V = %lf\n", invJ[2][0], b[3]/sixV);
	printf("[0][1]= %lf \t c2/6V = %lf\n", invJ[0][1], c[1]/sixV);
	printf("[1][1]= %lf \t c3/6V = %lf\n", invJ[1][1], c[2]/sixV);
	printf("[2][1]= %lf \t c4/6V = %lf\n", invJ[2][1], c[3]/sixV);
	printf("[0][2]= %lf \t d2/6V = %lf\n", invJ[0][2], d[1]/sixV);
	printf("[1][2]= %lf \t d3/6V = %lf\n", invJ[1][2], d[2]/sixV);
	printf("[2][2]= %lf \t d4/6V = %lf\n", invJ[2][2], d[3]/sixV);
	getchar(); */

	return Volume;
}
