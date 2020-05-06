#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct
{
	double x,y;
}NodeType;

typedef struct
{
	int Vertex[3];
}ElementType;

int main(int argc, char **argv)
{
	int nel, nnodes, I, index, aux1, aux2, aux3, aux4, aux5;
	double z;	
	NodeType *Node;
	ElementType *Element;
	char label[300], filename[300];
	char *aux;
	const char delimiter[3]=" \n";
	FILE *In, *Out;

	if (argc<3){
		printf("Use ./gmsh2Mesh <msh file> <problem>\n");
		//Exemplo ./gmsh2Mesh cavity.msh CAVITY
		exit(1);
	}

	if ((In=fopen(argv[1],"r"))==NULL){
		printf("File %s not found!\n", argv[1]);
		exit(1);
	}		
	
	//Reading header information
	do{
		fscanf(In,"%s",label);
	}while(strcmp(label, "$Nodes")!=0);

	fscanf(In,"%d",&nnodes);
	
	if ((Node = calloc(nnodes,sizeof(NodeType)))==NULL){
		printf("Memmory allocation error in Node!\n");
		exit(1);
	}

	//Reading nodes
	for (I=0; I<nnodes; I++)
		fscanf(In,"%d%lf%lf%lf", &index, &Node[I].x, &Node[I].y, &z);
	
	//Reading header information
	do{
		fscanf(In,"%s",label);
	}while(strcmp(label, "$Elements")!=0);

	fscanf(In,"%d\n",&nel);

	do{
		fgets(label, 300, In);
		aux = strtok(label, delimiter);
		aux = strtok(NULL, delimiter);
		nel--;
	}while (strcmp(aux, "2")!=0);

	nel++;
	
	if ((Element = calloc(nel,sizeof(ElementType)))==NULL){
		printf("Memmory allocation error in Element!\n");
		exit(1);
	}

	aux = strtok(NULL, delimiter);
	aux = strtok(NULL, delimiter);
	aux = strtok(NULL, delimiter);
	aux = strtok(NULL, delimiter);
	Element[0].Vertex[0] = atoi(aux);
	aux = strtok(NULL, delimiter);
	Element[0].Vertex[1] = atoi(aux);
	aux = strtok(NULL, delimiter);
	Element[0].Vertex[2] = atoi(aux);
	
	for(I=1; I<nel;I++)
		fscanf(In, "%d%d%d%d%d%d%d%d", &aux1, &aux2, &aux3, &aux4, &aux5, &Element[I].Vertex[0], &Element[I].Vertex[1], &Element[I].Vertex[2]);  	
	
	fclose(In);

	if (strcasecmp(argv[2],"hemker")==0){
		sprintf(filename,"%s_%d_%d.dat",argv[2],nnodes,nel);
		Out = fopen(filename,"w");
		fprintf(Out,"%d\n",nnodes);
		for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			if (fabs(Node[I].x*Node[I].x + Node[I].y*Node[I].y - 1.0)<=1e-15 || fabs(Node[I].x + 3.0)<=1e-15)
				fprintf(Out, "0\n");
			else
				fprintf(Out, "1\n");
		}
		fprintf(Out,"%d\n",nel);
		for(I=0;I<nel;I++)
			fprintf(Out, "%d\t%d\t%d\n", Element[I].Vertex[0]-1, Element[I].Vertex[1]-1, Element[I].Vertex[2]-1);
		
		fclose(Out);

	}
	else if (strcasecmp(argv[2],"explosion")==0){
		sprintf(filename,"%s_%d_%d.dat",argv[2],nnodes,nel);
		Out = fopen(filename,"w");
		fprintf(Out,"%d\n",nnodes);
		for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			if (fabs(Node[I].x)<=1e-15 || fabs(Node[I].y)<=1e-15 || fabs(Node[I].x-2)<=1e-15 || fabs(Node[I].y-2)<=1e-15)
				fprintf(Out, "0\t0\t0\t0\n");
			else
				fprintf(Out, "1\t1\t1\t1\n");
		}
		fprintf(Out,"%d\n",nel);
		for(I=0;I<nel;I++)
			fprintf(Out, "%d\t%d\t%d\t-1\n", Element[I].Vertex[0]-1, Element[I].Vertex[1]-1, Element[I].Vertex[2]-1);
		
		fclose(Out);

	}else if (strcasecmp(argv[2],"cavity")==0){
		sprintf(filename,"%s_%d_%d.dat",argv[2],nnodes,nel);
		Out = fopen(filename,"w");
		fprintf(Out,"%d\n",nnodes);
		for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			/*
			if(fabs(Node[I].y-1.)<=1e-15 && fabs(Node[I].x-0.)<=1e-15)
				fprintf(Out, "0\t0\t0\n");
			else
			*/			
			if (fabs(Node[I].y-1)<=1e-15 || fabs(Node[I].y-0)<=1e-15 || fabs(Node[I].x-1)<=1e-15 || fabs(Node[I].x-0)<=1e-15)
				fprintf(Out, "0\t0\t1\n");
			else 
				fprintf(Out, "1\t1\t1\n");
		}
		fprintf(Out,"%d\n",nel);
		for(I=0;I<nel;I++)
			fprintf(Out, "%d\t%d\t%d\t-1\n", Element[I].Vertex[0]-1, Element[I].Vertex[1]-1, Element[I].Vertex[2]-1);
		
		fclose(Out);

	}else if (strcasecmp(argv[2],"channel")==0){
		sprintf(filename,"%s_%d_%d.dat",argv[2],nnodes,nel);
		Out = fopen(filename,"w");
		fprintf(Out,"%d\n",nnodes);
		//dominio old		
		/*for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			if(fabs(Node[I].x-0.)<=1e-15 || fabs(Node[I].y-0.)<=1e-15 || fabs(Node[I].y-1.5)<=1e-15)
				fprintf(Out, "0\t0\t1\n");
			else if (fabs(Node[I].y-0.75)<=1e-15 && Node[I].x <= 7.5 )
				fprintf(Out, "0\t0\t1\n");
			else if (fabs(Node[I].x-7.5)<=1e-15 && Node[I].y <= 0.75 )
				fprintf(Out, "0\t0\t1\n");
			else if (fabs(Node[I].x-30.0)<=1e-15 )
				fprintf(Out, "1\t0\t1\n");
			else 
				fprintf(Out, "1\t1\t1\n");
		}*/
		// Dominio Erturk artigo 2008 com h_i = h = 0,75 
		for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			if(fabs(Node[I].x-0.)<=1e-15 || fabs(Node[I].y-0.)<=1e-15 || fabs(Node[I].y-1.5)<=1e-15)
				fprintf(Out, "0\t0\t1\n");
			else if (fabs(Node[I].y-0.75)<=1e-15 && Node[I].x <= 15.0 )
				fprintf(Out, "0\t0\t1\n");
			else if (fabs(Node[I].x-15)<=1e-15 && Node[I].y <= 0.75 )
				fprintf(Out, "0\t0\t1\n");
			else if (fabs(Node[I].x-240.0)<=1e-15 )
				fprintf(Out, "1\t0\t1\n");
			else 
				fprintf(Out, "1\t1\t1\n");
		}
		fprintf(Out,"%d\n",nel);
		for(I=0;I<nel;I++)
			fprintf(Out, "%d\t%d\t%d\t-1\n", Element[I].Vertex[0]-1, Element[I].Vertex[1]-1, Element[I].Vertex[2]-1);
		
		fclose(Out);

	}else if (strcasecmp(argv[2],"exata")==0){
		sprintf(filename,"%s_%d_%d.dat",argv[2],nnodes,nel);
		Out = fopen(filename,"w");
		fprintf(Out,"%d\n",nnodes);
		for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			if(fabs(Node[I].y-0.)<=1e-15 && fabs(Node[I].x-0.)<=1e-15)
				fprintf(Out, "0\t0\t0\n");
			else if (fabs(Node[I].y-1)<=1e-15 || fabs(Node[I].y-0)<=1e-15 || fabs(Node[I].x-1)<=1e-15 || fabs(Node[I].x-0)<=1e-15)
				fprintf(Out, "0\t0\t1\n");
			else 
				fprintf(Out, "1\t1\t1\n");
		}
		fprintf(Out,"%d\n",nel);
		for(I=0;I<nel;I++)
			fprintf(Out, "%d\t%d\t%d\t-1\n", Element[I].Vertex[0]-1, Element[I].Vertex[1]-1, Element[I].Vertex[2]-1);
		
		fclose(Out);

	}else if (strcasecmp(argv[2],"cylinder")==0){
		sprintf(filename,"%s_%d_%d.dat",argv[2],nnodes,nel);
		Out = fopen(filename,"w");
		fprintf(Out,"%d\n",nnodes);
		 //dominio old		
		/*for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			if(fabs(Node[I].x-0.) <= 1e-15 || fabs(Node[I].y-0.) <= 1e-15 || fabs(Node[I].y-9.) <= 1e-15 || fabs(((Node[I].x-4.5)*(Node[I].x-4.5) + (Node[I].y-4.5)*(Node[I].y-4.5)) - 0.25) <= 1e-6 )
				fprintf(Out, "0\t0\t1\n");
			else 
				fprintf(Out, "1\t1\t1\n");
		}*/
		//Dominio Turek/Volker 
		/*for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			if(fabs(Node[I].x-0.) <= 1e-15 || fabs(Node[I].y-0.) <= 1e-15 || fabs(Node[I].y-0.41) <= 1e-15 || fabs(((Node[I].x-0.2)*(Node[I].x-0.2) + (Node[I].y-0.2)*(Node[I].y-0.2)) - 0.0025) <= 1e-6 )
				fprintf(Out, "0\t0\t1\n");
			else 
				fprintf(Out, "1\t1\t1\n");
		}*/
		//Dominio Turek/Volker Jhon com fronteira preparadas para calculo de coef. Lift e Drag
		for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			if(fabs(Node[I].x-0.) <= 1e-15 || fabs(Node[I].y-0.) <= 1e-15 || fabs(Node[I].y-0.41) <= 1e-15 || fabs(Node[I].y-2.2) <= 1e-15)
				fprintf(Out, "-1\t-1\t1\n"); // (= 1) incognita, diferente c.c. -1 borda externa. -2 borda do cilindro 
			else if(fabs(((Node[I].x-0.2)*(Node[I].x-0.2) + (Node[I].y-0.2)*(Node[I].y-0.2)) - 0.0025) <= 1e-6 )
				fprintf(Out, "-2\t-2\t1\n"); // (= 1) incognita, diferente c.c. -1 borda externa. -2 borda do cilindro 
			else 
				fprintf(Out, "1\t1\t1\n");   // (= 1) incognita, diferente c.c. -1 borda externa. -2 borda do cilindro 
		}
		//Dominio Erturk
		/*for(I=0;I<nnodes;I++){
			fprintf(Out, "%lf\t%lf\t", Node[I].x, Node[I].y);
			if(fabs(Node[I].x+5.0) <= 1e-15 || fabs(Node[I].y-5.0) <= 1e-15 || fabs(((Node[I].x-0.0)*(Node[I].x-0.0) + (Node[I].y-0.0)*(Node[I].y-0.0)) - 0.01) <= 1e-15 )
				fprintf(Out, "0\t0\t1\n");
			else if(fabs(Node[I].y-0.0) <= 1e-15 || fabs(Node[I].x-20.0) <= 1e-15 )
				fprintf(Out, "1\t0\t1\n");
			else 
				fprintf(Out, "1\t1\t1\n");
		}
		*/
		fprintf(Out,"%d\n",nel);
		for(I=0;I<nel;I++)
			fprintf(Out, "%d\t%d\t%d\t-1\n", Element[I].Vertex[0]-1, Element[I].Vertex[1]-1, Element[I].Vertex[2]-1);
		
		fclose(Out); 
	}
	else{
		printf("Problem not defined!\n");		
		exit(1);
	}	
	return 0;
}


