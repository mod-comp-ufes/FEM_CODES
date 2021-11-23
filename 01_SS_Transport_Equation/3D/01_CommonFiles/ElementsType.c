/* Definimos os tipos de elemento de acordo com a fronteira que ele está.
 * 1: plano y = yInf
 * 2: plano x = xSup
 * 3: plano y = ySup
 * 4: plano x = xInf
 * 5: plano z = zSup
 * 6: plano z = zInf
 * 7: esquina dos planos y = 0, z = 0 e x = 0 vertice (xInf,yInf,zInf)
 * 8: esquina dos planos y = 0, z = 0 e x = 1 vertice (xSup,yInf,yInf)
 * 9: esquina dos planos z = 0, y = 1 e x = 1 vertice (1,1,0)
 * 10: esquina dos planos z = 0, y = 1 e x = 0 vertice (0,1,0)
 * 11: esquina dos planos y = 0, z = 1 e x = 0 vertice (0,0,1)
 * 12: esquina dos planos y = 0, z = 1 e x = 1 vertice (1,0,1)
 * 13: esquina dos planos z = 1, y = 1 e x = 1 vertice (1,1,1)
 * 14: esquina dos planos z = 1, y = 1 e x = 0 vertice (0,1,1)
 * 15: reta entre os planos y = 0 e z = 1
 * 16: reta entre os planos y = 0 e x = 1
 * 17: reta entre os planos y = 0 e z = 0
 * 18: reta entre os planos y = 0 e x = 0
 * 19: reta entre os planos x = 1 e z = 0
 * 20: reta entre os planos x = 1 e y = 1
 * 21: reta entre os planos x = 1 e z = 1
 * 22: reta entre os planos y = 1 e z = 0
 * 23: reta entre os planos y = 1 e x = 0
 * 24: reta entre os planos y = 1 e z = 1
 * 25: reta entre os planos z = 1 e x = 0
*/

#include "SSTransportEquation3D.h"

void ElementsType(ParametersType *Parameters, NodeType *Node, ElementType *Element)
{

	int I, J, p, cont1, cont2, cont3, cont4, cont5, cont6;
	int nel = Parameters->nel;
	double xInf = Parameters->xInf;
	double xSup = Parameters->xSup;
	double yInf = Parameters->yInf;
	double ySup = Parameters->ySup;
	double zInf = Parameters->zInf;
	double zSup = Parameters->zSup;
	
	for(I = 0; I < nel; I++){
		cont1 = 0; cont2 = 0; cont3 = 0; cont4 = 0; cont5 = 0; cont6 = 0;
		for(J = 0; J < NNOEL; J++){
			p = Element[I].Vertex[J];
			if( Node[p].y == yInf){ // plano 1: y = yInf 
				cont1 = cont1 + 1;
			}
			if(Node[p].x == xSup){ // plano 2: x = xSup
				cont2 = cont2 + 1;
			}
			if(Node[p].y == ySup){ // plano 3: y = ySup
				cont3 = cont3 + 1;
			}
			if(Node[p].x == xInf){ // plano 4: x = xInf
				cont4 = cont4 + 1;
			}
			if(Node[p].z == zSup){ // plano 5: z = zSup
				cont5 = cont5 + 1;
			}
			if(Node[p].z == zInf){ // plano 6: z = zInf
				cont6 = cont6 + 1;
			}
		}
		
		if(cont1 == 3){ // plano 1: y = 0 
			Element[I].Type = 1;
			if (cont1 == cont6){ // reta entre os planos y = 0 e z = 0
				Element[I].Type = 17;
				if(cont1 == cont4){ // vértice (0,0,0)
					Element[I].Type = 7;
				}
				if(cont1 == cont2){ // vértice (Lx,0,0)
					Element[I].Type = 8;
				}
			}else if (cont1 == cont5){ // reta entre os planos y = 0 e z = Lz
				Element[I].Type = 15;
				if(cont1 == cont4){ // vértice (0,0,Lz)
					Element[I].Type = 11;
				}
				if(cont1 == cont2){ // vértice (Lx,0,Lz)
					Element[I].Type = 12;
				}
			}else if (cont1 == cont2){ // reta entre os planos y = 0 e x = Lx
				Element[I].Type = 16;
			}else if (cont1 == cont4){ // reta entre os planos y = 0 e x = 0
				Element[I].Type = 18;
			}
			
		}else if(cont2 == 3){ // plano 2: x = Lx
			Element[I].Type = 2;
			if (cont2 == cont6){ // reta entre os planos x = Lx e z = 0
				Element[I].Type = 19;
			}else if (cont2 == cont5){ // reta entre os planos x = Lx e z = Lz
				Element[I].Type = 18;
			}
			
		}else if(cont3 == 3){ // plano 3: y = Ly
			Element[I].Type = 3;
			if (cont3 == cont6){ // reta entre os planos y = Ly e z = 0
				Element[I].Type = 22;
				if(cont3 == cont4){ // vértice (0,Ly,0)
					Element[I].Type = 10;
				}
				if(cont3 == cont2){ // vértice (Lz,Ly,0)
					Element[I].Type = 9;
				}
			}else if (cont3 == cont5){ // reta entre os planos y = Ly e z = Lz
				Element[I].Type = 24;
				if(cont3 == cont4){ // vértice (0,Ly,Lz)
					Element[I].Type = 14;
				}
				if(cont3 == cont2){ // vértice (Lx,Ly,Lz)
					Element[I].Type = 13;
				}
			}else if (cont3 == cont2){ // reta entre os planos y = Ly e x = Lx
				Element[I].Type = 20;
			}else if (cont3 == cont4){ // reta entre os planos y = Ly e x = 0
				Element[I].Type = 23;
			}
			
		}else if(cont4 == 3){ // plano 4: x = 0
			Element[I].Type = 4;
			if (cont4 == cont6){ // reta entre os planos x = 0 e z = 0
				Element[I].Type = 20;
			}else if (cont4 == cont5){ // reta entre os planos x = 0 e z = Lz
				Element[I].Type = 25;
			}
			
		}else if(cont5 == 3){ // plano 5: z = Lz
			Element[I].Type = 5;
			
		}else if(cont6 == 3){ // plano 6: z = 0
			Element[I].Type = 6;
			
		}else{ // não está em nenhum contorno
			Element[I].Type = -1;
		}
	}
	
}
