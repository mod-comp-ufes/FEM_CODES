# Malhas oriundas do GMSH
Para gerar as malhas oriundas do gmsh basta executar os passos:

- Passo 1:\
Descrever o arquivo de configuracao do gmsh (<problem>.geo)

-Passo 2:\
Gerar a malha como gmsh

- Passo 3:\
O passo 2 vai ter como resultado um arquivo <problem>.msh. Utilize esse arquivo da seguinte forma:
```bash
gmsh2Mesh <problem>.msh <PROBLEM>
```

Note que a definição de <PROBLEM> é apenas utilizada para geração do arquivo <problem>_<nNodes>_<nElements>.dat

Para adiconar novos problemas, basta implementar a função `setmark_according_boundary()` e ligá-la a uma das implementações existentes (gmsh_2_2_2D.c ou gmsh_4_1_2D), onde são adicionadas as marcas de cada nó, com **0: nó prescrito, 1: nó incógnita**. Para exemplos de ligação, vide Makefile. Também note que existem implementações diferentes de acordo com o tipo de saída do GMSH. gmsh_2_2_2D refere-se à saída 2.2 em problemas bidimensionais, enquanto gmsh_4_1_2D refere-se à saída 4.1 também em problemas bidimensionais.

```c
int setmark_according_boundary(NodeType *node, int *mark, int nNode)
/*
    Arguments:
        node (NodeType*). Nós lidos do arquivo .msh. Sua implementação é descrita no arquivo gmsh2Mesh.h.
        mark (int*). Marcas a serem adicionadas.
        nNode (int). Número de nós.

    Return:
        NDOF (int). Número de graus de liberdade.

    Obs. Valores das marcas: Note que estas são escritas da esquerda para direita, de acordo com o valor de NDOF.
    {{0,0,0,0}/*mark 0*/, {0,0,0,1}/*mark 1*/, {0,0,1,0}/*mark 2*/, {0,0,1,1}/*mark 3*/,
     {0,1,0,0}/*mark 4*/, {0,1,0,1}/*mark 5*/, {0,1,1,0}/*mark 6*/, {0,1,1,1}, /*mark 7*/
     {1,0,0,0}/*mark 8*/, {1,0,0,1}/*mark 9*/, {1,0,1,0}/*mark 10*/, {1,0,1,1}/*mark 11*/,
     {1,1,0,0}/*mark 12*/, {1,1,0,1}/*mark 13*/, {1,1,1,0}/*mark 14*/, {1,1,1,1}, /*mark 15*/}
*/
```