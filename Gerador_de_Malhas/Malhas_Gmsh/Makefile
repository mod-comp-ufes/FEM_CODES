EXE=gmsh2Mesh

hemker:
	gcc problems/hemker/setmark_according_boundary.c gmsh_2_2_2D.c -o $(EXE)

explosion:
	gcc problems/explosion/setmark_according_boundary.c gmsh_2_2_2D.c -o $(EXE)

cavity:
	gcc problems/cavity/setmark_according_boundary.c gmsh_2_2_2D.c -o $(EXE)

channel:
	gcc problems/channel/setmark_according_boundary.c gmsh_2_2_2D.c -o $(EXE)

exata:
	gcc problems/exata/setmark_according_boundary.c gmsh_2_2_2D.c -o $(EXE)

cylinder:
	gcc problems/cylinder/setmark_according_boundary.c gmsh_2_2_2D.c -o $(EXE)

validacao:
	gcc problems/validacao/setmark_according_boundary.c gmsh_4_1_2D.c -o $(EXE)


clean:
	rm $(EXE)