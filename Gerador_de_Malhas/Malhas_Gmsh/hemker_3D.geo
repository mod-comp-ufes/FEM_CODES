//Gmsh projeto criado em 08:02:2018 - Problema Cilindro 2D
cl = 1.0;
Point(1)  = {  0, 0, 0, cl};
Point(2)  = {  9, 0, 0, cl};
Point(3)  = { 20, 0, 0, cl};
Point(4)  = { 20, 9, 0, cl};
Point(5)  = {  9, 9, 0, cl};
Point(6)  = {  0, 9, 0, cl};
//Pontos do cilindro
//Raio = 0.5; 
Point(7)  = { 4.14644660, 4.14644660, 0, cl};
Point(8)  = { 4.85355339, 4.14644660, 0, cl};
Point(9)  = { 4.85355339, 4.85355339, 0, cl};
Point(10) = { 4.14644660, 4.85355339, 0, cl};
Point(11) = { 4.5, 4.5, 0, cl};

Line(1) = { 1, 2}; Transfinite Line {1} = 10 Using Progression 1;
Line(2) = { 2, 3}; Transfinite Line {2} = 10 Using Progression 1;
Line(3) = { 4, 5}; Transfinite Line {3} = 10 Using Progression 1;
Line(4) = { 5, 6}; Transfinite Line {4} = 10 Using Progression 1;
Line(5) = { 1, 6}; Transfinite Line {5} = 10 Using Progression 1;
Line(6) = { 2, 5}; Transfinite Line {6} = 10 Using Progression 1;
Line(7) = { 3, 4}; Transfinite Line {7} = 10 Using Progression 1;

// Linhas do cilindro
Circle(8)  = { 7, 11, 8}; Transfinite Line {8} = 10 Using Progression 1;
Circle(9)  = { 8, 11, 9}; Transfinite Line {9} = 10 Using Progression 1;
Circle(10) = { 9, 11, 10}; Transfinite Line {10} = 10 Using Progression 1;
Circle(11) = { 10, 11, 7}; Transfinite Line {11} = 10 Using Progression 1;

//Linhas dos Blocos
Line(12) = { 1, 7}; Transfinite Line {12} = 10 Using Progression 1;
Line(13) = { 2, 8}; Transfinite Line {13} = 10 Using Progression 1;
Line(14) = { 5, 9}; Transfinite Line {14} = 10 Using Progression 1;
Line(15) = { 6, 10}; Transfinite Line {15} = 10 Using Progression 1;

// Superficies
Line Loop(16) = {12, 8, -13, -1};
Plane Surface(17) = {16};
Line Loop(18) = {13, 9, -14, -6};
Plane Surface(19) = {18};
Line Loop(20) = {14, 10, -15, -4};
Plane Surface(21) = {20};
Line Loop(22) = {15, 11, -12, 5};
Plane Surface(23) = {22};
Line Loop(24) = {2, 7, 3, -6};
Plane Surface(25) = {24};

Transfinite Surface {17};
Transfinite Surface {19};
Transfinite Surface {21};
Transfinite Surface {23};
Transfinite Surface {25};

Recombine Surface {17};
Recombine Surface {19};
Recombine Surface {21};
Recombine Surface {23};
Recombine Surface {25};

Extrude {0, 0, 1} {
  Surface{17, 19, 21, 23, 25};
  Layers{1};
  Recombine;
}






