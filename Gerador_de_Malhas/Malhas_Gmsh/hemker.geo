//Gmsh projeto criado em 08:02:2018 - Problema Cilindro 2D
cl = 1.0;
/*
Na = 13; Ra = 1.00; 
Nb = 11; Rb = 1.00;
Nc = 19; Rc = 1.00;
Nd = 23; Rd = 1.00;
Nr = 21; Rr = 0.94;
*/

Na = 19; Ra = 1.00; 
Nb = 17; Rb = 5.00;
Nc = 29; Rc = 5.00;
Nd = 33; Rd = 1.00;
Nr = 25; Rr = 0.92;
  
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

Line(1) = { 1, 2}; Transfinite Line {1} = Na Using Progression Ra;
Line(2) = { 2, 3}; Transfinite Line {2} = Nd Using Progression Rd;
Line(3) = { 4, 5}; Transfinite Line {3} = Nd Using Progression Rd;
Line(4) = { 5, 6}; Transfinite Line {4} = Na Using Progression Ra;
Line(5) = { 1, 6}; Transfinite Line {5} = Nb Using Bump Rb; // Bump: refina progressivamente
Line(6) = { 2, 5}; Transfinite Line {6} = Nc Using Bump Rc;
Line(7) = { 3, 4}; Transfinite Line {7} = Nc Using Bump Rd;

// Linhas do cilindro
Circle(8)  = { 7, 11, 8}; Transfinite Line {8} = Na Using Progression 0.944;
Circle(9)  = { 8, 11, 9}; Transfinite Line {9} = Nc Using Progression 1.00;
Circle(10) = { 9, 11, 10}; Transfinite Line {10} = Na Using Progression 1.064;
Circle(11) = { 10, 11, 7}; Transfinite Line {11} = Nb Using Progression 1.00;

//Linhas dos Blocos
Line(12) = { 1, 7}; Transfinite Line {12} = Nr Using Progression Rr;
Line(13) = { 2, 8}; Transfinite Line {13} = Nr Using Progression Rr;
Line(14) = { 5, 9}; Transfinite Line {14} = Nr Using Progression Rr;
Line(15) = { 6, 10}; Transfinite Line {15} = Nr Using Progression Rr;

// Superficies
Line Loop(16) = {-12, 1, 13, -8};
Plane Surface(17) = {16};		
Line Loop(18) = {-13, 6, 14, -9};
Plane Surface(19) = {18};		
Line Loop(20) = {-14, 4, 15, -10};
Plane Surface(21) = {20};		
Line Loop(22) = {-15, -5, 12, -11};
Plane Surface(23) = {22};		
Line Loop(24) = {2, 7, 3, -6};
Plane Surface(25) = {24};		

Transfinite Surface {17} Left; 
Transfinite Surface {19} Alternate;
Transfinite Surface {21} Right;
Transfinite Surface {23} Alternate;
Transfinite Surface {25} Alternate;

/*
// Recombina os triangulos em quadrados
Recombine Surface {17};
Recombine Surface {19};
Recombine Surface {21};
Recombine Surface {23};
Recombine Surface {25};

Extrude {0, 0, 1} {
  Surface{17, 19, 21, 23, 25};
  Layers{1};
  Recombine;
};
*/

Mesh.Smoothing = 100;


