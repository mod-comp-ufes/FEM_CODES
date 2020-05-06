//Gmsh projeto criado em 08:02:2018 - Problema Cilindro 2D
cl = 1.0;
/*
Na = 13; Ra = 1.00; 
Nb = 11; Rb = 1.00;
Nc = 19; Rc = 1.00;
Nd = 23; Rd = 1.00;
Nr = 21; Rr = 0.94;
*/

Na = 89; Ra = 1.03; 
Nb = 89; Rb = 4.00;
Nc = 101; Rc = 1.01;
Nd = 105; Rd = 1.02; Rd2 = 0.98;
Nr = 75; Rr = 0.96;
  
Point(1)  = {  -5.0, 0.0, 0.0, cl};
Point(2)  = {  -0.1, 0.0, 0.0, cl};
Point(3)  = {   0.1, 0.0, 0.0, cl};
Point(4)  = {   5.0, 0.0, 0.0, cl};
Point(5)  = {  20.0, 0.0, 0.0, cl};
Point(6)  = {  20.0, 5.0, 0.0, cl};
Point(7)  = {   5.0, 5.0, 0.0, cl};
Point(8)  = {  -5.0, 5.0, 0.0, cl};
//Pontos do cilindro
//Raio = 0.1; 
Point(9)  = { -0.0707106781, 0.0707106781, 0, cl};
Point(10) = {  0.0707106781, 0.0707106781, 0, cl};
Point(11) = {          0.0,         0.0, 0, cl};

Line(1) = { 1, 2}; Transfinite Line {1} = Nr Using Progression Rr;//Using Progression Ra;
Line(2) = { 4, 3}; Transfinite Line {2} = Nr Using Progression Rr;
Line(3) = { 7, 4}; Transfinite Line {3} = Na Using Progression Rr;
Line(4) = { 7, 8}; Transfinite Line {4} = Na Using Progression 1.00;//Using Progression Ra;
Line(5) = { 8, 1}; Transfinite Line {5} = Na Using Progression Rr; // Bump: refina progressivamente
Line(6) = { 4, 5}; Transfinite Line {6} = Nc Using Progression Rc;
Line(7) = { 6, 5}; Transfinite Line {7} = Na Using Progression Rr;
Line(8) = { 7, 6}; Transfinite Line {8} = Nc Using Progression Rc;

// Linhas do cilindro
Circle(9)  = { 2, 11,  9}; Transfinite Line { 9} = Na Using Progression 1.00;
Circle(10) = { 9, 11, 10}; Transfinite Line {10} = Na Using Progression 1.00;//1.064;
Circle(11) = { 10, 11, 3}; Transfinite Line {11} = Na Using Progression 1.00;

//Linhas dos Blocos
Line(12) = { 8,  9}; Transfinite Line {12} = Nr Using Progression Rr;
Line(13) = { 7, 10}; Transfinite Line {13} = Nr Using Progression Rr;

// Superficies
Line Loop(14) = {  1,   9, -12,  5}; Plane Surface(15) = {14};		
Line Loop(16) = { -2,  -3, 13, 11}; Plane Surface(17) = {16};		
Line Loop(18) = { 10, -13,  4, 12}; Plane Surface(19) = {18};		
Line Loop(20) = {  6,  -7, -8,  3}; Plane Surface(21) = {20};		

Transfinite Surface {15} Alternate;
Transfinite Surface {17} Alternate;//Left; 
Transfinite Surface {19} Alternate;
Transfinite Surface {21} Alternate;//Right;

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



