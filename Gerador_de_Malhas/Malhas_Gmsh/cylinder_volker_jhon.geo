//Gmsh projeto criado em 08:02:2018 - Problema Cilindro 2D
cl = 1.0;
/*
Na = 13; Ra = 1.00; 
Nb = 11; Rb = 1.00;
Nc = 19; Rc = 1.00;
Nd = 23; Rd = 1.00;
Nr = 21; Rr = 0.94;

Na = 89; Ra = 4.00; 
Nb = 89; Rb = 4.00;
Nc = 89; Rc = 4.00;
Nd = 105; Rd = 1.02; Rd2 = 0.98;
Nr = 75; Rr = 0.96;
*/
Na = 45; Ra = 4.00; 
Nb = 45; Rb = 4.00;
Nc = 45; Rc = 4.00;
Nd = 85; Rd = 1.02; Rd2 = 0.98;
Nr = 45; Rr = 0.96;
//  
Point(1)  = {     0,    0, 0, cl};
Point(2)  = {  0.41,    0, 0, cl};
Point(3)  = {   2.2,    0, 0, cl};
Point(4)  = {   2.2, 0.41, 0, cl};
Point(5)  = {  0.41, 0.41, 0, cl};
Point(6)  = {     0, 0.41, 0, cl};
//Pontos do cilindro
//Raio = 0.1; 
Point(7)  = { 0.16464466, 0.16464466, 0, cl};
Point(8)  = { 0.23535534, 0.16464466, 0, cl};
Point(9)  = { 0.23535534, 0.23535534, 0, cl};
Point(10) = { 0.16464466, 0.23535534, 0, cl};
Point(11) = { 0.2, 0.2, 0, cl};

Line(1) = { 1, 2}; Transfinite Line {1} = Na Using Bump Ra;//Using Progression Ra;
Line(2) = { 2, 3}; Transfinite Line {2} = Nd Using Progression Rd;
Line(3) = { 4, 5}; Transfinite Line {3} = Nd Using Progression Rd2;
Line(4) = { 5, 6}; Transfinite Line {4} = Na Using Bump Ra;//Using Progression Ra;
Line(5) = { 1, 6}; Transfinite Line {5} = Nb Using Bump Rb; // Bump: refina progressivamente
Line(6) = { 2, 5}; Transfinite Line {6} = Nc Using Bump Rc;
Line(7) = { 3, 4}; Transfinite Line {7} = Nc Using Bump Rd;

// Linhas do cilindro
Circle(8)  = { 7, 11, 8}; Transfinite Line {8} = Na Using Progression 1.00; //0.944;
Circle(9)  = { 8, 11, 9}; Transfinite Line {9} = Nc Using Progression 1.00;
Circle(10) = { 9, 11, 10}; Transfinite Line {10} = Na Using Progression 1.00;//1.064;
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

Transfinite Surface {17} Alternate;//Left; 
Transfinite Surface {19} Alternate;
Transfinite Surface {21} Alternate;//Right;
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



