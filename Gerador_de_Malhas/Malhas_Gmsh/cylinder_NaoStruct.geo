//parametro de malha
lc = 0.00515;
lc2= 0.0025;
lc3= 0.002;
lc4= 0.01;

Point(1) = {    0,     0, 0, lc}; //{-3, -3, 0, lc};
Point(2) = {  2.2,     0, 0, lc4}; //{9, -3, 0, lc};
Point(3) = {  2.2,  0.41, 0, lc4}; //{9, 3, 0, lc};
Point(4) = {    0,  0.41, 0, lc}; //{-3, 3, 0, lc};
Point(5) = {  0.2,  0.25, 0, lc3}; //{0, 1, 0, lc};
Point(6) = {  0.2,  0.15, 0, lc3}; //{0, -1, 0, lc};
Point(7) = {  0.2,   0.2, 0, lc2}; //{0, 0, 0, lc};
Point(8) = {  0.15,   0.2, 0, lc2}; //{0, 0, 0, lc};
Point(9) = {  0.25,   0.2, 0, lc3}; //{0, 0, 0, lc};
Point(10) = { 0.35,   0.2, 0, lc3}; //{0, 0, 0, lc};
Point(11) = { 0.45,   0.2, 0, lc3}; //{0, 0, 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {5,7,8};
Circle(6) = {8,7,6};
Circle(7) = {6,7,9};
Circle(8) = {9,7,5};

Line Loop(7) = {1,2,3,4}; //retângulo

Line Loop(8) = {5,6,7,8};//buraco

Plane Surface(9) = {7,8}; 

Point{10,11} In Surface{9};

//Physical Point(1) = {1,2};

Dirichlet = 1;				//marca dos nós de fronteira
Neumann = 2;				//marca dos nós de fronteira
Physical Line(Dirichlet) = {2,5,6};	//nós de fronteira
Physical Line(Neumann) = {1,3,4}; 	//nós de fronteira
Physical Point(Dirichlet) = {1,4,5,6};	//nós de fronteira
Physical Point(Neumann) = {2,3};	//nós de fronteira
//Physical Surface("My fancy surface label") = {9}
Physical Surface("My surface") = {9} ;

//Color Grey50{ Surface{ 8 }; }
//Color Purple{ Surface{ 9 }; }
//Color Red{ Line{ 1:4 }; }
//Color Yellow{ Line{ 5:6 }; }

