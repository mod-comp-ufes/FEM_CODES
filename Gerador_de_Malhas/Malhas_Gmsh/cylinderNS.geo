//parametro de malha
lc = 0.01;
lc2= 0.005;
//
Point(1) = {    0,    0, 0, lc};
Point(2) = {  2.2,    0, 0, lc};
Point(3) = {  2.2, 0.41, 0, lc};
Point(4) = {    0, 0.41, 0, lc};
Point(5) = {  0.2,  0.2, 0, lc2};
Point(6) = {  0.2,  0.15, 0, lc2};
Point(7) = {  0.2,  0.25, 0, lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Circle(5) = {6,5,7};
Circle(6) = {7,5,6};

Line Loop(7) = {1,2,3,4}; //retângulo

Line Loop(8) = {5,6};//buraco

Plane Surface(9) = {7,8}; 

//Physical Point(1) = {1,2};

//Dirichlet = 1;				//marca dos nós de fronteira
//Neumann = 2;				//marca dos nós de fronteira
//Physical Line(Dirichlet) = {2,5,6};	//nós de fronteira
//Physical Line(Neumann) = {1,3,4}; 	//nós de fronteira
//Physical Point(Dirichlet) = {1,4,5,6};	//nós de fronteira
//Physical Point(Neumann) = {2,3};	//nós de fronteira
//Physical Surface("My fancy surface label") = {9}
Physical Surface("My surface") = {9} ;

//Color Grey50{ Surface{ 8 }; }
//Color Purple{ Surface{ 9 }; }
//Color Red{ Line{ 1:4 }; }
//Color Yellow{ Line{ 5:6 }; }

