n1=10;

Point(1) = {0, -1, 0, 0};
Point(2) = {10, -1, 0, 0};
Point(3) = {0, 0, 0, 0};
Point(4) = {-6, 0, 0, 0};
Point(5) = {-6, 2, 0, 0};
Point(6) = {10, 2, 0, 0};
Point(7) = {0, 2, 0, 0};
Point(8) = {-0, 1, 0, 0};
Point(9) = {-6, 1, 0, 0};
Line(1) = {2, 1};
Line(2) = {1, 3};
Line(3) = {3, 4};
Line(4) = {4, 9};
Line(5) = {9, 5};
Line(6) = {5, 7};
Line(7) = {7, 6};
Line(8) = {3, 8};
Line(9) = {8, 7};
Line(10) = {8, 9};
Point(10) = {10, 1, -0, 0};
Point(11) = {10, 0, 0, 0};
Line(11) = {2, 11};
Line(12) = {11, 10};
Line(13) = {10, 6};
Line(14) = {11, 3};
Line(15) = {10, 8};

Line Loop(27) = {7, -13, 15, 9};
Plane Surface(28) = {27};
Line Loop(29) = {6, -9, 10, 5};
Plane Surface(30) = {29};
Line Loop(31) = {10, -4, -3, 8};
Plane Surface(32) = {31};
Line Loop(33) = {12, 15, -8, -14};
Plane Surface(34) = {33};
Line Loop(35) = {14, -2, -1, 11};
Plane Surface(36) = {35};
// Making Mesh
Transfinite Line {5,9}=n1 Using Progression 0.8;
Transfinite Line {6,10}=n1 Using Progression 1;
Transfinite Line {9,13}=n1 Using Progression 0.8;
Transfinite Line {7,15}=n1 Using Progression 1;
Transfinite Line {4,8}=n1 Using Progression 1.2;
Transfinite Line {10,3}=n1 Using Progression 1;
Transfinite Line {12}=n1 Using Progression 1;
Transfinite Line {15,14}=n1 Using Progression 1;
Transfinite Line {2,11}=n1 Using Progression 1.2;
Transfinite Line {14,1}=n1 Using Progression 1.0;
Transfinite Surface "*";
Recombine Surface "*";

// Ruller to measure y+ for the given Reynolds number








