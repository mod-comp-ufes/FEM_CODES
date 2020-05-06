cl = 0.05;
Point(1) = {0, 0.75, 0, cl};
Point(2) = {7.5, 0.75, 0, cl};
Point(3) = {7.5, 1.5, 0, cl};
Point(4) = {0, 1.5, 0, cl};
Point(5) = {7.5, 0, 0, cl};
Point(6) = {30, 0, 0, cl};
Point(7) = {30, 1.5, 0, cl};
//
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 3};
//
Line Loop(15) = {1, 2, 3, 4};
Plane Surface(16) = {15};
Line Loop(17) = {5, 6, 7, 8, -2};
Plane Surface(18) = {17};
//
// Put n points total on combination of curves
//Transfinite Line{1,3} = 30;
//Transfinite Line{2,4,5,7,8,11,14} = 6;
//Transfinite Line{6,9} = 18;
//
// Put nsp points with a refinement toward the extremities on curve 2
//Transfinite Line{10,-12,13} = 70 Using Progression 1.00;
//
// Define the Surface as transfinite, by specifying the four corners of the
// transfinite interpolation
//
//Transfinite Surface{20} Left; 
//
//Transfinite Surface{16} Left;
//Transfinite Surface{22} Left;
//Transfinite Surface{18} = {5,6,11,3} Left;
//
// (Note that the list on the right hand side refers to points, not curves. When
// the surface has only 3 or 4 points on its boundary the list can be
// omitted. The way triangles are generated can be controlled by appending
// "Left", "Right" or "Alternate" after the list.)
//
// Recombine the triangles into quads
//Recombine Surface{30}; 
//
// Apply an elliptic smoother to the grid
//Mesh.Smoothing = 100;

//Physical Surface(1) = 1;

