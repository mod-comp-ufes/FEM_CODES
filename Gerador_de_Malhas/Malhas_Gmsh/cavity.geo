cl = 0.5;
nsp = 180;  //numero de sub particoes

Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {1, 1, 0, cl};
Point(4) = {0, 1, 0, cl};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};

// Put nsp points with a refinement toward the extremities on curve 2

Transfinite Line{1} = nsp Using Bump 0.05;
Transfinite Line{2} = nsp Using Bump 0.05;
Transfinite Line{3} = nsp Using Bump 0.05;
Transfinite Line{4} = nsp Using Bump 0.05;

// Define the Surface as transfinite, by specifying the four corners of the
// transfinite interpolation
Transfinite Surface{6} = {1,2,3,4} Alternate;

// (Note that the list on the right hand side refers to points, not curves. When
// the surface has only 3 or 4 points on its boundary the list can be
// omitted. The way triangles are generated can be controlled by appending
// "Left", "Right" or "Alternate" after the list.)

// Recombine the triangles into quads
//Recombine Surface{1}; 

// Apply an elliptic smoother to the grid
Mesh.Smoothing = 100;



