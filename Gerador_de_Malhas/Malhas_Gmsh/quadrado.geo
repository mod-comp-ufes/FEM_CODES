cl1 = 0.01;
cl2 = cl1*5.0;

Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};

Field[1] = Box;
Field[1].VIn = cl2;
Field[1].VOut = cl1s;
Field[1].XMin = 0.25;
Field[1].XMax = 0.75;
Field[1].YMin = 0.25;
Field[1].YMax = 0.75;

Background Field = 1;
