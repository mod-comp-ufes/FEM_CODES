cl = 0.2;

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

Field[1] = Attractor;
Field[1].NodesList = {1,2,3,4};
//Field[1].NNodesByEdge = 100;
//Field[1].EdgesList = {1,2,4};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = cl/30;
Field[2].LcMax = cl;
Field[2].DistMin = 0.3;
Field[2].DistMax = 0.5;

Background Field = 2;
