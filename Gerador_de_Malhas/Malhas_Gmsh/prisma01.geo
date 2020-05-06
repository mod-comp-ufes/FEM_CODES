Point(7) = {0.2, 0.2, 0, 1.0};
Point(8) = {0.2, 0.1, 0, 1.0};
Point(9) = {-0, 0.3, 0, 1.0};
Point(10) = {0.25, 0.2, 0, 1.0};
Point(11) = {0.3, 0.1, 0, 1.0};
Line(10) = {8, 11};
Line(11) = {11, 10};
Line(12) = {10, 7};
Line(13) = {7, 8};
Line Loop(14) = {13, 10, 11, 12};
Plane Surface(15) = {14};
Transfinite Line {10:13} = 10;
Transfinite Surface{15};

