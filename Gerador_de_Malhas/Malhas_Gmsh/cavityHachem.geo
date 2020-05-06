lr1 = 0.008264463;
lr2 = 0.011494253;
Point(1) = {0, 0, 0, lr1};
Point(2) = {0.5, 0, 0, lr1};
Point(3) = {1, 0, 0, lr1};
Point(4) = {1, 0.5, 0, lr1};
Point(5) = {1, 1, 0, lr1};
Point(6) = {0.5, 1, 0, lr1};
Point(7) = {0, 1, 0, lr1};
Point(8) = {0, 0.5, 0, lr1};
//
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
//
Line Loop(5) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(6) = {5};
//
i = 0.055;
Point(9) = {i, i, 0, lr2};
Point(10) = {0.5, i, 0, lr2};
Point(11) = {1.0-i, i, 0, lr2};
Point(12) = {1-i, 0.5, 0, lr2};
Point(13) = {1-i, 1-i, 0, lr2};
Point(14) = {0.5, 1-i, 0, lr2};
Point(15) = {i, 1-i, 0, lr2};
Point(16) = {i, 0.5, 0, lr2};
Point{ 9, 10, 11, 12, 13, 14, 15, 16} In Surface{6};

