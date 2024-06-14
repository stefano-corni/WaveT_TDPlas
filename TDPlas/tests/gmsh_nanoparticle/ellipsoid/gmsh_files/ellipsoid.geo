//+
Point(1) = {-10, 0, 0, 1.0};
//+
Point(2) = {+10, 0, 0, 1.0};
//+
Point(3) = {0, +5, 0, 1.0};
//+
Point(4) = {0, -5, 0, 1.0};
//+
Point(5) = {0, 0, +5, 1.0};
//+
Point(6) = {0, 0, 0, 1.0};
//+
Point(7) = {0, 0, -5, 1.0};
//+
Circle(1) = {5, 6, 4};
//+
Circle(2) = {4, 6, 7};
//+
Circle(3) = {7, 6, 3};
//+
Circle(4) = {3, 6, 5};
//+
Ellipse(5) = {2, 6, 1, 4};
//+
Ellipse(6) = {4, 6, 1, 1};
//+
Ellipse(7) = {1, 6, 2, 3};
//+
Ellipse(8) = {2, 6, 1, 3};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Line Loop(2) = {2, 3, -7, -6};
//+
Ellipse(9) = {5, 6, 1, 1};
//+
Ellipse(10) = {5, 6, 2, 2};
//+
Ellipse(11) = {7, 6, 1, 1};
//+
Ellipse(12) = {7, 6, 2, 2};
//+
Line Loop(3) = {11, -6, -5, -12};
//+
Line Loop(4) = {8, -7, -9, 10};
//+
Line Loop(5) = {1, 6, -9};
//+
Surface(1) = {5};
//+
Line Loop(6) = {10, 5, -1};
//+
Surface(2) = {6};
//+
Line Loop(7) = {6, -11, -2};
//+
Surface(3) = {7};
//+
Line Loop(8) = {12, 5, 2};
//+
Surface(4) = {8};
//+
Line Loop(9) = {9, 7, 4};
//+
Surface(5) = {9};
//+
Line Loop(10) = {8, 4, 10};
//+
Surface(6) = {10};
//+
Line Loop(11) = {7, -3, 11};
//+
Surface(7) = {11};
//+
Line Loop(12) = {3, -8, -12};
//+
Surface(8) = {12};
//+
Physical Surface("ellipsoid") = {6, 5, 7, 1, 3, 4, 8, 2};
