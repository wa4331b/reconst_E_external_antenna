lc = 0.01;
lx = 0.02;
ly = 0.02;
ly_dipole = 0.001;
lz = 0.075;

xOffset = 0.075;

Point(1) = {-lx/2.0,-ly/2.0,-lz/2.0, lc};
Point(2) = { lx/2.0,-ly/2.0,-lz/2.0, lc};
Point(3) = { lx/2.0, ly/2.0,-lz/2.0, lc};
Point(4) = {-lx/2.0, ly/2.0,-lz/2.0, lc};
Point(5) = {-lx/2.0,-ly/2.0, lz/2.0, lc};
Point(6) = { lx/2.0,-ly/2.0, lz/2.0, lc};
Point(7) = { lx/2.0, ly/2.0, lz/2.0, lc};
Point(8) = {-lx/2.0, ly/2.0, lz/2.0, lc};
Translate {xOffset, 0, 0} {
  Point{6, 7, 8, 5, 2, 3, 4, 1};
}

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};
//+
Line Loop(14) = {1, 2, 3, 4};
//+
Plane Surface(15) = {14};
//+
Line Loop(16) = {5, 6, 7, 8};
//+
Plane Surface(17) = {16};
//+
Line Loop(18) = {9, 5, -10, -1};
//+
Plane Surface(19) = {18};
//+
Line Loop(21) = {11, 7, -12, -3};
//+
Plane Surface(22) = {21};
//+
Line Loop(23) = {12, 8, -9, -4};
//+
Plane Surface(24) = {23};
//+
Line Loop(25) = {10, 6, -11, -2};
//+
Plane Surface(26) = {25};

lc = 0.01;
lx = 0.05;
x_center = 0;
y_center = 0;
z_center = 0;

p = 14;

Point(101) = {x_center, y_center, z_center, lc};
Point(102) = {x_center - lx, y_center, z_center, lc};
Point(103) = {x_center,y_center-lx,z_center,lc};
Point(104) = {x_center-lx*Sqrt(2)/2,y_center-lx*Sqrt(2)/2,z_center,lc};
Point(105) = {x_center+lx*Sqrt(2)/2,y_center-lx * Sqrt(2)/2,z_center,lc};
Point(106) = {x_center+lx,y_center,z_center,lc};

// the basis circles
Circle(101) = {102, 101, 104};
Circle(102) = {104, 101, 103};
Circle(103) = {103, 101, 105};
Circle(104) = {105, 101, 106};


Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{101};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{105};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{108};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{111};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{114};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{117};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{120};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{123};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{102};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{129};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{133};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{137};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{141};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{145};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{149};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{153};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{103};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{161};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{165};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{169};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{173};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{177};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{181};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{185};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{104};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{193};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{196};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{199};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{202};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{205};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{208};
}
Extrude {{1,0,0}, {0,0,0}, Pi/4} {
  Line{211};
}
//+
Characteristic Length {106} = 0.005;
