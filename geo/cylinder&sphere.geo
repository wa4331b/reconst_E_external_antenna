// parameters
lc = 0.02;
lc2 = 0.02;
lx = 0.01;
lz = 0.08;
x0 = 0.0;

// and now the geometry
Point(1) = {0, 0, 0, lc};
Point(2) = {-lx, 0, 0, lc};
Point(3) = {0, -lx, 0, lc};
Point(4) = {lx, 0, 0, lc};
Point(5) = {0, lx, 0, lc};
Rotate {{1, 0, 0}, {0, 0, 0}, Pi/2} {
  Point{1, 2, 3, 4, 5};
}
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Extrude {0,lz/2.0,0} {
  Line{2,1,4,3};
}
Extrude {0,-lz/2.0,0} {
  Line{1,4,3,2};
}
Line Loop(37) = {9,5,17,13};
Plane Surface(38) = {37};
Line Loop(39) = {21,33,29,25};
Plane Surface(40) = {39};

Translate{x0, 0, 0} {
	Line{1, 2, 3, 4};
	Surface{38, 40};
}

lc = 0.01;
lx = 0.05;
x_center = 0.02 + lx;
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
Characteristic Length {102} = 0.02;
