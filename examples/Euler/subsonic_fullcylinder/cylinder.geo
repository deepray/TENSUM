Ri=1;
Ro=20*Ri;

//ntheta=100;
//nr=100	;
//r=1.05;

Point(1) = {0,0,0};
Point(2) = {-Ri,0,0};
Point(3) = {0,Ri,0};
Point(4) = {Ri,0,0};
Point(5) = {0,-Ri,0};
Point(6) = {Ro,0,0};
Point(7) = {0,Ro,0};
Point(8) = {-Ro,0,0};
Point(9) = {0,-Ro,0};

Circle(1) = {2,1,3};
Line(2) = {2,8};
Circle(3) = {3,1,4};
Line(4) = {3,7};
Circle(5) = {4,1,5};
Line(6) = {4,6};
Circle(7) = {5,1,2};
Line(8) = {5,9};

Circle(9) = {8,1,7};
Circle(10) = {7,1,6};
Circle(11) = {6,1,9};
Circle(12) = {9,1,8};

Line Loop(1) = {1,4,-9,-2};
Line Loop(2) = {3,6,-10,-4};
Line Loop(3) = {5,8,-11,-6};
Line Loop(4) = {7,2,-12,-8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Field[1] = Attractor;
Field[1].NodesList = {1};

Field[2] = MathEval;
Field[2].F = Sprintf("(F1/8.0)^2 + %g", 1 / 10000.0);

Background Field = 2;

//Transfinite Line{1,3,5,7} = 100;
//Transfinite Line{9,10,11,12} = 10;
//Transfinite Line{2,4,6,8} = nr Using Progression 2;

//Transfinite Surface(1) = {2,3,7,8};
//Transfinite Surface(2) = {3,4,6,7};
//Transfinite Surface(3) = {4,5,9,6};
//Transfinite Surface(4) = {5,2,8,9};

Physical Surface(100000) = {1,2,3,4};

Physical Line(100001) = {1,3,5,7}; // cylinder
Physical Line(100002) = {9,10,11,12}; // farfield

