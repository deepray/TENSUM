// Transfinte mesh for naca0012
// Mesh.RecombineAll = 1;

lc = 5;

Li = 50.0; // distance of inflow boundary from origin
Lo = 50.0; // distance of outflow boundary from origin

n  = 65; // points on upper/lower surface of airfoil used to define airfoil
         // These points may not appear in the mesh.

m = 2*n-2; // total number of points on airfoil without repetition
             // LE and TE points are common to upper/lower surface

nle = n; // point number of LE = no. of points on upper surface
         // Point(1) is trailing edge

eps = 0.4;

// NACA0012 profile
// formula taken from http://turbmodels.larc.nasa.gov/naca0012_val.html
Macro NACA0012
   x2 = x * x;
   x3 = x * x2;
   x4 = x * x3;
   y = 0.594689181*(0.298222773*Sqrt(x) 
       - 0.127125232*x - 0.357907906*x2 + 0.291984971*x3 - 0.105174606*x4);
Return

// put points on upper surface of airfoil
For i In {1:n}
   theta = Pi * (i-1) / (n-1);
   x = 0.5 * (Cos(theta) + 1.0);
   Call NACA0012;
   Point(i) = {x, y, 0.0};
   xx[i] = x;
   yy[i] = y;
EndFor

// put points on lower surface of airfoil, use upper surface points and reflect
For i In {n+1:m}
   Point(i) = {xx[2*n-i], -yy[2*n-i], 0.0}; 
EndFor


Spline(1) = {1:n};
Spline(2) = {n:m, 1};

Point(1000) = {1.0, Li, 0.0,lc};
Point(1001) = {0.0, Li, 0.0,lc};
Point(1002) = {-Li, 0.0, 0.0,lc};
Point(1003) = {0.0, -Li, 0.0,lc};
Point(1004) = {1.0, -Li, 0.0,lc};

Point(1005) = {Lo, 0.0, 0.0,lc};
Point(1006) = {Lo, Li, 0.0,lc};
Point(1007) = {Lo, -Li, 0.0,lc};

Line(3) = {nle, 1002};
Line(4) = {1, 1000};
Line(5) = {1000, 1001};
Circle(6) = {1001, nle, 1002};
Circle(7) = {1002, nle, 1003};
Line(8) = {1003, 1004};
Line(9) = {1004, 1};
Compound Line(10) = {5,6};
Compound Line(11) = {7,8};
Line(12) = {1,1005};
Line(13) = {1005, 1006};
Line(14) = {1005, 1007};
Line(15) = {1006, 1000};
Line(16) = {1007, 1004};

Line Loop(1) = {-1,4,10,-3};
Plane Surface(1) = {1};

Line Loop(2) = {-2,3,11,9};
Plane Surface(2) = {2};

Line Loop(3) = {12,13,15,-4};
Plane Surface(3) = {3};

Line Loop(4) = {-9, -16, -14, -12};
Plane Surface(4) = {4};

Field[1] = Attractor;
Field[1].NodesList = {nle};


//
// LcMax -                         /------------------
//                               /
//                             /
//                           /
// LcMin -o----------------/
//        |                |       |
//     Attractor       DistMin   DistMax
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc / 500;
Field[2].LcMax = lc;
Field[2].DistMin = 0.05;
Field[2].DistMax = 0.1;

Field[3] = Attractor;
Field[3].NodesList = {1};


//
// LcMax -                         /------------------
//                               /
//                             /
//                           /
// LcMin -o----------------/
//        |                |       |
//     Attractor       DistMin   DistMax
Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = lc / 500;
Field[4].LcMax = lc;
Field[4].DistMin = 0.05;
Field[4].DistMax = 0.1;

Field[5] = Attractor;
Field[5].NNodesByEdge = 100;
Field[5].EdgesList = {1,2};


//
// LcMax -                         /------------------
//                               /
//                             /
//                           /
// LcMin -o----------------/
//        |                |       |
//     Attractor       DistMin   DistMax
Field[6] = Threshold;
Field[6].IField = 5;
Field[6].LcMin = lc / 500;
Field[6].LcMax = lc;
Field[6].DistMin = 0.1;
Field[6].DistMax = 40;

// Many other types of fields are available: see the reference manual for a
// complete list. You can also create fields directly in the graphical user
// interface by selecting Define->Fields in the Mesh module.

// Finally, let's use the minimum of all the fields as the background mesh field
Field[7] = Min;
Field[7].FieldsList = {2, 4, 6};
Background Field = 7;


// set physical labels for use in boundary condition
Physical Surface(100) = {1,2,3,4};
Physical Line(1) = {10,11,15,16}; // inflow and top/bottom
Physical Line(2) = {1};     // upper surface of airfoil
Physical Line(3) = {2};     // lower surface of airfoil
Physical Line(4) = {13,14}; // outflow
