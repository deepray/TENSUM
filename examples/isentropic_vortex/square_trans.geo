L = 10;// channel length
H = 10; // height

n1 = 100;  // normal to BL
n2 = 100;  // along plate

Point(1) = {-L/2,-H/2,0};
Point(2) = {L/2,-H/2,0};
Point(3) = {L/2,H/2,0};
Point(4) = {-L/2,H/2,0};
Point(5) = {-L/2,0,0};
Point(6) = {L/2,0,0};

Line(1) = {1,2};
Line(2) = {4,3};
Line(3) = {5,6};
Line(4) = {1,5};
Line(5) = {5,4};
Line(6) = {2,6};
Line(7) = {6,3};

Line Loop(1) = {1,6,-3,-4};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,6,5};

Line Loop(2) = {3,7,-2,-5};
Ruled Surface(2) = {2};
Transfinite Surface(2) = {5,6,3,4};

Transfinite Line{4,6} = n1;
Transfinite Line{5,7} = n1;
Transfinite Line{1,2,3} = n2;

Physical Surface(100000) = {1,2};

Physical Line(100001) = {4,5}; // inlet
Physical Line(100002) = {6,7}; // outlet
Physical Line(100003) = {1}; // bottom
Physical Line(100004) = {2}; // top

Periodic Line {1} = {2};
Periodic Line {4} = {6};
Periodic Line {5} = {7};
