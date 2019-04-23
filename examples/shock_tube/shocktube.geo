L = 1;// channel length
H = 0.4; // height

n1 = 101;  // normal to BL
n2 = 81;  // along plate

Point(1) = {0,0,0};
Point(2) = {L,0,0};
Point(3) = {L,H,0};
Point(4) = {0,H,0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1,2,3,4};
Ruled Surface(1) = {1};
Transfinite Surface(1) = {1,2,3,4};

Transfinite Line{1,-3} = n1;//n1 Using Bump .05;
Transfinite Line{2,-4} = n2;

Physical Surface(100000) = {1};

Physical Line(100001) = {4}; // inlet
Physical Line(100002) = {2}; // outlet
Physical Line(100003) = {1}; // bottom
Physical Line(100004) = {3}; // top

Periodic Line{1} = {-3};
