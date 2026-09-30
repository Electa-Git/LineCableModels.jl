If(!Exists(H)) H=.16; EndIf
Point(1)={0,0,0,H};
Point(2)={.2,0,0,H/4}; Point(3)={0,.2,0,H/4};
Point(4)={-.2,0,0,H/4}; Point(5)={0,-.2,0,H/4};
Point(6)={-3,-3,0,H}; Point(7)={0,-3,0,H}; Point(8)={3,-3,0,H};
Point(9)={3,3,0,H}; Point(10)={-3,3,0,H};
Circle(1)={2,1,3}; Circle(2)={3,1,4}; Circle(3)={4,1,5}; Circle(4)={5,1,2};
Line(5)={6,7}; Line(6)={7,8}; Line(7)={8,9}; Line(8)={9,10}; Line(9)={10,6};
Curve Loop(1)={1,2,3,4}; Curve Loop(2)={5,6,7,8,9};
Plane Surface(1)={1}; Plane Surface(2)={2,1};
Line(10)={7,5}; Curve{10} In Surface{2};
Physical Surface(3001)={1}; Physical Surface(9001)={1}; Physical Surface(1004)={2};
Physical Curve(4001)={1,2,3,4}; Physical Curve(2001)={5,6,7,8,9};
Physical Curve(7001)={10}; Physical Point(8001)={7};
Mesh.MshFileVersion=2.2;
