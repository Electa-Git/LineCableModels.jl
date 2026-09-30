If(!Exists(N)) N = 12; EndIf
Point(1) = {0,0,0,1/N}; Point(2) = {1,0,0,1/N};
Point(3) = {1,1,0,1/N}; Point(4) = {0,1,0,1/N};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Curve Loop(1) = {1,2,3,4}; Plane Surface(1) = {1};
Transfinite Curve {1:4} = N+1;
Transfinite Surface {1} = {1,2,3,4};
Physical Surface(1) = {1}; Physical Curve(2) = {1,2,3,4};
Mesh.MshFileVersion = 2.2;
