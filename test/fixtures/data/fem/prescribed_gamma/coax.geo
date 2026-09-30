If(!Exists(H)) H = .08; EndIf
Point(1) = {0,0,0,H};
Point(2) = {.2,0,0,H/4}; Point(3) = {0,.2,0,H/4};
Point(4) = {-.2,0,0,H/4}; Point(5) = {0,-.2,0,H/4};
Point(6) = {2,0,0,H}; Point(7) = {0,2,0,H};
Point(8) = {-2,0,0,H}; Point(9) = {0,-2,0,H};
Circle(1) = {2,1,3}; Circle(2) = {3,1,4};
Circle(3) = {4,1,5}; Circle(4) = {5,1,2};
Circle(5) = {6,1,7}; Circle(6) = {7,1,8};
Circle(7) = {8,1,9}; Circle(8) = {9,1,6};
Curve Loop(1) = {1,2,3,4}; Curve Loop(2) = {5,6,7,8};
Plane Surface(1) = {1}; Plane Surface(2) = {2,1};
Line(9) = {9,5}; Curve {9} In Surface {2};
Physical Surface(3001) = {1}; Physical Surface(9001) = {1};
Physical Surface(1002) = {2};
Physical Curve(4001) = {1,2,3,4}; Physical Curve(2001) = {5,6,7,8};
Physical Curve(7001) = {9}; Physical Point(8001) = {9};
Mesh.MshFileVersion = 2.2;
