// Piecewise edge field with a jump at y=0.37. Both oriented measurement
// curves are shared edges of the field mesh, including at the jump.
SetFactory("Built-in");
x()={0,.3,.7,1}; y()={0,.37,1}; orientation()={1,1,-1,1};
For j In {0:2}
  For i In {0:3}
    Point(1+i+4*j)={x(i),y(j),0,.15};
  EndFor
EndFor
For j In {0:2}
  For i In {0:2}
    Line(1+i+3*j)={1+i+4*j,2+i+4*j};
  EndFor
EndFor
For j In {0:1}
  For i In {0:3}
    If(orientation(i)>0)
      Line(10+i+4*j)={1+i+4*j,5+i+4*j};
    Else
      Line(10+i+4*j)={5+i+4*j,1+i+4*j};
    EndIf
  EndFor
EndFor
For j In {0:1}
  For i In {0:2}
    Curve Loop(1+i+3*j)={1+i+3*j,orientation(i+1)*(11+i+4*j),
      -(4+i+3*j),-orientation(i)*(10+i+4*j)};
    Plane Surface(1+i+3*j)={1+i+3*j};
  EndFor
EndFor
Physical Surface("Lower",100)={1,2,3};
Physical Surface("Upper",101)={4,5,6};
Physical Curve("Forward",301)={11,15};
Physical Curve("Reverse",302)={12,16};
Mesh.MshFileVersion=2.2;
