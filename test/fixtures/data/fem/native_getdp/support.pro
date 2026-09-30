DefineConstant[LinePoints=4];
Group {
  Lower=Region[100]; Upper=Region[101]; Domain=Region[{Lower,Upper}];
  Forward=Region[301]; Reverse=Region[302]; Paths=Region[{Forward,Reverse}];
}
Function {
  exact[Lower]=Complex[1,2]*Vector[0,2,0];
  exact[Upper]=Complex[1,2]*Vector[0,6,0];
}
Jacobian { {Name Volume; Case { {Region All; Jacobian Vol;} } }
           {Name LineMeasure; Case { {Region All; Jacobian Sur;} } } }
Integration {
  {Name I1; Case { {Type Gauss; Case { {GeoElement Triangle; NumberOfPoints 3;} } } } }
  {Name I2; Case { {Type Gauss; Case { {GeoElement Line; NumberOfPoints LinePoints;} } } } }
}
FunctionSpace { {Name H; Type Form1; BasisFunction {
  {Name s; NameOfCoef a; Function BF_Edge;
   Support Region[{Domain,Paths}]; Entity EdgesOf[All];}
} } }
Formulation { {Name Projection; Type FemEquation; Quantity {
  {Name e; Type Local; NameOfSpace H;}
} Equation {
  Integral { [Dof{e},{e}]; In Domain; Jacobian Volume; Integration I1; }
  Integral { [-exact[],{e}]; In Domain; Jacobian Volume; Integration I1; }
} } }
Resolution { {Name Solve; System {
  {Name A; NameOfFormulation Projection; Type ComplexValue; Frequency 1;}
} Operation { Generate[A]; Solve[A]; PostOperation[Integrals]; } } }
PostProcessing { {Name Fields; NameOfFormulation Projection; Quantity {
  {Name Circulation; Value { Integral {
    [{e}*Tangent[]]; In Paths; Jacobian LineMeasure; Integration I2;
  } } }
  {Name Upwards; Value { Integral {
    [CompY[{e}]]; In Paths; Jacobian LineMeasure; Integration I2;
  } } }
} } }
PostOperation { {Name Integrals; NameOfPostProcessing Fields; Operation {
  Print[Circulation[Forward], OnGlobal, Format Table, File "forward.txt"];
  Print[Circulation[Reverse], OnGlobal, Format Table, File "reverse.txt"];
  Print[Upwards[Reverse], OnGlobal, Format Table, File "upwards.txt"];
} } }
