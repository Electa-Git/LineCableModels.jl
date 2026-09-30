Group { Domain=Region[100]; Measurement=Region[101]; Loop=Region[102]; }
Function { exact[]=Complex[1,2]*Vector[1-Y[],X[],0]; }
Jacobian { {Name J; Case {
 {Region Region[{Measurement,Loop}]; Jacobian Sur;}
 {Region All; Jacobian Vol;}
} } }
Integration { {Name I; Case { {Type Gauss; Case {
  {GeoElement Triangle; NumberOfPoints 3;}
  {GeoElement Line; NumberOfPoints 2;}
} } } } }
FunctionSpace { {Name H; Type Form1; BasisFunction {
 {Name s; NameOfCoef a; Function BF_Edge; Support Region[{Domain,Measurement,Loop}]; Entity EdgesOf[All];}
} } }
Formulation { {Name Projection; Type FemEquation; Quantity {
 {Name e; Type Local; NameOfSpace H;}
} Equation {
 Integral { [Dof{e},{e}]; In Domain; Jacobian J; Integration I; }
 Integral { [-exact[],{e}]; In Domain; Jacobian J; Integration I; }
} } }
Resolution { {Name Solve; System { {Name A; NameOfFormulation Projection; Type ComplexValue; Frequency 1;} }
 Operation { Generate[A]; Solve[A]; SaveSolution[A]; PostOperation[Integrals]; }
} }
PostProcessing { {Name Fields; NameOfFormulation Projection; Quantity {
 {Name E; Value { Term { [{e}]; In Domain; Jacobian J; } } }
 {Name Circulation; Value { Integral { [{e}*Tangent[]]; In Region[{Measurement,Loop}]; Jacobian J; Integration I; } } }
} } }
PostOperation { {Name Integrals; NameOfPostProcessing Fields; Operation {
 Print[Circulation[Measurement], OnGlobal, Format Table, File "open.txt"];
 Print[Circulation[Loop], OnGlobal, Format Table, File "loop.txt"];
 Print[E, OnPoint{0.31,0.44,0}, Format Table, File "offmesh.txt"];
} } }
