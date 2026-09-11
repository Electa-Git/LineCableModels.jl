// One selector serves ONELAB and the headless CLI (-setnumber Physics 0|1).
DefineConstant[
  Physics = {0, Choices{0="quasi-tem", 1="quasi-fw"},
    Name "LineCableModels/FEM/physics", Label "Physics"}
];
// Immutable inputs for one frequency and its requested terminal excitations.
If(!Exists(ModelDataPath))
  Error("Pass the immutable input path with -setstring ModelDataPath");
EndIf
Include ModelDataPath;
eps0 = 8.8541878128e-12;
mu0 = 1.2566370614359173e-6;
UnitSource = 1.0;
UnitTransverseSource = 1.0;
If(!Exists(Val_Rint))
  Val_Rint = DomainRadius;
EndIf
If(!Exists(Val_Rext))
  Val_Rext = ShellOuterRadius;
EndIf
Group {
  // Keep the two primary physical regions explicit here. DOMAIN_INF is an
  // overlapping Gmsh inventory group; selecting it directly does not preserve
  // GetDP's region dispatch for the VolSphShell Jacobian.
  AirInfJacobian = Region[{AIR_INF}];
  EarthInfJacobian = Region[{EARTH_INF}];
  DomainInf = Region[{AirInfJacobian, EarthInfJacobian}];
}
Include "jacobian_integration.pro";
If(Physics == 0)
  Include "quasi-tem.pro";
ElseIf(Physics == 1)
  Include "quasi-full.pro";
Else
  Error("Physics must be 0 (quasi-tem) or 1 (quasi-fw)");
EndIf
