// LineCableModelsFEM immutable one-frequency/one-basis GetDP entry point.
// Julia publishes the immutable run-local data path in the shared ONELAB
// database before launching this client.
ModelDataPath = GetString["LineCableModels/FEM/model_data_path"];
Include ModelDataPath;
eps0 = 8.8541878128e-12;
mu0 = 1.2566370614359173e-6;
UnitSource = 1.0;
GammaQuasiTEMRe = 0.0;
GammaQuasiTEMIm = 1.0e-12;
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
Include "quasi_tem.pro";
