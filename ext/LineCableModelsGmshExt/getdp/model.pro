// LineCableModelsFEM immutable one-frequency/one-basis GetDP entry point.
// Julia publishes the immutable run-local data path in the shared ONELAB
// database before launching this client.
ModelDataPath = GetString["LineCableModels/FEM/model_data_path"];
Include ModelDataPath;
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
