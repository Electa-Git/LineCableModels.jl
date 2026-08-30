// LineCableModelsFEM single production GetDP entry point.
// Julia publishes the immutable run-local data path in the shared ONELAB
// database before launching this client.
ModelDataPath = GetString["LineCableModels/FEM/model_data_path"];
Include ModelDataPath;
Group {
  DomainInf = Region[{DOMAIN_INF}];
}
Include "jacobian_integration.pro";
Include "quasi_tem.pro";
