@testitem "Gmsh FEM / core package remains Gmsh-independent" tags=[
    :core_only,
    :extension
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
    @test_throws r"Load Gmsh" import_data(:msh, "saved.msh")
    @test_throws r"Load Gmsh" import_data(:pos, "saved.pos")
    @test_throws r"Load Gmsh" export_data(:onelab)
    @test LineCableModels.LineCableModelsFEM <:
          LineCableModels.AbstractFormulation
    @test supertype(LineCableModels.LineCableModelsFEM) === LineCableModels.AbstractFormulation

    formulation = LineCableModels.Formulation(
        :LineCableModelsFEM;
        options = (ideal_transposition = false,))
    execution = computation_options(LineCableModelsFEM, ComputationOptions())
    @test execution isa ComputationOptions
    @test formulation isa LineCableModels.LineCableModelsFEM
    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
end
