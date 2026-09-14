@testitem "Gmsh FEM / core package remains Gmsh-independent" tags=[
    :core_only,
    :extension
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
    @test LineCableModels.LineCableModelsFEM <:
          LineCableModels.AbstractFormulation
    @test supertype(LineCableModels.LineCableModelsFEM) === LineCableModels.AbstractFormulation

    formulation = LineCableModels.Formulation(
        :LineCableModelsFEM;
        options = (ideal_transposition = false,))
    execution = computation_options(LineCableModelsFEM, (ui=false,))
    @test execution isa ComputationOptions
    @test !execution.ui
    @test formulation isa LineCableModels.LineCableModelsFEM
    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
end
