@testitem "Gmsh FEM / core package remains Gmsh-independent" tags=[
    :core_only,
    :extension
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
    @test LineCableModels.LineCableModelsFEM <:
          LineCableModels.Engine.AbstractFormulationBackend
    @test LineCableModels.LineCableModelsFEMOptions <:
          LineCableModels.Engine.AbstractFormulationOptions

    formulation = LineCableModels.Formulation(
        :LineCableModelsFEM;
        options = (ideal_transposition = false,),
        fem_options = (ui = false,)
    )
    @test formulation isa LineCableModels.LineCableModelsFEM
    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
end
