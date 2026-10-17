@testitem "Gmsh FEM / core package remains Gmsh-independent" tags=[
    :core_only,
    :extension,
    :engine
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
    @test LineCableModels.LineCableModelsFEM <:
          LineCableModels.AbstractFormulation
    @test supertype(LineCableModels.LineCableModelsFEM) === LineCableModels.AbstractFormulation

    formulation = LineCableModels.Formulation(
        :fem;
        options = (ideal_transposition = false,))
    execution = computation_options(LineCableModelsFEM, ComputationOptions((ui=false,)))
    @test execution isa ComputationOptions
    @test !execution.data.ui
    @test formulation isa LineCableModels.LineCableModelsFEM
    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
end

@testitem "Gmsh FEM / nested option lifting without a solve" tags=[:unit, :engine] begin
    nested = LineCableModelsFEM(; options=(reduce_bundle=Grid((true, false)),))
    complete = LineCableModelsFEM(;
        options=Grid(((reduce_bundle=true,), (reduce_bundle=false,))))
    @test NamedTuple.(collect(nested)) == NamedTuple.(collect(complete))
end
