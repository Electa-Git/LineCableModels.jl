@testitem "Gmsh FEM / core package remains Gmsh-independent" tags=[
    :core_only,
    :extension
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
    @test_throws MethodError import_data(:msh, "saved.msh")
    @test_throws MethodError import_data(:pos, "saved.pos")
    @test_throws MethodError export_data(:onelab)
    @test !isdefined(LineCableModels, :LineCableModelsFEM)
    @test_throws MethodError Formulation(:LineCableModelsFEM)
    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt) === nothing
end
