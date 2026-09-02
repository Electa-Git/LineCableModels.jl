@testitem "Extensions / Makie boundary / unavailable graphics fail explicitly" tags=[
    :extension,
    :core_only
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsMakieExt) === nothing
    @test Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt) === nothing
    @test Base.get_extension(LineCableModels, :LineCableModelsGLMakieExt) === nothing
    @test Base.get_extension(LineCableModels, :LineCableModelsWGLMakieExt) === nothing
    @test isdefined(LineCableModels, :PlotBuilder)
    @test !isdefined(LineCableModels, :set_backend!)
    @test_throws ArgumentError LineCableModels.preview(nothing)
    @test_throws ArgumentError LineCableModels.show_material_scale()
    @test_throws ArgumentError LineCableModels.plot(nothing)
    @test_throws MethodError LineCableModels.materialcolors(:rho)
    @test_throws MethodError LineCableModels.export_svg(nothing)
    @test_throws MethodError LineCableModels.plotwindow(identity; title = "headless")
end
