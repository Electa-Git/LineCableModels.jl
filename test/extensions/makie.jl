@testitem "Extensions / Makie boundary / unavailable graphics fail explicitly" tags=[
    :extension,
    :core_only
] begin
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsMakieExt) === nothing
    @test Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt) === nothing
    @test Base.get_extension(LineCableModels, :LineCableModelsGLMakieExt) === nothing
    @test Base.get_extension(LineCableModels, :LineCableModelsWGLMakieExt) === nothing
    @test_throws ArgumentError LineCableModels.DataModel.preview(nothing)
    @test_throws ArgumentError LineCableModels.DataModel.show_material_scale()
    @test_throws ArgumentError LineCableModels.Engine.plot(nothing)
    @test_throws MethodError LineCableModels.PlotBuilder.export_svg(nothing)

    backends = LineCableModels.PlotBuilder
    @test backends.current_backend_symbol() === :none
    for backend in (:cairo, :gl, :wgl)
        @test !backends.backend_available(backend)
        @test !backends.backend_available(Val(backend))
        @test_throws ArgumentError backends.set_backend!(backend)
        @test_throws ArgumentError backends.set_backend!(Val(backend))
        @test backends.make_screen(backend, "headless") === nothing
        @test backends.make_screen(Val(backend), "headless") === nothing
    end
    @test_throws ArgumentError backends.backend_available(:unknown)
    @test_throws ArgumentError backends.backend_available(Val(:unknown))
    @test_throws ArgumentError backends.set_backend!(:unknown)
    @test_throws ArgumentError backends.set_backend!(Val(:unknown))
    @test_throws ArgumentError backends.ensure_backend!()
    @test_throws ArgumentError backends.set_backend!(:cairo)
    @test_throws ArgumentError backends.set_backend!(Val(:cairo))
    @test_throws ArgumentError backends.renderfig(nothing)
    @test backends.next_fignum() + 1 == backends.next_fignum()
end
