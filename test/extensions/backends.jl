@testitem "Extensions / backend dispatch / loaded extension activation" tags=[:extension] begin
    import LineCableModels
    using TOML

    packages=(
        (:cairo, :CairoMakie, :LineCableModelsCairoMakieExt),
        (:gl, :GLMakie, :LineCableModelsGLMakieExt),
        (:wgl, :WGLMakie, :LineCableModelsWGLMakieExt)
    )
    dependencies=get(TOML.parsefile(Base.active_project()), "deps", Dict{String, Any}())
    installed=filter(entry -> haskey(dependencies, String(entry[2])), packages)
    if isempty(installed)
        @test Base.get_extension(LineCableModels, :LineCableModelsMakieExt) === nothing
    else
        for (backend, package_name, extension_name) in installed
            Core.eval(@__MODULE__, Expr(:using, Expr(:., package_name)))
            extension=Base.get_extension(LineCableModels, extension_name)
            makie_extension=Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
            @test extension !== nothing
            @test makie_extension !== nothing
            @test Base.invokelatest(extension.activate!) === backend
            @test Base.invokelatest(makie_extension.current_backend_symbol) === backend
        end
    end
    @test !isdefined(LineCableModels, :set_backend!)
end

@testitem "Extensions / backend dispatch / explicit Cairo sugar" tags=[:visual] setup=[
    TestFixtures
] begin
    using CairoMakie
    using LineCableModels

    extension=Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    cairo_extension=Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt)
    @test extension !== nothing
    @test cairo_extension !== nothing
    @test cairo_extension.activate!() === :cairo
    @test extension.current_backend_symbol() === :cairo
    @test !isdefined(LineCableModels, :set_backend!)

    parameters=TestFixtures.two_conductor_results(
        frequencies = [50.0, 100.0, 500.0]
    )
    plots=Makie.plot(parameters; backend = :cairo, display_plot = false)
    @test plots isa Vector{UIPlot}
    @test length(plots) == 4
    @test all(plot -> plot.figure isa Makie.Figure, plots)
    @test_throws ArgumentError Makie.plot(
        parameters; backend = :gl, display_plot = false
    )
end
