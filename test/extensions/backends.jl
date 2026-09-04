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

    copper=Material(kind = :conductor, rho = 1.72e-8, eps_r = 1.0, mu_r = 1.0)
    insulation_material=Material(
        kind = :insulator, rho = Inf, eps_r = 2.3, mu_r = 1.0
    )
    perfect_insulator=extension._base_material_color(Inf)
    @test perfect_insulator.r > 0.8
    @test perfect_insulator.g > 0.8
    @test perfect_insulator.b > 0.8
    @test extension._base_material_color(1e10) != perfect_insulator
    strand_radius=0.5e-3
    core_radius=sqrt(7)*strand_radius
    core=stranded(
        copper;
        shape = Disk(strand_radius),
        layers = 1,
        n = 6,
        lay = nothing,
        compact = FillFactor(1),
        boundary = Disk(core_radius)
    )
    design=build(
        CableDesign,
        "outlined-compacted-strands",
        terminal(:phase, core),
        insulation(insulation_material; t = 1e-3)
    )
    polygons=extension._native_design_shapes(
        design, 0.0, 0.0; display_legend = false
    )
    @test length(polygons) == 8
    @test all(polygon -> polygon.stroke == (:black, 0.35), polygons[1:7])
    @test all(polygon -> polygon.width == 0.6, polygons[1:7])
    @test polygons[8].stroke === :transparent
    @test polygons[8].width == 0.0

    @test_throws ArgumentError Makie.plot(
        parameters; backend = :gl, display_plot = false
    )
end
