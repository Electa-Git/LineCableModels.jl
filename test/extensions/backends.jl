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
    silver=extension._conductor_color(1.72e-8)
    bronze=extension._conductor_color(2.5e-7)
    air=extension._dielectric_color(1.0)
    amber=extension._dielectric_color(10.0)
    @test silver.b > bronze.b
    @test bronze.r > bronze.b
    @test air.r > 0.95
    @test air.g > 0.95
    @test air.b > 0.95
    @test amber.r > amber.g > amber.b

    conductive_eps=Material(
        kind = :conductor, rho = copper.rho, eps_r = 8.0, mu_r = 1.0
    )
    lossy_insulation=Material(
        kind = :insulator, rho = 1.0e6, eps_r = 2.3, mu_r = 1.0
    )
    magnetic_copper=Material(
        kind = :conductor, rho = copper.rho, eps_r = 1.0, mu_r = 300.0
    )
    earth=(rho = 100.0, eps_r = 4.0, mu_r = 1.0)
    different_earth_eps=(rho = 100.0, eps_r = 20.0, mu_r = 1.0)
    @test extension._material_color(copper) ==
          extension._material_color(conductive_eps)
    @test extension._material_color(insulation_material) ==
          extension._material_color(lossy_insulation)
    @test extension._material_color(copper) !=
          extension._material_color(magnetic_copper)
    @test extension._material_color(earth) ==
          extension._material_color(different_earth_eps)
    strand_radius=0.5e-3
    core_radius=sqrt(7)*strand_radius
    core=stranded(
        copper;
        shape = Disk(strand_radius),
        lay = nothing,
        compact = true,
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
    inner_boundary=extension.RGB(102 / 255, 109 / 255, 118 / 255)
    @test all(polygon -> polygon.stroke == inner_boundary, polygons[1:7])
    @test all(polygon -> polygon.width == 0.5, polygons[1:7])
    @test polygons[8].stroke == inner_boundary
    @test polygons[8].width == 0.45

    enclosed_design=build(
        CableDesign,
        "outlined-enclosure",
        pipe(
            terminal(:phase, LineCableModels.core(copper; r = 0.5e-3));
            shape = Disk(2.0e-3),
            fill = insulation_material,
            wall = insulation(insulation_material; t = 0.5e-3)
        )
    )
    enclosed_polygons=extension._native_design_shapes(
        enclosed_design, 0.0, 0.0; display_legend = false
    )
    @test length(enclosed_polygons) == 3
    @test all(polygon -> polygon.stroke === :black, enclosed_polygons[2:3])
    @test all(polygon -> polygon.width == 0.8, enclosed_polygons[2:3])

    @test_throws ArgumentError Makie.plot(
        parameters; backend = :gl, display_plot = false
    )
end
