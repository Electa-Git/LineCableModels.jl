@testitem "Makie addons / preview grouping and late horizontal legend retain geometry" tags=[:visual] begin
    using CairoMakie

    copper = Material(kind=:conductor, rho=1.72e-8)
    dielectric = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    design = build(CableDesign, "preview-grouping",
        terminal(:core, core(copper; r=1e-3), insulation(dielectric; t=0.2e-3)))
    options = (; backend=:cairo, display_plot=false, controls=false, open_export=false)
    grouped = preview(design; options..., display_colorbars=false,
        legend_group=Dict(:core=>:metal, :insulation=>:dielectric),
        legend_labels=group -> uppercasefirst(String(group)))
    @test Set(values(grouped.addon_state.labels)) == Set(("Metal", "Dielectric"))
    @test isempty(grouped.colorbars)
    @test length(grouped.axes) == 1
    @test_throws ArgumentError preview(design; options..., legend_group=1)
    @test_throws ArgumentError preview(design; options..., legend_labels=1)

    collection = preview([design, design, design]; options..., layout=nothing)
    @test length(collection.axes) == 3
    @test length(collection.colorbars) == 3
    docked = preview(design; options..., display_legend=false, colorbar_position=:top)
    @test docked.legend === nothing
    geometry_plots = [copy(axis.scene.plots) for axis in docked.axes]
    figurelegend!(docked; position=:top, orientation=:horizontal, overflow=:show_all)
    Makie.colorbuffer(docked.figure)
    @test docked.legend.orientation[] == :horizontal
    bounds = docked.legend.layoutobservables.computedbbox[]
    scale_bounds = [bar.layoutobservables.computedbbox[] for bar in docked.colorbars]
    @test bounds.origin[1] + bounds.widths[1] <=
        minimum(box.origin[1] for box in scale_bounds) + 1
    figurelegend!(docked; position=:right, overflow=:show_all)
    figurelegend!(docked; position=:top, orientation=:horizontal, overflow=:show_all)
    Makie.colorbuffer(docked.figure)
    @test [axis.scene.plots for axis in docked.axes] == geometry_plots
    @test length(docked.colorbars) == 3
    @test all(isfinite, docked.legend.layoutobservables.computedbbox[].widths)

    system = build(LineCableSystem, design, (0.0, -0.1); connections=Dict(:core=>1))
    system_plot = preview(system; options..., display_colorbars=false)
    @test isempty(system_plot.colorbars)
    for zoom_factor in ("invalid", -1.0, Inf)
        @test_throws ArgumentError preview(system; options..., zoom_factor)
    end
end

@testitem "Makie addons / constant positive observations retain logarithmic limits" tags=[:visual] begin
    using CairoMakie

    frequency = [10.0, 100.0, 1000.0]
    impedance = fill(ComplexF64(1e-4, 1e-3), 1, 1, 3)
    admittance = fill(ComplexF64(0.0, 1e-8), 1, 1, 3)
    parameters = LineParameters(copy(impedance), copy(admittance), frequency)
    options = (; backend=:cairo, display_plot=false, controls=false,
        length_unit=:base, quantity_units=:base)
    plot = Makie.plot(parameters, (R, 1, 1, :); options..., yscale=:log10)
    Makie.colorbuffer(plot.figure)
    axis = only(plot.axes)
    @test axis.yscale[] === Makie.log10
    limits = axis.finallimits[]
    @test 0 < limits.origin[2] < 1e-4 < limits.origin[2] + limits.widths[2]
    @test all(isfinite, limits.widths)
    line = only(filter(item -> item isa Makie.Lines, axis.scene.plots))
    @test last.(line[1][]) ≈ fill(1e-4, 3)
    @test parameters.Z.values == impedance
    @test parameters.Y.values == admittance
end
