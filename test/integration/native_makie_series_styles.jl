@testitem "Makie addons / per-series markers and dashes reach plots and legends" tags=[:visual] begin
    using CairoMakie
    using LineCableModels

    impedance = reshape(ComplexF64[1 + 2im, 2 + 3im], 1, 1, 2)
    parameters = LineParameters(impedance, impedance .* 1e-6, [1.0, 10.0])
    sources = (; fem=parameters, proposed=parameters, xue=parameters)
    options = (backend=:cairo, display_plot=false, controls=false)
    page = LineCableModels.plot(sources, (R,); options...,
        series_labels=("FEM", "Proposed", "Xue"),
        series_attributes=((marker=:circle, markersize=8), (;), (linestyle=:dash,)))
    fem, markers = page.addon_state.groups[:result_1]
    xue = only(page.addon_state.groups[:result_3])
    @test fem isa Makie.Lines
    @test markers isa Makie.Scatter
    @test markers.marker[] == Makie.to_spritemarker(:circle)
    @test all(==(8), markers.markersize[])
    @test markers[1][] == fem[1][]
    fem.visible[] = false
    @test !markers.visible[]
    fem.visible[] = true
    @test markers.visible[]
    fem.color[] = :red
    @test markers.color[] == fem.color[]
    @test xue.linestyle[] == Makie.to_linestyle(:dash)
    entries = last(only(page.legend.entrygroups[]))
    @test [entry.label[] for entry in entries] == ["FEM", "Proposed", "Xue"]
    @test any(element -> element isa Makie.MarkerElement, first(entries).elements)
    @test !isempty(Makie.colorbuffer(page.figure))
    @test_throws ArgumentError LineCableModels.plot(sources, (R,); options...,
        series_attributes=((linestyle=:dash,),))
    @test_throws ArgumentError LineCableModels.plot(sources, (R,); options...,
        series_attributes=(unsupported_attribute=true,))
    mktempdir() do directory
        path = export_svg(page; path=joinpath(directory, "styles.svg"), open_file=false)
        @test filesize(path) > 0
        @test markers.marker[] == Makie.to_spritemarker(:circle)
        @test xue.linestyle[] == Makie.to_linestyle(:dash)
    end
end

@testitem "Makie addons / shared series attributes cover every plot family" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    using LineCableModels
    using Measurements

    options = (backend=:cairo, display_plot=false, controls=false, open_export=false)
    attributes = (color=:magenta,)
    parameters = TestFixtures.two_conductor_results()
    # A common override reaches every matrix facet on every Z/Y page.
    pages = LineCableModels.plot(parameters, (Z, Y); options...,
        series_attributes=attributes)
    @test length(pages) == 4
    @test all(page -> length(page.axes) == 4, pages)
    styled = UIPlot[pages...]
    push!(styled, Makie.plot(parameters.Z, frequencies(parameters), R;
        options..., series_attributes=attributes))
    push!(styled, Makie.plot(parameters.Y, frequencies(parameters), G;
        options..., series_attributes=attributes))

    # Table illustrations use the observation publication renderer.
    artifact = report(TableReportDefinition(((R, :, :, :),);
        illustration=true, plot_options=(; options..., series_attributes=attributes)),
        parameters)
    push!(styled, artifact.illustration)

    result = TestFixtures.cable_monte_carlo_result()
    for recipe in (Makie.hist, Makie.stairs, Makie.ecdfplot, Makie.lines, Makie.qqplot)
        push!(styled, recipe(result, R; options..., series_attributes=attributes))
    end

    copper = Material(kind=:conductor, rho=1.72e-8)
    design = build(CableDesign, "shared-style", terminal(:core, core(copper; r=0.01)))
    system = build(LineCableSystem, design, (0.0, -1.0); connections=Dict(:core=>1))
    for source in (design, [design, design], system)
        push!(styled, preview(source; options..., display_colorbars=false,
            series_attributes=attributes))
    end
    native = LineCableModels.plotwindow(; title="Shared style", options...,
        series_attributes=attributes) do canvas
        axis = Axis(canvas[1, 1])
        lines!(axis, [1.0, 2.0], [3.0, 4.0])
        scatter!(axis, [1.0, 2.0], [4.0, 3.0])
    end
    push!(styled, native)
    for page in styled
        @test !isempty(page.addon_state.groups)
        @test all(Makie.to_color(handle.color[]) == Makie.to_color(:magenta)
            for handles in values(page.addon_state.groups) for handle in handles)
    end

    # Adding markers retains measurement error bars and their visibility group.
    uncertain = complex.(measurement.(real.(Z(parameters)), 1e-6),
        measurement.(imag.(Z(parameters)), 1e-6))
    page = Makie.plot(SeriesImpedance(uncertain), frequencies(parameters), (R, 1, 1, :);
        options..., series_attributes=(marker=:circle, markersize=8, alpha=0.35))
    handles = only(values(page.addon_state.groups))
    @test count(handle -> handle isa Makie.Lines, handles) == 1
    @test count(handle -> handle isa Makie.Errorbars, handles) == 1
    @test count(handle -> handle isa Makie.Scatter, handles) == 1
    @test all(handle -> handle.alpha[] == 0.35,
        filter(handle -> handle isa Union{Makie.Lines, Makie.Scatter}, handles))
    @test !isempty(Makie.colorbuffer(page.figure))
end
