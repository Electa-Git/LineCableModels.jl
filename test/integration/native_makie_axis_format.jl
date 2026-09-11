@testitem "Makie addons / axis multipliers follow the displayed range" tags=[:visual] begin
    using CairoMakie

    frequency = [1.0, 10.0, 100.0]
    impedance = fill(1.0 + 2.0im, 1, 1, 3)
    susceptance = reshape([-1.0, 0.0, 1.0] .* 1e-14, 1, 1, 3)
    parameters = LineParameters(impedance, complex.(zero.(susceptance), susceptance), frequency)
    plot = Makie.plot(parameters, B; backend=:cairo, display_plot=false,
        controls=false, length_unit=:base, quantity_units=:base, clip=false)
    axis = only(plot.axes)
    Makie.colorbuffer(plot.figure)
    # Tiny values are not constant merely because their SI magnitude is small.
    @test axis.finallimits[].widths[2] < 1e-12
    @test axis.ytickformat[]([-1e-14, 0.0, 1e-14]) == ["-1", "0", "1"]

    ylims!(axis, -2e-8, 2e-8)
    Makie.colorbuffer(plot.figure)
    @test axis.ytickformat[]([-1e-8, 0.0, 1e-8]) == ["-1", "0", "1"]
    @test occursin("−8", repr(axis.ylabel[]))
    @test !occursin("−14", repr(axis.ylabel[]))
    axis.ylabel[] = "Caller response [S/m]"
    @test occursin("Caller response", repr(axis.ylabel[]))
    @test occursin("−8", repr(axis.ylabel[]))
    resize!(plot.figure, 1000, 600)
    Makie.colorbuffer(plot.figure)
    label = repr(axis.ylabel[])
    ticks = axis.ytickformat[]([-1e-8, 0.0, 1e-8])
    mktempdir() do directory
        path = export_svg(plot; path=joinpath(directory, "axis.svg"), open_file=false)
        @test isfile(path)
        @test repr(axis.ylabel[]) == label
        @test axis.ytickformat[]([-1e-8, 0.0, 1e-8]) == ticks
    end
    @test Y(parameters) == complex.(zero.(susceptance), susceptance)
    publication = observables(parameters, ((@observe B[1, 1, :]),);
        length_unit=:base, quantity_units=:base, clip=false)
    published_plot = Makie.plot(publication; backend=:cairo, display_plot=false,
        controls=false)
    published_axis = only(published_plot.axes)
    ylims!(published_axis, -2e-9, 2e-9)
    Makie.colorbuffer(published_plot.figure)
    @test published_axis.ytickformat[]([-1e-9, 0.0, 1e-9]) == ["-1", "0", "1"]
    @test occursin("−9", repr(published_axis.ylabel[]))
end

@testitem "Makie addons / plain mantissas and relative constant limits" tags=[:visual] begin
    using CairoMakie
    using Measurements: measurement
    extension = Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    for exponent in (-300, -14, -3, 0, 8, 300)
        format = extension._addon_linear_tickformat(exponent)
        @test format([-1.0, 0.0, 1.0] .* 10.0^exponent) == ["-1", "0", "1"]
        @test format([-0.0, 0.25, 0.5] .* 10.0^exponent) == ["0", "0.25", "0.5"]
    end
    @test extension._addon_linear_tickformat(-324)([0.0, nextfloat(0.0)]) ==
          ["0", "4.941"]
    format = extension._addon_linear_tickformat(8)
    ticks = [1e8, 1e8 + 1, 1e8 + 2]
    @test length(unique(format(ticks))) == 3
    @test all(!occursin(r"[eE]", label) for label in format(ticks))

    for (value, bounds) in (
        (0.0, (-1.0, 1.0)),
        (1e-14, (0.95e-14, 1.05e-14)),
        (-1e-14, (-1.05e-14, -0.95e-14)),
        (measurement(0.0, 1e-14), (-2e-14, 2e-14)),
        (measurement(1e-14, 2e-15), (6e-15, 1.4e-14)),
    )
        resistance = fill(value, 1, 1, 3)
        parameters = LineParameters(complex.(resistance, one.(resistance)),
            fill(1.0im, 1, 1, 3), [1.0, 10.0, 100.0])
        plot = Makie.plot(parameters, R; backend=:cairo, display_plot=false,
            controls=false, length_unit=:base, quantity_units=:base, clip=false)
        axis = only(plot.axes)
        Makie.colorbuffer(plot.figure)
        @test axis.finallimits[].origin[2] ≈ bounds[1]
        @test sum((axis.finallimits[].origin[2], axis.finallimits[].widths[2])) ≈ bounds[2]
    end
end

@testitem "Makie addons / scale controls and visibility keep multipliers synchronized" tags=[:visual] begin
    using CairoMakie
    frequency = [1.0, 10.0, 100.0]
    small = LineParameters(fill(1e-14 + 1.0im, 1, 1, 3), fill(1.0im, 1, 1, 3), frequency)
    large = LineParameters(fill(1e8 + 1.0im, 1, 1, 3), fill(1.0im, 1, 1, 3), frequency)
    plot = Makie.plot((; small, large), R; backend=:cairo, display_plot=false,
        controls=true, length_unit=:base, quantity_units=:base, clip=false)
    axis = only(plot.axes)
    Makie.colorbuffer(plot.figure)
    @test axis.ytickformat[]([0.0, 1e8]) == ["0", "1"]
    for item in filter(item -> item isa Makie.Lines, axis.scene.plots)
        last(first(item[1][])) > 1.0 && (item.visible[] = false)
    end
    Makie.colorbuffer(plot.figure)
    @test axis.finallimits[].widths[2] < 1e-12
    @test occursin("−14", repr(axis.ylabel[]))
    for _ in 1:2
        plot.controls[:ylog].active[] = true
        Makie.colorbuffer(plot.figure)
        @test axis.yscale[] === log10
        @test !occursin("× 10", repr(axis.ylabel[]))
        @test axis.ytickformat[] === Makie.automatic
        plot.controls[:ylog].active[] = false
        Makie.colorbuffer(plot.figure)
        @test axis.yscale[] === identity
        @test axis.ytickformat[]([1e-14]) == ["1"]
        @test length(findall("× 10", repr(axis.ylabel[]))) == 1
    end
end

@testitem "Makie addons / scientific formatting respects native customization" tags=[:visual] begin
    using CairoMakie
    custom = values -> fill("custom", length(values))
    plot = LineCableModels.plotwindow(; title="Native customization", backend=:cairo,
        display_plot=false, controls=false) do layout
        for (column, attributes) in enumerate((
            (; ytickformat=custom),
            (; yticks=([0.0, 1e8], ["low", "high"])),
            (;),
        ))
            axis = Axis(layout[1, column]; ylabel="Response", attributes...)
            lines!(axis, [1.0, 10.0], [1e7, 1e8])
        end
    end
    first_axis, labelled_axis, axis = plot.axes
    Makie.colorbuffer(plot.figure)
    @test first_axis.ytickformat[] === custom
    @test first_axis.ylabel[] == "Response"
    @test labelled_axis.yticks[] == ([0.0, 1e8], ["low", "high"])
    @test labelled_axis.ylabel[] == "Response"
    @test occursin("× 10", repr(axis.ylabel[]))
    axis.ytickformat[] = custom
    @test axis.ylabel[] == "Response"
    ylims!(axis, 1e6, 1e9)
    Makie.colorbuffer(plot.figure)
    @test axis.ytickformat[] === custom
    axis.ytickformat[] = Makie.automatic
    @test axis.ytickformat[]([0.0, 1e9]) == ["0", "1"]
    @test occursin("× 10", repr(axis.ylabel[]))
    for scale in (log10, sqrt)
        axis.yscale[] = scale
        Makie.colorbuffer(plot.figure)
        @test axis.ytickformat[] === Makie.automatic
        @test axis.ylabel[] == "Response"
    end
    axis.yscale[] = identity
    Makie.colorbuffer(plot.figure)
    @test axis.ytickformat[]([0.0, 1e9]) == ["0", "1"]
    @test length(findall("× 10", repr(axis.ylabel[]))) == 1
end

@testitem "Makie addons / previews and statistical plots share live multipliers" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    copper = Material(kind=:conductor, rho=1.7241e-8)
    design = @cable "axis-preview" begin
        @terminal :core begin
            core(copper; r=1e-3)
        end
    end
    result = TestFixtures.cable_monte_carlo_result()
    plots = (
        preview(design; backend=:cairo, display_plot=false, controls=false),
        Makie.hist(result, R; backend=:cairo, display_plot=false, controls=false),
        Makie.ecdfplot(result, R; backend=:cairo, display_plot=false, controls=false),
    )
    for plot in plots
        axis = only(plot.axes)
        xlims!(axis, -2e-12, 2e-12)
        ylims!(axis, -2e-12, 2e-12)
        Makie.colorbuffer(plot.figure)
        @test axis.xtickformat[]([-1e-12, 0.0, 1e-12]) == ["-1", "0", "1"]
        @test axis.ytickformat[]([-1e-12, 0.0, 1e-12]) == ["-1", "0", "1"]
        @test occursin("−12", repr(axis.xlabel[]))
        @test occursin("−12", repr(axis.ylabel[]))
    end
end

@testitem "Makie addons / native plot windows share axis formatting" tags=[:visual] begin
    using CairoMakie

    plot = LineCableModels.plotwindow(; title="Axis formatting", backend=:cairo,
        display_plot=false, controls=false) do layout
        axis = Axis(layout[1, 1]; xlabel="Time [s]", ylabel="Response")
        lines!(axis, [0.0, 1e-12, 2e-12], [-2e8, 0.0, 2e8])
    end
    axis = only(plot.axes)
    Makie.colorbuffer(plot.figure)
    @test axis.xtickformat[] isa Function
    @test axis.ytickformat[] isa Function
    if axis.xtickformat[] isa Function && axis.ytickformat[] isa Function
        @test axis.xtickformat[]([0.0, 1e-12, 2e-12]) == ["0", "1", "2"]
        @test axis.ytickformat[]([-2e8, 0.0, 2e8]) == ["-2", "0", "2"]
    end
end
