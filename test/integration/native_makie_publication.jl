@testitem "Makie addons / matrix selections preserve numerical observations" tags=[:visual] begin
    using CairoMakie
    using LinearAlgebra: diag

    frequency = [10.0, 100.0, 1000.0, 10_000.0]
    # Deliberately distinguish both off-diagonals and every frequency plane:
    # plotting must not assume reciprocity or silently reorder the selection.
    resistance = reshape(collect(1.0:16.0), 2, 2, 4) .* 1e-4
    inductance = reshape(collect(21.0:36.0), 2, 2, 4) .* 1e-7
    conductance = reshape(collect(41.0:56.0), 2, 2, 4) .* 1e-9
    capacitance = reshape(collect(61.0:76.0), 2, 2, 4) .* 1e-11
    angular = reshape(2pi .* frequency, 1, 1, :)
    impedance = complex.(resistance, angular .* inductance)
    admittance = complex.(conductance, angular .* capacitance)
    parameters = LineParameters(copy(impedance), copy(admittance), frequency)
    options = (; backend=:cairo, display_plot=false, controls=false,
        freq_unit=:kilo, length_unit=:base, quantity_units=:base, clip=false)
    samples = [4, 2, 1]

    @testset "standalone $selector" for (source, selector, expected) in (
        (SeriesImpedance(copy(impedance)), R, resistance),
        (SeriesImpedance(copy(impedance)), L, inductance),
        (ShuntAdmittance(copy(admittance)), G, conductance),
        (ShuntAdmittance(copy(admittance)), C, capacitance),
    )
        plot = Makie.plot(source, frequency,
            (selector, [2, 1], [2, 1], samples); options...)
        @test plot isa UIPlot
        @test Set(keys(plot.addon_state.panel_data)) ==
            Set(((1, 1), (1, 2), (2, 1), (2, 2)))
        for ((row, column), panel) in plot.addon_state.panel_data
            curve = only(filter(item -> item isa Makie.Lines, panel.axis.scene.plots))
            points = curve[1][]
            @test first.(points) ≈ frequency[samples] ./ 1000
            @test last.(points) ≈ expected[row, column, samples]
        end
        @test source.values == (source isa SeriesImpedance ? impedance : admittance)
    end

    @testset "diagonal $selector" for (selector, expected) in (
        (R, resistance), (L, inductance), (G, conductance), (C, capacitance),
    )
        plot = Makie.plot(parameters, (selector, diag, [2, 1], samples);
            frequencies=copy(frequency), options...)
        @test Set(keys(plot.addon_state.panel_data)) == Set(((1, 1), (2, 2)))
        for ((row, column), panel) in plot.addon_state.panel_data
            @test row == column
            curve = only(filter(item -> item isa Makie.Lines, panel.axis.scene.plots))
            @test first.(curve[1][]) ≈ frequency[samples] ./ 1000
            @test last.(curve[1][]) ≈ expected[row, row, samples]
        end
    end

    @test_throws ArgumentError Makie.plot(parameters, (R, 1, 1, :);
        frequencies=frequency .* 2, options...)
    @test_throws DimensionMismatch Makie.plot(parameters,
        ((R, 1, 1, 1:2), (R, 2, 2, 2:3)); options...)
    @test_throws BoundsError Makie.plot(parameters, (R, 3, 1, :); options...)
    @test_throws DimensionMismatch Makie.plot(SeriesImpedance(impedance),
        frequency[1:2], R; options...)
    @test_throws ArgumentError Makie.plot(ShuntAdmittance(admittance),
        [10.0, 100.0, Inf, 10_000.0], G; options...)
    @test_throws DomainError Makie.plot(SeriesImpedance(impedance),
        [0.0, 100.0, 1000.0, 10_000.0], L; options...)
    @test_throws DomainError Makie.plot(ShuntAdmittance(admittance),
        [0.0, 100.0, 1000.0, 10_000.0], C; options...)

    residue = fill(complex(eps(Float64) / 2, 1.0), 1, 1, 4)
    source = SeriesImpedance(residue)
    for clipped in (false, true)
        plot = Makie.plot(source, frequency, R; options..., clip=clipped)
        curve = only(filter(item -> item isa Makie.Lines, only(plot.axes).scene.plots))
        @test last.(curve[1][]) == fill(clipped ? 0.0 : eps(Float64) / 2, 4)
        @test real.(source.values) == fill(eps(Float64) / 2, 1, 1, 4)
    end
    @test Z(parameters) == impedance
    @test Y(parameters) == admittance
    @test frequencies(parameters) == frequency
end

@testitem "Makie addons / SVG exports preserve live state and existing files" tags=[:visual] begin
    using CairoMakie

    copper = Material(kind=:conductor, rho=1.7241e-8)
    xlpe = Material(kind=:insulator, rho=1e14, eps_r=3.5)
    design = @cable "publication-state" begin
        @terminal :core begin
            core(copper; r=1e-3)
            insulation(xlpe; t=0.3e-3)
        end
    end
    plot = preview(design; backend=:cairo, display_plot=false,
        controls=true, open_export=false, size=(750, 550))
    axis = only(plot.axes)
    overlay = lines!(axis, [-0.001, 0.001], [0.0, 0.0]; color=:red)
    annotation = text!(axis, 0, 0; text="caller annotation")
    axis.title[] = "Caller title"
    xlims!(axis, -0.002, 0.002)
    ylims!(axis, -0.002, 0.002)
    plot.figure.scene.backgroundcolor[] = Makie.to_color(:lightgray)
    # These are the caller-owned observables that publication export temporarily
    # changes, in addition to the live content that it must leave untouched.
    observables = Any[plot.figure.scene.backgroundcolor, axis.titlefont,
        axis.xlabelfont, axis.ylabelfont, axis.xticklabelfont, axis.yticklabelfont,
        plot.legend.labelfont, plot.legend.titlefont]
    for colorbar in plot.colorbars
        colorbar.labelfont[] = :italic
        colorbar.ticklabelfont[] = :bold
        append!(observables, (colorbar.labelfont, colorbar.ticklabelfont))
    end
    for role in (:regular, :italic, :bold)
        push!(observables, plot.figure.scene.theme[:fonts][role])
    end
    append!(observables, [control.blockscene.visible for control in values(plot.controls)])
    before = map(observable -> observable[], observables)
    rows = copy(plot.figure.layout.rowsizes)
    gap = plot.figure.layout.default_rowgap
    backend = Makie.current_backend()
    limits = axis.limits[]
    geometry = copy(overlay[1][])
    mktempdir() do directory
        for theme in (:default, :publication)
            path = joinpath(directory, "$theme.svg")
            @test export_svg(plot; path, theme, open_file=false) == path
            document = read(path, String)
            @test occursin("<svg", document)
            @test occursin("</svg>", document)
            @test map(observable -> observable[], observables) == before
            @test plot.figure.layout.rowsizes == rows
            @test plot.figure.layout.default_rowgap == gap
            @test Makie.current_backend() === backend
            @test axis.limits[] == limits
            @test axis.title[] == "Caller title"
            @test overlay[1][] == geometry
            @test annotation in axis.scene.plots
            @test_throws ArgumentError export_svg(plot; path, open_file=false)
            @test read(path, String) == document
        end
        @test_throws ArgumentError export_svg(plot;
            path=joinpath(directory, "invalid.png"), open_file=false)
        @test_throws ArgumentError export_svg(plot;
            path=joinpath(directory, "invalid.svg"), theme=:invalid, open_file=false)
        @test !isfile(joinpath(directory, "invalid.svg"))
        @test map(observable -> observable[], observables) == before
        @test plot.figure.layout.rowsizes == rows

        # A failed native save must also restore the live layout. Unix directory
        # permissions provide a real I/O failure without replacing Makie methods.
        if Sys.isunix() && ccall(:geteuid, Cuint, ()) != 0
            destination = mkdir(joinpath(directory, "read-only"))
            chmod(destination, 0o500)
            failure = try
                export_svg(plot; path=joinpath(destination, "failure.svg"),
                    theme=:publication, open_file=false)
                nothing
            catch exception
                exception
            finally
                chmod(destination, 0o700)
            end
            @test failure isa Exception
            @test !isfile(joinpath(destination, "failure.svg"))
            @test map(observable -> observable[], observables) == before
            @test plot.figure.layout.rowsizes == rows
            @test plot.figure.layout.default_rowgap == gap
            @test Makie.current_backend() === backend
            @test axis.limits[] == limits
        end

        # Automatic paths outside the checkout stay in the caller's directory.
        # An empty PATH exercises the honest unavailable-opener response without
        # invoking a desktop application or depending on one being installed.
        cd(directory) do
            plot.export_name = "  Test / Cable # 1  "
            plot.open_export = false
            first_path = export_svg(plot)
            @test dirname(first_path) == directory
            @test startswith(basename(first_path), "test_cable_1_")
            first_contents = read(first_path, String)
            second_path = export_svg(plot)
            @test second_path != first_path
            @test read(first_path, String) == first_contents
            if Sys.islinux()
                withenv("PATH" => "") do
                    path = joinpath(directory, "no-desktop.svg")
                    @test_logs (:info, r"automatic opening was unavailable") begin
                        export_svg(plot; path, open_file=true)
                    end
                    @test isfile(path)
                end
            end
        end
    end
end
