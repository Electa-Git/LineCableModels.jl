@testitem "Makie addons / construction-time panel legends preserve source ownership" tags=[:visual] begin
    using CairoMakie

    frequency = [1.0, 10.0, 100.0]
    impedance = reshape(ComplexF64[1 + 2im, 2 + 4im, 3 + 6im], 1, 1, :)
    admittance = impedance .* 1e-6
    reference = LineParameters(copy(impedance), copy(admittance), frequency)
    candidate = LineParameters(2impedance, 2admittance, frequency)
    sources = (; reference, candidate)
    requests = ((R, 1, 1, :), (X, 1, 1, :))
    options = (; backend=:cairo, display_plot=false, controls=false,
        layout=(1, 2), panel_titles=("Resistance", "Reactance"),
        series_labels=("baseline", "alternative"), legend_position=nothing,
        legend_overflow=:show_all, length_unit=:base, quantity_units=:base,
        freq_unit=:base, clip=false)

    plot = Makie.plot(sources, requests; options...,
        panel_legends=(1, 1) => :right)
    @test Set(keys(plot.panel_legends)) == Set(((1, 1),))
    @test plot.panel_legends[(1, 1)] isa Makie.Legend
    @test plot.addon_state.panel_legend_positions[(1, 1)] == :right
    @test plot.panel_legends[(1, 1)].orientation[] == :vertical

    local_labels = ("local baseline", "local alternative")
    inside = Makie.plot(sources, requests; options...,
        panel_legends=(
            (1, 1) => (position=:inside, anchor=(:left, :bottom),
                title="Local sources", legend_labels=local_labels,
                overflow=:show_all, labelsize=13),
            (1, 2) => false,
        ))
    legend = inside.panel_legends[(1, 1)]
    @test Set(keys(inside.panel_legends)) == Set(((1, 1),))
    @test legend.halign[] == :left
    @test legend.valign[] == :bottom
    @test legend.labelsize[] == 13
    @test first(only(legend.entrygroups[])) == "Local sources"
    @test Set(values(inside.addon_state.panel_data[(1, 1)].labels)) == Set(local_labels)
    @test [entry.label[] for entry in last(only(legend.entrygroups[]))] ==
        collect(local_labels)
    @test Set(values(inside.addon_state.labels)) == Set(("baseline", "alternative"))
    @test Set(values(inside.addon_state.panel_data[(1, 2)].labels)) ==
        Set(("baseline", "alternative"))

    moved = panellegend!(inside, (1, 1); position=:bottom, overflow=:show_all)
    @test moved === inside.panel_legends[(1, 1)]
    @test moved !== legend
    @test moved.orientation[] == :horizontal
    @test first(only(moved.entrygroups[])) == "Local sources"
    @test inside.addon_state.panel_legend_positions[(1, 1)] == :bottom
    @test panellegend!(inside, (1, 1); position=nothing) === nothing
    @test !haskey(inside.panel_legends, (1, 1))
    @test !haskey(inside.addon_state.panel_legend_positions, (1, 1))

    configured = Makie.plot(sources, requests; options...,
        panel_legends=Dict(
            (1, 1) => nothing,
            (1, 2) => (position=:top, overflow=:show_all,
                legend_labels=Dict("baseline"=>"renamed"), title="Reactance sources"),
        ))
    @test Set(keys(configured.panel_legends)) == Set(((1, 2),))
    @test configured.panel_legends[(1, 2)].orientation[] == :horizontal
    @test first(only(configured.panel_legends[(1, 2)].entrygroups[])) ==
        "Reactance sources"
    @test Set(values(configured.addon_state.panel_data[(1, 2)].labels)) ==
        Set(("renamed", "alternative"))
    @test [entry.label[] for entry in last(only(
        configured.panel_legends[(1, 2)].entrygroups[]))] == ["renamed", "alternative"]
    @test Set(values(configured.addon_state.labels)) == Set(("baseline", "alternative"))
    panellegend!(configured, (1, 2); position=:left, overflow=:show_all)
    @test configured.addon_state.panel_legend_positions[(1, 2)] == :left
    @test configured.panel_legends[(1, 2)].orientation[] == :vertical

    for rendered in (plot, inside, configured)
        Makie.colorbuffer(rendered.figure)
        @test [axis.title[] for axis in rendered.axes] == ["Resistance", "Reactance"]
        for (index, axis) in enumerate(rendered.axes)
            curves = filter(item -> item isa Makie.Lines, axis.scene.plots)
            @test length(curves) == 2
            expected = index == 1 ? real.(vec(impedance)) : imag.(vec(impedance))
            @test last.(curves[1][1][]) ≈ expected
            @test last.(curves[2][1][]) ≈ 2expected
            @test first.(curves[1][1][]) ≈ frequency
        end
    end
    @test Z(reference) == impedance
    @test Y(reference) == admittance
    @test Z(candidate) == 2impedance
    @test Y(candidate) == 2admittance

    for requested in (42, ((0, 1) => :right), ((true, 1) => :right),
            ((1, 1) => 42), ((1, 1) => (position=:inside, anchor=(:left,))))
        @test_throws ArgumentError Makie.plot(sources, requests; options...,
            panel_legends=requested)
    end
    @test_throws BoundsError Makie.plot(sources, requests; options...,
        panel_legends=(2, 1) => :right)
end
