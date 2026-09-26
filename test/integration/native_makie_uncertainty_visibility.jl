@testitem "Makie addons / coincident uncertainty intervals retain nested native glyphs" tags=[:visual] begin
    using CairoMakie, Measurements

    f = measurement.([1.0, 10.0, 100.0], [0.1, 0.2, 0.3])
    z = [measurement(i+2j+k, 0.1i+0.01j+0.001k) + im*measurement(i+j+k, 0.2)
         for i in 1:2, j in 1:2, k in eachindex(f)]
    source = LineParameters(z, z .* 1e-6, f)
    before = deepcopy((Z(source), Y(source), frequencies(source)))
    options = (backend = :cairo, display_plot = false, open_export = false,
        length_unit = :base, quantity_units = :base, freq_unit = :base, clip = false)
    pages = LineCableModels.plot(source, source, source; options..., ydata = (R,),
        layout = (1, 2), errorbar_sampling = :all,
        series_labels = ("first", "second", "third"), legend_position = :bottom)
    for page in pages
        for axis in page.axes
            curves = filter(p -> p isa Makie.Lines, axis.scene.plots)
            bars = filter(p -> p isa Makie.Errorbars, axis.scene.plots)
            @test length(curves)==3 && length(bars)==6
            @test all(isempty(p[1][]) for p in axis.scene.plots if p isa Makie.Scatter)
            @test all(line -> line.linewidth[] == 2, curves)
            @test all(line -> line[1][] == first(curves)[1][], curves)
            for direction in (:x, :y)
                intervals = filter(bar -> bar.direction[] === direction, bars)
                @test all(bar -> bar[1][] == first(intervals)[1][], intervals)
                @test length(first(intervals)[1][]) == length(f)
                # Wider caps AND stems expose earlier coincident intervals.
                # Exact styling constants are not part of the retained quantity data.
                caps = [bar.whiskerwidth[] for bar in intervals]
                stems = [bar.linewidth[] for bar in intervals]
                @test all(diff(caps) .< 0) && all(caps .> 0)
                @test all(diff(stems) .< 0) && all(stems .> 0)
            end
        end
        intervals = [p for axis in page.axes
                     for p in axis.scene.plots if p isa Makie.Errorbars]
        saved = [(copy(p[1][]), p.whiskerwidth[], p.linewidth[]) for p in intervals]
        entry = first(last(only(page.legend.entrygroups[])))
        Makie.toggle_visibility!(entry)
        Makie.toggle_visibility!(entry, true)
        for key in (:xlog, :ylog)
            page.controls[key].active[] = true
            @test !isempty(Makie.colorbuffer(page.figure))
            page.controls[key].active[] = false
        end
        resize!(page.figure, 1000, 700)
        mktempdir() do directory
            @test isfile(export_svg(page; path = joinpath(directory, "nested.svg"), open_file = false))
        end
        @test [(p[1][], p.whiskerwidth[], p.linewidth[]) for p in intervals] == saved
    end
    # Default nesting must not take ownership away from explicit native kwargs.
    custom = LineCableModels.plot(source, source; options..., ydata = ((R, 1, 1, :),),
        series_attributes = (
            (whiskerwidth = 15, linewidth = 2.5), (whiskerwidth = 7, linewidth = 0.8)))
    bars = filter(p -> p isa Makie.Errorbars, only(custom.axes).scene.plots)
    @test [bar.whiskerwidth[] for bar in bars] == [15, 15, 7, 7]
    @test [bar.linewidth[] for bar in bars] ≈ [2.5, 2.5, 0.8, 0.8]
    # Styling a caller's native plots must not reset unrelated attributes to
    # owned defaults merely because those primitives also happen to be bars.
    native = LineCableModels.plotwindow(;
        title = "Native bars", backend = :cairo, display_plot = false,
        controls = false, open_export = false, series_attributes = (color = :red,)) do canvas
        axis = Axis(canvas[1, 1])
        errorbars!(
            axis, [1.0, 2.0], [3.0, 4.0], [0.2, 0.3]; whiskerwidth = 15, linewidth = 2.5)
        errorbars!(
            axis, [1.0, 2.0], [3.0, 4.0], [0.2, 0.3]; whiskerwidth = 7, linewidth = 0.8)
    end
    bars = filter(p -> p isa Makie.Errorbars, only(native.axes).scene.plots)
    @test [bar.whiskerwidth[] for bar in bars] == [15, 7]
    @test [bar.linewidth[] for bar in bars] ≈ [2.5, 0.8]
    @test isequal(before, (Z(source), Y(source), frequencies(source)))
end

@testitem "Makie addons / staggered comparisons retain complete uncertainty bounds" tags=[:visual] begin
    using CairoMakie, Measurements

    f = collect(1.0:101.0)
    resistance = 1 .+ f ./ 100
    errors = fill(0.05, length(f))
    errors[1] = 3.0 # An endpoint interval, not a sparse interior glyph.
    z = reshape(complex.(measurement.(resistance, errors), 0.1), 1, 1, :)
    reference = LineParameters(z, z .* 1e-6, measurement.(f, 0.01))
    other_z = reshape(complex.(measurement.(1.05 .* resistance, 0.02), 0.1), 1, 1, :)
    candidate = LineParameters(other_z, other_z .* 1e-6, measurement.(f, 0.02))
    axes = (problems = [:one], formulations = [NamedTuple(Formulation())])
    result = ParametricResult(nothing, [candidate], axes, ComputationDetails((;)))
    before = deepcopy((Z(reference), Y(reference), frequencies(reference), Z(candidate)))
    options = (; backend = :cairo, display_plot = false, open_export = false,
        ydata = ((R, 1, 1, :),), reference, length_unit = :base, quantity_units = :base,
        freq_unit = :base, clip = false, fig_size = (1000, 650))
    page = LineCableModels.plot(result; options...)
    axis = only(page.axes)
    for (key, source) in zip(page.addon_state.order, (candidate, reference))
        group = page.addon_state.groups[key]
        line = only(filter(p -> p isa Makie.Lines, group))
        marker = only(filter(p -> p isa Makie.Scatter, group))
        bars = filter(p -> p isa Makie.Errorbars, group)
        @test length(line[1][]) == length(f)
        @test length(bars) == 2
        @test 0 < length(first(bars)[1][]) < length(f)
        @test first.(bars[1][1][]) == first.(bars[2][1][])
        for bar in bars
            indices = [findfirst(==(p[1]), f) for p in bar[1][]]
            expected = bar.direction[] === :x ? frequencies(source)[indices] :
                       R(source)[1, 1, indices]
            @test getindex.(bar[1][], 3) ≈ uncertainty.(expected)
            @test isdisjoint(first.(marker[1][]), first.(bar[1][]))
        end
        @test all(p -> p in line[1][], marker[1][])
        @test all(bar -> Makie.to_color(bar.color[]) == Makie.to_color(line.color[]), bars)
    end
    reference_group = page.addon_state.groups[last(page.addon_state.order)]
    reference_bar = only(filter(p -> p isa Makie.Errorbars && p.direction[] === :y, reference_group))
    @test first(f) ∉ first.(reference_bar[1][])
    @test axis.finallimits[].origin[2] <= resistance[1] - errors[1]
    @test sum((axis.finallimits[].origin[2], axis.finallimits[].widths[2])) >=
          resistance[1] + errors[1]
    # A hidden negative interval still rules out an ordinary logarithmic axis.
    page.controls[:ylog].active[] = true
    @test axis.yscale[] !== log10
    page.controls[:ylog].active[] = false
    entry = last(last(only(page.legend.entrygroups[])))
    Makie.toggle_visibility!(entry)
    @test all(!p.visible[] for p in reference_group)
    @test sum((axis.finallimits[].origin[2], axis.finallimits[].widths[2])) < 4.0
    Makie.toggle_visibility!(entry, true)
    for size in ((460, 400), (1200, 750))
        resize!(page.figure, size...)
        @test !isempty(Makie.colorbuffer(page.figure))
        @test axis.finallimits[].origin[2] <= resistance[1] - errors[1]
        @test sum((axis.finallimits[].origin[2], axis.finallimits[].widths[2])) >=
              resistance[1] + errors[1]
    end
    view = Makie.Rect2d(20.0, 1.2, 30.0, 0.5)
    axis.targetlimits[] = view
    resize!(page.figure, 900, 650)
    @test axis.targetlimits[] == view
    page.controls[:reset].clicks[] += 1
    mktempdir() do directory
        @test isfile(export_svg(page; path = joinpath(directory, "staggered.svg"), open_file = false))
    end

    for n in (1, 2)
        short = ParametricResult(nothing, [candidate[1:n]], axes, ComputationDetails((;)))
        plot = LineCableModels.plot(short; options..., reference = reference[1:n])
        for group in values(plot.addon_state.groups)
            marker = only(filter(p -> p isa Makie.Scatter, group))
            bars = filter(p -> p isa Makie.Errorbars, group)
            @test all(!isempty(bar[1][]) for bar in bars)
            @test all(isdisjoint(first.(marker[1][]), first.(bar[1][])) for bar in bars)
        end
        @test any(element -> element isa Makie.MarkerElement,
            first(last(only(plot.legend.entrygroups[]))).elements)
    end
    overridden = LineCableModels.plot(result; options..., errorbar_sampling = :all,
        series_attributes = (marker = :diamond, color = :magenta, whiskerwidth = 12))
    for group in values(overridden.addon_state.groups)
        marker = only(filter(p -> p isa Makie.Scatter, group))
        @test length(marker[1][]) == length(f)
        @test marker.marker[] == Makie.to_spritemarker(:diamond)
        @test all(p -> Makie.to_color(p.color[]) == Makie.to_color(:magenta), group)
        @test all(p.whiskerwidth[] == 12 for p in group if p isa Makie.Errorbars)
    end
    @test_throws r"errorbar_sampling" LineCableModels.plot(result; options..., errorbar_sampling = :invalid)
    @test isequal(before, (
        Z(reference), Y(reference), frequencies(reference), Z(candidate)))
end

@testitem "Makie addons / uncertainty legend actions preserve complete series" tags=[:visual] begin
    using CairoMakie, Measurements

    f = [1.0, 10.0, 100.0]
    sources = map((1.0, 10.0)) do offset
        z = [measurement(offset + i + 2j + k, 0.1offset + 0.01i + 0.001k) +
             im * measurement(2offset + i + j + k, 0.2offset)
             for i in 1:2, j in 1:2, k in eachindex(f)]
        LineParameters(z, z .* 1e-6, f)
    end
    page = LineCableModels.plot(sources...; ydata = (R,),
        backend = :cairo, display_plot = false, controls = true, open_export = false,
        series_labels = ("first", "second"), legend_position = :bottom,
        series_attributes = ((marker = :circle,), (marker = :utriangle,)),
        length_unit = :base, quantity_units = :base, clip = false)
    entries = last(only(page.legend.entrygroups[]))
    for entry in entries
        # Exercise the native legend action, not a direct assignment to the mean
        # line: the latter missed the exception in the error-bar child plots.
        Makie.toggle_visibility!(entry)
    end
    for axis in page.axes
        @test all(!plot.visible[]
        for plot in axis.scene.plots
        if plot isa Union{Makie.Lines, Makie.Errorbars, Makie.Scatter})
    end
    @test all(entry -> first(Makie.get_n_visible(entry)) == 0, entries)
    for entry in entries
        Makie.toggle_visibility!(entry, true)
    end
    for axis in page.axes
        @test all(plot.visible[]
        for plot in axis.scene.plots
        if plot isa Union{Makie.Lines, Makie.Errorbars, Makie.Scatter})
    end
    @test all(entry -> ==(Makie.get_n_visible(entry)...), entries)
end

@testitem "Makie addons / uncertainty coordinates and visibility survive every axis mode" tags=[:visual] begin
    using CairoMakie, Measurements

    f = [1.0, 10.0, 100.0]
    options = (
        backend = :cairo, display_plot = false, controls = true, open_export = false,
        length_unit = :base, quantity_units = :base, freq_unit = :base, clip = false)
    for (xspread, yspread, markers) in ((0, 0, false), (0, 1, false), (1, 0, true), (
        1, 1, true))
        frequency = measurement.(f, xspread .* 0.01f)
        sources = map((1.0, 10.0)) do offset
            z = [measurement(offset+i+2j+k, yspread*(0.1offset+0.01i+0.001k)) +
                 im*measurement(2offset+i+j+k, yspread*0.2offset)
                 for i in 1:2, j in 1:2, k in eachindex(f)]
            LineParameters(z, z .* 1e-6, frequency)
        end
        originals = map(p -> (copy(Z(p)), copy(Y(p)), copy(frequencies(p))), sources)
        attributes = markers ? ((marker = :circle,), (marker = :utriangle,)) : nothing
        page = LineCableModels.plot(
            sources...; options..., ydata = (R,), errorbar_sampling = :all,
            series_labels = ("first", "second"), series_attributes = attributes,
            legend_position = :bottom, panel_legends = (1, 1)=>:right)
        for (index, axis) in enumerate(page.axes)
            i, j = cld(index, 2), mod1(index, 2)
            curves = filter(p -> p isa Makie.Lines, axis.scene.plots)
            bars = filter(p -> p isa Makie.Errorbars, axis.scene.plots)
            @test length(bars) == 2*(xspread+yspread)
            for (source, curve) in zip(sources, curves)
                expected = observe(source, R)[i, j, :]
                @test first.(curve[1][]) ≈ f
                @test last.(curve[1][]) ≈ nominal.(expected)
                matching = filter(
                    bar -> Makie.to_color(bar.color[]) ==
                           Makie.to_color(curve.color[]), bars)
                @test length(matching) == xspread+yspread
                for bar in matching
                    @test getindex.(bar[1][], 1) ≈ f
                    @test getindex.(bar[1][], 2) ≈ nominal.(expected)
                    spread = bar.direction[] === :x ? uncertainty.(frequency) :
                             uncertainty.(expected)
                    @test getindex.(bar[1][], 3) ≈ spread
                    @test getindex.(bar[1][], 4) ≈ spread
                end
            end
        end
        # Panel actions are local; figure actions recover a genuinely mixed
        # panel state. No assumption about private dependency storage is needed.
        panel_entry = first(last(only(page.panel_legends[(1, 1)].entrygroups[])))
        Makie.toggle_visibility!(panel_entry)
        first_color = Makie.to_color(first(filter(p -> p isa Makie.Lines, first(page.axes).scene.plots)).color[])
        for (index, axis) in enumerate(page.axes), plot in axis.scene.plots

            plot isa Union{Makie.Lines, Makie.Errorbars, Makie.Scatter} || continue
            @test plot.visible[] ==
                  !(index == 1 && Makie.to_color(plot.color[]) == first_color)
        end
        entry = first(last(only(page.legend.entrygroups[])))
        visible, total = Makie.get_n_visible(entry)
        @test 0 < visible < total
        Makie.toggle_visibility!(entry, visible!=total)
        for cycle in 1:3
            Makie.toggle_visibility!(entry)
            for axis in page.axes, plot in axis.scene.plots

                plot isa Union{Makie.Lines, Makie.Errorbars, Makie.Scatter} || continue
                @test plot.visible[] == (Makie.to_color(plot.color[]) != first_color)
            end
            Makie.toggle_visibility!(entry, true)
        end
        for key in (:xlog, :ylog)
            page.controls[key].active[] = true
            @test !isempty(Makie.colorbuffer(page.figure))
            page.controls[key].active[] = false
        end
        for (source, original) in zip(sources, originals)
            @test isequal(Z(source), original[1])
            @test isequal(Y(source), original[2])
            @test isequal(frequencies(source), original[3])
        end
    end
end

@testitem "Makie addons / uncertainty legends survive recreation and publication" tags=[:visual] begin
    using CairoMakie, Measurements

    f = [1.0, 10.0, 100.0]
    z = [measurement(i+2j+k, 0.1i+0.01j+0.001k) + im*measurement(i+j+k, 0.2)
         for i in 1:3, j in 1:3, k in eachindex(f)]
    parameters = LineParameters(z, z .* 1e-6, measurement.(f, 0.01f))
    options = (
        backend = :cairo, display_plot = false, controls = true, open_export = false,
        length_unit = :base, quantity_units = :base, clip = false, fig_size = (1000, 700))
    pages = LineCableModels.plot(parameters; options..., ydata = (R,), layout = (2, 2),
        series_labels = ("uncertain",), series_attributes = (marker = :circle,), legend_position = :bottom)
    @test length.(getproperty.(pages, :axes)) == [4, 2, 2, 1]
    pages = copy(pages)
    push!(pages,
        Makie.plot(parameters.Z, frequencies(parameters), (R, 1, 1, :);
            options..., series_labels = ("uncertain",), legend_position = :bottom))
    # This exercises the publication renderer, not the matrix renderer used by
    # the benchmark inspector. Both must honour the same uncertainty visibility rule.
    artifact = report(
        TableReportDefinition(((R, 1, 1, :),); illustration = true,
            clip = false,
            plot_options = (backend = :cairo, display_plot = false,
                controls = true, open_export = false,
                fig_size = (1000, 700), series_labels = ("uncertain",), legend_position = :bottom)),
        parameters)
    push!(pages, artifact.illustration)
    single = Makie.plot(parameters.Y, frequencies(parameters), (G, 1, 1, :); options...)
    @test single.legend === nothing
    single_line = only(filter(p->p isa Makie.Lines, only(single.axes).scene.plots))
    single_line.visible[] = false
    @test all(!p.visible[] for p in only(single.axes).scene.plots if p isa Makie.Errorbars)
    single_line.visible[] = true
    @test all(p.visible[] for p in only(single.axes).scene.plots if p isa Makie.Errorbars)
    for page in pages
        for position in (:right, :inside, :bottom)
            legend = figurelegend!(page; position, max_fraction = 0.5)
            entries = last(only(legend.entrygroups[]))
            foreach(Makie.toggle_visibility!, entries)
            @test all(!plot.visible[] for axis in page.axes
            for plot in axis.scene.plots
            if plot isa Union{Makie.Lines, Makie.Errorbars, Makie.Scatter})
            # Recreating a hidden legend must not resurrect its series or lose
            # its visibility state. Glyph count is deliberately not an asserted requirement.
            previous_legend = legend
            legend = figurelegend!(page; position, max_fraction = 0.5, legend_labels = ("uncertain",))
            @test legend !== previous_legend
            @test all(entry -> first(Makie.get_n_visible(entry)) == 0,
                last(only(legend.entrygroups[])))
            foreach(entry -> Makie.toggle_visibility!(entry, true), last(only(legend.entrygroups[])))
            @test all(plot.visible[] for axis in page.axes
            for plot in axis.scene.plots
            if plot isa Union{Makie.Lines, Makie.Errorbars, Makie.Scatter})
        end
        @test !isempty(Makie.colorbuffer(page.figure))
    end
    page = last(pages)
    axis = only(page.axes)
    line = first(filter(p -> p isa Makie.Lines, axis.scene.plots))
    bar = first(filter(p -> p isa Makie.Errorbars, axis.scene.plots))
    bar.visible[] = false
    @test line.visible[] # Independent native edits are still allowed.
    figurelegend!(page; position = :bottom, max_fraction = 0.5, legend_labels = ("uncertain",))
    @test !bar.visible[] # A legend rebuild is not a series hide/show action.
    line.visible[] = false
    @test !bar.visible[]
    line.visible[] = true
    @test bar.visible[] # A subsequent series action restores its components.
    mktempdir() do directory
        path = export_svg(page; path = joinpath(directory, "uncertainty.svg"), open_file = false)
        @test filesize(path)>0
    end
end

@testitem "Makie addons / uncertainty legend shading follows owners after relayout" tags=[:visual] begin
    using CairoMakie, Measurements

    f = [1.0, 10.0, 100.0]
    z = reshape(complex.(measurement.([1.0, 2.0, 3.0], 0.4), 1.0), 1, 1, :)
    source = LineParameters(z, z .* 1e-6, f)
    page = LineCableModels.plot(
        source; ydata = (R,), backend = :cairo, display_plot = false,
        controls = false, open_export = false, series_labels = ("uncertain",),
        series_attributes = (marker = :circle,), legend_position = :bottom, fig_size = (
            900, 600))
    function legend_pixels(page)
        pixels = Makie.colorbuffer(page.figure)
        bounds = page.legend.layoutobservables.computedbbox[]
        viewport = page.figure.scene.viewport[]
        sx, sy = size(pixels, 2)/viewport.widths[1], size(pixels, 1)/viewport.widths[2]
        x1 = clamp(floor(Int, bounds.origin[1]*sx)+1, 1, size(pixels, 2))
        x2 = clamp(ceil(Int, (bounds.origin[1]+bounds.widths[1])*sx), x1, size(pixels, 2))
        y1 = clamp(floor(Int, (viewport.widths[2]-bounds.origin[2]-bounds.widths[2])*sy)+1, 1, size(pixels, 1))
        y2 = clamp(ceil(Int, (viewport.widths[2]-bounds.origin[2])*sy), y1, size(pixels, 1))
        return pixels[y1:y2, x1:x2]
    end
    shown = legend_pixels(page)
    Makie.toggle_visibility!(first(last(only(page.legend.entrygroups[]))))
    hidden = legend_pixels(page)
    @test shown != hidden
    figurelegend!(page; position = :bottom, max_fraction = 0.5, legend_labels = ("uncertain",))
    # Pixel comparison checks the actual shading listeners, not just the click
    # target list. Replacing a hidden legend must preserve its rendered state.
    @test legend_pixels(page) == hidden
    Makie.toggle_visibility!(first(last(only(page.legend.entrygroups[]))), true)
    @test legend_pixels(page) == shown
    figurelegend!(page; position = :bottom, max_fraction = 0.5, joinstyle = :round)
    @test !isempty(legend_pixels(page))
    axis = only(page.axes)
    view = Makie.Rect2d(2.0, 1.0, 30.0, 1.0)
    axis.targetlimits[] = view
    figurelegend!(page; position = :bottom, max_fraction = 0.5)
    # Reinitialising native legend shading is not a visibility change and
    # must not reset a caller's current zoom/pan.
    @test axis.targetlimits[] == view
end

@testitem "Makie addons / overflowing uncertainty entries retain their owners" tags=[:visual] begin
    using CairoMakie, Measurements
    f = [1.0, 10.0, 100.0]
    sources = Tuple(LineParameters(fill(measurement(Float64(i), 0.1)+im, 1, 1, 3),
                        fill(measurement(1e-6, 1e-8)+1e-6im, 1, 1, 3), f) for i in 1:24)
    legend_bounds = Observable(Rect2f(0, 0, 300, 100))
    page = LineCableModels.plot(
        sources...; ydata = (R,), backend = :cairo, display_plot = false,
        controls = false, open_export = false, fig_size = (620, 340), legend_position = :inside,
        legend_attributes = (bbox = legend_bounds,), legend_cap = 1.0,
        series_labels = Tuple("Uncertain result $i" for i in 1:24))
    Makie.colorbuffer(page.figure)
    entries = last(only(page.legend.entrygroups[]))
    @test last(entries).label[] == "(...)"
    resize!(page.figure, 1200, 950)
    legend_bounds[] = Rect2f(0, 0, 300, 900)
    Makie.colorbuffer(page.figure)
    entries = last(only(page.legend.entrygroups[]))
    @test last(entries).label[] == "Uncertain result 24"
    # This entry was absent from the compact legend. Normalising only visible
    # entries would allow its old derived-child click targets to reappear.
    Makie.toggle_visibility!(last(entries))
    axis = only(page.axes)
    @test !last(filter(p -> p isa Makie.Lines, axis.scene.plots)).visible[]
    @test !last(filter(p -> p isa Makie.Errorbars, axis.scene.plots)).visible[]
    @test first(filter(p -> p isa Makie.Lines, axis.scene.plots)).visible[]
    Makie.toggle_visibility!(last(entries), true)
    @test all(p.visible[]
    for p in axis.scene.plots if p isa Union{Makie.Lines, Makie.Errorbars})
end

@testitem "Makie / native data edits invalidate display support without source acquisition" tags=[:visual] begin
    using CairoMakie, Measurements
    frequency=collect(1.0:20.0)
    values=measurement.(fill(2.0, 20), [fill(0.1, 19); 20.0])
    raw=LineParameters(reshape(complex.(values, values), 1, 1, :), fill(1+2im, 1, 1, 20), frequency)
    observed=ObservedResult(raw, (R, X); clip = false, length_unit = :base)
    p=LineCableModels.plot(observed; ydata = (R,), errorbar_sampling = :staggered,
        backend = :cairo, display_plot = false, open_export = false)
    axis=only(p.axes)
    curve=only(filter(plot -> plot isa Makie.Lines, axis.scene.plots))
    bars=only(filter(plot -> plot isa Makie.Errorbars, axis.scene.plots))
    count=length(axis.scene.plots)
    saved=deepcopy(observed.quantities)
    native_intervals=copy(bars[1][])
    Makie.update!(curve, [1.0, 2.0, 3.0], [10.0, 20.0, 30.0])
    @test length(curve[1][])==3
    @test bars[1][]==native_intervals
    resetview!(p)
    @test sum((axis.targetlimits[].origin[2], axis.targetlimits[].widths[2]))>=30.0
    bars.visible[]=false
    resetview!(p)
    @test axis.targetlimits[].origin[2]>0
    for _ in 1:3
        Makie.resize!(p.figure, 850, 650)
        Makie.resize!(p.figure, 900, 700)
    end
    @test !bars.visible[]
    @test length(axis.scene.plots)==count
    @test isequal(observed.quantities, saved)
end
