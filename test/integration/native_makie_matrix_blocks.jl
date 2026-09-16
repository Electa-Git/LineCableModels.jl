@testitem "Makie addons / matrix blocks retain coordinates and residual axis dimensions" tags=[:visual] begin
    using CairoMakie
    using LineCableModels

    frequency = 10.0 .^ range(-1, 7; length=13)
    z = [complex(10.0^(i-j)*(1+log10(f+1)), i+j+f*1e-3)
        for i in 1:5, j in 1:5, f in frequency]
    parameters = LineParameters(z, z.*1e-6, frequency)
    options = (; backend=:cairo, display_plot=false, controls=true,
        open_export=false, length_unit=:base, fig_size=(1200,800))
    pages = LineCableModels.plot(parameters; ydata=(R,), blocks=(2,3), options...,
        panel_legends=(5,5)=>(position=:inside, overflow=:show_all))
    @test length(pages) == 6
    @test [length(page.axes) for page in pages] == [6,4,6,4,3,2]
    @test [page.addon_state.matrix_block.index for page in pages] ==
        [(1,1),(1,2),(2,1),(2,2),(3,1),(3,2)]
    seen = Tuple{Int,Int}[]
    for page in pages
        @test page.addon_state.matrix_block.dimensions == (2,3)
        @test count(block -> block isa Makie.Axis, page.figure.content) == length(page.axes)
        @test page.figure.scene.viewport[].widths[1] > page.figure.scene.viewport[].widths[2]
        for (coordinate, panel) in page.addon_state.panel_data
            push!(seen, coordinate)
            i,j = coordinate
            curve = only(filter(plot -> plot isa Makie.Lines, panel.axis.scene.plots))
            @test first.(curve[1][]) ≈ frequency
            @test last.(curve[1][]) ≈ real.(z[i,j,:])
        end
        @test !isempty(Makie.colorbuffer(page.figure))
    end
    @test allunique(seen)
    @test Set(seen) == Set(Iterators.product(1:5,1:5))
    @test all(isempty(page.panel_legends) for page in pages[1:5])
    @test haskey(last(pages).panel_legends, (5,5))
    @test all(axis -> axis.xlabelvisible[], last(pages).axes)
    dimensions = [Tuple(axis.layoutobservables.computedbbox[].widths)
        for page in pages for axis in page.axes]
    @test all(size -> all(isapprox.(size, first(dimensions); atol=1)), dimensions)
    mktempdir() do directory
        for page in pages[[1,end]]
            before = [axis.layoutobservables.computedbbox[] for axis in page.axes]
            path = export_svg(page; path=joinpath(directory,page.export_name*".svg"), open_file=false)
            @test filesize(path) > 0
            @test [axis.layoutobservables.computedbbox[] for axis in page.axes] == before
            @test !isempty(Makie.colorbuffer(page.figure))
        end
    end

    # Selections retain their original slots, not a packed/renumbered submatrix.
    sparse = LineCableModels.plot(parameters; ydata=((R,[1,5],[2,5],:),),
        blocks=(2,3), options...)
    @test [page.addon_state.matrix_block.index for page in sparse] == [(1,1),(1,2),(3,1),(3,2)]
    @test Set(coordinate for page in sparse for coordinate in keys(page.addon_state.panel_data)) ==
        Set(((1,2),(1,5),(5,2),(5,5)))
    small = LineParameters(z[1:2,1:2,:], z[1:2,1:2,:].*1e-6, frequency)
    full = LineCableModels.plot(small; ydata=(R,), options...)
    @test full isa UIPlot && length(full.axes) == 4
    @test !haskey(full.addon_state,:matrix_block)
    one = LineCableModels.plot(small; ydata=(R,), blocks=(3,4), title="Example", options...)
    @test one isa UIPlot && length(one.axes) == 4
    @test endswith(one.export_name,"Series resistance (1,1)")
    for blocks in ((0,2),(-1,2),(true,2),(2.0,2),(2,),[2,2])
        @test_throws ArgumentError LineCableModels.plot(small; ydata=(R,), blocks, options...)
    end
    @test_throws ArgumentError LineCableModels.plot(small; ydata=(R,), blocks=(2,2),
        layout=(1,2), options...)
    family = LineCableModels.plot(small; ydata=(R,G), blocks=(1,2), options...,
        figure_title=("One line", "Two\nlines", "One line", "One line"),
        series_labels=("Saved curve",), legend_position=:bottom,
        panel_legends=(2,2)=>(position=:right, overflow=:show_all))
    rectangles = [axis.layoutobservables.computedbbox[] for page in family for axis in page.axes]
    @test all(rectangle -> all(isapprox.(rectangle.widths,first(rectangles).widths;atol=1)),rectangles)
    mktempdir() do directory
        for (index,page) in enumerate(family), theme in (:default,:publication)
            before = [axis.layoutobservables.computedbbox[] for axis in page.axes]
            rendered = Any[]
            observer = on(page.figure.scene, Makie.events(page.figure).tick) do tick
                tick.state === Makie.OneTimeRenderTick &&
                    push!(rendered,[axis.layoutobservables.computedbbox[] for axis in page.axes])
            end
            export_svg(page;path=joinpath(directory,"block-$index-$theme.svg"),theme,open_file=false)
            off(observer)
            @test !isempty(rendered) && all(==(before),rendered)
            @test [axis.layoutobservables.computedbbox[] for axis in page.axes] == before
        end
    end
    @test Z(parameters) == z
end

@testitem "Makie addons / stable sparse comparison markers and measured legend grid" tags=[:visual] begin
    using CairoMakie
    using LineCableModels

    frequency = 10.0 .^ range(-1,7; length=101)
    z = reshape(complex.(1 .+ log10.(frequency.+1), frequency.*1e-4),1,1,:)
    reference = LineParameters(z,z.*1e-6,frequency)
    records = [NamedTuple(Formulation(earth_impedance=id))
        for id in (:xue2018,:default,:pollaczek1926,:saad1996,:wedepohl1973)]
    result = ParametricResult(nothing,fill(reference,5),
        (problems=[:one],formulations=records), ComputationDetails((;)))
    options = (; backend=:cairo,display_plot=false,controls=true,open_export=false,
        fig_size=(1000,650),length_unit=:base)
    page = LineCableModels.plot(result; ydata=(R,),reference,options...)
    axis = only(page.axes)
    curves = filter(plot -> plot isa Makie.Lines,axis.scene.plots)
    markers = [only(filter(plot -> plot isa Makie.Scatter,page.addon_state.groups[group]))
        for group in page.addon_state.order]
    @test length(curves) == length(markers) == 6
    @test all(curve -> curve.linestyle[] == Makie.to_linestyle(:solid),curves)
    @test all(curve -> curve[1][] == first(curves)[1][],curves)
    @test allunique([curve.color[] for curve in curves])
    @test first(markers).marker[] == Makie.to_spritemarker(:circle)
    @test markers[3].marker[] == Makie.to_spritemarker(:utriangle) # default is not the first candidate
    @test Makie.alpha(first(markers).color[]) == 0
    @test all(marker -> 1 <= length(marker[1][]) < length(frequency),markers)
    @test all(marker -> all(point -> point in first(curves)[1][],marker[1][]),markers)
    @test allunique([first(marker[1][]) for marker in markers])
    @test first(first(markers)[1][]) == first(first(curves)[1][])
    @test last(first(markers)[1][]) == last(first(curves)[1][])
    @test allunique(first(markers)[1][])
    extension = Base.get_extension(LineCableModels,:LineCableModelsMakieExt)
    other_reference = extension._addon_marker_coordinates(first(curves),(3,10);endpoints=true)
    @test first(other_reference[]) == first(first(curves)[1][])
    @test last(other_reference[]) == last(first(curves)[1][])
    @test allunique(other_reference[])
    @test maximum(findall(plot -> plot isa Makie.Lines,axis.scene.plots)) <
        minimum(findall(plot -> plot isa Makie.Scatter,axis.scene.plots))
    entries = last(only(page.legend.entrygroups[]))
    @test length(entries) == 6
    @test page.legend.nbanks[] > 1
    @test all(entry -> any(element -> element isa Makie.MarkerElement,entry.elements),entries)
    Makie.toggle_visibility!(first(entries))
    @test !first(curves).visible[] && !first(markers).visible[]
    Makie.toggle_visibility!(first(entries))
    @test first(curves).visible[] && first(markers).visible[]

    filtered = LineCableModels.plot(result; ydata=(R,),reference,formulations=[5,2],options...)
    filtered_curves = filter(plot -> plot isa Makie.Lines,only(filtered.axes).scene.plots)
    filtered_markers = [only(filter(plot -> plot isa Makie.Scatter,filtered.addon_state.groups[group]))
        for group in filtered.addon_state.order]
    @test [curve.color[] for curve in filtered_curves] == [curves[i].color[] for i in [1,6,3]]
    @test [marker.marker[] for marker in filtered_markers] == [markers[i].marker[] for i in [1,6,3]]
    plain = LineCableModels.plot(result; ydata=(R,),reference,
        series_attributes=(marker=nothing,),options...)
    @test !any(plot -> plot isa Makie.Scatter,only(plain.axes).scene.plots)
    short_reference = reference[1:2]
    short_result = ParametricResult(nothing,fill(short_reference,5),result.axes, ComputationDetails((;)))
    short = LineCableModels.plot(short_result; ydata=(R,),reference=short_reference,options...)
    @test all(group -> any(plot -> plot isa Makie.Scatter && !isempty(plot[1][]),group),
        values(short.addon_state.groups))
    @test last(only(short.axes).scene.plots).marker[] == Makie.to_spritemarker(:utriangle)
    short_markers = only(filter(plot -> plot isa Makie.Scatter,
        short.addon_state.groups[first(short.addon_state.order)]))
    @test length(short_markers[1][]) == 2
    moved = figurelegend!(page;position=:top,overflow=:show_all)
    Makie.toggle_visibility!(first(last(only(moved.entrygroups[]))))
    @test !first(curves).visible[] && !first(markers).visible[]
    Makie.toggle_visibility!(first(last(only(moved.entrygroups[]))))
    @test first(curves).visible[] && first(markers).visible[]

    # Reflow both entries and labels while retaining native click targets.
    long_labels = ["F$i · earth Z=" * repeat("long_identifier",12) * "; earth Y=default" for i in 1:6]
    long = LineCableModels.plot(result; ydata=(R,),reference,series_labels=long_labels,options...)
    @test !isempty(Makie.colorbuffer(long.figure))
    @test long.legend.nbanks[] == 1
    @test all(occursin('\n',entry.label[]) for entry in last(only(long.legend.entrygroups[])))
    @test long.legend.layoutobservables.autosize[][1] <= long.figure.scene.viewport[].widths[1]
    @test Set(values(long.addon_state.labels)) == Set(long_labels)
    for size in ((700,500),(1400,900),(1000,650))
        resize!(page.figure,size...)
        @test !isempty(Makie.colorbuffer(page.figure))
        @test page.legend.layoutobservables.autosize[][1] <= size[1]
        @test length(last(only(page.legend.entrygroups[]))) == 6
        @test last(first(markers)[1][]) == last(first(curves)[1][])
        @test allunique(first(markers)[1][])
    end
    mktempdir() do directory
        @test isfile(export_svg(page;path=joinpath(directory,"overlap.svg"),open_file=false))
    end
    @test Z(reference) == z
end
