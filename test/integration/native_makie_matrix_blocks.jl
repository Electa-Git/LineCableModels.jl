@testitem "Makie addons / matrix pages retain coordinates and residual axis dimensions" tags=[:visual] begin
    using CairoMakie
    using LineCableModels

    frequency = 10.0 .^ range(-1, 7; length=13)
    z = [complex(10.0^(i-j)*(1+log10(f+1)), i+j+f*1e-3)
        for i in 1:5, j in 1:5, f in frequency]
    parameters = LineParameters(z, z.*1e-6, frequency)
    options = (; backend=:cairo, display_plot=false, controls=true,
        open_export=false, length_unit=:base, fig_size=(1200,800))
    pages = LineCableModels.plot(parameters; ydata=(R,), layout=(2,3), options...,
        panel_legends=(5,5)=>(position=:inside, overflow=:show_all))
    @test length(pages) == 6
    @test [length(page.axes) for page in pages] == [6,4,6,4,3,2]
    @test [page.addon_state.panel_page.index for page in pages] ==
        [(1,1),(1,2),(2,1),(2,2),(3,1),(3,2)]
    seen = Tuple{Int,Int}[]
    for page in pages
        @test page.addon_state.nominal_capacity == (2,3)
        @test count(block -> block isa Makie.Axis, page.figure.content) == length(page.axes)
        @test all(>(0),page.figure.scene.viewport[].widths)
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
        layout=(2,3), options...)
    @test [page.addon_state.panel_page.index for page in sparse] == [(1,1),(1,2),(3,1),(3,2)]
    @test Set(coordinate for page in sparse for coordinate in keys(page.addon_state.panel_data)) ==
        Set(((1,2),(1,5),(5,2),(5,5)))
    small = LineParameters(z[1:2,1:2,:], z[1:2,1:2,:].*1e-6, frequency)
    full = LineCableModels.plot(small; ydata=(R,), options...)
    @test full isa UIPlot && length(full.axes) == 4
    @test full.addon_state.panel_page.dimensions==(2,2)
    one = LineCableModels.plot(small; ydata=(R,), layout=(3,4), title="Example", options...)
    @test one isa UIPlot && length(one.axes) == 4
    @test one.export_name=="Example"
    for layout in ((0,2),(-1,2),(true,2),(2.0,2),(2,),[2,2])
        @test_throws ArgumentError LineCableModels.plot(small; ydata=(R,), layout, options...)
    end
    @test_throws ArgumentError LineCableModels.plot(small; ydata=(R,), blocks=(2,2),
        layout=(1,2), options...)
    family = LineCableModels.plot(small; ydata=(R,G), layout=(1,2), options...,
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
            # Export reclaims chrome, so frame origins may move. Data-frame
            # dimensions remain equal within native pixel rounding.
            @test !isempty(rendered) && all(all(all(isapprox.(actual.widths,expected.widths;atol=1))
                for (actual,expected) in zip(rectangles,before)) for rectangles in rendered)
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
    records = [Formulation(earth_impedance=id)
        for id in (:xue2018,:default,:pollaczek1926,:saad1996,:wedepohl1973)]
    source_id = LineCableModels.Grammar.gridpoint_id().source_id
    completed = [LineCableModels.Engine.retain_gridpoint(reference,
        LineCableModels.Grammar.gridpoint_id(;source_id,formulation_index=index);
        fields=merge(LineCableModels.Engine.completed_formulation(formula),
            (inputs=(resistivity=100.,),coordinates=["core"])))
        for (index,formula) in enumerate(records)]
    result = ParametricResult(nothing,completed,
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
    @test last(markers).marker[] == Makie.to_spritemarker(:circle)
    @test Makie.to_color(last(curves).color[]) == Makie.to_color(:black)
    @test allunique([marker.marker[] for marker in markers[1:end-1]])
    @test Makie.alpha(last(markers).color[]) == 0
    @test all(marker -> 1 <= length(marker[1][]) < length(frequency),markers)
    @test all(marker -> all(point -> point in first(curves)[1][],marker[1][]),markers)
    @test allunique([first(marker[1][]) for marker in markers])
    @test first(last(markers)[1][]) == first(first(curves)[1][])
    @test last(last(markers)[1][]) == last(first(curves)[1][])
    @test allunique(last(markers)[1][])
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
    @test [curve.color[] for curve in filtered_curves] == [curves[i].color[] for i in [5,2,6]]
    @test [marker.marker[] for marker in filtered_markers] == [markers[i].marker[] for i in [5,2,6]]
    @test [marker[1][] for marker in filtered_markers] == [markers[i][1][] for i in [5,2,6]]
    plain = LineCableModels.plot(result; ydata=(R,),reference,
        series_attributes=(marker=nothing,),options...)
    @test !any(plot -> plot isa Makie.Scatter,only(plain.axes).scene.plots)
    short_reference = reference[1:2]
    short_result = ParametricResult(nothing,[point[1:2] for point in completed],result.axes, ComputationDetails((;)))
    short = LineCableModels.plot(short_result; ydata=(R,),reference=short_reference,options...)
    @test all(group -> any(plot -> plot isa Makie.Scatter && !isempty(plot[1][]),group),
        values(short.addon_state.groups))
    short_markers = only(filter(plot -> plot isa Makie.Scatter,
        short.addon_state.groups[last(short.addon_state.order)]))
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
        @test last(last(markers)[1][]) == last(first(curves)[1][])
        @test allunique(last(markers)[1][])
    end
    mktempdir() do directory
        @test isfile(export_svg(page;path=joinpath(directory,"overlap.svg"),open_file=false))
    end
    @test Z(reference) == z
end

@testitem "Makie addons / nominal capacity calibrates frames and compact residuals" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    f=[1.,10.,100.]
    z=[complex(i+j+k,i-j+k) for i in 1:3,j in 1:3,k in 1:3]
    raw=LineParameters(z,z.*1e-6,f)
    pages=LineCableModels.plot(raw;ydata=(R,L),layout=(2,2),fig_size=(1000,700),
        controls=false,display_plot=false,backend=:cairo,clip=false)
    @test length(pages)==8
    @test [length(p.axes) for p in pages]==[4,2,2,1,4,2,2,1]
    @test [p.addon_state.panel_page.dimensions for p in pages[1:4]]==[(2,2),(2,1),(1,2),(1,1)]
    frames=[Tuple(axis.layoutobservables.computedbbox[].widths) for p in pages for axis in p.axes]
    @test all(frame -> all(isapprox.(frame,first(frames);atol=1)),frames)
    sizes=[Tuple(p.figure.scene.viewport[].widths) for p in pages]
    @test all(abs.(first(sizes).- (1000,700)).<=1)
    ext=Base.get_extension(LineCableModels,:LineCableModelsMakieExt)
    for p in pages[1:2]
        before=Tuple(p.figure.scene.viewport[].widths)
        for _ in 1:3
            ext._addon_edit_presentation!(() -> nothing,p)
            @test Tuple(p.figure.scene.viewport[].widths)==before
        end
    end
    @test sizes[2][1]<sizes[1][1]
    @test sizes[3][2]<sizes[1][2]
    @test all(sizes[4].<sizes[1])
    @test all(count(block -> block isa Makie.Axis,p.figure.content)==length(p.axes) for p in pages)
    sibling=[(Tuple(p.figure.scene.viewport[].widths),[axis.targetlimits[] for axis in p.axes]) for p in pages[2:end]]
    p=first(pages)
    before=[Tuple(axis.layoutobservables.computedbbox[].widths) for axis in p.axes]
    views=[axis.targetlimits[] for axis in p.axes]
    figuretitle!(p,"A new title\nwith two lines")
    figurelegend!(p;position=:bottom)
    @test all(all(isapprox.(Tuple(axis.layoutobservables.computedbbox[].widths),size;atol=1)) for (axis,size) in zip(p.axes,before))
    @test [axis.targetlimits[] for axis in p.axes]==views
    @test [(Tuple(p.figure.scene.viewport[].widths),[axis.targetlimits[] for axis in p.axes]) for p in pages[2:end]]==sibling
    oldsize=Tuple(p.figure.scene.viewport[].widths)
    Makie.resize!(p.figure,oldsize[1]+150,oldsize[2]-100)
    @test length(p.axes)==4
    @test p.addon_state.panel_page.dimensions==(2,2)
    @test [axis.targetlimits[] for axis in p.axes]==views
    @test_throws ArgumentError LineCableModels.plot(raw;ydata=(R,),blocks=(2,2),display_plot=false)

    design=TestFixtures.coaxial_design()
    previews=preview(fill(design,4);layout=(1,2),size=(1000,600),
        display_plot=false,backend=:cairo,controls=false,panel_titles=("one","two","three","four"))
    @test length(previews)==2
    @test [length(p.axes) for p in previews]==[2,2]
    @test Set(keys(previews[1].addon_state.panel_data))==Set((1,2))
    @test Set(keys(previews[2].addon_state.panel_data))==Set((3,4))
    @test [axis.title[] for p in previews for axis in p.axes]==["one","two","three","four"]
    @test all(length(p.colorbars)==3 for p in previews)
    strip=preview(fill(design,4);layout=(1,4),size=(1200,600),
        display_plot=false,backend=:cairo,controls=false)
    @test strip isa UIPlot
    @test length(strip.axes)==4
    @test all(isapprox(axis.scene.viewport[].widths[1],axis.scene.viewport[].widths[2];atol=1) for axis in strip.axes)
    @test all(axis.scene.viewport[].widths[1]>60 for axis in strip.axes)
end

@testitem "Makie addons / selected matrix footprints preserve original coordinates" tags=[:visual] begin
    using CairoMakie
    using LineCableModels
    ext=Base.get_extension(LineCableModels,:LineCableModelsMakieExt)
    f=[1.,10.,100.]
    z=[complex(10i+j+k,i-j+k) for i in 1:4,j in 1:4,k in 1:3]
    raw=LineParameters(z,z.*1e-6,f)
    options=(;controls=false,display_plot=false,backend=:cairo,clip=false,fig_size=(1000,650))
    for coordinate in ((1,1),(2,2))
        i,j=coordinate
        pages=LineCableModels.plot(raw,(Z,i,j,:);options...,series_labels=("retained",),
            legend_position=:inside,legend_attributes=(halign=:right,valign=:bottom))
        @test length(pages)==2 # quantities remain separate families
        for p in pages
            @test p.addon_state.nominal_capacity==(1,1)
            @test p.addon_state.panel_page.dimensions==(1,1)
            @test p.addon_state.panel_page.origin==coordinate
            @test Set(keys(p.addon_state.page_cells))==Set([(1,1)])
            @test Set(keys(p.addon_state.panel_data))==Set([coordinate])
            axis=only(p.axes)
            frame=axis.scene.viewport[]
            legend=p.legend.layoutobservables.computedbbox[]
            @test all(legend.origin .>= frame.origin .- 1)
            @test all(legend.origin+legend.widths .<= frame.origin+frame.widths .+ 1)
        end
        curve=only(filter(x -> x isa Makie.Lines,only(first(pages).axes).scene.plots))
        @test last.(curve[1][])≈real.(z[i,j,:]).*1000
    end
    selected=LineCableModels.plot(raw,(R,2,2:3,:);options...)
    @test selected isa UIPlot
    @test selected.addon_state.nominal_capacity==(1,2)
    @test selected.addon_state.panel_page.origin==(2,2)
    @test selected.addon_state.panel_page.dimensions==(1,2)
    split=LineCableModels.plot(raw,(R,2,2:3,:);options...,layout=(1,2))
    @test [p.addon_state.panel_page.index for p in split]==[(2,1),(2,2)]
    @test all(p -> p.addon_state.panel_page.dimensions==(1,1),split)
    @test [p.addon_state.panel_page.origin for p in split]==[(2,2),(2,3)]
    sparse=LineCableModels.plot(raw,(R,[1,3],[1,3],:);options...)
    @test length(sparse.axes)==4
    @test sparse.addon_state.panel_page.dimensions==(3,3)
    @test length(sparse.addon_state.page_cells)==9 # internal holes are intentional
    geometry=only(ext._addon_matrix_pages([(1,1),(3,3)],(3,4),(3,4)))
    @test geometry.dimensions==(3,3)
    @test geometry.positions==((1,1),(3,3))
    explicit=LineCableModels.plot(raw,(R,2,2,:);options...,layout=(1,2))
    full=LineCableModels.plot(raw,(R,2,1:2,:);options...,layout=(1,2))
    @test explicit.addon_state.panel_page.index==(2,1)
    @test explicit.addon_state.panel_page.origin==(2,2)
    @test explicit.addon_state.panel_page.dimensions==(1,1)
    @test all(isapprox.(only(explicit.axes).scene.viewport[].widths,
        first(full.axes).scene.viewport[].widths;atol=1))
    @test explicit.figure.scene.viewport[].widths[1]<full.figure.scene.viewport[].widths[1]
    cells=copy(explicit.addon_state.page_cells)
    for group in values(explicit.addon_state.groups), primitive in group
        primitive.visible[]=false
    end
    @test explicit.addon_state.page_cells==cells
end
