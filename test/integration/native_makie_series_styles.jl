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
    xue = only(filter(p -> p isa Makie.Lines,page.addon_state.groups[:result_3]))
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

@testitem "Makie addons / large style sets retain original slots and palette" tags=[:visual] begin
    using CairoMakie
    ext=Base.get_extension(LineCableModels,:LineCableModelsMakieExt)
    prefix=copy(ext._ADDON_CURVE_COLORS)
    styles=ext._addon_comparison_styles(collect(1:324),fill(:candidate,324),324)
    @test length(styles)==324
    @test ext._ADDON_CURVE_COLORS[eachindex(prefix)]==prefix
    # Values recorded from the pre-change palette, including distant slots.
    expected=(1=>(0.0,0.36,0.68),
        2=>(0.85,0.5808333333333333,0.04250000000000004),
        3=>(0.85,0.04250000000000004,0.85),
        18=>(0.85,0.3789583333333334,0.04250000000000004),
        81=>(0.31499999999999995,0.7,0.4112499999999999),
        162=>(0.4604166666666667,0.38249999999999995,0.85),
        324=>(0.5716666666666667,0.31499999999999995,0.7))
    for (slot,rgb) in expected
        color=styles[slot].attributes.color
        @test [color.r,color.g,color.b]≈collect(rgb) rtol=1e-14 atol=0
        @test styles[slot].phase==(slot,324)
    end
    sparse=ext._addon_comparison_styles([324,0,2,324],
        [:candidate,:reference,:candidate,:candidate],324)
    @test sparse[[1,3,4]]==styles[[324,2,324]]
    @test sparse[2].attributes.color==ext.RGB(0.,0.,0.)
    @test sparse[2].attributes.marker===:circle
    @test sparse[2].hollow && sparse[2].endpoints
    @test sparse[2].phase==(1,1) && sparse[2].priority==1
    @test only(ext._addon_comparison_styles([0],[:reference],0))==sparse[2]
    @test isempty(ext._addon_comparison_styles(Int[],Symbol[],0))
    @test_throws ArgumentError ext._addon_comparison_styles([0],[:candidate],1)
    @test ext._addon_comparison_styles([2,1],[:candidate,:candidate],324)==styles[[2,1]]

    for attributes in (nothing,(color=:red,),
            Tuple((linewidth=i,) for i in 1:324),[(linewidth=i,) for i in 1:324])
        normalized=ext._series_attributes(attributes,324)
        @test length(normalized)==324
        for slot in (1,18,324)
            expected_attributes=attributes===nothing ? (;) : attributes isa NamedTuple ? attributes : attributes[slot]
            @test normalized[slot]==expected_attributes
        end
    end
    @test_throws ArgumentError ext._series_attributes(((color=:red,),),324)
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

@testitem "Makie addons / automatic difference labels and chromatic candidate prefix" tags=[:visual] begin
    using CairoMakie
    using LineCableModels.Engine: retain_gridpoint,completed_formulation
    using LineCableModels.Grammar: gridpoint_id,observation_labels
    ext=Base.get_extension(LineCableModels,:LineCableModelsMakieExt)
    options=(backend=:cairo,display_plot=false,controls=false,open_export=false)
    source_id=gridpoint_id().source_id
    analytical=Formulation(earth_impedance=:unified,earth_admittance=:unified)
    fem=Formulation(:LineCableModelsFEM)
    z=reshape(complex.(1.:12.,21.:32.),2,2,3)
    function point(formulation,index;rho=100.,problem=1)
        inputs=(rho=rho,radius=.0425,field_descriptions=(rho=(name="electrical resistivity",unit="Ω·m"),radius=(name="radius",unit="m")))
        retain_gridpoint(LineParameters(copy(z),z.*1e-6,[1.,10.,100.]),
            gridpoint_id(;source_id,problem_index=problem,formulation_index=index);
            fields=merge(completed_formulation(formulation),(inputs=inputs,)))
    end
    a=point(analytical,1)
    # Synthetic numerical fixture retains a FEM identity; no FEM solve is claimed.
    reference=point(fem,2)
    pages=LineCableModels.plot(a;ydata=(R,X,G,B),reference,length_unit=:base,options...)
    @test length(pages)==4
    candidates=LineCableModels.plot(a;ydata=(R,X,G,B),length_unit=:base,options...)
    native_lines(axis)=filter(item -> item isa Makie.Lines,axis.scene.plots)
    expected_reference=description(fem;compact=true)*" (reference)"
    for (page,solo,request) in zip(pages,candidates,(R,X,G,B))
        method=request in (R,X) ? analytical.methods.earth_impedance : analytical.methods.earth_admittance
        @test Set(values(page.addon_state.labels))==Set((description(method;compact=true),expected_reference))
        for (axis,solo_axis) in zip(page.axes,solo.axes)
            curves=native_lines(axis)
            @test Makie.to_color(curves[1].color[])==Makie.to_color(ext._addon_comparison_color(1))
            @test curves[1].color[]==only(native_lines(solo_axis)).color[]
            @test Makie.to_color(curves[2].color[])==Makie.to_color(:black)
        end
    end
    prefix=[ext._addon_comparison_color(i) for i in 1:2]
    palette=[ext._addon_comparison_color(i) for i in 1:12]
    @test palette[1:2]==prefix
    @test length(unique(palette))==12
    @test all(ext._addon_candidate_color,palette)
    labs=ext.Oklab.(palette)
    distances=[sqrt(ext._addon_color_distance(labs[i],labs[j])) for i in 1:12 for j in i+1:12]
    @test minimum(distances)>.08 # measured engineering separation for this finite prefix
    population=[point(analytical,i;rho=100.0i,problem=i) for i in 1:12]
    p=LineCableModels.plot(population;ydata=(@observe(R[1,1,:]),),length_unit=:base,options...)
    @test [Makie.to_color(c.color[]) for c in native_lines(only(p.axes))]==Makie.to_color.(palette)
    @test all(label -> occursin("Ω·m",label) && !occursin("radius",label),values(p.addon_state.labels))
    small=LineCableModels.plot(population[1:2];ydata=(@observe(R[1,1,:]),),length_unit=:base,options...)
    @test [c.color[] for c in native_lines(only(small.axes))]==[c.color[] for c in native_lines(only(p.axes))[1:2]]
    filtered=LineCableModels.plot(population;ydata=(R,X,G,B),reference,problem=[2,5],layout=(1,2),length_unit=:base,options...)
    @test length(filtered)==8
    for page in filtered, axis in page.axes
        @test [Makie.to_color(c.color[]) for c in native_lines(axis)]==Makie.to_color.([palette[2],palette[5],ext.RGB(0.,0.,0.)])
    end
    # Rebuilding native guides, changing views and exporting use detached data.
    observed=pages[1].addon_state.observed
    Z(a).=NaN;Z(reference).=NaN
    p=LineCableModels.plot(observed[1];reference=observed[2],ydata=(R,),options...,controls=true)
    labels=copy(p.addon_state.labels)
    colors=[c.color[] for c in native_lines(first(p.axes))]
    owner=first(native_lines(first(p.axes)))
    owner.visible[]=false
    figurelegend!(p;position=:right)
    figurelegend!(p;position=:bottom)
    @test !owner.visible[]
    @test p.addon_state.labels==labels
    owner.visible[]=true
    axisscale!(p,:y,:log10);resetview!(p)
    resize!(p.figure.scene,(1100,800))
    view=[axis.finallimits[] for axis in p.axes]
    mktempdir() do directory
        export_svg(p;path=joinpath(directory,"automatic.svg"),theme=:publication,open_file=false)
        @test [axis.finallimits[] for axis in p.axes]==view
    end
    @test [c.color[] for c in native_lines(first(p.axes))]==colors
    @test p.addon_state.labels==labels
end
