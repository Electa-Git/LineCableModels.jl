@testitem "PlotBuilder / row overlays separate gridpoints and retain physical columns" tags=[:visual] begin
    using CairoMakie
    import LineCableModels.Engine as E
    import LineCableModels.Grammar as G

    f=[10.0,20.0,40.0,80.0]
    voltage=[complex(i+2j+k/10, i-j+k/5) for i in 1:3,j in 1:3,k in eachindex(f)]
    current=[complex(2i+j+k/5, j-i-k/10) for i in 1:3,j in 1:3,k in eachindex(f)]
    modal=LineParameters(ModalDomain(ModalOperators(voltage,current),fill(.01+.02im,3,4)),
        SeriesImpedance(voltage),ShuntAdmittance(current.*1e-6),f,
        ComputationDetails((inputs=(system=(line_length=20.0,),),
            phase_coordinates=["core","screen","armor"])))
    source_id=G.gridpoint_id().source_id
    modal_points=[E.retain_gridpoint(modal,
        G.gridpoint_id(;source_id,problem_index=index,formulation_index=1);
        fields=(inputs=(system=(line_length=20.0index,),),)) for index in 1:2]
    segments=collect(Gridspace{PropagationParameters}((Grid(modal_points),)))
    requests=((Tv,abs),(Tv,angle),(Ti,abs),(Ti,angle))
    observed=observables.(segments,Ref(requests))
    options=(overlay=:rows,backend=:cairo,display_plot=false,controls=false,open_export=false)
    lines(axis)=filter(p -> p isa Makie.Lines,axis.scene.plots)
    kept=UIPlot[]

    result=LineCableModels.plot(observed;layout=(1,2),options...)
    append!(kept,result)
    @test length(result)==16 # four components × two points × two pages
    expected=(abs.(voltage),rad2deg.(angle.(voltage)),abs.(current),rad2deg.(angle.(current)))
    for (component,values) in enumerate(expected), point in 1:2
        pair=result[(component-1)*4+(point-1)*2+1:(component-1)*4+point*2]
        @test [length(page.axes) for page in pair]==[2,1]
        @test [page.addon_state.panel_page.index for page in pair]==[(1,1),(2,1)]
        @test [page.addon_state.panel_page.coordinates for page in pair]==
            [((point,1),(point,2)),((point,3),)]
        @test all(page -> page.addon_state.displayed_indices==[point],pair)
        @test all(page -> page.addon_state.nominal_capacity==(1,2),pair)
        labels=G.observation_labels(observed;request=requests[component],fallback="")
        @test all(page -> occursin(labels[point],page.export_name),pair)
        for (mode,axis) in enumerate(vcat(pair[1].axes,pair[2].axes))
            @test endswith(string(axis.title[]),"Mode $mode")
            curves=lines(axis)
            @test length(curves)==3
            @test [curve.label[] for curve in curves]==
                ["Conductor core","Conductor screen","Conductor armor"]
            for row in 1:3
                @test first.(curves[row][1][])≈f rtol=8eps(Float32)
                @test last.(curves[row][1][])≈values[row,mode,:] rtol=8eps(Float32)
                @test curves[row].color[]==lines(result[1].axes[1])[row].color[]
            end
        end
    end
    # Each page is rendered by the same native renderer, with mode titles and
    # conductor legends. Retained matrix arrays are unchanged.
    @test !isempty(Makie.colorbuffer(first(result).figure))
    @test Tv(segments[1])==voltage && Ti(segments[1])==current

    automatic=LineCableModels.plot(observed;ydata=((Tv,abs),),options...)
    append!(kept,automatic)
    @test length(automatic)==2
    @test all(page -> page.addon_state.nominal_capacity==(2,2),automatic)
    @test all(page -> length(page.axes)==3,automatic)
    scalar=LineCableModels.plot(first(observed);ydata=((Tv,abs),),options...)
    singleton=LineCableModels.plot([first(observed)];ydata=((Tv,abs),),options...)
    raw=LineCableModels.plot(first(segments);ydata=((Tv,abs),),options...)
    append!(kept,[scalar,singleton,raw])
    @test scalar.addon_state.panel_page==singleton.addon_state.panel_page==raw.addon_state.panel_page
    @test scalar.export_name==singleton.export_name==raw.export_name
    @test !occursin(r"(?:Point|Result) 1",scalar.export_name)
    @test [[curve[1][] for curve in lines(axis)] for axis in scalar.axes]==
          [[curve[1][] for curve in lines(axis)] for axis in singleton.axes]
    named=LineCableModels.plot((short=segments[1],long=segments[2]);
        ydata=((Tv,abs),),options...)
    append!(kept,named)
    @test length(named)==2
    @test all(page -> [curve.label[] for curve in lines(page.axes[1])]==
        ["Conductor core","Conductor screen","Conductor armor"],named)
    different_names=observables(PropagationParameters(E.retain_gridpoint(modal,G.gridpoint_id();
        fields=(phase_coordinates=["inner","shield","outer"],))),requests)
    labelled=LineCableModels.plot([first(observed),different_names];ydata=((Tv,abs),),options...)
    append!(kept,labelled)
    @test all(axis -> [curve.label[] for curve in lines(axis)]==
        ["Conductor inner","Conductor shield","Conductor outer"],labelled[2].axes)
    @test [curve.color[] for curve in lines(labelled[1].axes[1])]==
        [curve.color[] for curve in lines(labelled[2].axes[1])]

    selected=@observe (Tv,abs)[[3,1],[3,1],2:3]
    reordered=LineCableModels.plot(first(observed);ydata=(selected,),
        layout=(1,1),options...)
    append!(kept,reordered)
    @test [only(page.addon_state.panel_page.coordinates) for page in reordered]==[(1,3),(1,1)]
    for (page,mode) in zip(reordered,[3,1]), (curve,row) in zip(lines(only(page.axes)),[3,1])
        @test first.(curve[1][])≈f[2:3] rtol=8eps(Float32)
        @test last.(curve[1][])≈abs.(voltage[row,mode,2:3]) rtol=8eps(Float32)
        @test curve.color[]==lines(scalar.axes[1])[row].color[]
    end
    styled=LineCableModels.plot(first(observed);ydata=(selected,),
        series_labels=("Armor","Core"),series_attributes=((color=:red,),(color=:blue,)),options...)
    push!(kept,styled)
    @test all(axis -> [curve.label[] for curve in lines(axis)]==["Armor","Core"],styled.axes)
    @test all(axis -> [curve.color[] for curve in lines(axis)]==Makie.to_color.([:red,:blue]),styled.axes)
    one=@observe (Ti,abs)[2,2,:]
    single_entry=LineCableModels.plot(first(observed);ydata=(one,),options...)
    push!(kept,single_entry)
    @test only(lines(only(single_entry.axes))).label[]=="Conductor screen"
    ranged=@observe (Tv,abs)[1:2,2:3,:]
    range_page=LineCableModels.plot(first(observed);ydata=(ranged,),options...)
    push!(kept,range_page)
    @test range_page.addon_state.panel_page.coordinates==((1,2),(1,3))
    @test last.(lines(range_page.axes[2])[2][1][])≈abs.(voltage[2,3,:]) rtol=8eps(Float32)

    filtered=LineCableModels.plot(observed;ydata=((Tv,abs),),problem=2,
        reference=first(observed),layout=(2,2),
        panel_titles=Dict((3,2)=>"Reference mode two"),
        panel_legends=Dict((3,2)=>:inside),options...)
    append!(kept,filtered)
    @test length(filtered)==2
    @test filtered[1].addon_state.panel_page.coordinates==((2,1),(2,2),(2,3))
    @test filtered[2].addon_state.panel_page.coordinates==((3,1),(3,2),(3,3))
    @test occursin("reference",lowercase(filtered[2].export_name))
    @test filtered[2].axes[2].title[]=="Reference mode two"
    @test haskey(filtered[2].panel_legends,(3,2))
    @test !haskey(filtered[1].panel_legends,(3,2))
    @test lines(filtered[1].axes[1])[1].color[]==lines(filtered[2].axes[1])[1].color[]
    assumptions=(selections=(Z=((:series,:given,(;)),),Y=((:shunt,:given,(;)),)),
        modal=(effective=(method=:default,),))
    equivalent=[observables(PropagationParameters(E.retain_gridpoint(modal,
        G.gridpoint_id(;source_id,problem_index=1,formulation_index=index);
        fields=assumptions)),requests) for index in 1:2]
    @test length(G.observation_groups(equivalent;request=(Tv,abs)))==1
    separate=LineCableModels.plot(equivalent;ydata=((Tv,abs),),options...)
    append!(kept,separate)
    @test length(separate)==2
    @test [page.addon_state.displayed_indices for page in separate]==[[1],[2]]

    # Invalid dimensions and positional overrides fail before backend activation.
    mixed=observables(first(segments),(gamma,requests...))
    @test_throws r"matrix observations" LineCableModels.plot(mixed;
        ydata=(gamma,(Tv,abs)),options...,backend=:invalid)
    @test_throws r"no completed comparison retains band" LineCableModels.plot(first(observed);
        ydata=(selected,),band=(20.0,40.0),options...,backend=:invalid)
    @test_throws r"different overlaid dimensions" LineCableModels.plot(mixed;
        ydata=(selected,(Ti,abs)),series_labels=("one","two"),options...,backend=:invalid)
    @test_throws r"different overlaid dimensions" LineCableModels.plot(mixed;
        ydata=(selected,(Ti,abs)),series_attributes=((color=:red,),(color=:blue,)),
        options...,backend=:invalid)
    foreach(page -> empty!(page.figure),kept)
end

@testitem "PlotBuilder / row overlays retain nonmodal quantities and uncertainty" tags=[:visual] begin
    using CairoMakie, Measurements
    f=[10.0,20.0,40.0]
    z=[complex(measurement(i+2j+k/10,.05i),measurement(j-i,.02))
       for i in 1:2,j in 1:2,k in eachindex(f)]
    source=LineParameters(z,z.*1e-6,f)
    observed=observables(source,(R,X);length_unit=:base,clip=false)
    page=LineCableModels.plot(observed;ydata=(R,),overlay=:rows,
        backend=:cairo,display_plot=false,open_export=false,errorbar_sampling=:all)
    @test length(page.axes)==2
    for (column,axis) in enumerate(page.axes)
        curves=filter(p -> p isa Makie.Lines,axis.scene.plots)
        bars=filter(p -> p isa Makie.Errorbars,axis.scene.plots)
        @test length(curves)==2 && length(bars)==2
        @test endswith(string(axis.title[]),"Conductor $column")
        for row in 1:2
            @test curves[row].label[]=="Conductor $row"
            @test last.(curves[row][1][])≈[row+2column+k/10 for k in eachindex(f)] rtol=8eps(Float32)
            @test getindex.(bars[row][1][],3)≈fill(.05row,3) rtol=8eps(Float32)
        end
    end
    @test !isempty(Makie.colorbuffer(page.figure))
    empty!(page.figure)

    # A band selects saved comparison samples; it is not a new plotting-side
    # frequency filter. The same retained selection works with row overlays.
    reference=LineCableModels.Engine.retain_gridpoint(source,LineCableModels.Grammar.gridpoint_id())
    comparisons=LineCableModels.Engine.compare(reference,source,[R];bands=((20.0,40.0),))
    retained=observables(source,(R,X);comparisons,length_unit=:base,clip=false)
    band_page=LineCableModels.plot(retained;ydata=(R,),overlay=:rows,band=(20.0,40.0),
        backend=:cairo,display_plot=false,controls=false,open_export=false)
    @test length(band_page.axes)==2
    for (column,axis) in enumerate(band_page.axes)
        curves=filter(p -> p isa Makie.Lines,axis.scene.plots)
        for row in 1:2
            @test first.(curves[row][1][])≈f[2:3] rtol=8eps(Float32)
            @test last.(curves[row][1][])≈[row+2column+k/10 for k in 2:3] rtol=8eps(Float32)
        end
    end
    empty!(band_page.figure)
end
