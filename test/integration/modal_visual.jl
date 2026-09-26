@testitem "ModalAnalysis / retained finite quantities render with completed formulations" tags=[:visual] begin
    using CairoMakie
    using LineCableModels.Engine: retain_gridpoint, completed_formulation
    using LineCableModels.Grammar: gridpoint_id, observation_labels

    z=reshape(ComplexF64[1+im,2+im,3+im],1,1,3)
    y=reshape(ComplexF64[1e-6im,2e-6im,3e-6im],1,1,3)
    phase=retain_gridpoint(LineParameters(z,y,[10.0,50.0,100.0]),gridpoint_id();
        fields=merge(completed_formulation(Formulation()),
            (inputs=(system=(line_length=10.0,),),coordinates=["core"])))
    modal=compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation(:default))
    segment=PropagationParameters(modal)
    observed=ObservedResult(segment,((H,abs),(Tv,abs));complete_pairs=true)
    @test observed.gridpoint.source_gridpoint==details(segment).data.source_gridpoint
    pages=LineCableModels.plot(observed; ydata=((H,abs),(Tv,abs)),
        backend=:cairo,display_plot=false,controls=false,open_export=false)
    pages=pages isa AbstractVector ? pages : [pages]
    @test !isempty(pages)
    @test all(page -> !isempty(page.axes),pages)

    alternate=compute(ModalAnalysisProblem(phase),ModalAnalysisFormulation(:default;
        options=(iteration=(max_iterations=101,),)))
    primary=[ObservedResult(value,(R,X)) for value in (modal,alternate)]
    labels=observation_labels(primary;request=R)
    @test labels[1]!=labels[2]
    primary_pages=LineCableModels.plot(primary;ydata=(R,),backend=:cairo,
        display_plot=false,controls=false,open_export=false)
    primary_pages=primary_pages isa AbstractVector ? primary_pages : [primary_pages]
    @test any(page -> !isempty(page.axes),primary_pages)
end

@testitem "PlotBuilder / full transformation overlays retain 324 coefficients" tags=[:visual] begin
    using CairoMakie
    n=18
    f=10.0 .^ range(1,5;length=321)
    voltage=ComplexF64[complex(i/20,j/30+k/10) for i in 1:n,j in 1:n,k in eachindex(f)]
    current=ComplexF64[complex(j/25,i/35+k/10) for i in 1:n,j in 1:n,k in eachindex(f)]
    roots=fill(0.01+0.02im,n,length(f))
    modal=LineParameters(ModalDomain(ModalOperators(voltage,current),roots),
        SeriesImpedance(ones(ComplexF64,n,n,length(f))),
        ShuntAdmittance(ones(ComplexF64,n,n,length(f))),f,
        ComputationDetails((inputs=(system=(line_length=20.0,),),)))
    observed=observables(PropagationParameters(modal),((Tv,abs),(Tv,angle),(Ti,abs),(Ti,angle)))
    pages=LineCableModels.plot(observed;overlay=:coordinates,
        backend=:cairo,display_plot=false,controls=true,open_export=false,
        series_attributes=(linewidth=2.5,))
    @test length(pages)==4
    for (page, expected) in zip(pages,(abs.(voltage),rad2deg.(angle.(voltage)),
            abs.(current),rad2deg.(angle.(current))))
        @test length(page.axes)==1
        curves=filter(plot->plot isa Makie.Lines,only(page.axes).scene.plots)
        @test length(curves)==n*n
        entries=last(only(page.legend.entrygroups[]))
        @test 1<length(entries)<n*n
        @test last(entries).label[]=="(...)"
        @test page.legend.layoutobservables.autosize[][2]<=0.5only(page.axes).scene.viewport[].widths[2]+1
        @test only(page.axes).scene.viewport[].widths[2]>100
        for i in 1:n,j in 1:n
            curve=curves[(i-1)*n+j]
            @test curve.label[]=="Conductor $i, Mode $j"
            @test getindex.(curve[1][],1)≈f rtol=8eps(Float32)
            @test getindex.(curve[1][],2)≈expected[i,j,:] rtol=8eps(Float32)
            @test curve.linewidth[]==2.5
        end
        @test !isempty(Makie.colorbuffer(page.figure))
        empty!(page.figure)
    end
end

@testitem "PlotBuilder / multimode vector orientation and physical matrix coefficients" tags=[:visual] begin
    using CairoMakie
    using LinearAlgebra
    import LineCableModels.Engine as E
    import LineCableModels.Grammar as G

    f=[10.0,20.0,40.0]
    roots=ComplexF64[1+2im 3+4im 5+6im;2+3im 4+5im 6+7im].*1e-3
    voltage=cat([ComplexF64[1 .2; .3 1] for _ in f]...;dims=3)
    current=cat([ComplexF64[1 .1; .4 1] for _ in f]...;dims=3)
    z=zeros(ComplexF64,2,2,3)
    y=copy(z)
    for k in eachindex(f), mode in 1:2
        z[mode,mode,k]=roots[mode,k]
        y[mode,mode,k]=roots[mode,k]
    end
    modal=LineParameters(E.ModalDomain(ModalOperators(voltage,current),roots),
        SeriesImpedance(z),ShuntAdmittance(y),f,
        ComputationDetails((inputs=(system=(line_length=20.0,),),
            phase_coordinates=["core","sheath"])))
    first_segment=PropagationParameters(modal)
    second_segment=PropagationParameters(first_segment;line_length=30.0)
    first_point=ObservedResult(first_segment,(gamma,Tv,velocity))
    second_point=ObservedResult(second_segment,(gamma,Tv,velocity))
    options=(backend=:cairo,display_plot=false,controls=true,open_export=false)
    pages(value)=value isa AbstractVector ? value : [value]
    native_lines(axis)=filter(plot -> plot isa Makie.Lines,axis.scene.plots)

    automatic=pages(LineCableModels.plot(first_point;ydata=(gamma,),options...))
    @test length(automatic)==2
    @test all(page -> length(page.axes)==1 && length(native_lines(only(page.axes)))==2,automatic)
    @test all(page -> page.addon_state.panel_page.coordinates==(1,),automatic)
    @test "Mode 1" in values(first(automatic).addon_state.labels)
    @test !isempty(Makie.colorbuffer(first(automatic).figure))
    one_cell=pages(LineCableModels.plot(first_point;ydata=(gamma,),layout=(1,1),options...))
    @test all(page -> length(page.axes)==1 && length(native_lines(only(page.axes)))==2,
        one_cell)
    scalar=pages(LineCableModels.plot(first_point;ydata=(gamma,),options...))
    collection=pages(LineCableModels.plot([first_point];ydata=(gamma,),options...))
    @test [page.addon_state.panel_page.coordinates for page in scalar]==
        [page.addon_state.panel_page.coordinates for page in collection]
    raw_alpha=LineCableModels.plot(first_segment;ydata=(alpha,),options...)
    @test length(raw_alpha.axes)==1 && length(native_lines(only(raw_alpha.axes)))==2
    singleton_alpha=LineCableModels.plot((first_segment,);ydata=(alpha,),options...)
    @test singleton_alpha isa UIPlot
    @test [q.values for q in only(singleton_alpha.addon_state.observed).quantities]==
        [q.values for q in only(raw_alpha.addon_state.observed).quantities]
    @test length(singleton_alpha.axes)==length(raw_alpha.axes) &&
        length(native_lines(only(singleton_alpha.axes)))==length(native_lines(only(raw_alpha.axes)))
    retained_alpha=LineCableModels.plot(ObservedResult(first_segment);ydata=alpha,options...)
    @test length(retained_alpha.axes)==1 && length(native_lines(only(retained_alpha.axes)))==2
    real_only=LineCableModels.plot(first_segment;ydata=((gamma,real),),options...)
    @test real_only isa UIPlot && length(native_lines(only(real_only.axes)))==2

    explicit=pages(LineCableModels.plot(first_point;ydata=(gamma,),layout=(2,1),options...))
    @test length(explicit)==2
    @test all(page -> length(page.axes)==2 && all(axis -> length(native_lines(axis))==1,page.axes),explicit)
    @test all(page -> page.addon_state.panel_page.coordinates==(1,2),explicit)
    multi=pages(LineCableModels.plot([first_point,second_point];ydata=(gamma,),options...))
    @test length(multi)==2
    @test all(page -> length(page.axes)==2 && all(axis -> length(native_lines(axis))==2,page.axes),multi)
    duplicates=pages(LineCableModels.plot([first_point,first_point];ydata=(gamma,),options...))
    @test length(duplicates)==2
    @test all(page -> length(page.axes)==2 &&
        all(axis -> length(native_lines(axis))==2,page.axes),duplicates)
    @test length(G.observation_groups([first_point,first_point];request=(gamma,real)))==2
    source_id=G.gridpoint_id().source_id
    assumptions=(selections=(Z=((:series,:given,(;)),),Y=((:shunt,:given,(;)),)),
        modal=(effective=(method=:default,),))
    equivalent=[ObservedResult(PropagationParameters(E.retain_gridpoint(modal,
        G.gridpoint_id(;source_id,problem_index=1,formulation_index=index);
        fields=assumptions)),(gamma,)) for index in 1:2]
    @test length(G.observation_groups(equivalent;request=(gamma,real)))==1
    grouped_auto=LineCableModels.plot(equivalent;ydata=(gamma,real),options...)
    @test length(grouped_auto.axes)==2
    @test all(axis -> length(native_lines(axis))==1,grouped_auto.axes)
    @test grouped_auto.addon_state.displayed_indices==[1]
    grouped_coordinates=LineCableModels.plot(equivalent;ydata=(gamma,real),
        overlay=:coordinates,options...)
    @test length(grouped_coordinates.axes)==2
    @test all(axis -> length(native_lines(axis))==2,grouped_coordinates.axes)
    @test grouped_coordinates.addon_state.displayed_indices==[1,2]
    small=pages(LineCableModels.plot([first_point,second_point];ydata=(gamma,),
        layout=(1,1),options...))
    @test length(small)==4
    @test [only(page.addon_state.panel_page.coordinates) for page in small]==[1,2,1,2]
    with_reference=pages(LineCableModels.plot(first_point;ydata=(gamma,),
        reference=second_point,options...))
    @test length(with_reference)==2
    @test all(page -> length(page.axes)==2 && all(axis -> length(native_lines(axis))==2,page.axes),with_reference)
    @test all(page -> any(occursin("reference",string(label)) for label in values(page.addon_state.labels)),with_reference)

    coordinate_pages=pages(LineCableModels.plot([first_point,second_point];
        ydata=(gamma,),overlay=:coordinates,layout=(1,1),options...))
    @test length(coordinate_pages)==4
    @test [only(page.addon_state.panel_page.coordinates) for page in coordinate_pages]==[1,2,1,2]
    @test all(page -> length(native_lines(only(page.axes)))==2,coordinate_pages)
    @test native_lines(only(coordinate_pages[1].axes))[1].color[]==
        native_lines(only(coordinate_pages[2].axes))[1].color[]
    second_mode=@observe (gamma,real)[2,:]
    selected_mode=LineCableModels.plot(first_point;ydata=(second_mode,),
        series_attributes=((color=:red,),),options...)
    selected_line=only(native_lines(only(selected_mode.axes)))
    @test selected_line.label[]=="Mode 2"
    @test last.(selected_line[1][])≈[2,4,6]
    @test selected_line.color[]==Makie.to_color(:red)
    reordered_modes=@observe (gamma,real)[[2,1],:]
    reordered=LineCableModels.plot(first_point;ydata=(reordered_modes,),
        series_attributes=((color=:red,),(color=:blue,)),options...)
    reordered_lines=native_lines(only(reordered.axes))
    @test [line.label[] for line in reordered_lines]==["Mode 2","Mode 1"]
    @test [last.(line[1][]) for line in reordered_lines]≈[[2,4,6],[1,3,5]]
    @test [line.color[] for line in reordered_lines]==collect(Makie.to_color.((:red,:blue)))
    reordered_default=LineCableModels.plot(first_point;ydata=(reordered_modes,),options...)
    @test [line.color[] for line in native_lines(only(reordered_default.axes))]==
        reverse([line.color[] for line in native_lines(only(automatic[1].axes))])
    relabelled=LineCableModels.plot(first_point;ydata=(reordered_modes,),
        series_labels=("Selected 2","Selected 1"),options...)
    @test [line.label[] for line in native_lines(only(relabelled.axes))]==
        ["Selected 2","Selected 1"]
    coordinate_reference=pages(LineCableModels.plot(first_point;ydata=(gamma,),
        reference=second_point,overlay=:coordinates,layout=(1,2),options...))
    @test length(coordinate_reference)==2
    @test all(page -> length(page.axes)==2 && all(axis -> length(native_lines(axis))==2,page.axes),coordinate_reference)
    @test all(page -> native_lines(page.axes[1])[1].color[]==
        native_lines(page.axes[2])[1].color[],coordinate_reference)
    distinct=[ObservedResult(PropagationParameters(E.retain_gridpoint(modal,
        G.gridpoint_id(;source_id,problem_index=index,formulation_index=1))),
        (gamma,)) for index in 1:3]
    stable_reference=LineCableModels.plot(distinct;ydata=(gamma,real),
        overlay=:coordinates,reference=second_point,problem=1,
        panel_legends=Dict(4=>:inside),options...)
    @test stable_reference.addon_state.panel_page.coordinates==(1,4)
    @test length(stable_reference.axes)==2
    filtered_modal=E.retain_gridpoint(modal,G.gridpoint_id(;problem_index=2))
    filtered_point=ObservedResult(PropagationParameters(filtered_modal),(gamma,))
    filtered=pages(LineCableModels.plot([first_point,filtered_point];ydata=(gamma,),
        problem=2,options...))
    @test length(filtered)==2
    @test all(page -> length(page.axes)==1 && length(native_lines(only(page.axes)))==2,filtered)
    gridpoint_pages=pages(LineCableModels.plot(first_point;ydata=(gamma,),
        overlay=:gridpoints,layout=(1,1),options...))
    @test length(gridpoint_pages)==4
    @test all(page -> length(native_lines(only(page.axes)))==1,gridpoint_pages)

    matrix_pages=pages(LineCableModels.plot(first_point;ydata=(Tv,),options...))
    @test length(matrix_pages)==2
    @test all(page -> length(page.axes)==4,matrix_pages)
    @test all(page -> Set(page.addon_state.panel_page.coordinates)==
        Set(((1,1),(1,2),(2,1),(2,2))),matrix_pages)
    @test all(axis -> startswith(string(axis.title[]),"Re(Tv)") &&
        occursin("Re(Tv)",string(axis.ylabel[])),
        first(matrix_pages).axes)
    @test any(page -> any(axis -> any(curve -> any(x -> !iszero(x),last.(curve[1][])),native_lines(axis)),page.axes),matrix_pages)
    matrix_coordinates=pages(LineCableModels.plot([first_point,second_point];
        ydata=((Tv,real),),overlay=:coordinates,layout=(1,2),options...))
    @test length(matrix_coordinates)==1
    @test length(only(matrix_coordinates).axes)==2
    @test all(axis -> length(native_lines(axis))==4,only(matrix_coordinates).axes)
    @test all(axis -> startswith(axis.title[],"Re(Tv)"),only(matrix_coordinates).axes)
    @test all(label -> occursin("Conductor",string(label)) && occursin("Mode",string(label)),
        values(only(matrix_coordinates).addon_state.labels))
    matrix_coefficient=@observe (Tv,real)[1,2,:]
    selected_matrix=LineCableModels.plot(first_point;ydata=(matrix_coefficient,),
        overlay=:coordinates,series_attributes=((color=:green,),),options...)
    selected_matrix_line=only(native_lines(only(selected_matrix.axes)))
    @test occursin("core",selected_matrix_line.label[])
    @test occursin("Mode 2",selected_matrix_line.label[])
    @test last.(selected_matrix_line[1][])≈fill(0.2,3)
    @test selected_matrix_line.color[]==Makie.to_color(:green)
    mixed=pages(LineCableModels.plot(first_point;ydata=(gamma,(Tv,real)),options...))
    @test [page.addon_state.nominal_capacity for page in mixed]==[(1,1),(1,1),(2,2)]
    mixed_explicit=pages(LineCableModels.plot(first_point;
        ydata=(gamma,(Tv,real)),layout=(2,2),options...))
    @test all(page -> page.addon_state.nominal_capacity==(2,2),mixed_explicit)
    @test_throws ArgumentError LineCableModels.plot(first_point;
        ydata=(gamma,(Tv,real)),series_labels=("Mode 1","Mode 2"),options...)
    @test_throws ArgumentError LineCableModels.plot(first_point;
        ydata=(gamma,(Tv,real)),series_attributes=((color=:red,),(color=:blue,)),options...)
    styled=LineCableModels.plot(first_point;ydata=(gamma,real),
        overlay=:coordinates,series_labels=("Mode A","Mode B"),
        series_attributes=((color=:red,),(linestyle=:dash,)),options...)
    @test length(native_lines(only(styled.axes)))==2
    @test native_lines(only(styled.axes))[1].color[]==Makie.to_color(:red)
    @test native_lines(only(styled.axes))[2].linestyle[]==Makie.to_linestyle(:dash)
    paneltitle!(styled,1,"Modes")
    @test only(styled.axes).title[]=="Modes"
    resetview!(styled)
    resize!(styled.figure,900,500)
    @test !isempty(Makie.colorbuffer(styled.figure))
    mktempdir() do directory
        path=export_svg(styled;path=joinpath(directory,"modal-overlay.svg"),open_file=false)
        @test filesize(path)>0
    end
end
