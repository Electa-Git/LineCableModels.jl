@testitem "Makie addons / observed UQ overlays retain uncertainty and reference styles" tags=[:visual] begin
    using CairoMakie,Measurements,Statistics
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    using LineCableModels.Engine: retain_gridpoint
    using LineCableModels.Grammar: gridpoint_id
    f=collect(range(2.,8.;length=7))
    omega=reshape(2pi.*f,1,1,:)
    include(joinpath(pkgdir(LineCableModels),"test/support/scenarios.jl"))
    parts=NamedTuple{(:R,:L,:C,:G)}(Tuple([CurrentScenarios.channel_value(Val(q),i,j,k)
        for i in 1:2,j in 1:2,k in eachindex(f)] for q in (:R,:L,:C,:G)))
    summaries=map(a -> map(x -> SampleSummary([.9x,x,1.1x]),a),parts)
    core=LineParameters(complex.(parts.R,omega.*parts.L),complex.(parts.G,omega.*parts.C),f)
    mc=MonteCarloResult(MonteCarlo(Formulation();trials=3,seed=1),
        [LineCableModels.materialize(core,summaries)],[summaries],nothing,nothing,UInt64(1),UInt64[2],[3])
    measured=map(a -> measurement.(1.05 .* a,0.03 .* a),parts)
    source_id=gridpoint_id().source_id
    candidates=[retain_gridpoint(LineParameters(factor.*complex.(measured.R,omega.*measured.L),
        factor.*complex.(measured.G,omega.*measured.C),f),gridpoint_id(;source_id,problem_index=index);
        fields=(inputs=(temperature=20index,),)) for (index,factor) in enumerate((1.,1.2))]
    lep=LinearErrorResult(LinearError(Formulation()),candidates)
    requests=(R,X,G,B,(statistics,R,mean),(statistics,R,std))
    artifact=report(BenchmarkTableDefinition(((statistics,R,mean),(statistics,R,std));bands=(:all,)),
        (reference=mc,candidate=lep,context=(id=:benchmark_uq_title_probe,));requests,
        observation_options=(length_unit=:base,))
    options=(backend=:cairo,display_plot=false,open_export=false)
    page=LineCableModels.plot(artifact;ydata=((R,2,1,:),),problem=2,options...,
        axis=(limits=((2.,8.),(.0005,.003)),),linewidth=3)
    axis=only(page.axes)
    @test axis.xlabelvisible[] && axis.xticklabelsvisible[]
    @test page.export_name=="benchmark_uq_title_probe — Series resistance"
    curves=filter(p -> p isa Makie.Lines,axis.scene.plots)
    bars=filter(p -> p isa Makie.Errorbars,axis.scene.plots)
    @test length(curves)==length(bars)==2
    @test Makie.to_color(last(curves).color[])==Makie.to_color(:black)
    @test last(bars).whiskerwidth[] > first(bars).whiskerwidth[] > 0
    for (curve,bar,source) in zip(curves,bars,(lep[2],mc[1]))
        @test last.(curve[1][])≈nominal.(R(source)[2,1,:])
        indices=[findfirst(==(point[1]),f) for point in bar[1][]]
        @test getindex.(bar[1][],2)≈nominal.(R(source)[2,1,indices])
        @test getindex.(bar[1][],3)≈uncertainty.(R(source)[2,1,indices])
        group=only(filter(handles -> curve in handles,collect(values(page.addon_state.groups))))
        marker=only(filter(p -> p isa Makie.Scatter,group))
        @test isdisjoint(first.(marker[1][]),first.(bar[1][]))
        @test curve.linewidth[]==bar.linewidth[]==3
    end
    page.controls[:ylog].active[]=true
    @test !isempty(Makie.colorbuffer(page.figure))
    @test axis.yscale[]===log10
    @test all(x -> x isa AbstractString,axis.yaxis.ticklabels[])
    full=LineCableModels.plot(artifact;ydata=((R,2,1,:),),problem=2,options...,errorbar_sampling=:all)
    for key in full.addon_state.order
        group=full.addon_state.groups[key]
        @test length(only(filter(p -> p isa Makie.Errorbars,group))[1][])==length(f)
        @test isempty(only(filter(p -> p isa Makie.Scatter,group))[1][])
    end
    # Persisted observations retain the same uncertainty curves without sources.
    mktempdir() do directory
        path=LineCableModels.save(artifact,joinpath(directory,"observed.jls"))
        restored=import_data(:observed,path)
        for point in candidates;Z(point).=NaN;end
        reloaded=report(BenchmarkTableDefinition(),restored.observed;reference=restored.reference)
        other=LineCableModels.plot(reloaded;ydata=((R,2,1,:),),problem=2,options...)
        other_curves=filter(p -> p isa Makie.Lines,only(other.axes).scene.plots)
        for (left,right) in zip(curves,other_curves)
            @test left[1][]==right[1][]
        end
    end
end

@testitem "Makie addons / detached UQ statistics retain matrix coordinates and bands" tags=[:visual] begin
    using CairoMakie,Statistics,Measurements,DataFrames
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f=10. .^ range(-1,7;length=13)
    omega=reshape(2pi.*f,1,1,:)
    include(joinpath(pkgdir(LineCableModels),"test/support/scenarios.jl"))
    parts=NamedTuple{(:R,:L,:C,:G)}(Tuple([CurrentScenarios.channel_value(Val(q),i,j,k)
        for i in 1:3,j in 1:3,k in eachindex(f)] for q in (:R,:L,:C,:G)))
    stats=map(a -> map(x -> SampleSummary([.9x,1.1x]),a),parts)
    core=LineParameters(complex.(parts.R,omega.*parts.L),complex.(parts.G,omega.*parts.C),f)
    reference=MonteCarloResult(MonteCarlo(Formulation();trials=2,seed=7),
        [LineCableModels.materialize(core,stats)],[stats],nothing,nothing,UInt64(7),UInt64[8],[2])
    measured=map(a -> measurement.(a,sqrt(2) * 0.1 .* a),parts)
    candidate=LinearErrorResult(LinearError(Formulation()),[LineParameters(
        complex.(measured.R,omega.*measured.L),complex.(measured.G,omega.*measured.C),measurement.(f,0.))])
    requests=((statistics,R,mean),(statistics,R,std),(statistics,B,mean))
    artifact=report(BenchmarkTableDefinition(requests;bands=(:all,:wide)),
        (reference,candidate,context=(id=:benchmark_uq_statistics,));observation_options=(length_unit=:base,))
    options=(backend=:cairo,display_plot=false,controls=false,open_export=false)
    pages=LineCableModels.plot(artifact;ydata=requests,layout=(2,2),options...,fig_size=(1100,750))
    @test length(pages)==12
    @test [length(page.axes) for page in pages]==repeat([4,2,2,1],3)
    @test first(pages).export_name=="benchmark_uq_statistics — Series resistance · mean"
    @test occursin("std",pages[5].export_name)
    for (index,page) in enumerate(pages)
        @test !isempty(Makie.colorbuffer(page.figure))
        request=requests[cld(index,4)]
        expected=observe(reference,statistics,request[2],request[3],1)
        for ((i,j),panel) in pairs(page.addon_state.panel_data)
            curves=filter(p -> p isa Makie.Lines,panel.axis.scene.plots)
            @test length(curves)==2
            @test first(curves)[1][]≈last(curves)[1][]
            @test last.(last(curves)[1][])≈expected[i,j,:]
            group=panel.groups[last(page.addon_state.order)]
            markers=only(filter(p -> p isa Makie.Scatter,group))
            @test last(markers[1][])[1]≈last(f)
        end
    end
    @test eltype(frequencies(only(candidate)))<:Measurement
    @test nrow(artifact.tables.statistics)>0
    selection=(statistics,R,mean,[3,1],[2],2:2:12)
    subset=LineCableModels.plot(artifact;ydata=(selection,),options...)
    @test length(subset.axes)==2
    for panel in values(subset.addon_state.panel_data),curve in filter(p -> p isa Makie.Lines,panel.axis.scene.plots)
        @test first.(curve[1][])≈f[2:2:12]
    end
    @test_throws ArgumentError LineCableModels.plot(artifact;ydata=(R,),options...)
    @test any(label -> occursin("empirical",label),values(first(pages).addon_state.labels))
    @test any(label -> occursin("first_order",label),values(first(pages).addon_state.labels))
    selected=LineCableModels.plot(artifact;ydata=requests[1:2],layout=(2,2),band=:wide,options...)
    @test length(selected)==8
    for page in selected,panel in values(page.addon_state.panel_data),curve in filter(p -> p isa Makie.Lines,panel.axis.scene.plots)
        @test all(>(1e6),first.(curve[1][]))
        @test last(curve[1][])[1]≈last(f)
    end
end
