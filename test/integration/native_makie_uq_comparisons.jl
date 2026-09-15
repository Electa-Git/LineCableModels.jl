@testitem "Makie addons / ordinary UQ overlays delegate uncertainty to its owner" tags=[:visual] begin
    using CairoMakie, Measurements, Statistics
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f = collect(range(2.0,8.0;length=7))
    omega = reshape(2pi.*f,1,1,:)
    include(joinpath(pkgdir(LineCableModels),"test/support/scenarios.jl"))
    parts=NamedTuple{(:R,:L,:C,:G)}(Tuple([CurrentScenarios.channel_value(Val(q),i,j,k)
        for i in 1:2,j in 1:2,k in eachindex(f)] for q in (:R,:L,:C,:G)))
    summaries = map(a -> map(x -> SampleSummary([0.9x,x,1.1x]),a),parts)
    core = LineParameters(complex.(parts.R,omega.*parts.L),complex.(parts.G,omega.*parts.C),f)
    second_core = LineParameters(1.2.*Z(core),1.2.*Y(core),f)
    second_summaries = map(a -> map(x -> SampleSummary(1.2.*[0.9x,x,1.1x]),a),parts)
    mc = MonteCarloResult(MonteCarlo(Formulation();trials=3,seed=1),
        [core,second_core],[summaries,second_summaries],nothing,nothing,UInt64(1),UInt64[1,2],[3,3])
    measured = map(a -> measurement.(1.05.*a,0.03.*a),parts)
    candidate = LineParameters(complex.(measured.R,omega.*measured.L),
        complex.(measured.G,omega.*measured.C),f)
    second_candidate = LineParameters(1.2.*Z(candidate),1.2.*Y(candidate),f)
    lep = LinearErrorResult(LinearError(Formulation()),[candidate,second_candidate])
    definition = BenchmarkTableDefinition(((statistics,R,mean),(statistics,R,std));
        bands=(:all,),pairing=((1,1),(2,2)))
    metadata = (port_order=["a","b"],)
    artifact = report(definition,(reference=(result=mc,metadata=metadata),
        candidate=(result=lep,metadata=metadata),
        context=(id=:benchmark_uq_title_probe,case_id=:uq_title_probe,collection=:test)))
    before = deepcopy((Z(candidate),Y(candidate),frequencies(candidate)))
    page = LineCableModels.plot(artifact;ydata=((R,2,1,:),),problem=2,
        backend=:cairo,display_plot=false,open_export=false,length_unit=:base,
        axis=(limits=((2.0,8.0),(.0005,.0030)),),linewidth=3)
    axis = only(page.axes)
    lines = filter(p -> p isa Makie.Lines,axis.scene.plots)
    bars = filter(p -> p isa Makie.Errorbars,axis.scene.plots)
    @test length(lines)==length(bars)==2
    @test page.export_name == "benchmark_uq_title_probe — Series resistance"
    # Data markers must not compete with the uncertainty glyphs. Native
    # Errorbars still own their cap geometry (which can itself use Scatter).
    @test !any(p -> p isa Makie.Scatter,axis.scene.plots)
    @test bars[1].whiskerwidth[] > bars[2].whiskerwidth[] > 0
    for (line,bar,source) in zip(lines,bars,(uncertain(mc,2),uncertain(lep,2)))
        @test last.(line[1][]) ≈ nominal.(R(source)[2,1,:])
        # Errorbars encode x,y,negative error,positive error in native Point4.
        @test getindex.(bar[1][],2) ≈ nominal.(R(source)[2,1,:])
        @test getindex.(bar[1][],3) ≈ uncertainty.(R(source)[2,1,:])
        @test line.linewidth[] == 3
        @test bar.linewidth[] == 3 # Explicit native style wins over nested defaults.
    end
    @test last.(lines[1][1][]) != last.(lines[2][1][])
    @test getindex.(bars[1][1][],3) != getindex.(bars[2][1][],3)
    page.controls[:ylog].active[]=true
    Makie.colorbuffer(page.figure)
    @test axis.yscale[] === log10
    @test all(x -> x isa AbstractString,axis.yaxis.ticklabels[])

    # A new result owner must enter the real comparison renderer without adding
    # a PlotBuilder type/name branch. Reuse retained comparisons, not new physics.
    struct PlotOwnerProbe{T,F} <: AbstractUncertaintyResult{T}
        values::Vector{T}
        formulation::F
    end
    Base.NamedTuple(x::PlotOwnerProbe) = (values=x.values,formulation=NamedTuple(x.formulation),details=(;))
    Base.length(x::PlotOwnerProbe) = length(x.values)
    Base.getindex(x::PlotOwnerProbe,i::Integer) = x.values[i]
    LineCableModels.uncertain(x::PlotOwnerProbe) = x.values
    probe = PlotOwnerProbe([candidate,second_candidate],LinearError(Formulation()))
    retained = merge(artifact.published,(candidate=merge(artifact.published.candidate,(result=probe,)),))
    # Use the owned saved-report path; the test result provides native cores via
    # the existing result/uncertain contracts, not a bespoke recipe.
    other = LineCableModels.plot(retained;ydata=((R,2,1,:),),problem=2,
        backend=:cairo,display_plot=false,open_export=false,length_unit=:base)
    other_lines = filter(p -> p isa Makie.Lines,only(other.axes).scene.plots)
    other_bars = filter(p -> p isa Makie.Errorbars,only(other.axes).scene.plots)
    @test !any(p -> p isa Makie.Scatter,only(other.axes).scene.plots)
    @test other_bars[1].linewidth[] > other_bars[2].linewidth[] > 0
    @test last.(last(other_lines)[1][]) ≈ nominal.(R(second_candidate)[2,1,:])
    @test isequal(before,(Z(candidate),Y(candidate),frequencies(candidate)))
    @test_throws ArgumentError LineCableModels.plot(artifact;problem=2,
        ydata=(R,(statistics,R,std)),backend=:cairo,display_plot=false)
end

@testitem "Makie addons / UQ comparisons use native statistics and matrix pages" tags=[:visual] begin
    using CairoMakie, Statistics, Measurements, DataFrames
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f=10.0 .^ range(-1, 7; length = 13)
    include(joinpath(pkgdir(LineCableModels),"test/support/scenarios.jl"))
    values=NamedTuple{(:R,:L,:C,:G)}(Tuple([CurrentScenarios.channel_value(Val(q),i,j,k)
        for i in 1:3,j in 1:3,k in eachindex(f)] for q in (:R,:L,:C,:G)))
    r=values.R
    omega=reshape(2pi .* f, 1, 1, :)
    stats=map(values) do array
        map(value -> SampleSummary([0.9value, 1.1value]), array)
    end
    core=LineParameters(
        values.R .+ im .* omega .* values.L, values.G .+
                                             im .* omega .* values.C, f)
    reference=MonteCarloResult(
        MonteCarlo(Formulation(); trials = 2, seed = 7), [core], [stats], nothing, nothing,
        UInt64(7), UInt64[8], [2])
    measured=map(values) do array
        measurement.(array, sqrt(2)*0.1 .* array)
    end
    candidate=LinearErrorResult(LinearError(Formulation()),
        [LineParameters(
            measured.R .+ im .* omega .* measured.L, measured.G .+
                                                     im .* omega .* measured.C, measurement.(f,0.0))])
    requests=((statistics, R, mean), (statistics, R, std), (statistics, B, mean))
    metadata=(port_order = ["a", "b", "c"],
        formulation = NamedTuple(LinearError(Formulation())), axes = nothing)
    artifact=report(BenchmarkTableDefinition(requests; bands = (:all,)),
        (reference = (result = reference, metadata = metadata),
            candidate = (result = candidate, metadata = metadata),
            context = (id=:benchmark_uq_statistics,case_id=:uq_statistics,collection=:test)))
    pages=LineCableModels.plot(
        artifact; backend = :cairo, ydata = requests, blocks = (2, 2),
        display_plot = false, controls = true, open_export = false, length_unit = :base, clip = false,
        fig_size = (1100, 750))
    @test length(pages)==12
    @test first(pages).export_name == "benchmark_uq_statistics — Series resistance · mean (1,1)"
    @test pages[5].export_name == "benchmark_uq_statistics — Series resistance · std (1,1)"
    @test [length(page.axes) for page in pages]==repeat([4, 2, 2, 1], 3)
    for (page_index,page) in enumerate(pages)
        @test page.figure.scene.viewport[].widths[1]>page.figure.scene.viewport[].widths[2]
        @test !isempty(Makie.colorbuffer(page.figure))
        request=requests[cld(page_index,4)]
        expected=observe(reference,statistics,request[2],request[3],1)
        for ((i,j),panel) in pairs(page.addon_state.panel_data)
            lines=filter(plot -> plot isa Makie.Lines, panel.axis.scene.plots)
            @test length(lines)==2
            @test first(lines)[1][]≈last(lines)[1][]
            # Matching curves alone would miss a shared mean/std or coordinate swap.
            @test last.(first(lines)[1][])≈expected[i,j,:]
            markers=filter(plot -> plot isa Makie.Scatter, panel.axis.scene.plots)
            @test last(first(markers)[1][])[1]≈last(f)
        end
    end
    @test occursin("mean", first(pages).export_name)
    @test occursin("std", pages[5].export_name)
    @test observe(candidate, statistics, R, std, 1)≈sqrt(2)*0.1 .* r
    # Saved LEP axes can be Measurement-typed even when frequency is exact.
    # Publication resolves its display cutoff without converting/mutating it.
    @test frequencies(only(candidate))==measurement.(f,0.0)
    @test eltype(frequencies(only(candidate)))<:Measurement
    @test size(artifact.table.statistics, 1)>0

    selected_request=@observe (statistics,R,mean)[1,[3,1],[2],2:2:12]
    subset=LineCableModels.plot(artifact; backend=:cairo,ydata=(selected_request,),
        display_plot=false,controls=false,open_export=false,length_unit=:base,
        clip=false,fig_size=(1100,750))
    @test length(subset.axes)==2
    for panel in Base.values(subset.addon_state.panel_data)
        for line in filter(plot -> plot isa Makie.Lines,panel.axis.scene.plots)
            @test [point[1] for point in line[1][]]≈f[2:2:12]
        end
    end
    @test_throws ArgumentError LineCableModels.plot(artifact; backend=:cairo,
        ydata=((statistics,R,mean,2,Colon(),Colon(),Colon()),),
        display_plot=false,controls=false,open_export=false)
    # Check the real composed output, not a retired private projection helper.
    @test Set(Base.values(first(pages).addon_state.labels))==Set(("Reference · Monte Carlo","LEP"))
    @test all(label -> !occursin("earth Y",label),Base.values(first(pages).addon_state.labels))

    # Composite children remain visible even when all candidates share them.
    # Distinct mean/std data above also guard this real wrapped path against a
    # shared coordinate or statistic swap hidden by agreeing legend strings.
    inner=Formulation(internal_impedance=(inner=:default,outer=:default,transfer=:default),
        earth_impedance=(air=:default,earth=:pollaczek1926,mixed=:default),
        earth_admittance=(air=:default,earth=:default,mixed=:default))
    composite_reference=MonteCarloResult(MonteCarlo(inner;trials=2,seed=7),
        [core],[stats],nothing,nothing,UInt64(7),UInt64[8],[2])
    composite_candidate=LinearErrorResult(LinearError(inner),collect(candidate))
    composite_report=report(BenchmarkTableDefinition(requests;bands=(:all,)),
        (reference=(result=composite_reference,metadata=(port_order=["a","b","c"],)),
            candidate=(result=composite_candidate,metadata=(port_order=["a","b","c"],))))
    composite_pages=LineCableModels.plot(composite_report;backend=:cairo,ydata=requests,
        display_plot=false,controls=false,open_export=false,length_unit=:base,clip=false)
    for (page,request) in zip(composite_pages,requests)
        suffix=request[2]===R ?
            "internal Z(inner)=default; internal Z(outer)=default; internal Z(transfer)=default; earth Z(air)=default; earth Z(earth)=Pollaczek; earth Z(mixed)=default" :
            "earth Y(air)=default; earth Y(earth)=default; earth Y(mixed)=default"
        names=[page.addon_state.labels[group] for group in page.addon_state.order]
        @test names==["Reference · Monte Carlo; "*suffix,"LEP; "*suffix]
        expected=observe(reference,statistics,request[2],request[3],1)
        for ((i,j),panel) in pairs(page.addon_state.panel_data)
            curves=filter(item->item isa Makie.Lines,panel.axis.scene.plots)
            @test last.(first(curves)[1][])≈expected[i,j,:]
        end
    end

    # A numerical publication with sufficient coordinates is not converted to
    # fake LineParameters or independent Measurements to reuse this renderer.
    retained=observables(reference, ((statistics, R, mean, 1), (statistics, R, std, 1));
        length_unit = :base, clip = false)
    context=merge(metadata, (basis = :pul, domain = :PhaseDomain, frequencies = f))
    detached_report=report(
        BenchmarkTableDefinition(((statistics, R, mean), (statistics, R, std)); bands = (
            :all, :wide)),
        (reference = (result = retained, metadata = context),
            candidate = (result = retained, metadata = context)))
    selected=LineCableModels.plot(
        detached_report; backend = :cairo, ydata = ((statistics,R,mean),(statistics,R,std)), blocks = (2, 2), band = :wide,
        display_plot = false, controls = false, open_export = false, length_unit = :base, fig_size = (
            1100, 750))
    @test length(selected)==8
    @test_throws r"uncertainty-bearing core" LineCableModels.plot(detached_report;
        backend=:cairo,ydata=(R,),display_plot=false)
    for page in selected
        @test !isempty(Makie.colorbuffer(page.figure))
        for panel in Base.values(page.addon_state.panel_data)
            for line in filter(plot -> plot isa Makie.Lines, panel.axis.scene.plots)
                @test all(point -> point[1]>1e6, line[1][])
                @test last(line[1][])[1]≈last(f)
            end
        end
    end
end
