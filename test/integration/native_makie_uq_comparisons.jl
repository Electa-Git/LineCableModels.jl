@testitem "Makie addons / UQ comparisons use native statistics and matrix pages" tags=[:visual] begin
    using CairoMakie, Statistics, Measurements, DataFrames
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f=10.0 .^ range(-1, 7; length = 13)
    r=[(i+j)*1e-3*(1+log10(1+x)) for i in 1:3, j in 1:3, x in f]
    omega=reshape(2pi .* f, 1, 1, :)
    values=(R = r, L = r .* 1e-4, C = r .* 1e-8, G = r .* 1e-5)
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
            candidate = (result = candidate, metadata = metadata)))
    pages=LineCableModels.plot(
        artifact; backend = :cairo, ydata = requests, blocks = (2, 2),
        display_plot = false, controls = true, open_export = false, length_unit = :base, clip = false,
        fig_size = (1100, 750))
    @test length(pages)==12
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
        earth_impedance=(air=:default,earth=:Pollaczek1926,mixed=:default),
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
            "internal Z(inner)=default; internal Z(outer)=default; internal Z(transfer)=default; earth Z(air)=default; earth Z(earth)=Pollaczek1926; earth Z(mixed)=default" :
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
    historical=report(
        BenchmarkTableDefinition(((statistics, R, mean), (statistics, R, std)); bands = (
            :all, :wide)),
        (reference = (result = retained, metadata = context),
            candidate = (result = retained, metadata = context)))
    selected=LineCableModels.plot(
        historical; backend = :cairo, ydata = (R,), blocks = (2, 2), band = :wide,
        display_plot = false, controls = false, open_export = false, length_unit = :base, fig_size = (
            1100, 750))
    @test length(selected)==8
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
