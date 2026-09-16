@testitem "UQ / scalar products preserve populations, units and eligibility" tags=[:unit] begin
    using Statistics, Measurements, DataFrames, Random
    using LineCableModels.Engine: compare, absolute_error, relative_error
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    UQ=LineCableModels.UQ
    f=[0.1, 50.0, 100.0, 1000.0, 1e6, 1e7]
    angular=reshape(2pi .* f, 1, 1, :)
    rng=Xoshiro(12)
    draws=randn(rng, 2, 2, length(f), 17)
    trials=(R = 1.0 .+ 0.1draws, L = 1e-3 .+ 1e-4draws,
        C = 1e-9 .+ 1e-10draws, G = 1e-6 .+ 1e-7draws)
    summaries=map(trials) do array
        [SampleSummary(vec(array[i, j, k, :])) for i in 1:2, j in 1:2, k in eachindex(f)]
    end
    means=map(array -> mean.(array), summaries)
    core=LineParameters(
        means.R .+ im .* angular .* means.L, means.G .+ im .* angular .* means.C, f;
        details = ComputationDetails(;coordinates = ["a", "b"],))
    core=LineCableModels.materialize(core,summaries)
    histograms=map(trials) do array
        [HistogramDensity(vec(array[i, j, k, :]); bins = 3)
         for i in 1:2, j in 1:2, k in eachindex(f)]
    end
    mc=MonteCarloResult(
        MonteCarlo(Formulation(); trials = 17, seed = 12,
            return_samples = true, return_histograms = true),
        [core], [summaries], [trials], [histograms], UInt64(12), UInt64[13], [17])
    expected=dropdims(std(trials.R .+ im .* reshape(angular, 1, 1, :, 1) .* trials.L; dims = 4); dims = 4)
    @test observe(mc, statistics, Z, std, 1)≈expected
    for stat in (mean, std, minimum, median, maximum,
        Base.Fix2(quantile, 0.05), Base.Fix2(quantile, 0.95))
        @test observe(
            mc, statistics, X, stat, 1)≈angular .* observe(mc, statistics, L, stat, 1)
        @test observe(
            mc, statistics, B, stat, 1)≈angular .* observe(mc, statistics, C, stat, 1)
    end
    pub=observables(mc, ((samples, R, 1), (LineCableModels.histograms, R, 1));
        clip = false, length_unit = :base)
    @test keys(pub.columns)==(:samples_R, :histograms_R)
    @test length(pub.metadata.observation_columns)==2
    @test only(mc)===core
    @test eachindex(mc)==Base.OneTo(1)
    restored=LineCableModels.ImportExport.deserialize_value(LineCableModels.ImportExport.serialize_value(mc))
    @test observe(restored, samples, R, 1)==trials.R
    @test NamedTuple(observe(
        restored, LineCableModels.histograms, R, 1, 1, 1, 1))==NamedTuple(histograms.R[1])
    @test observe(restored, statistics, C, Base.Fix2(quantile, 0.95),
        1)==observe(mc, statistics, C, Base.Fix2(quantile, 0.95), 1)
    @test UQ.confidence(restored, 1).mean_standard_error.R≈std.(summaries.R) ./ sqrt(17)
    @test UQ.confidence(restored, 1).spread_estimated
    retained=observables(mc,((statistics,R,mean,1),);
        frequency_unit=:kilo,length_unit=:kilo,clip=false)
    metadata=(basis=:pul,domain=:PhaseDomain,frequencies=f,port_order=["a","b"],
        formulation=NamedTuple(Formulation()),axes=nothing)
    detached_report=report(BenchmarkTableDefinition(((statistics,R,mean),);bands=(:all,)),
        (reference=(result=retained,metadata),candidate=(result=mc,metadata)))
    @test all(<(1e-14),skipmissing(detached_report.table.terms.absolute_rms))
    # A saved table may use kHz and ohms/km while its native peer uses Hz and
    # ohms/m. Equal column names must not silently combine different units.
    combined=detached_report.table.statistics
    contracts=LineCableModels.ReportBuilder.observation_columns(combined)
    @test LineCableModels.Units.label(contracts.R.unit)=="Ω/m"
    @test LineCableModels.Units.label(contracts.frequency.unit)=="Hz"
    for role in (:reference,:candidate)
        rows=filter(row -> row.role===role && row.statistic===:mean,combined)
        @test Set(rows.frequency)==Set(f)
        for row in eachrow(rows)
            k=only(findall(==(row.frequency),f))
            @test row.R≈means.R[row.row,row.column,k]
        end
    end
    @test retained.columns.frequency≈f[repeat(1:length(f);inner=4)] ./ 1000

    # Numerator prefixes matter too: the default statistics publication uses
    # mH and µF. Report values, unit labels and coordinate-specific mean errors
    # must agree without changing the retained moments.
    physical_report=report(BenchmarkTableDefinition(
        ((statistics,L,std),(statistics,C,mean));bands=(:all,)),(reference=mc,candidate=mc))
    physical=physical_report.table.statistics
    physical_contracts=LineCableModels.ReportBuilder.observation_columns(physical)
    @test LineCableModels.Units.label(physical_contracts.L.unit)=="H/m"
    @test LineCableModels.Units.label(physical_contracts.C.unit)=="F/m"
    for row in eachrow(filter(row -> row.statistic in (:mean,:std),physical))
        k=only(findall(==(row.frequency),f))
        stat=row.statistic===:mean ? mean : std
        @test row.L≈stat(summaries.L[row.row,row.column,k])
        @test row.C≈stat(summaries.C[row.row,row.column,k])
    end

    two=MonteCarloResult(
        mc.formulation, [core, core], [summaries, summaries], nothing, nothing,
        UInt64(12), UInt64[13, 14], [17, 17])
    @test_throws ArgumentError compare(two, two, (statistics, R, mean))
    @test_throws ArgumentError compare(two, two, (statistics, R, mean); pairing = (
        (1, 1), (2, 1)))
    errors=compare(two, two, (statistics, R, mean); pairing = ((1, 2), (2, 1)))
    @test length(errors)==2 &&
          all(error -> all(iszero, observe(error, absolute_error)), errors)
    report_two=report(
        BenchmarkTableDefinition(((statistics, R, mean), (statistics, R, std)); pairing = (
            (1, 2), (2, 1))),
        (reference = two, candidate = two))
    @test Set(report_two.table.terms.problem_index)==Set((1, 2))
    @test Set(report_two.table.overview.coverage.candidate_point)==Set((1,2))
    @test Set(zip(report_two.table.overview.coverage.point,
        report_two.table.overview.coverage.reference_point))==Set(((1,2),(2,1)))
    @test length(report_two.table.features)==4
    @test length(report_two.table.sampling.point)==4
    # Equal nominal endpoints and standard uncertainties can carry different
    # correlations. Coordinate identity and the original quantities survive
    # presentation-row deduplication through the full reporting stages.
    correlated_frequency=measurement(50.,.25)
    independent_frequency=measurement(50.,.25)
    uncertain_sources=map((correlated_frequency,independent_frequency)) do frequency
        LineParameters(fill(1.0+2im,1,1,1),fill(3e-6+4e-6im,1,1,1),[frequency];
            details=ComputationDetails(;coordinates=["core"],))
    end
    uncertain_pair=LinearErrorResult(LinearError(Formulation()),collect(uncertain_sources))
    correlated_report=report(BenchmarkTableDefinition(((statistics,R,mean),);
        bands=(:all,),pairing=((1,1),(2,2))),
        (reference=uncertain_pair,candidate=uncertain_pair))
    coverage=correlated_report.table.overview.coverage
    @test coverage.point==[1,2]
    @test uncertainty(coverage.first_Hz[1]-correlated_frequency)==0
    @test uncertainty(coverage.first_Hz[2]-independent_frequency)==0
    @test uncertainty(coverage.first_Hz[1]-coverage.first_Hz[2])≈sqrt(2)*.25
    # Multiple configurations stay explicit; selecting a comparison configuration
    # must not mislabel a source-owned MC sampling point or hide whole-call timings.
    for mime in (MIME"text/plain"(),MIME"text/html"())
        text=sprint(show,mime,report_two)
        @test occursin("configuration 1",text) && occursin("configuration 2",text)
        @test occursin("reference configuration 1",text) && occursin("reference configuration 2",text)
        @test occursin("R · mean",text)
        @test occursin("R · std (propagated uncertainty)",text)
        selected=sprint((io,value) -> show(io,mime,value;problem=2),report_two)
        @test occursin("R · mean · configuration 2",selected)
        @test !occursin("R · mean · configuration 1",selected)
        @test occursin("reference configuration 1",selected)
        @test occursin("MC sampling workload",text)
        @test occursin("CDF precision",text)
        @test !occursin("mean_standard_error",text)
    end
    percentiles=report(
        BenchmarkTableDefinition((
            (statistics, R, median), (statistics, R, Base.Fix2(quantile, 0.05)))),
        (reference = mc, candidate = mc))
    @test length(percentiles.table.features)==2
    @test Set(percentiles.table.terms.statistic)==Set((:median, :q05))
    @test Set(percentiles.table.statistics.statistic)==Set((
        :mean, :std, :min, :q05, :median, :q95, :max))
    @test Set(names(percentiles.table.statistics))==Set((
        "point", "frequency", "row", "column", "statistic",
        "R", "trials", "point_seed", "role", "estimator"))

    for T in (Float32, Float64, BigFloat), quantity in (R, B)

        frequency=T.(f)
        a=fill(T(1), 2, 2, length(f));
        b=copy(a)
        b[1, 2, 2]=T(1e-30)
        error=compare(a, b, quantity; frequencies = frequency, result_basis = :pul)
        @test ismissing(observe(error, relative_error)[1, 2])
        @test observe(error, absolute_error)[1, 2]>0
        @test eltype(observe(error, absolute_error))===Union{Missing, T}
        @test all(ismissing,
            observe(
                compare(a, b, quantity; frequencies = frequency,
                    result_basis = :pul, band = (2e7, 3e7)),
                relative_error))
        @test all(iszero,
            skipmissing(observe(
                compare(a, b, quantity; frequencies = frequency,
                    result_basis = :pul, band = (1e7, 1e7)),
                relative_error)))
    end
    # A physically nonzero LEP uncertainty remains visible and compared even
    # when its corresponding mean is zero; the two are separate requests.
    uncertain_core=LineParameters(
        complex.(measurement.(zeros(2, 2, 6), fill(1e-4, 2, 2, 6))),
        zeros(ComplexF64, 2, 2, 6), f)
    spread=LinearErrorResult(LinearError(Formulation()), [uncertain_core])
    @test all(iszero, observe(only(compare(spread, spread, (statistics, R, std))), relative_error))
    @test all(ismissing, observe(only(compare(spread, spread, (statistics, R, mean))), relative_error))
    lep=LinearErrorResult(LinearError(Formulation()), [core])
    @test_throws ArgumentError report(
        BenchmarkTableDefinition(((statistics, R, median),)), (
            reference = mc, candidate = lep))
    @test_throws ArgumentError compare(mc, lep, (statistics, R, Base.Fix2(quantile, 0.1)))
    @test_throws ArgumentError compare(mc, lep, (statistics, Z, median))

    @test @inferred(observe(mc, statistics, R, mean, 1)) == means.R
    @test @inferred(observe(mc, statistics, B, std, 1)) ≈ angular .* std.(summaries.C)
    @test @inferred(observe(mc, statistics, Z, std, 1)) ≈ expected
    # Mean/std publication must not visit or copy the optional trial cloud.
    large_samples=map(array -> repeat(array; outer = (1, 1, 1, 100)), trials)
    large_summaries=map(summaries) do array
        map(
            summary -> SampleSummary(summary.mean, summary.std, summary.min, summary.q05,
                summary.median, summary.q95, summary.max, 1700),
            array)
    end
    large=MonteCarloResult(
        MonteCarlo(Formulation(); trials = 1700, seed = 12, return_samples = true),
        [core], [large_summaries], [large_samples], nothing, UInt64(12), UInt64[13], [1700])
    definition=BenchmarkTableDefinition(((statistics, R, mean), (statistics, R, std)); bands = (:all,))
    report(definition, (reference = mc, candidate = mc))
    report(definition, (reference = large, candidate = large))
    small_bytes=@allocated report(definition, (reference = mc, candidate = mc))
    large_bytes=@allocated report(definition, (reference = large, candidate = large))
    @test large_bytes <= small_bytes+16_384
end
