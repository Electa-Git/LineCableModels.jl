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
    retained=ObservedResult(mc,1,((samples,R),(LineCableModels.histograms,R,1,1,1));
        length_unit=:base)
    @test length(retained.quantities)==2
    @test first(retained.quantities).values==trials.R
    @test only(last(retained.quantities).coordinates.rows)==1
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
    retained=ObservedResult(mc,1,((statistics,R,mean),);frequency_unit=:kilo,length_unit=:kilo)
    @test first(retained.quantities).coordinates.frequencies≈f./1000
    @test observe(retained,statistics,R,mean)≈1000means.R
    physical_report=report(BenchmarkTableDefinition(((statistics,L,std),(statistics,C,mean))),
        (reference=mc,candidate=mc);observation_options=(length_unit=:base,quantity_units=:base))
    @test length(physical_report.tables.features)==2
    @test observe(only(physical_report.observed),statistics,L,std)≈std.(summaries.L)
    @test observe(only(physical_report.observed),statistics,C,mean)≈mean.(summaries.C)
    @test all(iszero,skipmissing(physical_report.tables.terms.relative_rms_percent))

    two=MonteCarloResult(
        mc.formulation, [core, core], [summaries, summaries], nothing, nothing,
        UInt64(12), UInt64[13, 14], [17, 17])
    @test_throws ArgumentError compare(two, two, (statistics, R, mean))
    @test_throws ArgumentError compare(two, two, (statistics, R, mean); pairing = (
        (1, 1), (2, 1)))
    errors=compare(two, two, (statistics, R, mean); pairing = ((1, 2), (2, 1)))
    @test length(errors)==2 &&
          all(error -> all(iszero, observe(error, absolute_error)), errors)
    @test_throws ArgumentError report(BenchmarkTableDefinition(((statistics,R,mean),);
        pairing=((1,2),(2,1))),(reference=two,candidate=two))
    percentiles=report(BenchmarkTableDefinition(((statistics,R,median),(statistics,R,Base.Fix2(quantile,.05)))),
        (reference=mc,candidate=mc))
    @test length(percentiles.tables.features)==2
    @test Set(percentiles.tables.terms.statistic)==Set((:median,Symbol("quantile_0.05")))
    @test Set(percentiles.tables.statistics.statistic)==Set((:median,Symbol("quantile_0.05")))

    for T in (Float32, Float64, BigFloat), quantity in (R, B)

        frequency=T.(f)
        a=fill(T(1), 2, 2, length(f));
        b=copy(a)
        b[1, 2, 2]=T(1e-30)
        error=compare(a, b, quantity; frequencies = frequency, result_basis = :pul)
        @test ismissing(observe(error, relative_error)[1, 2])
        @test ismissing(observe(error,absolute_error)[1,2])
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
