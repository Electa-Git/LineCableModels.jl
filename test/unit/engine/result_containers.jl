@testitem "Engine / result containers / numeric and selector behavior" tags=[:unit] setup=[
    UseEngineSupport, TestNumerics, TestFixtures] begin
    using DataFrames

    design=TestFixtures.coaxial_design()
    constants=CableConstants(design)
    @test constants isa CableConstants{Float64}
    @test all(>(0), constants.R)
    @test all(>(0), constants.L)
    @test all(>(0), constants.C)
    @test all(>=(0), constants.G)
    @test basis(constants) === :pul
    constants_display=sprint(show, constants)
    @test occursin("CableConstants(assemblies=1", constants_display)
    @test resistance(constants) === constants.R
    @test inductance(constants) === constants.L
    @test capacitance(constants) === constants.C
    @test conductance(constants) === constants.G
    @test observe(constants, R) === constants.R
    @test observables(typeof(constants)) == (R, L, C, G)
    constants_observables=observables(
        constants,
        (R, L, C, G)
    )
    @test constants_observables isa ObservedResult
    @test observe(constants_observables,R)≈1_000constants.R
    constants_table=LineCableModels.ReportBuilder.tabulate(constants_observables,R)
    @test names(constants_table)==["assembly","value"]
    @test only(constants_table.assembly)==1
    @test constants_table.value≈1_000constants.R
    native_constants_table=DataFrame(constants)
    @test names(native_constants_table) == ["core", "R", "L", "C", "G"]
    @test native_constants_table.R == constants.R

    frequency=[50.0, 100.0, 200.0]
    angular=reshape(2π .* frequency, 1, 1, :)
    resistance_values=reshape(collect(1.0:12.0), 2, 2, 3) .* 1.0e-4
    inductance_values=reshape(collect(13.0:24.0), 2, 2, 3) .* 1.0e-7
    conductance_values=reshape(collect(25.0:36.0), 2, 2, 3) .* 1.0e-8
    capacitance_values=reshape(collect(37.0:48.0), 2, 2, 3) .* 1.0e-10
    impedance=complex.(resistance_values, inductance_values .* angular)
    admittance=complex.(conductance_values, capacitance_values .* angular)
    parameters=LineParameters(
        impedance,
        admittance,
        frequency;
        basis = :total
    )

    @test basis(parameters) === :total
    @test frequencies(parameters) == frequency
    @test nconductors(parameters) == 2
    @test nfrequencies(parameters) == 3
    @test Z(parameters) == impedance
    @test Y(parameters) == admittance
    @test Z(parameters, 1, 2) == impedance[1, 2, :]
    @test Z(parameters, 1, 2, 2) == impedance[1, 2, 2]
    @test Y(parameters, 2, 1, 1:2) == admittance[2, 1, 1:2]
    @test R(parameters, 1, 2, :) == resistance_values[1, 2, :]
    @test X(parameters, 1, 2, 2) == imag(impedance[1, 2, 2])
    @test L(parameters, 1, 2, 2:3) ≈ inductance_values[1, 2, 2:3]
    @test G(parameters, 2, 1) == conductance_values[2, 1, :]
    @test B(parameters, 2, 1, 2) == imag(admittance[2, 1, 2])
    @test C(parameters, 2, 1, 1:2) ≈ capacitance_values[2, 1, 1:2]
    @test series_impedance(parameters) === parameters.Z.values
    @test shunt_admittance(parameters) === parameters.Y.values
    @test resistance(parameters, 1, 1, 1) == R(parameters, 1, 1, 1)
    @test reactance(parameters, 1, 1, 1) == X(parameters, 1, 1, 1)
    @test inductance(parameters, 1, 1, 1) ≈ L(parameters, 1, 1, 1)
    @test conductance(parameters, 1, 1, 1) == G(parameters, 1, 1, 1)
    @test susceptance(parameters, 1, 1, 1) == B(parameters, 1, 1, 1)
    @test capacitance(parameters, 1, 1, 1) ≈ C(parameters, 1, 1, 1)
    @test L(parameters) ≈ inductance_values
    @test C(parameters) ≈ capacitance_values
    @test observe(parameters, Z) === parameters.Z.values
    @test observe(parameters, Y) === parameters.Y.values
    parameter_observables=observables(parameters,((Z,abs,1,2,:),(Z,angle,1,2,:));
        atol=(R=0.,X=0.,G=0.,B=0.))
    @test parameter_observables isa ObservedResult
    @test first(parameter_observables.quantities).coordinates.frequencies==frequency
    @test observe(parameter_observables,Z,abs)≈abs.(impedance[1,2,:])
    @test observe(parameter_observables,Z,angle)≈rad2deg.(angle.(impedance[1,2,:]))
    @test first(parameter_observables.quantities).coordinates.frequencies !== parameters.f

    series=SeriesImpedance(impedance; basis = :total)
    shunt=ShuntAdmittance(admittance; basis = :total)
    @test series_impedance(series) === series.values
    @test shunt_admittance(shunt) === shunt.values
    @test Z(series) == impedance
    @test Y(shunt) == admittance
    @test resistance(series, 1, 2, 2) == resistance_values[1, 2, 2]
    @test reactance(series, 1, 2, 2) == imag(impedance[1, 2, 2])
    @test conductance(shunt, 2, 1, 2) == conductance_values[2, 1, 2]
    @test susceptance(shunt, 2, 1, 2) == imag(admittance[2, 1, 2])
    for container in (series, shunt)
        @test size(container) == (2, 2, 3)
        @test size(container, 3) == 3
        @test axes(container) == (Base.OneTo(2), Base.OneTo(2), Base.OneTo(3))
        @test ndims(typeof(container)) == 3
        @test eltype(typeof(container)) == ComplexF64
        @test Base.IndexStyle(typeof(container)) == IndexCartesian()
        @test container[1, 2, 2] == container.values[1, 2, 2]
        @test basis(typeof(container)) === :total
        @test basis(container) === :total
    end
    @test occursin("2×2×3", sprint(show, series))
    @test occursin("unit=Ω", sprint(show, series))
    @test occursin("2×2×3", sprint(show, shunt))
    @test occursin("unit=S", sprint(show, shunt))

    reconstructed=LineParameters(series, shunt, frequency)
    @test reconstructed.Z === series
    @test reconstructed.Y === shunt
    @test occursin("LineParameters(phase domain; 2×2×3, basis=:total)", sprint(show, reconstructed))
    detailed=sprint(show, MIME"text/plain"(), reconstructed)
    @test occursin("LineParameters · phase domain", detailed)
    @test occursin("Z  2×2×3 · Ω", detailed)
    @test occursin("Y  2×2×3 · S", detailed)

    promoted=LineParameters(
        SeriesImpedance(ComplexF32.(impedance); basis = :total),
        shunt,
        frequency
    )
    @test eltype(promoted) === ComplexF64
    @test_throws ArgumentError SeriesImpedance(impedance; basis = :invalid)
    @test_throws ArgumentError ShuntAdmittance(admittance; basis = :invalid)
    @test_throws ArgumentError SeriesImpedance{ComplexF64, 1}(impedance)
    @test_throws ArgumentError ShuntAdmittance{ComplexF64, 1}(admittance)
    for retired_basis in (:per_length, :per_lenght, :per_unit_length)
        @test_throws ArgumentError SeriesImpedance(impedance; basis = retired_basis)
        @test_throws ArgumentError ShuntAdmittance(admittance; basis = retired_basis)
        @test_throws ArgumentError LineParameters(
            impedance,
            admittance,
            frequency;
            basis = retired_basis
        )
    end
    @test_throws ArgumentError LineParameters(
        SeriesImpedance(impedance; basis = :total),
        ShuntAdmittance(admittance; basis = :pul),
        frequency
    )

    selected=parameters[2:3]
    @test basis(selected) === :total
    @test domain(selected) === domain(parameters)
    @test frequencies(selected) == frequency[2:3]
    @test selected.Z.values == impedance[:, :, 2:3]
    @test frequencies(parameters[2]) == [100.0]
    @test frequencies(parameters[[3, 1]]) == [200.0, 50.0]
    @test frequencies(parameters[:]) == frequency
    @test_throws BoundsError parameters[4]

    requests=(
        @observe(R[:, :, :]),
        @observe(L[:, :, :]),
        @observe(G[:, :, :]),
        @observe(C[:, :, :])
    )
    retained=ObservedResult(parameters,requests;atol=(R=0.,L=0.,G=0.,C=0.))
    tables=LineCableModels.ReportBuilder.tabulate(retained)
    parameter_table=tables.Z.R
    @test names(parameter_table)==["frequency","[1,1]","[1,2]","[2,1]","[2,2]"]
    @test nrow(parameter_table)==3
    @test collect(parameter_table[1,2:5])==vec(resistance_values[:,:,1])[[1,3,2,4]]
    @test nrow(DataFrame(retained))==48
    @test DataFrames.metadata(parameter_table,"basis")===basis(parameters)
    @test LineCableModels.Units.label(LineCableModels.ReportBuilder.observation_columns(parameter_table)[Symbol("[1,1]")].unit)=="Ω"
    subset_table=LineCableModels.ReportBuilder.tabulate(retained,(R,2,1,2:3))
    @test names(subset_table)==["frequency","[2,1]"]
    @test subset_table.frequency==frequency[2:3]
    @test subset_table[!,2]==resistance_values[2,1,2:3]
    @test_throws Exception DataFrame(parameters)
    @test_throws DimensionMismatch ObservedResult(parameters,((R,1,1,:),(L,2,1,:));atol=(R=0.,L=0.,G=0.,C=0.))

    zero_frequency=LineParameters(
        impedance[:, :, 1:1],
        admittance[:, :, 1:1],
        [0.0]
    )
    @test R(zero_frequency, 1, 1, 1) == real(impedance[1, 1, 1])
    @test X(zero_frequency, 1, 1, 1) == imag(impedance[1, 1, 1])
    @test_throws DomainError L(zero_frequency)
    @test_throws DomainError L(zero_frequency, 1, 1)
    @test_throws DomainError C(zero_frequency, 1, 1, 1)
    @test_throws DomainError C(zero_frequency)
    dc=ObservedResult(zero_frequency,(R,L,G,C))
    @test all(ismissing,observe(dc,L))
    @test all(ismissing,observe(dc,C))
    @test observe(dc,R)==1000real.(impedance[:,:,1:1])

    @test parentmodule(which(DataFrame, (typeof(series),))) !==
          LineCableModels.ReportBuilder
    @test parentmodule(which(DataFrame, (typeof(shunt),))) !==
          LineCableModels.ReportBuilder

    @test_throws DimensionMismatch LineParameters(
        zeros(ComplexF64, 2, 3, 1),
        zeros(ComplexF64, 2, 3, 1),
        [50.0]
    )
    @test_throws DimensionMismatch LineParameters(
        zeros(ComplexF64, 2, 2, 2),
        zeros(ComplexF64, 2, 2, 2),
        [50.0]
    )
    @test_throws ArgumentError LineParameters(
        zeros(ComplexF64, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1),
        [Inf]
    )
end

@testitem "UQ / result products / statistical invariants" tags=[:unit] setup=[
    UseEngineSupport, TestNumerics, TestFixtures] begin
    using Measurements
    using Distributions
    using Random
    using Statistics

    values=[1.0, 2.0, 3.0, 4.0, 5.0]
    summary=SampleSummary(values)
    @test summary == SampleSummary(3.0, sqrt(2.5), 1.0, 1.2, 3.0, 4.8, 5.0, 5)
    @test SampleSummary([1, 2, 3]).mean === 2.0
    @test_throws ArgumentError SampleSummary(Float64[])
    @test_throws ArgumentError SampleSummary([1.0, Inf])
    @test_throws ArgumentError SampleSummary(2.0, -1.0, 1.0, 1.1, 2.0, 2.9, 3.0, 3)
    @test_throws ArgumentError SampleSummary(2.0, 1.0, 1.0, 2.1, 2.0, 2.9, 3.0, 3)
    @test_throws ArgumentError SampleSummary(2.0, 1.0, 1.0, 1.1, 2.0, 2.9, 3.0, 0)

    density=HistogramDensity([1.0, 3.0, 5.0], [0.25, 0.25])
    @test HistogramDensity([0, 1], [2]) isa HistogramDensity{Float64}
    @test pdf(density, 2.0) == 0.25
    @test density(2.0) == 0.25
    @test pdf(density, 6.0) == 0.0
    @test cdf(density, 0.0) == 0.0
    @test cdf(density, 3.0) == 0.5
    @test cdf(density, 6.0) == 1.0
    @test UQ.cumulative_probability(density, 0.0) == 0.0
    @test UQ.cumulative_probability(density, 3.0) == 0.5
    @test UQ.cumulative_probability(density, 6.0) == 1.0
    @test quantile(density, 0.25) == 2.0
    pairs=UQ.quantile_pairs(density, [4.0, 1.0, 3.0, 2.0])
    @test pairs.sample == [1.0, 2.0, 3.0, 4.0]
    @test pairs.model == [1.5, 2.5, 3.5, 4.5]
    @test pairs.reference == (1.0, 4.5)
    @test minimum(density) == 1.0
    @test maximum(density) == 5.0
    @test isfinite(logpdf(density, 2.0))
    @test isfinite(mean(density))
    @test std(density) > 0
    @test !isempty(modes(density))
    @test isfinite(rand(MersenneTwister(91), density))
    @test_throws DomainError quantile(density, -0.1)
    @test_throws ArgumentError UQ.quantile_pairs(density, Float64[])
    @test_throws ArgumentError UQ.quantile_pairs(density, [1.0, Inf])
    @test_throws ArgumentError HistogramDensity([0.0, 1.0], Float64[])
    @test_throws ArgumentError HistogramDensity([0.0, 1.0, 2.0], [1.0])
    @test_throws ArgumentError HistogramDensity([0.0, 1.0], [-1.0])
    @test_throws ArgumentError HistogramDensity([0.0, 1.0], [0.0])
    @test_throws ArgumentError HistogramDensity([0.0, 0.0], [1.0])
    @test_throws ArgumentError HistogramDensity([0.0, Inf], [1.0])

    # Explicit, distinguishable channel inputs control publication and retention;
    # no historical model calculation supplies these expectations.
    raw_samples=(R=reshape([2.0, 3.0, 5.0, 8.0], 1, :),
        L=reshape([11.0, 13.0, 17.0, 19.0].*1e-6, 1, :),
        C=reshape([23.0, 29.0, 31.0, 37.0].*1e-10, 1, :),
        G=reshape([41.0, 43.0, 47.0, 53.0].*1e-9, 1, :))
    retained_samples=vec(raw_samples.R)
    retained_mean=sum(retained_samples)/length(retained_samples)
    retained_std=sqrt(sum(abs2, retained_samples .- retained_mean)/
        (length(retained_samples)-1))
    expected_edges=[2.0, 5.0, 8.0]
    expected_density=[2/4/3, 2/4/3]
    complete=MonteCarloResult(
        MonteCarlo(Formulation(); trials=4, seed=2027,
            return_samples=true, return_histograms=true),
        [LineCableModels.materialize(CableConstants(mean(raw_samples.R), mean(raw_samples.L),
            mean(raw_samples.C), mean(raw_samples.G)),map(x->[SampleSummary(vec(x))],raw_samples))],
        [map(x->[SampleSummary(vec(x))], raw_samples)],
        [raw_samples], [map(x->[HistogramDensity(vec(x); bins=2)], raw_samples)],
        UInt64(2027), UInt64[2039], [4])
    @test UQ.root_seed(complete) == UInt64(2027)
    @test UQ.point_seed(complete, 1) == UInt64(2039)
    @test UQ.trial_count(complete, 1) == 4
    @test UQ.confidence(complete) == 0.95
    @test UQ.cdf_tolerance(complete) == 0.02
    @test UQ.sampling_distribution(complete) === :normal
    @test_throws BoundsError UQ.point_seed(complete, 2)
    @test_throws BoundsError UQ.trial_count(complete, 2)
    sample_only=MonteCarloResult(
        complete.formulation,
        complete.values,
        complete.stats,
        complete.sample_values,
        nothing,
        complete.root_seed,
        complete.point_seeds,
        complete.trial_counts
    )
    histogram_only=MonteCarloResult(
        complete.formulation,
        complete.values,
        complete.stats,
        nothing,
        complete.histogram_values,
        complete.root_seed,
        complete.point_seeds,
        complete.trial_counts
    )
    summaries_only=MonteCarloResult(
        complete.formulation,
        complete.values,
        complete.stats,
        nothing,
        nothing,
        complete.root_seed,
        complete.point_seeds,
        complete.trial_counts
    )
    for source in (complete,sample_only,histogram_only,summaries_only)
        @test only(observables(source)) isa ObservedResult
    end

    @test all(in(observables(typeof(complete))), (
        R, L, C, G,
        (statistics, R), (statistics, L), (statistics, C), (statistics, G),
        (samples, R), (samples, L), (samples, C), (samples, G),
        (histograms, R), (histograms, L), (histograms, C), (histograms, G)
    ))
    @test nominal(only(@inferred(observe(complete, R, 1)))) == retained_mean
    @test only(@inferred(observe(complete, statistics, R, mean, 1))) == retained_mean
    @test only(@inferred(observe(complete, statistics, R, std, 1))) ≈ retained_std
    retained_histogram=observe(complete, histograms, R, 1, 1)
    @test observe(complete, samples, R, 1, 1, :) == retained_samples
    @test retained_histogram.edges == expected_edges
    @test retained_histogram.density == expected_density
    @test_throws BoundsError observe(complete, samples, R, 2, 1, :)
    @test_throws BoundsError observe(complete, samples, R, 1, 2, :)
    @test_throws ArgumentError observe(histogram_only, samples, R, 1, 1, :)
    @test_throws ArgumentError observe(sample_only, histograms, R, 1, 1)
    @test_throws DimensionMismatch MonteCarloResult(
        complete.formulation,
        complete.values,
        [(R = summary, L = summary, C = summary, G = summary)],
        complete.sample_values,
        complete.histogram_values,
        complete.root_seed,
        complete.point_seeds,
        complete.trial_counts
    )
    @test_throws DimensionMismatch MonteCarloResult(
        complete.formulation,
        complete.values,
        complete.stats,
        [(R = [1.0], L = [1.0], C = [1.0], G = [1.0])],
        complete.histogram_values,
        complete.root_seed,
        complete.point_seeds,
        complete.trial_counts
    )

    retained=ObservedResult(complete,1,((statistics,R,mean,1),(samples,R,1,:),(histograms,R,1));
        units=(:milli,:milli,:milli))
    @test retained.quantities[1].quantity==quantity(R)
    @test retained.quantities[1].values==retained_mean*1e6
    @test retained.quantities[2].values==retained_samples.*1e6
    @test retained.quantities[3].distribution.edges==expected_edges.*1e6
    @test retained.quantities[3].values.density≈expected_density./1e6 rtol=4eps(Float64)

    summary_product=only(statistics(complete))
    sample_product=only(samples(complete))
    histogram_product=only(histograms(complete))
    @test keys(summary_product) == (:R, :L, :C, :G)
    @test keys(sample_product) == (:R, :L, :C, :G)
    @test keys(histogram_product) == (:R, :L, :C, :G)
    @test !applicable(observe, summary_product, R)
    @test !applicable(observe, sample_product, R, :)
    @test !applicable(observe, histogram_product, R)
    summary_observed=ObservedResult(complete,1,((statistics,R,1),);quantity_units=(R=:milli,))
    histogram_observed=ObservedResult(complete,1,((histograms,R,1),);units=(:milli,))
    @test length(summary_observed.quantities)==7
    summary_factor=scale_factor(R,basis(complete),first(summary_observed.quantities).unit)
    @test observe(summary_observed,statistics,R,mean)==summary_product.R[1].mean*summary_factor
    @test observe(summary_observed,statistics,R,std)==summary_product.R[1].std*abs(summary_factor)
    @test only(histogram_observed.quantities).distribution.edges==histogram_product.R[1].edges.*summary_factor
    @test only(histogram_observed.quantities).values.density≈histogram_product.R[1].density./summary_factor rtol=4eps(Float64)

    frequency=[50.0, 100.0]
    impedance=fill(1.0e-4+2.0e-4im, 2, 2, 2)
    admittance=fill(3.0e-8+4.0e-8im, 2, 2, 2)
    parameters=LineParameters(impedance, admittance, frequency)
    storage=LineCableModels.UQ._sample_storage(parameters, 2)
    LineCableModels.UQ._record_sample!(storage, parameters, 1, frequency)
    LineCableModels.UQ._record_sample!(storage, parameters, 2, frequency)
    @test @allocated(LineCableModels.UQ._record_sample!(storage, parameters, 2, frequency)) ==
          0
    @test storage.R[1, 1, 1, :] == fill(1.0e-4, 2)
    @test storage.L[1, 1, 1, :] ==
          fill(2.0e-4 / (2π * frequency[1]), 2)

    line_summary=fill(SampleSummary([1.0, 2.0]), size(impedance))
    line_histogram=fill(density, size(impedance))
    line_statistics=(R = line_summary, L = line_summary, C = line_summary, G = line_summary)
    line_histograms=(R = line_histogram, L = line_histogram,
        C = line_histogram, G = line_histogram)
    line_result=MonteCarloResult(
        complete.formulation,
        [LineCableModels.materialize(parameters,line_statistics)],
        [line_statistics],
        [storage],
        [line_histograms],
        complete.root_seed,
        complete.point_seeds,
        [2]
    )
    @test observe(line_result, samples, R, 1, 1, 1, 1, :) == fill(1.0e-4, 2)
    @test observe(line_result, frequencies, 1, :) == frequency
    line_observed=ObservedResult(line_result,1)
    @test first(line_observed.quantities).coordinates.frequencies==frequency
    @test first(line_observed.quantities).coordinates.frequency_unit==LineCableModels.Units.units(:base,:hertz)
    malformed_samples=(R = storage.R[:, :, 1:1, :], L = storage.L,
        C = storage.C, G = storage.G)
    @test_throws DimensionMismatch MonteCarloResult(
        complete.formulation,
        line_result.values,
        [line_statistics],
        [malformed_samples],
        [line_histograms],
        complete.root_seed,
        complete.point_seeds,
        [2]
    )
end
