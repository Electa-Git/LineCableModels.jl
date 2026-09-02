@testitem "Engine / result containers / numeric and selector behavior" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using DataFrames

    library=CablesLibrary()
    load!(
        library;
        file_name = joinpath(
            pkgdir(LineCableModels),
            "test",
            "fixtures",
            "data",
            "mv_cable_design.json"
        )
    )
    design=first(values(library.data))
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
    @test constants_observables isa LineCableModels.Grammar.ObservationPublication
    @test first(constants_observables).values ≈ 1_000constants.R
    @test keys(first(constants_observables)) == (:values, :quantity, :unit)
    constants_table=DataFrame(constants_observables)
    @test names(constants_table) == ["core", "R", "L", "C", "G"]
    @test only(constants_table.core) === :core
    @test constants_table.R ≈ 1_000constants.R
    @test_throws Exception DataFrame(constants)

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
    parameter_observables=observables(
        parameters,
        (
            (frequencies, Colon()),
            (R, 1, 2, Colon()),
            (Z, abs, 1, 2, Colon()),
            (Z, angle, 1, 2, Colon())
        )
    )
    @test parameter_observables isa LineCableModels.Grammar.ObservationPublication
    @test parameter_observables[1].values == frequency
    @test parameter_observables[2].values == resistance_values[1, 2, :]
    @test parameter_observables[4].values ≈
          rad2deg.(angle.(impedance[1, 2, :]))
    @test parameter_observables[1].values !== parameters.f

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
    parameter_table=DataFrame(observables(parameters, requests))
    @test parameter_table isa DataFrame
    @test names(parameter_table) == ["frequency", "row", "column", "R", "L", "G", "C"]
    @test nrow(parameter_table) == 12
    @test parameter_table[1:4, [:frequency, :row, :column]] ==
          DataFrame(
        frequency = fill(50.0, 4),
        row = [1, 1, 2, 2],
        column = [1, 2, 1, 2]
    )
    @test parameter_table.R[1:4] == vec(resistance_values[:, :, 1])[[1, 3, 2, 4]]
    observed_columns=LineCableModels.ReportBuilder.observation_columns(parameter_table)
    @test keys(observed_columns) == (:frequency, :R, :L, :G, :C)
    @test LineCableModels.Units.label(observed_columns.R.unit) == "Ω"
    @test LineCableModels.Units.label(observed_columns.L.unit) == "mH"
    @test LineCableModels.Units.label(observed_columns.G.unit) == "S"
    @test LineCableModels.Units.label(observed_columns.C.unit) == "μF"
    @test DataFrames.metadata(parameter_table, "basis") === basis(parameters)
    subset_table=DataFrame(observables(parameters, (@observe(R[2, 1, 2:3]),)))
    @test subset_table.row == [2, 2]
    @test subset_table.column == [1, 1]
    @test subset_table.frequency == frequency[2:3]
    @test subset_table.R == resistance_values[2, 1, 2:3]
    transformed_table=DataFrame(observables(
        parameters, (
            @observe((Z, abs)[:, :, :]),
            @observe((Y, angle)[:, :, :])
        )))
    @test names(transformed_table) ==
          ["frequency", "row", "column", "|Z|", "∠Y"]
    @test transformed_table[!, Symbol("|Z|")][1:4] ≈
          vec(abs.(impedance[:, :, 1]))[[1, 3, 2, 4]]
    @test transformed_table[!, Symbol("∠Y")][1:4] ≈
          rad2deg.(vec(angle.(admittance[:, :, 1]))[[1, 3, 2, 4]])
    @test_throws Exception DataFrame(parameters)
    @test_throws DimensionMismatch observables(parameters, (
        @observe(R[1, 1, :]),
        @observe(L[2, 1, :])
    ))

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
    @test DataFrame(observables(zero_frequency, (
        @observe(R[:, :, :]),
        @observe(G[:, :, :])
    ))) isa DataFrame
    @test_throws DomainError observables(zero_frequency,
        (
            @observe(R[:, :, :]),
            @observe(L[:, :, :]),
            @observe(G[:, :, :]),
            @observe(C[:, :, :])
        ))

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
    EngineTestSupport, UseEngineSupport, TestNumerics, TestFixtures] begin
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

    complete=TestFixtures.cable_monte_carlo_result()
    @test UQ.root_seed(complete) == UInt64(1)
    @test UQ.point_seed(complete, 1) == UInt64(1)
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
    @test_throws MethodError observables(complete)
    @test_throws MethodError observables(sample_only)
    @test_throws MethodError observables(histogram_only)
    @test_throws MethodError observables(summaries_only)

    @test observables(typeof(complete)) == (
        R, L, C, G,
        (statistics, R), (statistics, L), (statistics, C), (statistics, G),
        (samples, R), (samples, L), (samples, C), (samples, G),
        (histograms, R), (histograms, L), (histograms, C), (histograms, G)
    )
    @test only(@inferred(observe(complete, R, 1))) == 2.5
    @test only(@inferred(observe(complete, statistics, R, mean, 1))) == 2.5
    @test only(@inferred(observe(complete, statistics, R, std, 1))) ≈ sqrt(5 / 3)
    retained_samples=[1.0, 2.0, 3.0, 4.0]
    retained_histogram=observe(complete, histograms, R, 1, 1)
    @test observe(complete, samples, R, 1, 1, :) == retained_samples
    @test retained_histogram.edges == density.edges
    @test retained_histogram.density == density.density
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

    result_publication=observables(
        complete,
        (
            (statistics, R, mean, 1, 1),
            (samples, R, 1, 1, Colon()),
            (histograms, R, 1, 1)
        );
        units = (:milli, :milli, :milli)
    )
    @test result_publication[1].quantity == quantity(R)
    @test result_publication[1].values == 2.5e6
    @test result_publication[2].values == retained_samples .* 1.0e6
    @test result_publication[3].values.edges == density.edges .* 1.0e6
    @test result_publication[3].values.density ==
          density.density ./ 1.0e6

    summary_product=only(statistics(complete))
    sample_product=only(samples(complete))
    histogram_product=only(histograms(complete))
    @test keys(summary_product) == (:R, :L, :C, :G)
    @test keys(sample_product) == (:R, :L, :C, :G)
    @test keys(histogram_product) == (:R, :L, :C, :G)
    @test !applicable(observe, summary_product, R)
    @test !applicable(observe, sample_product, R, :)
    @test !applicable(observe, histogram_product, R)
    summary_publication=observables(
        complete,
        ((statistics, R, 1, 1),);
        units = (:milli,)
    )|>only
    histogram_publication=observables(
        complete,
        ((histograms, R, 1, 1),);
        units = (:milli,)
    )|>only
    summary_factor=scale_factor(R, basis(complete), summary_publication.unit)
    @test summary_publication.values.mean == summary_product.R[1].mean * summary_factor
    @test summary_publication.values.std == summary_product.R[1].std * abs(summary_factor)
    @test histogram_publication.values.edges ==
          histogram_product.R[1].edges .* summary_factor
    @test histogram_publication.values.density ==
          histogram_product.R[1].density ./ summary_factor

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
        [parameters],
        [line_statistics],
        [storage],
        [line_histograms],
        complete.root_seed,
        complete.point_seeds,
        [2]
    )
    @test observe(line_result, samples, R, 1, 1, 1, 1, :) == fill(1.0e-4, 2)
    @test observe(line_result, frequencies, 1, :) == frequency
    line_frequency=observables(
        line_result,
        ((frequencies, 1, Colon()),)
    )|>only
    @test line_frequency.values == frequency
    @test line_frequency.quantity == quantity(frequencies)
    @test line_frequency.unit == LineCableModels.Units.units(:base, :hertz)
    malformed_samples=(R = storage.R[:, :, 1:1, :], L = storage.L,
        C = storage.C, G = storage.G)
    @test_throws DimensionMismatch MonteCarloResult(
        complete.formulation,
        [parameters],
        [line_statistics],
        [malformed_samples],
        [line_histograms],
        complete.root_seed,
        complete.point_seeds,
        [2]
    )
end
