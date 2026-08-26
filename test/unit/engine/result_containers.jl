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
    constants=compute(CableConstantsProblem(design), Formulation())
    @test constants isa CableConstants{Float64}
    @test constants.R ≈ 2.7567652874268654e-5 rtol=2eps()
    @test constants.L ≈ 2.8718381083175005e-7 rtol=2eps()
    @test constants.C ≈ 4.1335723330313053e-10 rtol=2eps()
    @test basis(constants) === :pul
    constants_display=sprint(show, constants)
    @test occursin("CableConstants(R=", constants_display)
    @test occursin("Ω/m", constants_display)
    @test occursin("H/m", constants_display)
    @test occursin("F/m", constants_display)
    @test resistance(constants) === constants.R
    @test inductance(constants) === constants.L
    @test capacitance(constants) === constants.C
    @test observe(constants, R) === constants.R
    @test observables(typeof(constants)) == (R, L, C)
    constants_observables=observables(
        constants,
        (resistance = R, inductance = L, capacitance = C)
    )
    @test keys(constants_observables) == (:resistance, :inductance, :capacitance)
    @test constants_observables.resistance.values ≈ 1_000constants.R
    @test keys(constants_observables.resistance) == (:values, :quantity, :unit)
    @test DataFrame(constants).value == [constants.R, constants.L, constants.C]
    @test_throws ArgumentError compute(
        CableConstantsProblem(design),
        Formulation();
        options = (output_basis = :total,)
    )

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
            frequency = (frequencies, Colon()),
            resistance = (R, 1, 2, Colon()),
            impedance_magnitude = (Z, abs, 1, 2, Colon()),
            impedance_angle = (Z, angle, 1, 2, Colon())
        )
    )
    @test keys(parameter_observables) ==
          (:frequency, :resistance, :impedance_magnitude, :impedance_angle)
    @test parameter_observables.frequency.values == frequency
    @test parameter_observables.resistance.values == resistance_values[1, 2, :]
    @test parameter_observables.impedance_angle.values ≈
          rad2deg.(angle.(impedance[1, 2, :]))
    @test parameter_observables.frequency.values !== parameters.f

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
    @test occursin("2×2×3 [Ω]", sprint(show, series))
    @test occursin("2×2×3 [S]", sprint(show, shunt))

    reconstructed=LineParameters(series, shunt, frequency)
    @test reconstructed.Z === series
    @test reconstructed.Y === shunt
    @test occursin("LineParameters{ComplexF64} 2×2×3 [total]", sprint(show, reconstructed))
    detailed=sprint(show, MIME"text/plain"(), reconstructed)
    @test occursin("domain: PhaseDomain, frequencies: 3", detailed)
    @test occursin("Z [Ω]", detailed)
    @test occursin("Y [S]", detailed)

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

    parameter_table=DataFrame(parameters)
    @test parameter_table isa DataFrame
    @test names(parameter_table) == [
        "family", "row", "column", "frequency", "quantity", "value", "unit"
    ]
    @test nrow(parameter_table) == 48
    @test parameter_table[1:4, [:family, :row, :column, :frequency, :quantity]] ==
          DataFrame(
        family = fill(:series, 4),
        row = ones(Int, 4),
        column = ones(Int, 4),
        frequency = [50.0, 50.0, 100.0, 100.0],
        quantity = [:real, :imag, :real, :imag]
    )
    units_metadata=DataFrames.metadata(parameter_table, "units")
    @test units_metadata[(:series, :real)] == "Ω"
    @test units_metadata[(:series, :imag)] == "Ω"
    @test units_metadata[(:shunt, :real)] == "S"
    @test units_metadata[(:shunt, :imag)] == "S"
    @test DataFrames.metadata(parameter_table, "basis") === basis(parameters)
    headings_metadata=DataFrames.metadata(parameter_table, "headings")
    @test headings_metadata[(:series, :real)] == "Series resistance [Ω]"
    @test headings_metadata[(:shunt, :imag)] == "Shunt susceptance [S]"
    requests_metadata=DataFrames.metadata(parameter_table, "requests")
    @test requests_metadata[(:series, :real)] === R
    @test requests_metadata[(:series, :imag)] === X
    @test requests_metadata[(:shunt, :real)] === G
    @test requests_metadata[(:shunt, :imag)] === B
    parameter_rlgc=DataFrame(parameters, (R, L, G, C))
    rlgc_units=DataFrames.metadata(parameter_rlgc, "units")
    @test rlgc_units[(:series, :R)] == "Ω"
    @test rlgc_units[(:series, :L)] == "mH"
    @test rlgc_units[(:shunt, :C)] == "μF"
    @test_throws ArgumentError DataFrame(parameters; tol = -1.0)

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
    @test DataFrame(zero_frequency) isa DataFrame
    @test_throws DomainError DataFrame(zero_frequency, (R, L, G, C))

    standalone_series=DataFrame(series; freqs = frequency)
    standalone_shunt=DataFrame(shunt, (G, C); freqs = frequency)
    @test standalone_series isa DataFrame
    @test standalone_shunt isa DataFrame
    @test unique(standalone_series.family) == [:series]
    @test unique(standalone_shunt.family) == [:shunt]
    @test unique(standalone_series.frequency) == frequency
    @test unique(standalone_shunt.frequency) == frequency
    @test unique(standalone_series.quantity) == [:real, :imag]
    @test unique(standalone_shunt.quantity) == [:G, :C]
    @test_throws ArgumentError DataFrame(series)
    @test_throws ArgumentError DataFrame(shunt)
    @test LineCableModels.ReportBuilder.clip(1.0 + 2.0im, 1.0) ==
          1.0 + 2.0im
    @test_throws DimensionMismatch DataFrame(series; freqs = [50.0])
    @test_throws ArgumentError DataFrame(series; freqs = [50.0, Inf, 200.0])
    @test_throws ArgumentError DataFrame(series; tol = Inf)

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
    @test quantile(density, 0.25) == 2.0
    @test minimum(density) == 1.0
    @test maximum(density) == 5.0
    @test isfinite(logpdf(density, 2.0))
    @test isfinite(mean(density))
    @test std(density) > 0
    @test !isempty(modes(density))
    @test isfinite(rand(MersenneTwister(91), density))
    @test_throws DomainError quantile(density, -0.1)
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
        R, L, C,
        (statistics, R), (statistics, L), (statistics, C),
        (samples, R), (samples, L), (samples, C),
        (histograms, R), (histograms, L), (histograms, C)
    )
    @test @inferred(observe(complete, R, 1)) == 2.5
    @test @inferred(observe(complete, statistics, R, mean, 1)) == 2.5
    @test @inferred(observe(complete, statistics, R, std, 1)) ≈ sqrt(5 / 3)
    retained_samples=[1.0, 2.0, 3.0, 4.0]
    retained_histogram=observe(complete, histograms, R, 1)
    @test observe(complete, samples, R, 1, :) == retained_samples
    @test retained_histogram.edges == density.edges
    @test retained_histogram.density == density.density
    @test_throws BoundsError observe(complete, samples, R, 2, :)
    @test_throws BoundsError observe(complete, samples, R, 1, 5)
    @test_throws ArgumentError observe(histogram_only, samples, R, 1, :)
    @test_throws ArgumentError observe(sample_only, histograms, R, 1)

    result_publication=observables(
        complete,
        (
            mean_resistance = (statistics, R, mean, 1),
            sample_resistance = (samples, R, 1, Colon()),
            resistance_density = (histograms, R, 1)
        );
        units = (
            mean_resistance = :milli,
            sample_resistance = :milli,
            resistance_density = :milli
        )
    )
    @test result_publication.mean_resistance.quantity == quantity(R)
    @test result_publication.mean_resistance.values == 2.5e6
    @test result_publication.sample_resistance.values == retained_samples .* 1.0e6
    @test result_publication.resistance_density.values.edges == density.edges .* 1.0e6
    @test result_publication.resistance_density.values.density ==
          density.density ./ 1.0e6

    summary_product=only(statistics(complete))
    sample_product=only(samples(complete))
    histogram_product=only(histograms(complete))
    @test observe(summary_product, R) === summary_product.R
    @test observe(summary_product, R, mean) == summary_product.R.mean
    @test observe(summary_product, R, std) == summary_product.R.std
    @test observe(sample_product, R, :) == sample_product.R[:]
    @test observe(histogram_product, R) === histogram_product.R
    published=observables(summary_product, (mean_resistance = (R, mean),))
    @test published.mean_resistance.values == summary_product.R.mean * 1_000
    @test keys(published.mean_resistance) == (:values, :quantity, :unit)
    summary_publication=observables(
        summary_product,
        (resistance = R,);
        units = (resistance = :milli,)
    ).resistance
    histogram_publication=observables(
        histogram_product,
        (resistance = R,);
        units = (resistance = :milli,)
    ).resistance
    summary_factor=scale_factor(R, basis(summary_product), summary_publication.unit)
    @test summary_publication.values.mean == summary_product.R.mean * summary_factor
    @test summary_publication.values.std == summary_product.R.std * abs(summary_factor)
    @test histogram_publication.values.edges ==
          histogram_product.R.edges .* summary_factor
    @test histogram_publication.values.density ==
          histogram_product.R.density ./ summary_factor

    frequency=[50.0, 100.0]
    impedance=fill(1.0e-4+2.0e-4im, 2, 2, 2)
    admittance=fill(3.0e-8+4.0e-8im, 2, 2, 2)
    parameters=LineParameters(impedance, admittance, frequency)
    storage=LineCableModels.UQ._sample_storage(parameters, 2)
    LineCableModels.UQ._record_sample!(storage, parameters, 1, frequency)
    LineCableModels.UQ._record_sample!(storage, parameters, 2, frequency)
    @test @allocated(LineCableModels.UQ._record_sample!(storage, parameters, 2, frequency)) ==
          0
    @test observe(storage, R, 1, 1, 1, :) == fill(1.0e-4, 2)
    @test observe(storage, L, 1, 1, 1, :) ==
          fill(2.0e-4 / (2π * frequency[1]), 2)
end
