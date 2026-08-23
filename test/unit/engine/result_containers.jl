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
    @test basis(constants) === :per_length
    @test resistance(constants) === constants.R
    @test inductance(constants) === constants.L
    @test capacitance(constants) === constants.C
    constants_observables=observables(constants)
    @test constants_observables isa NamedTuple
    @test keys(constants_observables) == (:resistance, :inductance, :capacitance)
    @test values(constants_observables) == (constants.R, constants.L, constants.C)
    @test !ismutabletype(typeof(constants_observables))
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
    @test series_impedance(parameters) === parameters.Z
    @test shunt_admittance(parameters) === parameters.Y
    @test resistance(parameters, 1, 1, 1) == R(parameters, 1, 1, 1)
    @test reactance(parameters, 1, 1, 1) == X(parameters, 1, 1, 1)
    @test inductance(parameters, 1, 1, 1) ≈ L(parameters, 1, 1, 1)
    @test conductance(parameters, 1, 1, 1) == G(parameters, 1, 1, 1)
    @test susceptance(parameters, 1, 1, 1) == B(parameters, 1, 1, 1)
    @test capacitance(parameters, 1, 1, 1) ≈ C(parameters, 1, 1, 1)
    @test L(parameters) ≈ inductance_values
    @test C(parameters) ≈ capacitance_values
    parameter_observables=observables(parameters)
    @test parameter_observables isa NamedTuple
    @test keys(parameter_observables) ==
          (:frequency, :series_impedance, :shunt_admittance)
    @test parameter_observables.frequency == parameters.f
    @test parameter_observables.series_impedance == parameters.Z
    @test parameter_observables.shunt_admittance == parameters.Y
    @test parameter_observables.frequency !== parameters.f
    @test parameter_observables.series_impedance.values !== parameters.Z.values
    @test parameter_observables.shunt_admittance.values !== parameters.Y.values
    @test !ismutabletype(typeof(parameter_observables))

    series=SeriesImpedance(impedance; basis = :total)
    shunt=ShuntAdmittance(admittance; basis = :total)
    @test series_impedance(series) === series
    @test shunt_admittance(shunt) === shunt
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
    @test_throws ArgumentError LineParameters(
        SeriesImpedance(impedance; basis = :total),
        ShuntAdmittance(admittance; basis = :per_length),
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

    series_frames, shunt_frames=DataFrame(parameters)
    @test DataFrames.metadata(series_frames[1, 1], "units")[:real] == "Ω"
    @test DataFrames.metadata(series_frames[1, 1], "units")[:imag] == "Ω"
    @test DataFrames.metadata(shunt_frames[1, 1], "units")[:real] == "S"
    @test DataFrames.metadata(shunt_frames[1, 1], "units")[:imag] == "S"
    series_rl, shunt_gc=DataFrame(parameters, (R, L, G, C))
    @test DataFrames.metadata(series_rl[1, 1], "units")[:R] == "Ω"
    @test DataFrames.metadata(series_rl[1, 1], "units")[:L] == "mH"
    @test DataFrames.metadata(shunt_gc[1, 1], "units")[:C] == "μF"
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
    @test DataFrame(zero_frequency) isa Tuple
    @test_throws DomainError DataFrame(zero_frequency, (R, L, G, C))

    standalone_series=DataFrame(series)
    standalone_shunt=DataFrame(shunt, (G, C); freqs = frequency)
    @test standalone_series[1, 1].frequency == [1.0, 2.0, 3.0]
    @test standalone_shunt[1, 1].frequency == frequency
    @test names(standalone_series[1, 1]) == ["frequency", "real", "imag"]
    @test names(standalone_shunt[1, 1]) == ["frequency", "G", "C"]
    @test Engine._clip_field(1.0 + 2.0im, 1.0) == 1.0 + 2.0im
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
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
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
end
