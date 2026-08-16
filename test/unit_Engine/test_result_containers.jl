@testitem "Primitive results retain numeric and selector behavior" setup = [defaults] begin
    using DataFrames

    library = CablesLibrary()
    load!(
        library;
        file_name=joinpath(pkgdir(LineCableModels), "test", "cable_test.json"),
    )
    design = first(values(library.data))
    constants = compute!(CableConstantsProblem(design), Formulation())
    @test constants isa CableConstants{Float64}
    @test constants.R ≈ 2.7567652874268654e-5 rtol=2eps()
    @test constants.L ≈ 2.8718381083175005e-7 rtol=2eps()
    @test constants.C ≈ 4.1335723330313053e-10 rtol=2eps()
    @test basis(constants) === :per_length
    @test resistance(constants) === constants.R
    @test inductance(constants) === constants.L
    @test capacitance(constants) === constants.C
    @test DataFrame(constants).value == [constants.R, constants.L, constants.C]
    @test_throws ArgumentError compute!(
        CableConstantsProblem(design),
        Formulation();
        options=(output_basis=:total,),
    )

    frequency = [50.0, 100.0, 200.0]
    angular = reshape(2π .* frequency, 1, 1, :)
    resistance_values = reshape(collect(1.0:12.0), 2, 2, 3) .* 1.0e-4
    inductance_values = reshape(collect(13.0:24.0), 2, 2, 3) .* 1.0e-7
    conductance_values = reshape(collect(25.0:36.0), 2, 2, 3) .* 1.0e-8
    capacitance_values = reshape(collect(37.0:48.0), 2, 2, 3) .* 1.0e-10
    impedance = complex.(resistance_values, inductance_values .* angular)
    admittance = complex.(conductance_values, capacitance_values .* angular)
    parameters = LineParameters(
        impedance,
        admittance,
        frequency;
        basis=:total,
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

    selected = parameters[2:3]
    @test basis(selected) === :total
    @test domain(selected) === domain(parameters)
    @test frequencies(selected) == frequency[2:3]
    @test selected.Z.values == impedance[:, :, 2:3]

    series_frames, shunt_frames = DataFrame(parameters)
    @test DataFrames.metadata(series_frames[1, 1], "units")[:real] == "Ω"
    @test DataFrames.metadata(series_frames[1, 1], "units")[:imag] == "Ω"
    @test DataFrames.metadata(shunt_frames[1, 1], "units")[:real] == "S"
    @test DataFrames.metadata(shunt_frames[1, 1], "units")[:imag] == "S"
    series_rl, shunt_gc = DataFrame(parameters, (R, L, G, C))
    @test DataFrames.metadata(series_rl[1, 1], "units")[:R] == "Ω"
    @test DataFrames.metadata(series_rl[1, 1], "units")[:L] == "mH"
    @test DataFrames.metadata(shunt_gc[1, 1], "units")[:C] == "μF"
    @test_throws ArgumentError DataFrame(parameters; tol=-1.0)

    zero_frequency = LineParameters(
        impedance[:, :, 1:1],
        admittance[:, :, 1:1],
        [0.0],
    )
    @test R(zero_frequency, 1, 1, 1) == real(impedance[1, 1, 1])
    @test X(zero_frequency, 1, 1, 1) == imag(impedance[1, 1, 1])
    @test_throws DomainError L(zero_frequency)
    @test_throws DomainError C(zero_frequency, 1, 1, 1)
    @test DataFrame(zero_frequency) isa Tuple
    @test_throws DomainError DataFrame(zero_frequency, (R, L, G, C))

    @test_throws DimensionMismatch LineParameters(
        zeros(ComplexF64, 2, 3, 1),
        zeros(ComplexF64, 2, 3, 1),
        [50.0],
    )
    @test_throws DimensionMismatch LineParameters(
        zeros(ComplexF64, 2, 2, 2),
        zeros(ComplexF64, 2, 2, 2),
        [50.0],
    )
    @test_throws ArgumentError LineParameters(
        zeros(ComplexF64, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1),
        [Inf],
    )
end

@testitem "Generic Monte Carlo result products validate their invariants" setup = [defaults] begin
    using Distributions
    using Random
    using Statistics

    values = [1.0, 2.0, 3.0, 4.0, 5.0]
    summary = SampleSummary(values)
    @test summary == SampleSummary(3.0, sqrt(2.5), 1.0, 1.2, 3.0, 4.8, 5.0)
    @test SampleSummary([1, 2, 3]).mean === 2.0
    @test_throws ArgumentError SampleSummary(Float64[])
    @test_throws ArgumentError SampleSummary([1.0, Inf])
    @test_throws ArgumentError SampleSummary(2.0, -1.0, 1.0, 1.1, 2.0, 2.9, 3.0)
    @test_throws ArgumentError SampleSummary(2.0, 1.0, 1.0, 2.1, 2.0, 2.9, 3.0)

    density = HistogramPDF([1.0, 3.0, 5.0], [0.25, 0.25])
    @test HistogramPDF([0, 1], [2]) isa HistogramPDF{Float64}
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
    @test_throws ArgumentError HistogramPDF([0.0, 1.0], Float64[])
    @test_throws ArgumentError HistogramPDF([0.0, 1.0, 2.0], [1.0])
    @test_throws ArgumentError HistogramPDF([0.0, 1.0], [-1.0])
    @test_throws ArgumentError HistogramPDF([0.0, 1.0], [0.0])
    @test_throws ArgumentError HistogramPDF([0.0, 0.0], [1.0])
    @test_throws ArgumentError HistogramPDF([0.0, Inf], [1.0])
end
