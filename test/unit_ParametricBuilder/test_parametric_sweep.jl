@testitem "ParametricSweep: typed deterministic results and accessors" setup = [defaults] begin
    frequency = [50.0, 100.0, 200.0]
    omega = reshape(2π .* frequency, 1, 1, :)
    resistance_values = reshape(collect(1.0:12.0), 2, 2, 3) .* 1.0e-4
    inductance_values = fill(2.0e-7, 2, 2, 3)
    conductance_values = fill(3.0e-9, 2, 2, 3)
    capacitance_values = fill(4.0e-10, 2, 2, 3)

    first_result = LineParameters(
        complex.(resistance_values, inductance_values .* omega),
        complex.(conductance_values, capacitance_values .* omega),
        frequency
    )
    second_result = LineParameters(
        complex.(2resistance_values, 2inductance_values .* omega),
        complex.(2conductance_values, 2capacitance_values .* omega),
        frequency
    )
    sweep = @inferred ParametricSweep(
        [(temperature = 20.0,), (temperature = 90.0,)],
        [first_result, second_result]
    )

    @test sweep isa AbstractVector{typeof(first_result)}
    @test @inferred(length(sweep)) == 2
    @test @inferred(ncases(sweep)) == 2
    @test size(sweep) == (2,)
    @test @inferred(sweep[1]) === first_result
    @test collect(sweep) == [first_result, second_result]
    @test cases(sweep) == [(temperature = 20.0,), (temperature = 90.0,)]
    @test results(sweep) == [first_result, second_result]

    selected = @inferred sweep[2:2]
    @test selected isa ParametricSweep
    @test cases(selected) == [(temperature = 90.0,)]
    @test only(selected) === second_result
    @test sweep[:] !== sweep
    @test results(sweep[:]) == results(sweep)
    @test copy(sweep) !== sweep
    @test results(copy(sweep)) == results(sweep)

    @test @inferred(basis(sweep)) === :per_length
    @test @inferred(domain(sweep)) === PhaseDomain
    @test @inferred(frequencies(sweep)) == frequency
    @test @inferred(nconductors(sweep)) == 2
    @test @inferred(nfrequencies(sweep)) == 3
    @test @inferred(basis(sweep, 2)) === :per_length
    @test @inferred(frequencies(sweep, 2)) == frequency

    @test @inferred(Z(sweep)) == [Z(first_result), Z(second_result)]
    @test @inferred(Y(sweep, 1, 2, 1)) == Y(first_result, 2, 1)
    @test @inferred(R(sweep, 2, 1, 2)) == R(second_result, 1, 2)
    @test X(sweep, 1, 2, 2, 2:3) == X(first_result, 2, 2, 2:3)
    @test L(sweep, 2, 1, 1, :) ≈ L(second_result, 1, 1, :)
    @test G(sweep, 1, 2, 1, 2) == G(first_result, 2, 1, 2)
    @test B(sweep, 2, 1, 2) == B(second_result, 1, 2)
    @test C(sweep, 1, 2, 2, 1:2) ≈ C(first_result, 2, 2, 1:2)
    @test series_impedance(sweep) == [
        series_impedance(first_result),
        series_impedance(second_result)
    ]
    @test series_impedance(sweep, 2) === series_impedance(second_result)
    @test shunt_admittance(sweep, 1) === shunt_admittance(first_result)
    @test capacitance(sweep, 2, 1, 1) ≈ C(second_result, 1, 1)

    different_frequency = LineParameters(
        Array(Z(second_result)),
        Array(Y(second_result)),
        [60.0, 120.0, 240.0]
    )
    varying_axes = ParametricSweep([:first, :second], [first_result, different_frequency])
    @test frequencies(varying_axes, 2) == [60.0, 120.0, 240.0]
    @test_throws ArgumentError frequencies(varying_axes)

    empty_sweep = @inferred ParametricSweep(Symbol[], typeof(first_result)[])
    @test isempty(empty_sweep)
    @test size(empty_sweep) == (0,)
    @test isempty(cases(empty_sweep))
    @test isempty(results(empty_sweep))
    @test isempty(sweep[1:0])
    @test_throws ArgumentError basis(empty_sweep)
    @test_throws ArgumentError ParametricSweep(Any[], Any[])
    @test_throws DimensionMismatch ParametricSweep([:one], [first_result, second_result])
    @test_throws ArgumentError ParametricSweep(Any[:one, 2], [first_result, second_result])
    total_result = LineParameters(
        Array(Z(first_result)),
        Array(Y(first_result)),
        frequency;
        basis = :total
    )
    @test_throws ArgumentError ParametricSweep(
        [:per_length, :total],
        LineParameters[first_result, total_result]
    )
end
