@testitem "Engine / observations / complex, diagonal, and frequency selections" tags=[:unit] begin
    using LinearAlgebra: diag
    frequency = [0.1, 50.0, 1e6]
    order = reshape(collect(1.0:12.0), 2, 2, 3)
    resistance = order .* 1e-4
    inductance = order .* 1e-7
    conductance = order .* 1e-9
    capacitance = order .* 1e-10
    omega = reshape(2pi .* frequency, 1, 1, :)
    impedance = resistance .+ im .* omega .* inductance
    admittance = conductance .+ im .* omega .* capacitance

    for result_basis in (:pul, :total)
        line = LineParameters(copy(impedance), copy(admittance), frequency; basis=result_basis)
        series, shunt = line.Z, line.Y
        @test basis(typeof(line)) === result_basis
        @test basis(typeof(series)) === result_basis
        @test basis(typeof(shunt)) === result_basis
        for (selector, standalone, values) in (
                (Z, series, impedance), (R, series, resistance),
                (X, series, omega .* inductance), (L, series, inductance),
                (Y, shunt, admittance), (G, shunt, conductance),
                (B, shunt, omega .* capacitance), (C, shunt, capacitance))
            @test (selector, diag) in observables(typeof(line))
            @test (selector, diag) in observables(typeof(standalone))
            expected = [values[row, row, sample] for row in 1:2, sample in 1:3]
            @test @inferred(observe(line, selector, diag)) ≈ expected
            @test @inferred(observe(line, selector, diag, [2, 1], [3, 1])) ≈ expected[[2, 1], [3, 1]]
            @test @inferred(observe(line, selector, diag, 2, 3)) ≈ expected[2, 3]
            @test observe(line, selector, 2, 1) ≈ values[2, 1, :]
            arguments = selector in (L, C) ? (frequency,) : ()
            @test observe(standalone, selector, diag, arguments...) ≈ expected
            @test observe(standalone, selector, diag, arguments..., [2, 1], [3, 1]) ≈ expected[[2, 1], [3, 1]]
        end
        for (selector, standalone, values) in ((Z, series, impedance), (Y, shunt, admittance)),
                transform in (abs, angle)
            expected = transform.(values)
            @test @inferred(observe(line, selector, transform)) ≈ expected
            @test @inferred(observe(standalone, selector, transform, 2, 1, [3, 1])) ≈ expected[2, 1, [3, 1]]
        end
        @test C(line, 2, 1) ≈ capacitance[2, 1, :]
        @test Z(line) == impedance
        @test Y(line) == admittance
        @test frequencies(line) == frequency
    end
end
