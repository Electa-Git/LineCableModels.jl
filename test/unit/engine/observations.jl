@testitem "Engine / observations / native selectors and detached publication" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport
] begin
    const U=LineCableModels.Units

    frequency=[50.0, 100.0]
    resistance=reshape([1.0, 2.0], 1, 1, 2) .* 1.0e-4
    inductance=reshape([3.0, 4.0], 1, 1, 2) .* 1.0e-7
    conductance=reshape([5.0, 6.0], 1, 1, 2) .* 1.0e-8
    capacitance=reshape([7.0, 8.0], 1, 1, 2) .* 1.0e-10
    angular=reshape(2π .* frequency, 1, 1, :)
    impedance=complex.(resistance, inductance .* angular)
    admittance=complex.(conductance, capacitance .* angular)
    parameters=LineParameters(impedance, admittance, frequency)

    @test observables(typeof(parameters)) == (
        frequencies, Z, Y, R, X, L, G, B, C,
        (Z, abs), (Z, angle), (Y, abs), (Y, angle)
    )
    @test @inferred(observe(parameters, Z, 1, 1, 1)) === impedance[1, 1, 1]
    @test @inferred(observe(parameters, R, 1, 1, 1)) === resistance[1, 1, 1]
    @test @inferred(observe(parameters, L, 1, 1, 1)) ≈ inductance[1, 1, 1]
    @test @observe(R[1, 1, :]) == (R, 1, 1, Colon())
    @test @observe((Z, abs)[:, :, :]) ==
          (Z, abs, Colon(), Colon(), Colon())
    @test @observe((Z, angle)[:, :, :]) ==
          (Z, angle, Colon(), Colon(), Colon())
    @test @observe(parameters, Z[1, 1, :]) == impedance[1, 1, :]
    @test @observe(parameters, (Z, abs)[1, 1, :]) == abs.(impedance[1, 1, :])
    @test @observe(parameters, (Z, angle)[1, 1, :]) == angle.(impedance[1, 1, :])
    i, j, samples=1, 1, 1:2
    @test @observe(parameters, L[i, j, samples]) ≈ inductance[1, 1, :]
    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@observe parameters Z[1, 1])
    )
    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@observe (Z, abs, angle)[1, 1, :])
    )
    observe(parameters, Z, 1, 1, 1)
    observe(parameters, R, 1, 1, 1)
    @test @allocated(observe(parameters, Z, 1, 1, 1)) == 0
    @test @allocated(observe(parameters, R, 1, 1, 1)) == 0
    @test Z(parameters, 1, 1, 1) === observe(parameters, Z, 1, 1, 1)
    @test L(parameters, 1, 1, :) ≈ observe(parameters, L, 1, 1, :)

    requests=(
        frequency = (frequencies, Colon()),
        resistance = (R, 1, 1, Colon()),
        phase = (Z, angle, 1, 1, Colon())
    )
    target=U.units(:milli, :ohm; per = (:kilo, :meter))
    published=observables(parameters, requests; units = (resistance = target,))
    @test keys(published) == keys(requests)
    for payload in values(published)
        @test keys(payload) == (:values, :quantity, :unit)
    end
    @test published.resistance.values ≈ 1.0e6 .* resistance[1, 1, :]
    @test published.resistance.unit == target
    @test published.phase.values ≈ rad2deg.(angle.(impedance[1, 1, :]))
    @test published.phase.quantity isa
          U.Quantity{(:series_impedance, :phase_angle)}
    repeated_quantity=observables(parameters,
        (
            first_resistance = (R, 1, 1, Colon()),
            second_resistance = (R, 1, 1, 1)
        ))
    @test keys(repeated_quantity) == (:first_resistance, :second_resistance)
    prefixed=observables(parameters, (resistance = (R, 1, 1, Colon()),);
        units = (resistance = :micro,))
    @test prefixed.resistance.unit ==
          U.units(:micro, :ohm; per = (:kilo, :meter))
    @test LineCableModels.Grammar.validate_observables(
        parameters,
        requests,
        (resistance = target,)
    ) == (frequencies, R, (Z, angle))
    @test Base.ispublic(LineCableModels.Grammar, :validate_observables)
    @test !isdefined(LineCableModels, :validate_observables)

    published.frequency.values[1]=0.0
    published.resistance.values[1]=0.0
    @test frequencies(parameters)[1] == 50.0
    @test R(parameters, 1, 1, 1) == resistance[1, 1, 1]

    @test_throws MethodError observables(parameters)
    @test_throws MethodError observables(parameters, Dict(:resistance => R))
    @test_throws ArgumentError observables(parameters, (invalid = identity,))
    @test_throws ArgumentError observables(
        parameters,
        (resistance = R,);
        units = (unknown = target,)
    )
    @test_throws ArgumentError observables(
        parameters,
        (resistance = R,);
        units = (resistance = U.units(:base, :farad),)
    )

    series=SeriesImpedance(impedance)
    shunt=ShuntAdmittance(admittance)
    @test observe(series, Z) === series.values
    @test observe(series, Z, abs, 1, 1, :) == abs.(impedance[1, 1, :])
    @test observe(series, L, frequency) ≈ inductance
    @test_throws DimensionMismatch observe(series, L, [50.0])
    @test_throws DomainError observe(series, L, [0.0, 100.0])
    @test observe(shunt, Y) === shunt.values
    @test observe(shunt, Y, angle, 1, 1, :) == angle.(admittance[1, 1, :])
    @test observe(shunt, C, frequency) ≈ capacitance
    @test_throws DimensionMismatch observe(shunt, C, [50.0])
    @test_throws DomainError observe(shunt, C, [50.0, 0.0])
end
