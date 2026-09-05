@testitem "Engine / observations / native selectors and detached publication" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport
] begin
    using LinearAlgebra: diag

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
        (Z, abs), (Z, angle), (Y, abs), (Y, angle),
        (Z, diag), (Y, diag), (R, diag), (X, diag),
        (L, diag), (G, diag), (B, diag), (C, diag)
    )
    @test @inferred(observe(parameters, Z, 1, 1, 1)) === impedance[1, 1, 1]
    @test @inferred(observe(parameters, R, 1, 1, 1)) === resistance[1, 1, 1]
    @test @inferred(observe(parameters, L, 1, 1, 1)) ≈ inductance[1, 1, 1]
    @test @observe(R[1, 1, :]) == (R, 1, 1, Colon())
    @test @observe((Z, abs)[:, :, :]) ==
          (Z, abs, Colon(), Colon(), Colon())
    @test @observe((Z, angle)[:, :, :]) ==
          (Z, angle, Colon(), Colon(), Colon())
    @test LineCableModels.Grammar.request_identity(@observe(R[:, :, :])) === R
    @test LineCableModels.Grammar.request_indices(@observe(R[:, :, :])) ==
          (Colon(), Colon(), Colon())
    @test LineCableModels.Grammar.request_identity(@observe((Z, abs)[:, :, :])) ==
          (Z, abs)
    @test @observe(parameters, Z[1, 1, :]) == impedance[1, 1, :]
    @test @observe(parameters, (Z, abs)[1, 1, :]) == abs.(impedance[1, 1, :])
    @test @observe(parameters, (Z, angle)[1, 1, :]) == angle.(impedance[1, 1, :])
    i, j, samples=1, 1, 1:2
    @test @observe(parameters, L[i, j, samples]) ≈ inductance[1, 1, :]
    @test @observe(parameters, Z[1, 1]) == impedance[1, 1, :]
    @test @observe((Z, diag)[:, :]) == (Z, diag, Colon(), Colon())
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
        (frequencies, Colon()),
        (R, 1, 1, Colon()),
        (Z, angle, 1, 1, Colon())
    )
    target=U.units(:milli, :ohm; per = (:kilo, :meter))
    published=observables(parameters, requests; units = (nothing, target, nothing))
    @test @inferred(observables(parameters, ((R, 1, 1, Colon()),))) isa
          LineCableModels.Grammar.ObservationPublication
    @test published isa LineCableModels.Grammar.ObservationPublication
    for payload in published
        @test keys(payload) == (:values, :quantity, :unit)
    end
    @test published[2].values ≈ 1.0e6 .* resistance[1, 1, :]
    @test published[2].unit == target
    @test published[3].values ≈ rad2deg.(angle.(impedance[1, 1, :]))
    @test published[3].quantity isa
          U.Quantity{(:series_impedance, :phase_angle)}
    @test_throws ArgumentError observables(parameters,
        (
            (R, 1, 1, Colon()),
            (R, 1, 1, 1)
        ))
    prefixed=observables(parameters, ((R, 1, 1, Colon()),);
        units = (:micro,))
    @test only(prefixed).unit ==
          U.units(:micro, :ohm; per = (:kilo, :meter))
    @test LineCableModels.Grammar.validate_observables(
        parameters,
        requests,
        (nothing, target, nothing)
    ) == (frequencies, R, (Z, angle))
    @test Base.ispublic(LineCableModels.Grammar, :validate_observables)
    @test !isdefined(LineCableModels, :validate_observables)

    published[1].values[1]=0.0
    published[2].values[1]=0.0
    @test frequencies(parameters)[1] == 50.0
    @test R(parameters, 1, 1, 1) == resistance[1, 1, 1]

    @test_throws MethodError observables(parameters)
    @test_throws MethodError observables(parameters, Dict(:resistance => R))
    @test_throws MethodError observables(parameters, (invalid = identity,))
    @test_throws DimensionMismatch observables(
        parameters,
        (R,);
        units = (target, target)
    )
    @test_throws ArgumentError observables(
        parameters,
        (R,);
        units = (U.units(:base, :farad),)
    )
    @test LineCableModels.Grammar.detach(eps(Float64) / 2, 1.0, true) == 0.0
    @test LineCableModels.Grammar.detach(eps(Float64) / 2, 1.0, false) > 0

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
