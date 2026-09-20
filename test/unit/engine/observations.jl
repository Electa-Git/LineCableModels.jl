@testitem "Engine / observations / native selectors and detached publication" tags=[:unit] setup=[
    UseEngineSupport
] begin
    using LinearAlgebra: diag
    using Statistics: mean

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
    # Three selectors carry owner-defined statistical requests. The macro
    # preserves them; whether a result supports the request belongs to dispatch.
    @test @observe((statistics, R, mean)[1, 1]) == (statistics, R, mean, 1, 1)
    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@observe (Z, abs, angle, identity)[1, 1, :])
    )
    observe(parameters, Z, 1, 1, 1)
    observe(parameters, R, 1, 1, 1)
    @test @allocated(observe(parameters, Z, 1, 1, 1)) == 0
    @test @allocated(observe(parameters, R, 1, 1, 1)) == 0
    @test Z(parameters, 1, 1, 1) === observe(parameters, Z, 1, 1, 1)
    @test L(parameters, 1, 1, :) ≈ observe(parameters, L, 1, 1, :)

    requests=((R,1,1,:),(X,1,1,:))
    target=U.units(:milli,:ohm;per=(:kilo,:meter))
    observed=observables(parameters,requests;quantity_units=(R=target,))
    @test observed isa ObservedResult
    @test observe(observed,R)≈1.0e6.*resistance[1,1,:]
    @test first(observed.quantities).unit==target
    polar=ObservedResult(parameters,((Z,abs,1,1,:),(Z,angle,1,1,:)))
    @test observe(polar,Z,angle)≈rad2deg.(angle.(impedance[1,1,:]))
    prefixed=ObservedResult(parameters,requests;quantity_units=(R=:micro,))
    @test first(prefixed.quantities).unit==U.units(:micro,:ohm;per=(:kilo,:meter))
    @test LineCableModels.Grammar.validate_observables(parameters,requests)==(R,X)
    @test Base.ispublic(LineCableModels.Grammar,:validate_observables)
    @test !isdefined(LineCableModels,:validate_observables)
    first(observed.quantities).values[1]=0
    first(observed.quantities).coordinates.frequencies[1]=0
    @test frequencies(parameters)[1]==50
    @test R(parameters,1,1,1)==resistance[1,1,1]
    @test observables(parameters) isa ObservedResult
    @test_throws MethodError observables(parameters,Dict(:resistance=>R))
    @test_throws MethodError observables(parameters,(invalid=identity,))
    @test_throws DimensionMismatch observables(parameters,requests;units=(target,))
    @test_throws ArgumentError observables(parameters,requests;quantity_units=(R=U.units(:base,:farad),))
    @test LineCableModels.Grammar.detach(eps(Float64)/2,1.0)==eps(Float64)/2

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

@testitem "Engine / observations / detached records and explicit tables" tags=[:unit] begin
    using DataFrames, Tables
    using LineCableModels.ReportBuilder: tabulate
    const GR=LineCableModels.Grammar
    constants=CableConstants([:a,:b],[1e-4,2e-4],[3e-7,4e-7],[5e-10,6e-10],[7e-9,8e-9],50.)
    observed=observables(constants;length_unit=:base,quantity_units=:base)
    @test !Tables.istable(typeof(observed))
    @test !hasproperty(observed,:columns)
    @test length(observed.quantities)==4
    @test tabulate(observed,R).value==constants.R
    @test tabulate(observed,L).value==constants.L
    @test tabulate(observed,R).assembly==[1,2]
    table=tabulate(observed,R)
    table.value[1]=99
    @test constants.R==[1e-4,2e-4]
    @test observe(observed,R)==constants.R
    @test size(DataFrame(observed),1)==8
    @test sprint(show,MIME"text/plain"(),observed)==sprint(show,observed)
    @test occursin("4 quantities",sprint(show,observed))
    @test_throws MethodError iterate(observed)
    nullable=Union{Missing,Float64}[missing,eps()/4,Inf,-Inf,NaN,1.]
    @test isequal(GR.detach(nullable,2.0),[missing,eps()/2,Inf,-Inf,NaN,2.])
end
