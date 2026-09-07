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

@testitem "Engine / publications / native Tables and collection interfaces" tags=[:unit] begin
    using DataFrames
    const GR = LineCableModels.Grammar
    const Tables = GR.Tables
    constants = CableConstants([:a, :b], [1e-4, 2e-4], [3e-7, 4e-7],
        [5e-10, 6e-10], [7e-9, 8e-9], 50.0)
    published = observables(constants, (R, L, C, G);
        length_unit=:base, quantity_units=:base)

    @test Tables.istable(typeof(published))
    @test Tables.columnaccess(typeof(published))
    @test Tables.columnnames(published) == (:core, :R, :L, :C, :G)
    schema = @inferred Tables.schema(published)
    @test schema.names == (:core, :R, :L, :C, :G)
    @test schema.types == (Symbol, Float64, Float64, Float64, Float64)
    for (index, name) in enumerate(schema.names)
        column = Tables.getcolumn(published, name)
        @test Tables.getcolumn(published, index) === column
        @test column == getproperty(constants, name === :core ? :cores : name)
    end
    @test Tables.rowtable(published) == collect(constants)
    @test DataFrame(published) == DataFrame(Tables.columntable(published))
    @test length(published) == 4 # Observations, not table rows.
    @test firstindex(published) == 1
    @test lastindex(published) == 4
    @test first(published) === published[begin]
    @test last(published) === published[end]
    @test Tuple(published) == (published[1], published[2], published[3], published[4])
    @test Base.tail(published) == (published[2], published[3], published[4])
    @test occursin("2 rows", sprint(summary, published))
    @test occursin("2 rows × 5 columns", sprint(show, published))
    @test sprint(show, MIME"text/plain"(), published) == sprint(show, published)
    @test_throws BoundsError published[5]

    # Native table consumers can mutate their detached columns, not the source.
    Tables.getcolumn(published, :R)[1] = 99.0
    Tables.getcolumn(published, :core)[1] = :changed
    @test constants.R == [1e-4, 2e-4]
    @test constants.cores == [:a, :b]
    @test published.metadata.basis === :pul
    @test published.metadata.row_order == (:core, :R, :L, :C, :G)
    @test keys(published.metadata.observation_columns) == (:R, :L, :C, :G)

    empty_table = GR.ObservationPublication((), (;),
        (basis=:pul, row_order=(), observation_columns=(;)))
    @test isempty(empty_table)
    @test isempty(Tables.columnnames(empty_table))
    @test occursin("0 rows", sprint(summary, empty_table))
    @test occursin("0 rows × 0 columns", sprint(show, empty_table))
    @test_throws DimensionMismatch GR.ObservationPublication((),
        (x=[1.0], y=[1.0, 2.0]), empty_table.metadata)
    @test_throws ArgumentError GR.ObservationPublication((), (;), (basis=:pul,))

    # Display clipping must not discard phase or silently convert missing data.
    complex_values = ComplexF64[eps()/2 + eps()/4 * im, 2 - 3im]
    @test GR.detach(complex_values, 2.0, true) == 2 .* complex_values
    nullable = Union{Missing, Float64}[missing, eps()/4, Inf, -Inf, NaN, 1.0]
    detached = GR.detach(nullable, 2.0, true)
    @test isequal(detached, [missing, 0.0, Inf, -Inf, NaN, 2.0])
    @test !ismissing(nullable[2]) && nullable[2] > 0
end
