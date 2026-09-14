@testitem "Engine / observation resolution is physical and shared with RMS" tags=[:unit] begin
    using LinearAlgebra: diag
    using LineCableModels.Engine: compare
    f = [1.0, 1e3, 1e7]
    z = fill(1.0 + im, 2, 2, 3)
    y = fill(1e-13 + 1e-18im, 2, 2, 3)
    y[2, 2, :] .= 1e-4 .+ 2π .* f .* 1e-9im
    source = LineParameters(PhaseDomain, z, y, f)
    original = deepcopy((z, y, f))
    request = ((G, 1, 1, :),)
    a = observables(source, request; length_unit=:base)
    b = observables(source, request; length_unit=:kilo)
    @test all(iszero, only(a).values)
    @test only(b).values == 1000 .* only(a).values
    @test only(observables(source, request; clip=false, length_unit=:base)).values == real.(y[1, 1, :])
    @test all(iszero, only(observables(source, ((B, 1, 1, :),))).values)
    @test all(iszero, only(observables(source, ((Y, 1, 1, :),))).values)
    @test all(ismissing, only(observables(source, ((Y, angle, 1, 1, :),))).values)
    @test only(observables(source, request; atol=0, length_unit=:base)).values == real.(y[1, 1, :])
    @test only(observables(source, ((G, :, :, :),); length_unit=:base)).values[1, 1, :] == only(a).values
    @test only(observables(source, ((G, diag, :, :),); length_unit=:base)).values[1, :] == only(a).values
    @test only(observables(source, ((G, 1, 1, 2:3),); length_unit=:base)).values == only(a).values[2:3]
    @test a.metadata.observation_columns.G.resolution.kind === :declared_floor
    @test a.metadata.observation_columns.G.resolution.unresolved_count == 3
    @test_throws ArgumentError observables(source, ((R, 1, 1, :), (G, 1, 1, :)); atol=1e-12)
    @test_throws ArgumentError observables(source, request; atol=(oops=0.0,))
    for quantity in (R, X, L, G, B, C, Z, Y)
        error = compare(source, source, quantity)
        @test error.details.resolution.revision == LineCableModels.Engine.OBSERVABLE_RESOLUTION_REVISION
        @test error.details.resolution.unit == LineCableModels.Units.native_unit(quantity, :pul)
    end
    @test compare(source, source, B).details.atol ≈ 2π .* f .* compare(source, source, C).details.atol
    @test compare(source, source, X).details.atol ≈ 2π .* f .* compare(source, source, L).details.atol
    @test (z, y, f) == original
end

@testitem "Engine / resolution boundaries preserve signals, types and raw absolute errors" tags=[:unit] begin
    using LineCableModels.Engine: compare
    using LinearAlgebra: norm
    for T in (Float32, Float64, BigFloat)
        cutoff = T(1e-12)
        f = T[1, 1e3, 1e7]
        g = reshape(T[cutoff/2, cutoff, 2cutoff], 1, 1, :)
        z = fill(complex(one(T), one(T)), 1, 1, 3)
        y = complex.(g, zero(T))
        a = LineParameters(PhaseDomain, z, y, f)
        b = LineParameters(PhaseDomain, z, 2y, f)
        original = deepcopy((z, y, f))
        request = ((G, 1, 1, :),)
        base = only(observables(a, request; length_unit=:base, atol=(G=cutoff,))).values
        kilo = only(observables(a, request; length_unit=:kilo, atol=(G=cutoff,))).values
        @test eltype(base) === T
        @test base == T[0, 0, 2cutoff]
        @test kilo ≈ 1000base
        for normalization in (:reference_rms, :pointwise), (left, right) in ((a, b), (b, a))
            error = compare(left, right, G; normalization, atol=cutoff)
            @test ismissing(only(error.relative))
            @test only(error.absolute) ≈ norm(vec(g)) / sqrt(T(3))
            @test Base.nonmissingtype(eltype(error.absolute)) === T
            @test only(error.details.unresolved_samples).reference > 0
            @test only(error.details.unresolved_samples).candidate > 0
            @test occursin("no samples were omitted", only(error.details.normalization_reason))
        end
        @test !ismissing(only(compare(a, b, G; band=:wide, atol=cutoff).relative))
        @test only(observables(a, ((G, 1, 1, 3:3),); length_unit=:base, atol=cutoff)).values == base[3:3]
        total = LineParameters(PhaseDomain, 10z, 10y, f; basis=:total)
        @test only(observables(total, request; atol=10cutoff)).values ≈ 10base
        @test compare(a, b, G; atol=0).relative ≈ fill(one(T), 1, 1)
        @test (z, y, f) == original
        # Stable norm evaluation must not square tiny Float32 data into zero.
        small = fill(T(1e-30), 1, 1, 3)
        @test only(compare(small, 2small).absolute) > 0
        @test only(compare(small, 2small).relative) ≈ one(T)
    end
    f = [1.0, 10.0, 100.0]
    z = fill(1.0 + im, 2, 2, 3)
    y = fill(2e-12 + 1e-6im, 2, 2, 3)
    a = LineParameters(PhaseDomain, z, copy(y), f)
    b = LineParameters(PhaseDomain, z, copy(y), f)
    b.Y.values[2, 2, :] .= 1e20
    request = ((G, 1, 1, :),)
    @test only(observables(a, request)).values == only(observables(b, request)).values
    @test all(>(0), only(observables(a, request)).values)
    standalone = ShuntAdmittance(fill(1e-18im, 2, 2, 3))
    raw = observables(standalone, ((Y, 1, 1, :),); length_unit=:base)
    @test only(raw).values == fill(1e-18im, 3)
    @test raw.metadata.observation_columns.Y.resolution.kind === :unassessed
    @test all(iszero, only(observables(standalone, ((Y, 1, 1, :),); frequencies=f)).values)
    @test all(iszero, only(observables(standalone, ((Y, 1, 1, :),); atol=1e-12)).values)
    @test_throws ArgumentError observables(a, request; frequencies=[1., 2., 3.])
    @test_throws DimensionMismatch observables(standalone, request; frequencies=[1.])
    for invalid in (-1, Inf, NaN, true, "bad", (G=-1.0,))
        @test_throws ArgumentError observables(a, request; atol=invalid)
    end
    a.Y.values[1, 1, 1] = NaN + im
    @test isnan(first(only(observables(a, request)).values))
    @test_throws ArgumentError compare(a, b, G)
    zerosource = LineParameters(PhaseDomain, z, y, [0., 10., 100.])
    @test_throws DomainError observables(zerosource, ((C, 1, 1, 1),))
    @test all(isfinite, only(observables(zerosource, ((C, 1, 1, 2:3),))).values)
    # Different operand precisions must use the same respective decisions as
    # their own publications, including a value between rounded cutoffs.
    single = LineParameters(PhaseDomain, ones(ComplexF32,1,1,3),
        fill(ComplexF32(1e-6),1,1,3), f)
    double = LineParameters(PhaseDomain, ones(ComplexF64,1,1,3),
        fill(complex((Float64(Float32(1e-12)) + 1e-12)/2),1,1,3), f)
    @test all(iszero, only(observables(double, ((G,1,1,:),))).values)
    @test ismissing(only(compare(single, double, G).relative))
    @test ismissing(only(compare(double, single, G).relative))
    @test_throws ArgumentError observables(single, ((G,1,1,:),); atol=1e100)
end

@testitem "UQ / detachment retains small sample spread and ordering" tags=[:unit] begin
    summary = LineCableModels.UQ.SampleSummary([1e-18, 2e-18, 3e-18])
    detached = LineCableModels.Grammar.detach(summary, 1000.0, true)
    @test detached.std == 1000summary.std > 0
    @test detached.min <= detached.q05 <= detached.median <= detached.q95 <= detached.max
    @test detached.mean == 1000summary.mean
    @test detached.n == summary.n
end

@testitem "Engine / resolution does not erase physical uncertainty" tags=[:extension] begin
    using Measurements
    a = measurement(1e-18, 2e-19)
    b = 2a
    detached = LineCableModels.Grammar.detach([a, b], 1000.0, true)
    @test Measurements.value(detached[1]) == 1000 * Measurements.value(a)
    @test Measurements.uncertainty(detached[1]) == 1000 * Measurements.uncertainty(a)
    @test iszero(Measurements.uncertainty(detached[2] - 2detached[1]))
    @test iszero(Measurements.uncertainty(detached[1] - 1000a))
    source = LineParameters(PhaseDomain, fill(complex(a, b), 1, 1, 2),
        fill(complex(a, b), 1, 1, 2), [1.0, 2.0])
    result = only(observables(source, ((G, 1, 1, :),); length_unit=:base)).values
    @test all(iszero, Measurements.value.(result))
    @test all(iszero, Measurements.uncertainty.(result))
    raw = only(observables(source, ((G, 1, 1, :),); length_unit=:base, clip=false)).values
    @test all(value -> iszero(Measurements.uncertainty(value - a)), raw)
    @test all(value -> Measurements.value(value) == Measurements.value(a), raw)
end

@testitem "Engine / nominal and spread resolution are independent projections" tags=[:extension] begin
    using Measurements
    using LinearAlgebra: diag
    using LineCableModels.Grammar: observation_resolution
    using LineCableModels.Engine: compare
    for T in (Float32, Float64, BigFloat)
        cutoff = T(1e-12)
        values = measurement.(T[cutoff/2, cutoff/2, 2cutoff, 2cutoff],
            T[cutoff/2, 2cutoff, cutoff/2, 2cutoff])
        y = reshape(complex.(values, zero.(values)), 1, 1, :)
        f = T[1, 10, 100, 1000]
        source = LineParameters(PhaseDomain, copy(y), y, f)
        saved = deepcopy((Z(source), Y(source), source.f))
        before = compare(source, source, G)
        for length_unit in (:base, :kilo)
            factor = length_unit === :base ? T(1) : T(1000)
            publication = observables(source, ((G, 1, 1, :),); length_unit, atol=cutoff)
            result = only(publication).values
            @test eltype(result) === eltype(values)
            @test Measurements.value.(result) ≈ factor .* T[0, 0, 2cutoff, 2cutoff]
            @test Measurements.uncertainty.(result) ≈ factor .* T[0, 2cutoff, 0, 2cutoff]
            @test iszero(Measurements.uncertainty(result[2] - factor * values[2]))
            @test iszero(Measurements.uncertainty(result[4] - factor * values[4]))
            metadata = publication.metadata.observation_columns.G.resolution
            @test metadata.unresolved_count == metadata.uncertainty_unresolved_count == 2
            raw = only(observables(source, ((G, 1, 1, :),); length_unit, clip=false)).values
            @test all(iszero, Measurements.uncertainty.(raw .- factor .* values))
            @test Measurements.value.(raw) == factor .* Measurements.value.(values)
        end
        clean = only(observables(source, ((G, 1, 1, :),); length_unit=:base, atol=cutoff)).values
        diagonal = only(observables(source, ((G, diag, :, :),); length_unit=:base, atol=cutoff)).values
        @test vec(diagonal) == clean
        @test only(observables(source, ((G, 1, 1, 2),); length_unit=:base, atol=cutoff)).values == clean[2]
        @test only(observables(source, ((G, 1, 1, :),); length_unit=:base, atol=0)).values == values
        after = compare(source, source, G)
        @test isequal(before.absolute, after.absolute) && isequal(before.relative, after.relative)
        @test isequal((Z(source), Y(source), source.f), saved)
    end
    # A complex uncertainty has two components; the generic Number uncertainty
    # fallback must not classify it as deterministic. Phase uses nominal magnitude.
    value = complex(measurement(1e-18, 2e-12), measurement(1e-18, 3e-12))
    resolution = observation_resolution(value, Y; atol=1e-12)
    @test resolution.unresolved && !resolution.uncertainty_unresolved
    @test observation_resolution(value, Y; atol=4e-12).uncertainty_unresolved
    source = LineParameters(fill(value, 1, 1, 1), fill(value, 1, 1, 1), [1.0])
    @test ismissing(only(only(observables(source, ((Y, angle, 1, 1, :),); atol=1e-12)).values))
    clean = only(only(observables(source, ((Y, 1, 1, :),); atol=1e-12, length_unit=:base)).values)
    @test nominal(clean) == 0
    @test uncertainty(real(clean) - real(value)) == uncertainty(imag(clean) - imag(value)) == 0
    for value in (NaN, Inf, missing)
        nonfinite_resolution = observation_resolution(Union{Missing,Float64}[value], G)
        @test !only(nonfinite_resolution.unresolved) && !only(nonfinite_resolution.uncertainty_unresolved)
    end
    unassessed = observation_resolution(nothing, identity)
    @test unassessed.unresolved === unassessed.uncertainty_unresolved === nothing
end
