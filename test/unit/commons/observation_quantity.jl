@testitem "Commons / fallback observations preserve inferred values and scalar counts" tags=[:unit, :commons] begin
    const G = LineCableModels.Commons
    struct ObservationSource{T}
        values::T
    end
    G.observe(source::ObservationSource, ::typeof(R)) = source.values
    G.observe(source::ObservationSource, ::typeof(R), index::Int) = source.values[index]
    G.basis(::ObservationSource) = :pul

    for value in (Float32(1), 1.0, big"1.0", missing, [1.0, 2.0])
        source = ObservationSource(value)
        product = @inferred G.observation_quantity(source, R)
        @test (@inferred G._validate_observed_quantity(product)) === nothing
        @test product.coordinates.extent == (value isa AbstractArray ? size(value) : ())
        @test typeof(product.values) === typeof(value)
        value isa AbstractArray && @test product.values !== value
        expected = G.detach(value, LineCableModels.Units.scale_factor(
            LineCableModels.Units.native_unit(product.quantity, :pul), product.unit))
        @test isequal(product.values, expected)
        malformed = merge(product, (coordinates=merge(product.coordinates, (extent=(3,),)),))
        @test_throws DimensionMismatch G._validate_observed_quantity(malformed)
    end
    # Native map may narrow a union element type based on the actual values.
    # Preserve missing values without requiring one concrete array return type.
    mixed = G.observation_quantity(ObservationSource([missing, 2.0]), R)
    @test ismissing(first(mixed.values))
    @test mixed.values[2] == G.observation_quantity(ObservationSource(2.0), R).values
    @test (@inferred G._validate_observed_quantity(mixed)) === nothing
    source = ObservationSource([1.0, 2.0])
    selected = @inferred G.observation_quantity(source, (R, 2))
    scalar = @inferred G.observation_quantity(ObservationSource(2.0), R)
    @test selected.values == scalar.values
    @test selected.coordinates.indices == (2,)
    @test (@inferred G._validate_observed_quantity(selected)) === nothing
    for request in (:R, [R], nothing)
        @test_throws ArgumentError G.observation_quantity(source, request)
    end
    @test_throws MethodError G.observation_quantity(source, ())
    @test_throws BoundsError G.observation_quantity(source, (R, 3))
end

@testitem "Commons / detached thresholds retain precision and uncertainty dependencies" tags=[:unit, :measurements] begin
    using Measurements
    const G = LineCableModels.Commons
    shared = measurement(0.5, 0.01)
    cutoffs = (real=Float32[0.25, 0.5], imaginary=(value=shared, absent=nothing),
        unavailable=missing, complex=ComplexF32(0.25, 0.5))
    detached = @inferred G.detach(cutoffs, 1000.0)
    @test detached.real == Float32[250, 500]
    @test eltype(detached.real) === Float32
    @test detached.real !== cutoffs.real
    @test detached.complex === ComplexF32(250, 500)
    @test detached.imaginary.absent === nothing
    @test detached.unavailable === missing
    @test uncertainty(detached.imaginary.value - 1000shared) == 0
    detached.real[1] = 0
    @test cutoffs.real[1] == 0.25f0
end

@testitem "Commons / observation equality retains uncertainty identity and collection shape" tags=[:unit, :measurements] begin
    using Measurements
    const G = LineCableModels.Commons
    shared = measurement(1.0, 0.1)
    independent = measurement(1.0, 0.1)
    @test isequal(shared, independent)
    @test (@inferred G._same_observed_values(shared, shared))
    @test !(@inferred G._same_observed_values(shared, independent))
    for wrap in (identity, x -> [x], x -> (x, missing), x -> (value=x, absent=nothing),
            x -> Union{Missing,typeof(x)}[missing, x])
        @test G._same_observed_values(wrap(shared), wrap(shared))
        @test !G._same_observed_values(wrap(shared), wrap(independent))
    end
    @test G._same_observed_values([NaN, -0.0], [NaN, -0.0])
    @test !G._same_observed_values([0.0], [-0.0])
    @test !G._same_observed_values([1.0, 2.0], reshape([1.0, 2.0], 1, 2))
    @test !G._same_observed_values((x=shared,), (y=shared,))
    @test !G._same_observed_values((shared,), (shared, shared))
end
