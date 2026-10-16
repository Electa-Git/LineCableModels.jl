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

@testitem "Commons / detached thresholds retain precision and uncertainty dependencies" tags=[:unit, :commons] begin
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

@testitem "Commons / frequency and time sample coordinates" tags=[:unit, :commons] begin
    import LineCableModels.Commons as G
    import LineCableModels.Units as U
    product=G.observation_quantity(TransientResponse([0.,1.,2.],ones(2,3,3),
        zeros(2,3,3),["a","b","c"],ComputationDetails()),voltage)
    matrix=product.coordinates
    @test G._validate_observed_quantity(product)===nothing
    @test length(matrix.labels)==2 && length(matrix.column_labels)==3
    @test_throws DimensionMismatch G._validate_observed_quantity(merge(product,
        (coordinates=merge(matrix,(labels=["sending"],)),)))
    @test_throws DimensionMismatch G._validate_observed_quantity(merge(product,
        (coordinates=merge(matrix,(column_labels=["a","b"],)),)))
    for kind in (:matrix,:diagonal,:vector)
        c=kind===:matrix ? matrix : kind===:diagonal ?
            (kind,indices=(),rows=[1,2],columns=[1,2],samples=[1,2,3],
                extent=(2,2,3),labels=matrix.labels,domain=:LineEnd) :
            (kind,indices=(),axis=:end,axis_label="End",positions=[1,2],samples=[1,2,3],
                extent=(2,3),labels=matrix.labels,domain=:LineEnd)
        c=Base.structdiff(c,(times=nothing,time_unit=nothing))
        values=kind===:matrix ? ones(2,3,3) : ones(2,3)
        record=merge(product,(values,coordinates=c))
        @test_throws ArgumentError G._validate_observed_quantity(record)
        for (axis,unit_key,unit,points) in ((:frequencies,:frequency_unit,:hertz,[0.,1.,2.]),
                (:times,:time_unit,:second,[-1.,0.,1.]))
            coordinates=merge(c,NamedTuple{(axis,unit_key)}((points,U.units(:base,unit))))
            valid=merge(record,(;coordinates))
            @test G._validate_observed_quantity(valid)===nothing
            @test G._validate_observed_quantity(merge(valid,
                (coordinates=merge(coordinates,NamedTuple{(axis,)}((nothing,))),)))===nothing
            @test_throws DimensionMismatch G._validate_observed_quantity(merge(valid,
                (coordinates=merge(coordinates,NamedTuple{(axis,)}(([0.,1.],))),)))
            @test_throws ArgumentError G._validate_observed_quantity(merge(valid,
                (coordinates=merge(coordinates,NamedTuple{(axis,)}(([0.,Inf,2.],))),)))
            @test_throws ArgumentError G._validate_observed_quantity(merge(valid,
                (coordinates=merge(coordinates,NamedTuple{(unit_key,)}((U.units(:base,:meter),))),)))
            both=merge(coordinates,(times=[0.,1.,2.],time_unit=U.units(:base,:second),
                frequencies=[0.,1.,2.],frequency_unit=U.units(:base,:hertz)))
            @test_throws ArgumentError G._validate_observed_quantity(merge(valid,(coordinates=both,)))
        end
    end
end
