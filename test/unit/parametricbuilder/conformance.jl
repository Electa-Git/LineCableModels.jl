@testitem "ParametricBuilder / Gridspace / composition conformance" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    mutable struct CountingGrid{V <: Tuple}<:PB.AbstractGrid
        vals::V
        iterate_calls::Int
        length_calls::Int
    end
    CountingGrid(values::Tuple)=CountingGrid(values, 0, 0)
    function Base.iterate(grid::CountingGrid, state...)
        grid.iterate_calls+=1
        return iterate(grid.vals, state...)
    end
    function Base.length(grid::CountingGrid)
        grid.length_calls+=1
        return length(grid.vals)
    end

    product_left=PB.Grid((1, 2, 3))
    product_right=PB.Grid((10, 20))
    product=PB.Gridspace{Tuple}(tuple, (product_left, product_right))
    @test collect(product) == [
        (1, 10), (2, 10), (3, 10),
        (1, 20), (2, 20), (3, 20)
    ]

    reused=PB.Grid((1, 2))
    reused_product=PB.Gridspace{Tuple}(tuple, (reused, reused))
    separate_product=PB.Gridspace{Tuple}(
        tuple,
        (PB.Grid((1, 2)), PB.Grid((1, 2)))
    )
    @test collect(reused_product) == [(1, 1), (2, 1), (1, 2), (2, 2)]
    @test collect(reused_product) == collect(separate_product)
    @test collect(PB.Gridspace{Tuple}(
        tuple,
        (reused, reused);
        combine = :zip
    )) == [(1, 1), (2, 2)]

    nested_zip=PB.Gridspace{Tuple}(
        tuple,
        (
            PB.Gridspace{Tuple}(
                tuple,
                (PB.Grid((1, 2, 3)), PB.Grid((10, 20, 30)));
                combine = :zip
            ),
            PB.Grid(:fixed)
        );
        combine = :zip
    )
    @test collect(nested_zip) == [
        ((1, 10), :fixed),
        ((2, 20), :fixed),
        ((3, 30), :fixed)
    ]

    left=CountingGrid(Tuple(1:2_000))
    right=CountingGrid(Tuple(2_001:4_000))
    counted_product=PB.Gridspace{Tuple}(tuple, (left, right))
    @test eltype(counted_product) === Any
    @test Base.IteratorEltype(typeof(counted_product)) isa Base.EltypeUnknown
    @test length(counted_product) == 4_000_000
    @test left.iterate_calls == 0
    @test right.iterate_calls == 0

    zipped=PB.Gridspace{Tuple}(tuple, (left, right); combine = :zip)
    left.iterate_calls=0
    right.iterate_calls=0
    @test length(zipped) == 2_000
    @test left.iterate_calls == 0
    @test right.iterate_calls == 0
    @test length(collect(zipped)) == 2_000
    @test 2_000 <= left.iterate_calls <= 2_001
    @test 2_000 <= right.iterate_calls <= 2_001

    materializations=Ref(0)
    build=(values...)->(materializations[]+=1; values)
    @test_throws DimensionMismatch PB.Gridspace{Tuple}(
        build,
        (PB.Grid((1, 2)), PB.Grid((10, 20, 30)));
        combine = :zip
    )
    @test materializations[] == 0

    @test_throws MethodError PB.Grid((1, 2); key = :shared)
    @test_throws ArgumentError PB.Gridspace{Tuple}(tuple, ((1, 2),))
    @test !applicable(getindex, product, 1)
    @test rand(Random.Xoshiro(0x1234), product) in collect(product)
end

@testitem "ParametricBuilder / Gridspace / structural realization" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    using Measurements
    import LineCableModels.ParametricBuilder as PB

    struct DuplicateValue end
    (::DuplicateValue)(value)=(value, value)

    duplicated=PB.Gridspace{Tuple}(
        DuplicateValue(),
        (PB.Grid(10.0, PB.AbsoluteError(0.5)),)
    )
    duplicated_point=first(PB.points(duplicated))
    propagated=PB.materialize(duplicated_point)
    @test propagated[1] === propagated[2]
    @test Measurements.cov(propagated[1], propagated[2]) ==
          Measurements.uncertainty(propagated[1])^2

    rng=MersenneTwister(42)
    @test all(1:32) do _
        draw=LineCableModels.realize(rng, duplicated_point, :normal)
        draw[1] === draw[2]
    end

    independent=PB.Gridspace{Tuple}(
        tuple,
        (
            PB.Grid(10.0, PB.AbsoluteError(0.5)),
            PB.Grid(10.0, PB.AbsoluteError(0.5))
        );
        combine = :zip
    )
    independent_values=PB.materialize(first(PB.points(independent)))
    @test independent_values[1] !== independent_values[2]
    @test iszero(Measurements.cov(independent_values...))

    reused=PB.Grid(10.0, PB.AbsoluteError(0.5))
    reused_values=PB.materialize(first(PB.points(PB.Gridspace{Tuple}(
        tuple,
        (reused, reused);
        combine = :zip
    ))))
    @test reused_values[1] !== reused_values[2]
    @test iszero(Measurements.cov(reused_values...))
end

@testitem "ParametricBuilder / Gridspace / supplied uncertainty realization" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    events=Symbol[]
    child=PB.Gridspace{Tuple}(
        (left, right)->begin
            push!(events, :child)
            (left, right)
        end,
        (
            PB.Grid(10.0, PB.AbsoluteError(0.5)),
            PB.Grid(20.0, PB.AbsoluteError(0.0))
        )
    )
    parent=PB.Gridspace{Tuple}(
        (nested, value, tag)->begin
            push!(events, :parent)
            (nested, value, tag)
        end,
        (
            child,
            PB.Grid(30.0, PB.AbsoluteError(3.0)),
            PB.Grid((:tag,))
        )
    )
    point=first(LineCableModels.points(parent))
    descriptors=@inferred LineCableModels.uncertainties(point)
    @test nominal.(descriptors) == (10.0, 20.0, 30.0)
    @test PB.uncertainty.(descriptors) == (0.5, 0.0, 3.0)

    arguments=@inferred LineCableModels.realize_arguments(
        point,
        (11.0, 20.0, 33.0)
    )
    @test events == [:child]
    supplied=@inferred LineCableModels.realize(point, arguments)
    @test events == [:child, :parent]
    @test supplied == ((11.0, 20.0), 33.0, :tag)

    empty!(events)
    @test_throws DimensionMismatch LineCableModels.realize_arguments(
        point,
        (11.0, 20.0)
    )
    @test isempty(events)
    @test LineCableModels.realize(
        Xoshiro(0x1234),
        point,
        :normal
    ) == (
        (10.541923815394892, 20.0),
        32.2168218083917,
        :tag
    )
    @test LineCableModels.realize(
        Xoshiro(0x1234),
        point,
        (rng, mean, sigma)->mean + sigma
    ) == ((10.5, 20.0), 33.0, :tag)
    comparison=LineCableModels.realize_arguments(
        point,
        (10.5, 20.0, 33.0)
    )
    @test LineCableModels.realize(point, comparison) ==
          ((10.5, 20.0), 33.0, :tag)

    duplicated=PB.Gridspace{Tuple}(
        value->(value, value),
        (PB.Grid(1.0, PB.AbsoluteError(0.1)),)
    )
    @test length(LineCableModels.uncertainties(
        first(LineCableModels.points(duplicated))
    )) == 1

    reused=PB.Grid(1.0, PB.AbsoluteError(0.1))
    repeated=PB.Gridspace{Tuple}(tuple, (reused, reused); combine = :zip)
    repeated_descriptors=LineCableModels.uncertainties(
        first(LineCableModels.points(repeated))
    )
    @test length(repeated_descriptors) == 2
    @test repeated_descriptors[1] == repeated_descriptors[2]

    opaque=(PB.UncertainValue(2.0, 0.2), :opaque)
    atomic=PB.Gridspace{Tuple}(tuple, (PB.Grid((opaque,)),))
    @test isempty(LineCableModels.uncertainties(
        first(LineCableModels.points(atomic))
    ))

    struct PositiveTarget<:AbstractProblemDefinition
        value::Float64
    end
    function LineCableModels.validate(value::PositiveTarget)
        value.value > 0 || throw(DomainError(value.value, "value must be positive"))
        return value
    end
    positive=PB.Gridspace{PositiveTarget}(
        PositiveTarget,
        (PB.Grid(1.0, PB.AbsoluteError(0.1)),)
    )
    positive_point=first(LineCableModels.points(positive))
    invalid_arguments=LineCableModels.realize_arguments(positive_point, (-1.0,))
    @test_throws DomainError LineCableModels.realize(positive_point, invalid_arguments)
end

@testitem "ParametricBuilder / Gridspace / inference and allocation contracts" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    struct ScalarTarget
        value::Float64
    end
    struct BuildScalarTarget end
    (::BuildScalarTarget)(value)=ScalarTarget(value)

    function deterministic_sum(space, repetitions)
        total=0.0
        for _ in 1:repetitions, value in space

            total+=value.value
        end
        return total
    end

    function realization_sum(rng, point, repetitions)
        total=0.0
        for _ in 1:repetitions
            total+=LineCableModels.realize(rng, point, :normal).value
        end
        return total
    end

    function supplied_sum(point, repetitions)
        total=0.0
        for _ in 1:repetitions
            arguments=LineCableModels.realize_arguments(point, (2.0,))
            total+=LineCableModels.realize(point, arguments).value
        end
        return total
    end

    deterministic=PB.Gridspace{ScalarTarget}(
        BuildScalarTarget(),
        (PB.Grid((1.0, 2.0, 3.0)),)
    )
    deterministic_point=first(PB.points(deterministic))
    @test eltype(deterministic) === ScalarTarget
    @test Base.IteratorEltype(typeof(deterministic)) isa Base.HasEltype
    @test @inferred(first(deterministic)) == ScalarTarget(1.0)
    @test @inferred(PB.materialize(deterministic_point)) == ScalarTarget(1.0)

    uncertain=PB.Gridspace{ScalarTarget}(
        BuildScalarTarget(),
        (PB.Grid(1.0, PB.AbsoluteError(0.1)),)
    )
    uncertain_point=first(PB.points(uncertain))
    @test eltype(uncertain) === Any
    @test Base.IteratorEltype(typeof(uncertain)) isa Base.EltypeUnknown
    rng=Random.Xoshiro(0x1234)
    @test @inferred(LineCableModels.realize(rng, uncertain_point, :normal)) isa ScalarTarget
    @test @inferred(LineCableModels.uncertainties(uncertain_point)) ==
          (PB.UncertainValue(1.0, 0.1),)
    supplied_arguments=@inferred LineCableModels.realize_arguments(
        uncertain_point,
        (2.0,)
    )
    @test @inferred(LineCableModels.realize(
        uncertain_point,
        supplied_arguments
    )) == ScalarTarget(2.0)

    deterministic_sum(deterministic, 1)
    realization_sum(rng, uncertain_point, 1)
    supplied_sum(uncertain_point, 1)
    @test @allocated(deterministic_sum(deterministic, 10_000)) == 0
    rng=Random.Xoshiro(0x1234)
    @test @allocated(realization_sum(rng, uncertain_point, 10_000)) == 0
    @test @allocated(supplied_sum(uncertain_point, 10_000)) == 0
end

@testitem "ParametricBuilder / Gridspace / removed traversal machinery" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    import LineCableModels.ParametricBuilder as PB

    for name in (
        :AutomaticGridKey,
        :NamedGridKey,
        :ConstantAxis,
        :GridBinding,
        :AxisSelection,
        :ResolvedGridValue,
        :Configuration,
        :configurations,
        :configuration_manifest,
        :_same_grid_key,
        :_axis_bindings,
        :_compatible_bindings,
        :_merged_bindings,
        :_gridspace_axis,
        :_AbstractDefinition,
        :_ManifestState,
        :_manifest_tree,
        :_coupling_manifest
    )
        @test !isdefined(PB, name)
        @test !isdefined(LineCableModels, name)
    end
    @test :Gridpoint ∉ names(PB)
    @test :points ∉ names(PB)
    @test :realize ∉ names(PB)
    @test :uncertainties ∉ names(PB)
    @test !isdefined(PB, :realize)
    @test !isdefined(PB, :realize_arguments)
    @test !isdefined(PB, :uncertainties)
    @test !Base.isexported(LineCableModels, :realize)
    @test !Base.isexported(LineCableModels, :realize_arguments)
    @test !Base.isexported(LineCableModels, :uncertainties)
    @test Base.ispublic(LineCableModels, :realize)
    @test Base.ispublic(LineCableModels, :realize_arguments)
    @test Base.ispublic(LineCableModels, :uncertainties)
end
