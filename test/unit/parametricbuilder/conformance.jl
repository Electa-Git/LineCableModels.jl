@testitem "ParametricBuilder / Gridspace / composition conformance" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    mutable struct CountingGrid{V <: Tuple}<:PB.AbstractGrid
        vals::V
        iterate_calls::Int
        length_calls::Int
    end
    CountingGrid(values::Tuple) = CountingGrid(values, 0, 0)
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
    (::DuplicateValue)(value) = (value, value)

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
        draw=PB.realize(rng, duplicated_point, :normal)
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

@testitem "ParametricBuilder / Gridspace / inference and allocation contracts" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    struct ScalarTarget
        value::Float64
    end
    struct BuildScalarTarget end
    (::BuildScalarTarget)(value) = ScalarTarget(value)

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
            total+=PB.realize(rng, point, :normal).value
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
    @test @inferred(PB.realize(rng, uncertain_point, :normal)) isa ScalarTarget

    deterministic_sum(deterministic, 1)
    realization_sum(rng, uncertain_point, 1)
    @test @allocated(deterministic_sum(deterministic, 10_000)) == 0
    rng=Random.Xoshiro(0x1234)
    @test @allocated(realization_sum(rng, uncertain_point, 10_000)) == 0
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
end
