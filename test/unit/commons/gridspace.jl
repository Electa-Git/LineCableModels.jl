@testitem "Grid / finite sources" tags=[:unit, :commons] setup=[
    UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    deterministic=PB.Grid((1, 3, 2))
    @test PB.Grid(deterministic) === deterministic
    @test collect(PB.Grid(1:3)) == [1, 2, 3]
    @test deterministic[2] == 3
    @test extrema(deterministic) == (1, 3)
    @test size(deterministic) == (3,)

    relative=PB.Grid((-10.0, 20.0), (5.0, 10.0))
    absolute=PB.Grid((10.0, 20.0), PB.AbsoluteError((0.5, 2.0)))
    @test length(relative) == 4
    @test relative[4] isa PB.UncertainValue
    @test extrema(relative) == (-11.0, 22.0)
    @test extrema(absolute) == (8.0, 22.0)
    @test_throws BoundsError relative[0]
    @test_throws BoundsError absolute[5]

    zero_sigma=PB.UncertainValue(4, 0)
    @test rand(MersenneTwister(1), zero_sigma) === 4.0
    uncertain=PB.UncertainValue(10.0, 2.0)
    @test isfinite(rand(MersenneTwister(10), uncertain; distribution = :normal))
    @test abs(rand(MersenneTwister(10), uncertain; distribution = :uniform) - 10.0) <=
          2sqrt(3.0)
    @test rand(
        MersenneTwister(10), uncertain;
        distribution = (_rng, nominal, sigma)->nominal+3sigma
    ) == 16.0
    @test_throws ArgumentError rand(
        MersenneTwister(1), uncertain; distribution = :cauchy)
    @test isfinite(rand(MersenneTwister(1), relative))
    @test rand(MersenneTwister(1), PB.Grid(3.0)) == 3.0
    @test rand(MersenneTwister(1), PB.Grid((:a, :b, :c))) in (:a, :b, :c)
    @test rand(PB.Grid(:fixed)) === :fixed

    @test_throws ArgumentError PB.Grid(1.0, -1.0)
    @test_throws ArgumentError PB.Grid(Inf, 1.0)
    @test_throws ArgumentError PB.Grid(1.0, PB.AbsoluteError(-1.0))
    @test_throws ArgumentError PB.Grid(:symbol, 1.0)
    @test_throws ArgumentError PB.Grid(:symbol, PB.AbsoluteError(1.0))
end

@testitem "Gridspace / product and zip" tags=[:unit, :commons] setup=[
    UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    product_space=PB.Gridspace{Tuple}(
        tuple,
        (PB.Grid((1, 2, 3)), PB.Grid((10, 20)))
    )
    @test (@inferred collect(product_space)) == [
        (1, 10), (2, 10), (3, 10),
        (1, 20), (2, 20), (3, 20)
    ]
    @test length(product_space) == 6
    @test size(product_space) == (6,)
    @test eltype(product_space) === Tuple{Int, Int}
    @test Base.IteratorEltype(typeof(product_space)) isa Base.HasEltype
    @test isconcretetype(eltype(product_space))
    @test PB.Grid(product_space) === product_space

    zip_space=PB.Gridspace{Tuple}(
        tuple,
        (PB.Grid((1, 2, 3)), PB.Grid((10, 20, 30)), PB.Grid(:fixed));
        combine = :zip
    )
    @test collect(zip_space) == [
        (1, 10, :fixed),
        (2, 20, :fixed),
        (3, 30, :fixed)
    ]
    @test length(zip_space) == 3
    @test_throws DimensionMismatch PB.Gridspace{Tuple}(
        tuple,
        (PB.Grid((1, 2)), PB.Grid((10, 20, 30)));
        combine = :zip
    )

    child=PB.Gridspace{Tuple}(
        tuple,
        (PB.Grid((1, 2)), PB.Grid((10, 20)));
        combine = :zip
    )
    parent=PB.Gridspace{Tuple}(
        tuple,
        (child, PB.Grid((:a, :b)))
    )
    @test collect(parent) == [
        ((1, 10), :a),
        ((2, 20), :a),
        ((1, 10), :b),
        ((2, 20), :b)
    ]
    @test collect(PB.Gridspace{Tuple}(
        tuple,
        (child, PB.Grid((:a, :b)));
        combine = :zip
    )) == [((1, 10), :a), ((2, 20), :b)]

    zero_argument_product = PB.Gridspace{Tuple}(tuple, ())
    @test (@inferred LineCableModels._combinations(zero_argument_product)) == ((),)
    @test (@inferred collect(zero_argument_product)) == [()]
    zero_argument_zip = PB.Gridspace{Tuple}(tuple, (); combine=:zip)
    @test (@inferred LineCableModels._combinations(zero_argument_zip)) == ((),)
    @test (@inferred collect(zero_argument_zip)) == [()]
    @test (@inferred collect(zip_space)) == [(1, 10, :fixed), (2, 20, :fixed), (3, 30, :fixed)]
    @test (@inferred collect(child)) == [(1, 10), (2, 20)]
    @test isempty(PB.Gridspace{Tuple}(tuple, (PB.Grid(()),)))
    empty_space=PB.Gridspace{Tuple}(tuple, (PB.Grid(()),))
    @test eltype(empty_space) === Any
    @test Base.IteratorEltype(typeof(empty_space)) isa Base.EltypeUnknown
    @test_throws ArgumentError PB.Gridspace{Tuple}(tuple, (1, 2))
    @test_throws ArgumentError PB.Gridspace{String}(
        identity, (PB.Grid((1, 2)),)
    )
    @test_throws ArgumentError PB.Gridspace{Tuple}(
        tuple, (PB.Grid(1),); combine = :outer)
    @test !applicable(getindex, product_space, 1)
    @test rand(MersenneTwister(1), product_space) in collect(product_space)
    @test rand(product_space) in collect(product_space)
    zipped_draws=Set(
        rand(MersenneTwister(seed), zip_space) for seed in 1:20
    )
    @test zipped_draws ⊆ Set(collect(zip_space))
    @test_throws ArgumentError rand(
        MersenneTwister(1),
        PB.Gridspace{Tuple}(tuple, (PB.Grid(()),))
    )
end

@testitem "Gridspace / composition conformance" tags=[:unit, :commons] setup=[
    UseEngineSupport] begin
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

    @test_throws ArgumentError PB.Gridspace{Tuple}(tuple, ((1, 2),))
    @test !applicable(getindex, product, 1)
    @test rand(Random.Xoshiro(0x1234), product) in collect(product)
end

@testitem "Gridspace / inference and allocation checks" tags=[:unit, :commons] setup=[
    UseEngineSupport] begin
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

    deterministic_sum(deterministic, 1)
    realization_sum(rng, uncertain_point, 1)
    @test @allocated(deterministic_sum(deterministic, 10_000)) == 0
    rng=Random.Xoshiro(0x1234)
    @test @allocated(realization_sum(rng, uncertain_point, 10_000)) == 0
end
