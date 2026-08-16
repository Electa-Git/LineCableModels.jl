@testitem "ParametricBuilder / Gridspace / explicit axes composition coupling and macros" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Random
    using LineCableModels.ParametricBuilder

    @test collect(Grid((1, 2, 3))) == [1, 2, 3]
    @test collect(Grid(1:3)) == [1, 2, 3]

    relative=collect(Grid((10.0, 20.0), (5.0, 10.0)))
    @test length(relative) == 4
    @test nominal(relative[1]) == 10.0
    @test standard_uncertainty(relative[1]) == 0.5

    absolute=only(Grid(10.0, AbsoluteError(0.25)))
    @test nominal(absolute) == 10.0
    @test standard_uncertainty(absolute) == 0.25

    product_space=Gridspace{Tuple}(
        tuple,
        (Grid((1, 2, 3)), Grid((10, 20))),
        (:left, :right)
    )
    @test collect(product_space) == [
        (1, 10), (2, 10), (3, 10),
        (1, 20), (2, 20), (3, 20)
    ]

    direct_constant=Gridspace{Tuple}(
        tuple,
        (Grid((1, 2)), :fixed, [10, 20]),
        (:value, :tag, :atomic_collection)
    )
    @test collect(direct_constant) == [
        (1, :fixed, [10, 20]),
        (2, :fixed, [10, 20])
    ]

    zip_space=Gridspace{Tuple}(
        tuple,
        (Grid((1, 2, 3)), Grid((10, 20, 30)), Grid(:fixed)),
        (:left, :right, :tag);
        combine = :zip
    )
    @test collect(zip_space) == [
        (1, 10, :fixed),
        (2, 20, :fixed),
        (3, 30, :fixed)
    ]
    bad_zip=Gridspace{Tuple}(
        tuple,
        (Grid((1, 2)), Grid((10, 20, 30))),
        (:left, :right);
        combine = :zip
    )
    @test_throws DimensionMismatch length(bad_zip)

    child=Gridspace{Tuple}(
        tuple,
        (Grid((1, 2)), Grid((10, 20))),
        (:x, :y);
        combine = :zip
    )
    parent=Gridspace{Tuple}(
        tuple,
        (child, Grid((:a, :b))),
        (:child, :tag)
    )
    @test collect(parent) == [
        ((1, 10), :a),
        ((2, 20), :a),
        ((1, 10), :b),
        ((2, 20), :b)
    ]
    zipped_parent=Gridspace{Tuple}(
        tuple,
        (child, Grid((:a, :b))),
        (:child, :tag);
        combine = :zip
    )
    @test collect(zipped_parent) == [
        ((1, 10), :a),
        ((2, 20), :b)
    ]

    shared=Grid((1, 2, 3))
    coupled=Gridspace{Tuple}(
        tuple,
        (shared, shared),
        (:left, :right)
    )
    @test collect(coupled) == [(1, 1), (2, 2), (3, 3)]
    @test length(coupled) == 3

    independent=Gridspace{Tuple}(
        tuple,
        (Grid((1, 2)), Grid((1, 2))),
        (:left, :right)
    )
    @test length(independent) == 4

    uncertain_space=Gridspace{Tuple}(
        tuple,
        (Grid(10.0, 5.0), Grid(:constant)),
        (:value, :tag)
    )
    configuration=only(configurations(uncertain_space))
    @test has_uncertainty(configuration)
    @test configuration_manifest(configuration).value.style == :relative
    direct=materialize(configuration)
    @test direct[1] isa Measurement
    draw=rand(MersenneTwister(42), configuration)
    @test draw[1] isa Float64
    @test draw[2] == :constant

    shared_uncertainty=Grid(10.0, 5.0)
    shared_space=Gridspace{Tuple}(
        tuple,
        (shared_uncertainty, shared_uncertainty),
        (:left, :right)
    )
    shared_configuration=only(configurations(shared_space))
    shared_draw=rand(MersenneTwister(91), shared_configuration)
    @test shared_draw[1] === shared_draw[2]
    shared_direct=materialize(shared_configuration)
    @test shared_direct[1] === shared_direct[2]
    @test Measurements.cov(shared_direct[1], shared_direct[2]) ==
          Measurements.uncertainty(shared_direct[1])^2

    @gridspace @relax struct MacroVault{T <: Real}
        value::T
        label::Symbol=:default
    end

    @relax @gridspace struct ReverseMacroVault{T <: Real}
        value::T
        label::Symbol=:default
    end

    @test collect(MacroVault(; value = Grid((1.0, 2.0)), label = :ok)) == [
        MacroVault(1.0, :ok),
        MacroVault(2.0, :ok)
    ]
    @test collect(ReverseMacroVault(; value = Grid((3.0, 4.0)))) == [
        ReverseMacroVault(3.0, :default),
        ReverseMacroVault(4.0, :default)
    ]
    @test MacroVault(Float32(1), :stable) isa MacroVault{Float32}

    abstract type AbstractMacroVault{T} end
    @gridspace @relax struct ConcreteMacroVault{T <: Real} <: AbstractMacroVault{T}
        value::T
    end
    concrete=ConcreteMacroVault(1.0)
    @test convert(AbstractMacroVault{Float32}, concrete) isa
          ConcreteMacroVault{Float32}

    @gridspace struct AtomicCollections{T}
        payload::T
    end
    matrix=[1.0 2.0; 3.0 4.0]
    atomic=only(AtomicCollections(; payload = matrix))
    @test atomic.payload === matrix
end

@testitem "ParametricBuilder / Grid / indexing, extrema, and local sampling" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    deterministic=PB.Grid((1, 3, 2))
    @test PB.Grid(deterministic) === deterministic
    @test_throws ArgumentError PB.Grid(deterministic; key = :replacement)
    @test deterministic[2] == 3
    @test extrema(deterministic) == (1, 3)
    @test size(deterministic) == (3,)
    @test eltype(typeof(deterministic)) === Int

    relative=PB.Grid((-10.0, 20.0), (5.0, 10.0); key = :shared)
    absolute=PB.Grid((10.0, 20.0), PB.AbsoluteError((0.5, 2.0)))
    @test size(relative) == (4,)
    @test relative[4] isa PB.UncertainValue
    @test PB.uncertainty_style(relative[1]) isa PB.RelativeUncertainty
    @test PB.uncertainty_style(absolute[1]) isa PB.AbsoluteUncertainty
    @test extrema(relative) == (-11.0, 22.0)
    @test extrema(absolute) == (8.0, 22.0)
    @test_throws BoundsError relative[0]
    @test_throws BoundsError absolute[5]

    first_key=PB.AutomaticGridKey(Ref(nothing))
    alias_key=PB.AutomaticGridKey(first_key.token)
    other_key=PB.AutomaticGridKey(Ref(nothing))
    @test first_key == alias_key
    @test isequal(first_key, alias_key)
    @test first_key != other_key
    @test hash(first_key) == hash(alias_key)

    zero_sigma=PB.UncertainValue(4, 0, PB.AbsoluteUncertainty())
    @test rand(MersenneTwister(1), zero_sigma) === 4.0
    @test rand(zero_sigma) === 4.0
    uncertain=PB.UncertainValue(10.0, 2.0, PB.AbsoluteUncertainty())
    normal=rand(MersenneTwister(10), uncertain; distribution = :normal)
    uniform=rand(MersenneTwister(10), uncertain; distribution = :uniform)
    custom=rand(
        MersenneTwister(10),
        uncertain;
        distribution = (_rng, nominal, sigma)->nominal+3sigma
    )
    @test isfinite(normal)
    @test abs(uniform - 10.0) ≤ 2sqrt(3.0)
    @test custom == 16.0
    @test_throws ArgumentError rand(MersenneTwister(1), uncertain; distribution = :cauchy)
    @test_throws ArgumentError rand(MersenneTwister(1), uncertain; distribution = missing)
    @test_throws ArgumentError rand(MersenneTwister(1), relative)
    @test rand(PB.Grid(3.0)) == 3.0
end

@testitem "ParametricBuilder / Gridspace / interface and failure boundaries" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    space=PB.Gridspace{Pair}((PB.Grid((1, 2)), :fixed))
    @test PB.target_type(typeof(space)) === Pair
    @test PB.target_type(space) === Pair
    @test PB.Grid(space) === space
    @test_throws ArgumentError PB.Grid(space; key = :coupled)
    @test size(space) == (2,)
    @test space[2] == (2 => :fixed)
    @test eltype(typeof(space)) === Pair

    configuration=first(PB.configurations(space))
    @test PB.target_type(typeof(configuration)) === Pair
    @test PB.target_type(configuration) === Pair
    @test !PB.has_uncertainty(configuration)
    @test !PB.has_uncertainty(space)
    @test !PB.has_uncertainty(:constant)

    axis=PB.ConstantAxis(:only)
    @test length(axis) == 1
    @test axis[1] === :only
    @test collect(axis) == [:only]
    @test_throws BoundsError axis[2]

    singleton=PB.Gridspace{Tuple}(tuple, (PB.Grid(3.0),), (:value,))
    @test rand(MersenneTwister(2), singleton) == (3.0,)
    @test rand(singleton) == (3.0,)
    @test_throws ArgumentError rand(MersenneTwister(2), space)
    empty_space=PB.Gridspace{Tuple}(tuple, (PB.Grid(()),), (:value,))
    @test isempty(collect(empty_space))
    @test_throws ArgumentError rand(MersenneTwister(2), empty_space)

    struct MissingGridspaceSpec<:PB.AbstractSpec{Int} end
    missing_spec=MissingGridspaceSpec()
    @test_throws MethodError PB.gridspace(missing_spec)

    @test PB.recast(Float32, [1.0, 2.0]) == Float32[1, 2]
    @test PB.recast(Float32, (1.0, :fixed)) == (Float32(1), :fixed)
    @test PB.recast(Float32, (value = 1.0, label = :fixed)) ==
          (value = Float32(1), label = :fixed)
    @test PB.recast(Float32, :fixed) === :fixed
    @test PB._promoted_numeric_type((Float32[1], (2.0,), :fixed)) === Float64
    @test_throws ArgumentError PB._promoted_numeric_type((:left, :right))

    escaped=Expr(:escape, Expr(:call, :identity, :x))
    @test PB._strip_escapes(escaped) == Expr(:call, :identity, :x)
    @test_throws ArgumentError PB._get_struct_node(:(begin
    end))
    @test_throws ArgumentError PB._get_struct_node(:(begin
        struct A
            x
        end
        struct B
            y
        end
    end))
    @test_throws ArgumentError PB._parse_fields(:(begin
        println("no field")
    end))
    malformed=Expr(:struct, false, Expr(:call, :Bad), Expr(:block, :x))
    @test_throws ArgumentError PB._extract_struct_name(malformed)

    @relax struct RelaxedCollections{T <: Real}
        scalar::T
        vector::Vector{T}
        tuple::Tuple{T, T}
        metadata::NamedTuple{(:value,), Tuple{T}}
        label::Symbol
    end
    relaxed=RelaxedCollections(
        Float32(1),
        [2.0, 3.0],
        (Float32(4), 5.0),
        (value = Float32(6),),
        :ok
    )
    @test relaxed isa RelaxedCollections{Float64}
    @test eltype(relaxed) === Float64
    @test eltype(typeof(relaxed)) === Float64
    relaxed32=convert(RelaxedCollections{Float32}, relaxed)
    @test relaxed32 isa RelaxedCollections{Float32}
    @test PB.recast(Float32, relaxed) isa RelaxedCollections{Float32}
end
