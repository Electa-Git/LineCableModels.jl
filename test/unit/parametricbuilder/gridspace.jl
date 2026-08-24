@testitem "ParametricBuilder / Grid / finite sources" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    deterministic=PB.Grid((1, 3, 2))
    @test PB.Grid(deterministic) === deterministic
    @test collect(PB.Grid(1:3)) == [1, 2, 3]
    @test deterministic[2] == 3
    @test extrema(deterministic) == (1, 3)
    @test size(deterministic) == (3,)
    @test fieldnames(typeof(deterministic)) == (:vals,)

    relative=PB.Grid((-10.0, 20.0), (5.0, 10.0))
    absolute=PB.Grid((10.0, 20.0), PB.AbsoluteError((0.5, 2.0)))
    @test length(relative) == 4
    @test relative[4] isa PB.UncertainValue
    @test fieldnames(typeof(relative)) == (:vals, :rel_err)
    @test fieldnames(typeof(absolute)) == (:vals, :abs_err)
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
    @test_throws ArgumentError rand(MersenneTwister(1), relative)
    @test rand(MersenneTwister(1), PB.Grid(3.0)) == 3.0

    @test_throws ArgumentError PB.Grid(1.0, -1.0)
    @test_throws ArgumentError PB.Grid(Inf, 1.0)
    @test_throws ArgumentError PB.Grid(1.0, PB.AbsoluteError(-1.0))
end

@testitem "ParametricBuilder / Gridspace / product and zip" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    import LineCableModels.ParametricBuilder as PB

    product_space=PB.Gridspace{Tuple}(
        tuple,
        (PB.Grid((1, 2, 3)), PB.Grid((10, 20)))
    )
    @test collect(product_space) == [
        (1, 10), (2, 10), (3, 10),
        (1, 20), (2, 20), (3, 20)
    ]
    @test length(product_space) == 6
    @test size(product_space) == (6,)
    @test PB.target_type(product_space) === Tuple
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

    @test collect(PB.Gridspace{Tuple}(tuple, ())) == [()]
    @test isempty(PB.Gridspace{Tuple}(tuple, (PB.Grid(()),)))
    @test_throws ArgumentError PB.Gridspace{Tuple}(tuple, (1, 2))
    @test_throws ArgumentError PB.Gridspace{Tuple}(
        tuple, (PB.Grid(1),); combine = :outer)
    @test !applicable(getindex, product_space, 1)
    @test_throws MethodError rand(product_space)
end

@testitem "ParametricBuilder / Gridspace / recursive point resolution" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    child=PB.Gridspace{Tuple}(tuple, (PB.Grid(10.0, 5.0), PB.Grid(:tag)))
    parent=PB.Gridspace{Tuple}(tuple, (child, PB.Grid((1, 2))))
    point=first(PB.points(parent))
    @test PB.has_uncertainty(parent)
    @test PB.has_uncertainty(point)
    @test PB.materialize(point)[1][1] isa Measurement
    draw=PB.realize(MersenneTwister(42), point, :normal)
    @test draw[1][1] isa Float64
    @test draw[1][2] == :tag
    @test draw[2] == 1

    struct MissingGridspaceDefinition end
    @test_throws MethodError PB.Gridspace(MissingGridspaceDefinition())
end

@testitem "ParametricBuilder / macros / strict and lifted construction" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    import LineCableModels.ParametricBuilder as PB

    PB.@gridspace PB.@relax struct MacroVault{T <: Real}
        value::T
        label::Symbol=:default
    end

    PB.@relax PB.@gridspace struct ReverseMacroVault{T <: Real}
        value::T
        label::Symbol=:default
    end

    @test collect(MacroVault(; value = PB.Grid((1.0, 2.0)), label = :ok)) == [
        MacroVault(1.0, :ok),
        MacroVault(2.0, :ok)
    ]
    @test only(MacroVault(; value = 1.0)) == MacroVault(1.0, :default)
    @test collect(ReverseMacroVault(; value = PB.Grid((3.0, 4.0)))) == [
        ReverseMacroVault(3.0, :default),
        ReverseMacroVault(4.0, :default)
    ]

    PB.@gridspace struct AtomicCollections{T}
        payload::T
        label::Symbol=:default
    end
    matrix=[1.0 2.0; 3.0 4.0]
    @test only(AtomicCollections(; payload = matrix)).payload === matrix
    varied=AtomicCollections(; payload = matrix, label = PB.Grid((:a, :b)))
    @test all(value -> value.payload === matrix, varied)
end
