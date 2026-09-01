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

@testitem "ParametricBuilder / Gridspace / product and zip" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
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

    @test collect(PB.Gridspace{Tuple}(tuple, ())) == [()]
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

@testitem "ParametricBuilder / Material / invariant class and scalar precision" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    material=@inferred PB.Material(kind = :conductor, rho = Float32(1.7e-8))
    @test material isa LineCableModels.Materials.Material{Float32}
    @test material.kind === :conductor

    space=PB.Material(
        kind = PB.Grid((:conductor, :semicon)),
        rho = Float32(1.7e-8)
    )
    @test space isa PB.Gridspace{LineCableModels.Materials.Material}
    @test isconcretetype(typeof(space.build))
    @test eltype(space) === LineCableModels.Materials.Material{Float32}
    @test Base.IteratorEltype(typeof(space)) isa Base.HasEltype
    @test @inferred(first(space)) isa LineCableModels.Materials.Material{Float32}
    @test getproperty.(collect(space), :kind) == [:conductor, :semicon]
    @test all(value -> eltype(value) === Float32, space)
    sampled=@inferred rand(MersenneTwister(8), space)
    @test sampled.kind in (:conductor, :semicon)
    @test eltype(sampled) === Float32

    resistivity_space=PB.Material(
        kind = :conductor,
        rho = PB.Grid((1.0, 100.0)),
        eps_r = 1.0,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.0
    )
    @test eltype(resistivity_space) ===
          LineCableModels.Materials.Material{Float64}
    @test @inferred(first(resistivity_space)) isa
          LineCableModels.Materials.Material{Float64}

    mixed_space=PB.Material(
        kind = :conductor,
        rho = PB.Grid((1.0, 100.0)),
        eps_r = PB.Grid(1.0, 5.0),
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.0
    )
    @test PB.has_uncertainty(mixed_space)
    @test eltype(mixed_space) === Any
    @test Base.IteratorEltype(typeof(mixed_space)) isa Base.EltypeUnknown
    @test rand(MersenneTwister(9), mixed_space) isa
          LineCableModels.Materials.Material{Float64}

    malformed_error=try
        PB.Material(
            1.0;
            kind = :insulator,
            rho = 1.97e14,
            eps_r = PB.Grid((2.3,)),
            mu_r = 1.0,
            T0 = 20.0,
            alpha = 0.0
        )
    catch error
        error
    end
    @test malformed_error isa ArgumentError
    @test occursin(
        "Grid(values, relative_error)", sprint(showerror, malformed_error)
    )
end

@testitem "ParametricBuilder / Gridspace / recursive point resolution" tags=[:unit] setup=[
    EngineTestSupport, UseEngineSupport, TestNumerics] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    child=PB.Gridspace{Tuple}(tuple, (PB.Grid(10.0, 5.0), PB.Grid(:tag)))
    parent=PB.Gridspace{Tuple}(tuple, (child, PB.Grid((1, 2))))
    point=first(PB.points(parent))
    @test point isa LineCableModels.Gridpoint{Tuple}
    @test PB.has_uncertainty(parent)
    @test eltype(parent) === Any
    @test Base.IteratorEltype(typeof(parent)) isa Base.EltypeUnknown
    @test PB.has_uncertainty(point)
    @test PB.materialize(point)[1][1] isa Measurement
    arguments=PB.realize_arguments(MersenneTwister(42), point, :normal)
    draw=PB.realize(point, arguments)
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

    PB.@relax struct RelaxedMixedFields{T <: Real}
        value::T
        count::Int
        enabled::Bool
        label::Symbol
    end

    @test collect(MacroVault(; value = PB.Grid((1.0, 2.0)), label = :ok)) == [
        MacroVault(1.0, :ok),
        MacroVault(2.0, :ok)
    ]
    scalar=MacroVault(; value = 1.0)
    @test scalar isa MacroVault
    @test scalar == MacroVault(1.0, :default)
    @test collect(ReverseMacroVault(; value = PB.Grid((3.0, 4.0)))) == [
        ReverseMacroVault(3.0, :default),
        ReverseMacroVault(4.0, :default)
    ]
    mixed=RelaxedMixedFields(Float32(1), 3, true, :fixed)
    @test eltype(mixed) === Float32
    @test mixed.count === 3
    @test mixed.enabled === true
    @test mixed.label === :fixed
    converted=convert(RelaxedMixedFields{Float64}, mixed)
    @test converted.value === 1.0
    @test converted.count === 3
    @test converted.enabled === true
    @test converted.label === :fixed

    PB.@gridspace struct AtomicCollections{T}
        payload::T
        label::Symbol=:default
    end
    matrix=[1.0 2.0; 3.0 4.0]
    @test AtomicCollections(; payload = matrix).payload === matrix
    varied=AtomicCollections(; payload = matrix, label = PB.Grid((:a, :b)))
    @test all(value -> value.payload === matrix, varied)
end
