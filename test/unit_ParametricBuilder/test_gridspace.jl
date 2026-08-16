@testitem "Gridspace: explicit axes, local composition, coupling, and macros" setup = [defaults] begin
    using Random
    using LineCableModels.ParametricBuilder

    @test collect(Grid((1, 2, 3))) == [1, 2, 3]
    @test collect(Grid(1:3)) == [1, 2, 3]

    relative = collect(Grid((10.0, 20.0), (5.0, 10.0)))
    @test length(relative) == 4
    @test nominal(relative[1]) == 10.0
    @test standard_uncertainty(relative[1]) == 0.5

    absolute = only(Grid(10.0, AbsoluteError(0.25)))
    @test nominal(absolute) == 10.0
    @test standard_uncertainty(absolute) == 0.25

    product_space = Gridspace{Tuple}(
        tuple,
        (Grid((1, 2, 3)), Grid((10, 20))),
        (:left, :right),
    )
    @test collect(product_space) == [
        (1, 10), (2, 10), (3, 10),
        (1, 20), (2, 20), (3, 20),
    ]

    direct_constant = Gridspace{Tuple}(
        tuple,
        (Grid((1, 2)), :fixed, [10, 20]),
        (:value, :tag, :atomic_collection),
    )
    @test collect(direct_constant) == [
        (1, :fixed, [10, 20]),
        (2, :fixed, [10, 20]),
    ]

    zip_space = Gridspace{Tuple}(
        tuple,
        (Grid((1, 2, 3)), Grid((10, 20, 30)), Grid(:fixed)),
        (:left, :right, :tag);
        combine=:zip,
    )
    @test collect(zip_space) == [
        (1, 10, :fixed),
        (2, 20, :fixed),
        (3, 30, :fixed),
    ]
    bad_zip = Gridspace{Tuple}(
        tuple,
        (Grid((1, 2)), Grid((10, 20, 30))),
        (:left, :right);
        combine=:zip,
    )
    @test_throws DimensionMismatch length(bad_zip)

    child = Gridspace{Tuple}(
        tuple,
        (Grid((1, 2)), Grid((10, 20))),
        (:x, :y);
        combine=:zip,
    )
    parent = Gridspace{Tuple}(
        tuple,
        (child, Grid((:a, :b))),
        (:child, :tag),
    )
    @test collect(parent) == [
        ((1, 10), :a),
        ((2, 20), :a),
        ((1, 10), :b),
        ((2, 20), :b),
    ]
    zipped_parent = Gridspace{Tuple}(
        tuple,
        (child, Grid((:a, :b))),
        (:child, :tag);
        combine=:zip,
    )
    @test collect(zipped_parent) == [
        ((1, 10), :a),
        ((2, 20), :b),
    ]

    shared = Grid((1, 2, 3))
    coupled = Gridspace{Tuple}(
        tuple,
        (shared, shared),
        (:left, :right),
    )
    @test collect(coupled) == [(1, 1), (2, 2), (3, 3)]
    @test length(coupled) == 3

    independent = Gridspace{Tuple}(
        tuple,
        (Grid((1, 2)), Grid((1, 2))),
        (:left, :right),
    )
    @test length(independent) == 4

    uncertain_space = Gridspace{Tuple}(
        tuple,
        (Grid(10.0, 5.0), Grid(:constant)),
        (:value, :tag),
    )
    configuration = only(configurations(uncertain_space))
    @test has_uncertainty(configuration)
    @test configuration_manifest(configuration).value.style == :relative
    direct = materialize(configuration)
    @test direct[1] isa Measurement
    draw = rand(MersenneTwister(42), configuration)
    @test draw[1] isa Float64
    @test draw[2] == :constant

    shared_uncertainty = Grid(10.0, 5.0)
    shared_space = Gridspace{Tuple}(
        tuple,
        (shared_uncertainty, shared_uncertainty),
        (:left, :right),
    )
    shared_configuration = only(configurations(shared_space))
    shared_draw = rand(MersenneTwister(91), shared_configuration)
    @test shared_draw[1] === shared_draw[2]
    shared_direct = materialize(shared_configuration)
    @test shared_direct[1] === shared_direct[2]
    @test Measurements.cov(shared_direct[1], shared_direct[2]) ==
          Measurements.uncertainty(shared_direct[1])^2

    @gridspace @relax struct MacroVault{T<:Real}
        value::T
        label::Symbol=:default
    end

    @relax @gridspace struct ReverseMacroVault{T<:Real}
        value::T
        label::Symbol=:default
    end

    @test collect(MacroVault(; value=Grid((1.0, 2.0)), label=:ok)) == [
        MacroVault(1.0, :ok),
        MacroVault(2.0, :ok),
    ]
    @test collect(ReverseMacroVault(; value=Grid((3.0, 4.0)))) == [
        ReverseMacroVault(3.0, :default),
        ReverseMacroVault(4.0, :default),
    ]
    @test MacroVault(Float32(1), :stable) isa MacroVault{Float32}

    abstract type AbstractMacroVault{T} end
    @gridspace @relax struct ConcreteMacroVault{T<:Real} <: AbstractMacroVault{T}
        value::T
    end
    concrete = ConcreteMacroVault(1.0)
    @test convert(AbstractMacroVault{Float32}, concrete) isa
          ConcreteMacroVault{Float32}

    @gridspace struct AtomicCollections{T}
        payload::T
    end
    matrix = [1.0 2.0; 3.0 4.0]
    atomic = only(AtomicCollections(; payload=matrix))
    @test atomic.payload === matrix
end
