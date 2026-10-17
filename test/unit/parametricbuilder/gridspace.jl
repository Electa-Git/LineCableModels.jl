@testitem "ParametricBuilder / joint input identity and bounded covariance" tags=[:unit, :extension, :measurements] begin
    using Measurements, Random, Statistics
    source = Grid(1.,10.)
    joint = Gridspace{NamedTuple{(:inner,:outer)}}(s->(inner=0.01s,outer=0.012s),(source,))
    p = only(joint)
    @test Measurements.cov(p.inner,p.outer) ≈ 0.01*0.01*0.012
    @test uncertainty(p.outer-p.inner) ≈ 0.0002
    a,b = only(Gridspace{Tuple}(tuple,(source,source)))
    @test iszero(Measurements.cov(a,b))
    calls = Ref(0)
    sampler = (rng,mu,sigma)->(calls[]+=1;mu+sigma*(2rand(rng)-1))
    rand(Xoshiro(2),joint;distribution=sampler)
    @test calls[] == 1
    calls[] = 0
    rand(Xoshiro(2),Gridspace{Tuple}(tuple,(source,source));distribution=sampler)
    @test calls[] == 2
    for combine in (:product,:zip)
        nested = Gridspace{Tuple}(tuple,(joint,Grid((:a,:b)));combine)
        @test length(nested) == 2
        @test last.(collect(nested)) == [:a,:b]
        @test all(v->uncertainty(first(v).outer/first(v).inner)<1e-14,nested)
    end
    rng = Xoshiro(42)
    draws = [rand(rng,joint;distribution=:uniform) for _ in 1:10_000]
    scales = getproperty.(draws,:inner)./0.01
    @test all(s->1-sqrt(3)*0.1 <= s <= 1+sqrt(3)*0.1,scales)
    @test all(p->p.outer>p.inner>0,draws)
    @test mean(scales) ≈ 1 atol=0.004
    @test std(scales) ≈ 0.1 atol=0.003
    @test mean(scales.^2) ≈ 1.01 atol=0.008
    @test Statistics.cov(getproperty.(draws,:inner),getproperty.(draws,:outer)) ≈
        0.01*0.01*0.012 rtol=0.06
end

@testitem "ParametricBuilder / Material / invariant class and scalar precision" tags=[:unit, :parametric] setup=[
    UseEngineSupport] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    # A literal class infers the concrete material type through constant propagation.
    material=@inferred (() -> PB.Material(kind = :conductor, rho = Float32(1.7e-8)))()
    @test material isa LineCableModels.Materials.Material{Float32, :conductor}
    @test material.kind === :conductor

    space=PB.Material(
        kind = PB.Grid((:conductor, :semicon)),
        rho = Float32(1.7e-8)
    )
    # The class is a type parameter. A grid over classes yields materials of two types.
    @test space isa PB.Gridspace{LineCableModels.Materials.Material}
    @test eltype(space) === LineCableModels.Materials.Material
    @test Base.IteratorEltype(typeof(space)) isa Base.EltypeUnknown
    @test first(space) isa LineCableModels.Materials.Material{Float32}
    @test getproperty.(collect(space), :kind) == [:conductor, :semicon]
    @test all(value -> eltype(value) === Float32, space)
    sampled=rand(MersenneTwister(8), space)
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
    # A fixed class keeps each material concrete.
    @test eltype(resistivity_space) === LineCableModels.Materials.Material
    @test @inferred(first(resistivity_space)) isa
          LineCableModels.Materials.Material{Float64, :conductor}
    @test @inferred(rand(MersenneTwister(8), resistivity_space)).kind === :conductor
    @test eltype(@inferred(collect(resistivity_space))) ===
          LineCableModels.Materials.Material{Float64, :conductor}

    mixed_space=PB.Material(
        kind = :conductor,
        rho = PB.Grid((1.0, 100.0)),
        eps_r = PB.Grid(1.0, 5.0),
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.0
    )
    @test PB.has_uncertainty(mixed_space)
    @test eltype(mixed_space) === LineCableModels.Materials.Material
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

@testitem "ParametricBuilder / Gridspace / recursive point resolution" tags=[:unit, :measurements] setup=[
    UseEngineSupport, TestNumerics] begin
    using Random
    import LineCableModels.ParametricBuilder as PB

    child=PB.Gridspace{Tuple}(tuple, (PB.Grid(10.0, 5.0), PB.Grid(:tag)))
    parent=PB.Gridspace{Tuple}(tuple, (child, PB.Grid((1, 2))))
    point=first(PB.points(parent))
    @test point isa LineCableModels.Gridpoint{Tuple}
    @test PB.has_uncertainty(parent)
    @test eltype(parent) === Tuple
    @test Base.IteratorEltype(typeof(parent)) isa Base.EltypeUnknown
    @test PB.has_uncertainty(point)
    @test PB.materialize(point)[1][1] isa Measurement
    arguments=LineCableModels.realize_arguments(MersenneTwister(42), point, :normal)
    draw=LineCableModels.realize(point, arguments)
    @test draw[1][1] isa Float64
    @test draw[1][2] == :tag
    @test draw[2] == 1

    struct MissingGridspaceDefinition end
    @test_throws MethodError PB.Gridspace(MissingGridspaceDefinition())
end

@testitem "ParametricBuilder / parameterize / scalar, Grid and atomic inputs of a user lift" tags=[:unit, :parametric] setup=[
    UseEngineSupport] begin
    struct LabeledValue{T <: Real}
        value::T
        label::Symbol
    end

    function LabeledValue(
            value::Union{Real, AbstractGrid, Gridspace},
            label::Union{Symbol, AbstractGrid, Gridspace} = :default;
            combine::Symbol = :product
    )
        return LineCableModels.parameterize(
            LabeledValue, LabeledValue, (value, label); combine)
    end

    @test collect(LabeledValue(Grid((1.0, 2.0)), :ok)) == [
        LabeledValue(1.0, :ok),
        LabeledValue(2.0, :ok)
    ]
    scalar=LabeledValue(1.0)
    @test scalar isa LabeledValue
    @test scalar == LabeledValue(1.0, :default)

    struct AtomicCollections{T <: AbstractMatrix}
        payload::T
        label::Symbol
    end

    function AtomicCollections(
            payload::Union{AbstractMatrix, AbstractGrid, Gridspace},
            label::Union{Symbol, AbstractGrid, Gridspace} = :default;
            combine::Symbol = :product
    )
        return LineCableModels.parameterize(
            AtomicCollections, AtomicCollections, (payload, label); combine)
    end

    matrix=[1.0 2.0; 3.0 4.0]
    @test AtomicCollections(matrix).payload === matrix
    varied=AtomicCollections(matrix, Grid((:a, :b)))
    @test all(value -> value.payload === matrix, varied)
end
