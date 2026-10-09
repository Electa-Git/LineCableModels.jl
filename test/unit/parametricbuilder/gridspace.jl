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

    material=@inferred PB.Material(kind = :conductor, rho = Float32(1.7e-8))
    @test material isa LineCableModels.Materials.Material{Float32}
    @test material.kind === :conductor

    space=PB.Material(
        kind = PB.Grid((:conductor, :semicon)),
        rho = Float32(1.7e-8)
    )
    @test space isa PB.Gridspace{LineCableModels.Materials.Material}
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

@testitem "ParametricBuilder / Gridspace / recursive point resolution" tags=[:unit, :measurements] setup=[
    UseEngineSupport, TestNumerics] begin
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
    arguments=LineCableModels.realize_arguments(MersenneTwister(42), point, :normal)
    draw=LineCableModels.realize(point, arguments)
    @test draw[1][1] isa Float64
    @test draw[1][2] == :tag
    @test draw[2] == 1

    struct MissingGridspaceDefinition end
    @test_throws MethodError PB.Gridspace(MissingGridspaceDefinition())
end

@testitem "ParametricBuilder / macros / strict and lifted construction" tags=[:unit, :parametric] setup=[
    UseEngineSupport] begin
    import LineCableModels.ParametricBuilder as PB

    PB.@gridspace struct MacroVault{T <: Real}
        value::T
        label::Symbol=:default
    end

    @test collect(MacroVault(; value = PB.Grid((1.0, 2.0)), label = :ok)) == [
        MacroVault(1.0, :ok),
        MacroVault(2.0, :ok)
    ]
    scalar=MacroVault(; value = 1.0)
    @test scalar isa MacroVault
    @test scalar == MacroVault(1.0, :default)
    PB.@gridspace struct AtomicCollections{T}
        payload::T
        label::Symbol=:default
    end
    matrix=[1.0 2.0; 3.0 4.0]
    @test AtomicCollections(; payload = matrix).payload === matrix
    varied=AtomicCollections(; payload = matrix, label = PB.Grid((:a, :b)))
    @test all(value -> value.payload === matrix, varied)
end
