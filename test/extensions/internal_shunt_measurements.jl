@testitem "Measurements / internal shunt / correlated implicit sensitivities" tags=[:extension] begin
    using Measurements, LinearAlgebra, Random, Statistics
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels),"test","support","internal_shunt.jl"))
    epsilon = measurement(2.5,0.0025)
    design = internal_shunt_test_design(;epsilon,count=4)
    domain,_ = internal_shunt_test_domain(design)
    values = E._shunt_values(domain,Formulation().methods)
    level = (wire=16,order=8,quadrature=64,modes=128)
    propagated = E.internal_shunt_response(values,domain;level)
    @test all(isfinite,Measurements.value.(propagated.C))
    @test all(isfinite,Measurements.uncertainty.(propagated.C))
    @test maximum(Measurements.uncertainty.(propagated.C)) > 0
    @test abs(Measurements.cor(vec(propagated.C))[1,2]) ≈ 1 rtol=1e-8
    deterministic(eps) = begin
        d,_ = internal_shunt_test_domain(internal_shunt_test_design(epsilon=eps,count=4))
        E.internal_shunt_response(E._shunt_values(d,Formulation().methods),d;level).C
    end
    central = (deterministic(2.5001)-deterministic(2.4999))/0.0002
    @test Measurements.derivative.(propagated.C,Ref(epsilon)) ≈ central rtol=2e-3
    @test Measurements.uncertainty.(propagated.C) ≈ abs.(central)*0.0025 rtol=2e-3
    # Equal nominal values must not cause sharing across independent inputs.
    independent,_ = internal_shunt_test_domain(internal_shunt_test_design(
        epsilon=measurement(2.5,0.0025),count=4))
    @test !E._shunt_domain_equal(domain,independent)
    rng = MersenneTwister(314)
    samples = [deterministic(2.5+0.0025randn(rng))[1,1] for _ in 1:24]
    # Three standard errors of the sampled standard deviation; no requirement
    # that a finite MC sample exactly equals the linearized variance.
    sigma = Measurements.uncertainty(propagated.C[1,1])
    @test abs(std(samples)-sigma) < 3sigma/sqrt(2(length(samples)-1))
    @test abs(mean(samples)-Measurements.value(propagated.C[1,1])) < 4sigma/sqrt(length(samples))
    radius = measurement(0.2e-3,0.1e-6)
    geometry_domain,_ = internal_shunt_test_domain(internal_shunt_test_design(;radius,count=4))
    geometry_values = E._shunt_values(geometry_domain,Formulation().methods)
    geometry_result = E.internal_shunt_response(geometry_values,geometry_domain;level)
    function varying_radius(r)
        d,_ = internal_shunt_test_domain(internal_shunt_test_design(radius=r,count=4))
        E.internal_shunt_response(E._shunt_values(d,Formulation().methods),d;level).C
    end
    geometric_central = (varying_radius(0.2e-3+1e-8)-varying_radius(0.2e-3-1e-8))/2e-8
    @test Measurements.derivative.(geometry_result.C,Ref(radius)) ≈ geometric_central rtol=5e-3
    @test maximum(Measurements.uncertainty.(geometry_result.C)) > 0
    centre = measurement(0.0,1e-6)
    other_centre = measurement(0.0,1e-6)
    @test E._shunt_same(centre,centre,0.01)
    @test !E._shunt_same(centre,other_centre,0.01)
    @test !E._shunt_same(centre,0.0,0.01)
    wire_angle = measurement(0.31,1e-4)
    angle_domain,_ = internal_shunt_test_domain(internal_shunt_test_design(;wire_angle,count=4))
    angle_values = E._shunt_values(angle_domain,Formulation().methods)
    angle_result = E.internal_shunt_response(angle_values,angle_domain;level)
    function varying_angle(angle)
        d,_ = internal_shunt_test_domain(internal_shunt_test_design(wire_angle=angle,count=4))
        E.internal_shunt_response(E._shunt_values(d,Formulation().methods),d;level).C
    end
    angular_central = (varying_angle(0.31001)-varying_angle(0.30999))/2e-5
    @test Measurements.derivative.(angle_result.C,Ref(wire_angle)) ≈ angular_central rtol=5e-3
    @test maximum(Measurements.uncertainty.(angle_result.C)) > 0
end

@testitem "Measurements / internal shunt / public parametric studies" tags=[:extension,:integration] begin
    using Measurements, Statistics
    include(joinpath(pkgdir(LineCableModels),"test","support","internal_shunt.jl"))
    space = Gridspace{CableConstantsProblem}(epsilon->CableConstantsProblem(
        internal_shunt_test_design(;epsilon,tapes=false,count=4)),
        (Grid(2.5,AbsoluteError(0.0025)),))
    problem = ParametricProblem(space)
    inner = CableConstantsFormulation(shunt_model=formula(:boundary;
        options=(resolution=(wire=16,order=8,quadrature=64,modes=128),)))
    linear = compute(problem,LinearError(inner))
    sampled = compute(problem,MonteCarlo(inner;trials=12,seed=314,
        distribution=:normal,return_samples=true,retain_details=true))
    @test sampled.trial_counts == [12]
    @test all(isempty,sampled.details.data.failures)
    @test only(sampled.values).details.data.shunt_model.effective === :boundary
    @test only(linear.values).details.data.shunt_model.effective === :boundary
    IE=LineCableModels.ImportExport
    restored=IE.deserialize_value(IE.serialize_value(linear))
    original=only(linear.values)
    saved=only(restored.values)
    @test Measurements.value.(saved.C) == Measurements.value.(original.C)
    @test Measurements.uncertainty.(saved.C) ≈ Measurements.uncertainty.(original.C)
    @test details(saved).data.shunt_model.effective === :boundary
    # One shared uncertain input must remain one source across distinct points.
    duplicated=LineCableModels.UQ.LinearErrorResult(LinearError(inner),[original,original], ComputationDetails((;)))
    copied=IE.deserialize_value(IE.serialize_value(duplicated))
    @test Measurements.uncertainty(copied.values[1].C[1]-copied.values[2].C[1]) == 0
    strict = CableConstantsFormulation(shunt_model=formula(:boundary;
        options=(resolution=(wire=100_000,),)))
    @test_throws BoundarySolveError compute(problem,MonteCarlo(strict;
        trials=2,seed=314,on_error=:retry,retain_details=true))
    # Blueprint sharing must not erase independent Measurements sources.
    independent = CableConstantsProblem(internal_shunt_test_design(
        epsilon=measurement(2.5,0.0025),tapes=false,count=4))
    original_problem = first(space)
    blueprints = only(LineCableModels.Engine.flatten(LineCableModelsCoaxial(),
        [original_problem.design, independent.design], eltype(original_problem), [inner]))
    @test blueprints[1].shunt[1].C !== blueprints[2].shunt[1].C
    @test Measurements.uncertainty(blueprints[1].shunt[1].C[1,1] -
        blueprints[2].shunt[1].C[1,1]) > 0
    cap = only(only(linear.values).C)
    @test Measurements.uncertainty(cap) > 0
    @test abs(only(only(sampled.values).C)-Measurements.value(cap)) <
        4Measurements.uncertainty(cap)/sqrt(12)
end
