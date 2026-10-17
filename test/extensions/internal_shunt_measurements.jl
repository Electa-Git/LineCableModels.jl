@testitem "Measurements / internal shunt / wire contacts preserve fixed topology" tags=[:extension, :measurements] begin
    using Measurements
    SM = LineCableModels.Engine.ShuntModel
    radius = measurement(0.2, 1e-5)
    one_layer = [(ri = 1.0, ro = 3.0, epsilon = 1.0)]
    uncertain = [(x = measurement(1.8, 0.0), y = zero(radius), r = radius, coat = radius,
            epsilon = zero(radius)+1, terminal = 2, layer = 1),
        (x = measurement(2.2, 0.0), y = zero(radius), r = radius, coat = radius,
            epsilon = zero(radius)+1, terminal = 2, layer = 1)]
    @test_throws SM.BoundarySolveError SM.validate(
        uncertain, SM.DielectricGroup, (cable_id = "uncertain-contact",), one_layer, 3.0)
    scale = measurement(1.0, 1e-3)
    coordinated = [(x = 1.8scale, y = zero(scale), r = 0.2scale, coat = 0.2scale,
            epsilon = zero(scale)+1, terminal = 2, layer = 1),
        (x = 2.2scale, y = zero(scale), r = 0.2scale, coat = 0.2scale,
            epsilon = zero(scale)+1, terminal = 2, layer = 1)]
    @test SM.validate(coordinated, SM.DielectricGroup,
        (cable_id = "fixed-contact",), one_layer, 3.0) === coordinated
end

@testitem "Measurements / internal shunt / correlated implicit sensitivities" tags=[:extension, :measurements, :slow] begin
    using Measurements, LinearAlgebra, Random, Statistics
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    epsilon = measurement(2.5, 0.0025)
    design = internal_shunt_test_design(; epsilon, count = 4)
    group, _ = internal_shunt_test_group(design)
    level = (wire = 16, order = 8, quadrature = 64, modes = 128)
    boundary = E.ShuntModel.Formula(:boundary; options = (resolution = level,))
    propagated = E.ShuntModel.internal_capacitance(boundary, group)
    @test all(isfinite, Measurements.value.(propagated.C))
    @test all(isfinite, Measurements.uncertainty.(propagated.C))
    @test maximum(Measurements.uncertainty.(propagated.C)) > 0
    @test abs(Measurements.cor(vec(propagated.C))[1, 2]) ≈ 1 rtol=1e-8
    deterministic(eps) = begin
        d,
        _ = internal_shunt_test_group(internal_shunt_test_design(epsilon = eps, count = 4))
        E.ShuntModel.internal_capacitance(boundary, d).C
    end
    central = (deterministic(2.5001)-deterministic(2.4999))/0.0002
    @test Measurements.derivative.(propagated.C, Ref(epsilon)) ≈ central rtol=2e-3
    @test Measurements.uncertainty.(propagated.C) ≈ abs.(central)*0.0025 rtol=2e-3
    # Equal nominal values must not cause sharing across independent inputs.
    independent,
    _ = internal_shunt_test_group(internal_shunt_test_design(
        epsilon = measurement(2.5, 0.0025), count = 4))
    @test !all(pair -> first(pair) === last(pair), zip(Vector(group), Vector(independent)))
    rng = MersenneTwister(314)
    samples = [deterministic(2.5+0.0025randn(rng))[1, 1] for _ in 1:24]
    # 3 standard errors of the sampled standard deviation. No requirement
    # that a finite MC sample exactly equals the linearized variance.
    sigma = Measurements.uncertainty(propagated.C[1, 1])
    @test abs(std(samples)-sigma) < 3sigma/sqrt(2(length(samples)-1))
    @test abs(mean(samples)-Measurements.value(propagated.C[1, 1])) <
          4sigma/sqrt(length(samples))
    radius = measurement(0.2e-3, 0.1e-6)
    geometry_group,
    _ = internal_shunt_test_group(internal_shunt_test_design(; radius, count = 4))
    geometry_result = E.ShuntModel.internal_capacitance(boundary, geometry_group)
    function varying_radius(r)
        d, _ = internal_shunt_test_group(internal_shunt_test_design(radius = r, count = 4))
        E.ShuntModel.internal_capacitance(boundary, d).C
    end
    geometric_central = (varying_radius(0.2e-3+1e-8)-varying_radius(0.2e-3-1e-8))/2e-8
    @test Measurements.derivative.(geometry_result.C, Ref(radius)) ≈ geometric_central rtol=5e-3
    @test maximum(Measurements.uncertainty.(geometry_result.C)) > 0
    center = measurement(0.0, 1e-6)
    other_center = measurement(0.0, 1e-6)
    @test E.ShuntModel._shunt_same(center, center, 0.01)
    @test !E.ShuntModel._shunt_same(center, other_center, 0.01)
    @test !E.ShuntModel._shunt_same(center, 0.0, 0.01)
    wire_angle = measurement(0.31, 1e-4)
    angle_group,
    _ = internal_shunt_test_group(internal_shunt_test_design(; wire_angle, count = 4))
    angle_result = E.ShuntModel.internal_capacitance(boundary, angle_group)
    function varying_angle(angle)
        d,
        _ = internal_shunt_test_group(internal_shunt_test_design(wire_angle = angle, count = 4))
        E.ShuntModel.internal_capacitance(boundary, d).C
    end
    angular_central = (varying_angle(0.31001)-varying_angle(0.30999))/2e-5
    @test Measurements.derivative.(angle_result.C, Ref(wire_angle)) ≈ angular_central rtol=5e-3
    @test maximum(Measurements.uncertainty.(angle_result.C)) > 0
end

@testitem "Measurements / internal shunt / public parametric studies" tags=[
    :extension, :integration,
    :measurements,
    :slow] begin
    using Measurements, Statistics
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    space = Gridspace{CableConstantsProblem}(
        epsilon->CableConstantsProblem(
            internal_shunt_test_design(; epsilon, tapes = false, count = 4)),
        (Grid(2.5, AbsoluteError(0.0025)),))
    problem = ParametricProblem(space)
    inner = CableConstantsFormulation(shunt_model = formula(:boundary;
        options = (resolution = (wire = 16, order = 8, quadrature = 64, modes = 128),)))
    linear = compute(problem, LinearError(inner))
    sampled = compute(problem,
        MonteCarlo(inner; trials = 12, seed = 314,
            distribution = :normal, return_samples = true, retain_details = true))
    @test sampled.trial_counts == [12]
    @test all(isempty, sampled.details.data.failures)
    @test only(sampled.values).details.data.shunt_model.requested === :boundary
    @test only(linear.values).details.data.shunt_model.requested === :boundary
    IE=LineCableModels.ImportExport
    restored=IE.deserialize_value(IE.serialize_value(linear))
    original=only(linear.values)
    saved=only(restored.values)
    @test Measurements.value.(saved.C) == Measurements.value.(original.C)
    @test Measurements.uncertainty.(saved.C) ≈ Measurements.uncertainty.(original.C)
    @test details(saved).data.shunt_model.requested === :boundary
    # One shared uncertain input must remain one source across distinct points.
    duplicated=LineCableModels.UQ.LinearErrorResult(
        LinearError(inner), [original, original], ComputationDetails((;)))
    copied=IE.deserialize_value(IE.serialize_value(duplicated))
    @test Measurements.uncertainty(copied.values[1].C[1]-copied.values[2].C[1]) == 0
    strict = CableConstantsFormulation(shunt_model = formula(:boundary;
        options = (resolution = (wire = 16, order = 8, quadrature = 64, modes = 128),)),
        insulation_admittance = :lossy)
    @test_throws BoundarySolveError compute(
        problem, MonteCarlo(strict;
            trials = 2, seed = 314, on_error = :resample, retain_details = true))
    # Sharing solved groups must not erase independent Measurements sources.
    independent = CableConstantsProblem(internal_shunt_test_design(
        epsilon = measurement(2.5, 0.0025), tapes = false, count = 4))
    original_problem = first(space)
    designs = [original_problem.design, independent.design]
    blueprints = [LineCableModels.DataModel.flatten(design,
                      eltype(original_problem)) for design in designs]
    blocks = only(LineCableModels.Engine.internal_potential_coefficients(designs, blueprints,
        [inner])).blocks
    @test blocks[1].C !== blocks[2].C
    @test Measurements.uncertainty(blocks[1].C[1, 1] - blocks[2].C[1, 1]) > 0
    cap = only(only(linear.values).C)
    @test Measurements.uncertainty(cap) > 0
    @test abs(only(only(sampled.values).C)-Measurements.value(cap)) <
          4Measurements.uncertainty(cap)/sqrt(12)
end
