@testitem "Engine / internal shunt / explicit boundary workflow and reuse" tags=[:integration] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels), "test", "support", "internal_shunt.jl"))
    design = internal_shunt_test_design(count = 4)
    system = build(LineCableSystem, [design, design], [(-0.02, -1.0), (0.02, -1.0)];
        connections = [Dict(:inner=>1, :middle=>0, :reference=>0),
            Dict(:inner=>2, :middle=>0, :reference=>0)])
    problem = LineParametersProblem(
        system; earth_props = homogeneous(rho = 100.0), frequencies = [0.1, 50.0, 1e7])
    physical = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
    boundary = formula(
        :boundary; options = (resolution = (
            wire = 32, order = 16, quadrature = 128, modes = 256),))
    first_formula = Formulation(shunt_model = boundary; options = physical)
    other_formula = Formulation(shunt_model = boundary,
        earth_impedance = formula(:unified; options=(Γ=1e-4im,)); options = physical)
    reduced_formula = Formulation(shunt_model = boundary; options = (reduce_bundle = true,
            kron_reduction = true, ideal_transposition = false))
    results = compute(problem, [first_formula, other_formula, reduced_formula]; options = (trace = true,))
    first_result = first(results)
    shunt = details(first_result).data.shunt_model
    @test shunt.effective === :boundary
    @test shunt.solves == 1
    @test length(shunt.domains) == 2
    @test details(results[2]).data.trace.Pin == details(first_result).data.trace.Pin
    @test typeof(results[1]) == typeof(results[2])
    trace = details(first_result).data.trace
    reduced = E.reduce_primitive_matrices(trace.Z, trace.P,
        problem.system.connection_order, reduced_formula.options)
    @test observe(results[3], Z) ≈ reduced.Z
    for k in eachindex(problem.frequencies)
        @test observe(results[3], Y)[:, :, k] ≈
              (2pi*im*problem.frequencies[k]) .* inv(reduced.P[:, :, k]) rtol=1e-10
    end
    @test details(results[3]).data.shunt_model.solves == 1
    constants = @inferred CableConstants(design; frequency = 50.0,
        formulation = CableConstantsFormulation(shunt_model = boundary))
    local_domain, blueprint = internal_shunt_test_domain(design)
    blueprint = E.flatten(LineCableModelsCoaxial(), design,
        CableConstantsFormulation(shunt_model = boundary))
    @test constants.C[1] ≈ blueprint.shunt[1].C[1, 1] rtol=1e-10
    @test constants.G[1] == 0
    lossy = Formulation(insulation_admittance = :lossy; options = physical)
    mixed = compute(problem, [first_formula, lossy])
    @test length(mixed) == 2
    @test typeof(mixed[1]) == typeof(mixed[2])
    @test details(mixed[2]).data.shunt_model.effective === :coaxial
    @test details(mixed[2]).data.shunt_model.solves == 0
    reference_temperature = Formulation(shunt_model = boundary,
        temperature_dependence = nothing; options = physical)
    blueprints = E.flatten(LineCableModelsCoaxial(), system.designs, Float64,
        [first_formula, other_formula, lossy, reference_temperature])
    @test blueprints[1] === blueprints[2]
    @test blueprints[1] === blueprints[4]
    @test blueprints[1][1].shunt[1].C === blueprints[1][2].shunt[1].C
    @test all(bp -> isempty(bp.shunt), blueprints[3])
    input = E.lineinput(problem, blueprints[1])
    @test input.cable.shunt[1].terminals == 1:3
    @test input.cable.shunt[2].terminals == 4:6
    @test input.cable.shunt[2].C === blueprints[1][2].shunt[1].C
    execution = computation_options(LineCableModelsCoaxial, ComputationOptions((;)))
    @test observe(
        E._compute(
            LineCableModelsCoaxial(), problem, other_formula, execution, input), Y) ==
          observe(results[2], Y)
    # The same physical API supports independent translated assemblies inside
    # one design; their local operators must be scattered into distinct ports.
    left = internal_shunt_test_design(count = 4, suffix = "_left")
    right = internal_shunt_test_design(count = 4, suffix = "_right")
    pair = build(CableDesign, "two-local-domains",
        assembly(at(left.origin, -0.01, 0.0), at(right.origin, 0.01, 0.0)))
    pair_blueprint = E.flatten(LineCableModelsCoaxial(), pair)
    domains = E.ShuntModel.internal_shunt_domains([pair], [pair_blueprint])
    @test length(domains) == 2
    @test getproperty.(domains, :terminals) == [1:3, 4:6]
    @test CableConstants(pair; formulation = CableConstantsFormulation(shunt_model = boundary)).C ≈
          fill(constants.C[1], 2) rtol=1e-8
end
