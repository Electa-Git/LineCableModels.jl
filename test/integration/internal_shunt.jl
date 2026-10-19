@testitem "Engine / internal shunt / explicit boundary workflow and reuse" tags=[:integration, :parametric, :slow] begin
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
    results = compute(problem, Grid((first_formula, other_formula, reduced_formula));
        options = (trace = true,))
    first_result = first(results)
    shunt = details(first_result).data.shunt_model
    @test shunt.requested === :boundary
    @test shunt.solves == 1
    @test length(shunt.domains) == 2
    @test details(results[2]).data.trace.Pin == details(first_result).data.trace.Pin
    @test typeof(results[1]) == typeof(results[2])
    trace = details(first_result).data.trace
    options = reduced_formula.options.data
    plan = LineCableModels.Commons.ReductionPlan(problem.system.connection_order;
        options.reduce_bundle, options.kron_reduction, options.ideal_transposition)
    buffers = LineCableModels.Commons.ReductionBuffers{eltype(trace.Z)}(plan)
    for k in eachindex(problem.frequencies)
        reduced_Z = similar(trace.Z, length(plan.keep), length(plan.keep))
        reduced_Y = similar(reduced_Z)
        LineCableModels.Commons.reduce_line_matrices!(reduced_Z, reduced_Y, trace.Z[:, :, k],
            trace.P[:, :, k], 2pi*im*problem.frequencies[k], plan, buffers)
        @test observe(results[3], Z)[:, :, k] ≈ reduced_Z
        @test observe(results[3], Y)[:, :, k] ≈ reduced_Y rtol=1e-10
    end
    @test details(results[3]).data.shunt_model.solves == 1
    constants = @inferred CableConstants(design; frequency = 50.0,
        formulation = CableConstantsFormulation(shunt_model = boundary))
    local_group, blueprint = internal_shunt_test_group(design)
    response = only(E.internal_potential_coefficients([design], [blueprint],
        [CableConstantsFormulation(shunt_model = boundary)]))
    @test constants.C[1] ≈ response.blocks[1].C[1, 1] rtol=1e-10
    @test constants.G[1] == 0
    lossy = Formulation(insulation_admittance = :lossy; options = physical)
    mixed = compute(problem, Grid((first_formula, lossy)))
    @test length(mixed) == 2
    @test typeof(mixed[1]) == typeof(mixed[2])
    @test details(mixed[2]).data.shunt_model.requested === :equivalent
    @test details(mixed[2]).data.shunt_model.solves == 0
    reference_temperature = Formulation(shunt_model = boundary,
        temperature_dependence = nothing; options = physical)
    blueprints = [LineCableModels.DataModel.flatten(design, Float64)
                  for design in system.designs]
    shunts = E.internal_potential_coefficients(system.designs, blueprints,
        [first_formula, other_formula, lossy, reference_temperature])
    # Formulations whose groups have the same physical inputs share one solve.
    @test all(index -> shunts[index].blocks[1].C === shunts[1].blocks[1].C, (2, 4))
    @test all(index -> shunts[index].details.solves == 1, (1, 2, 4))
    @test isempty(shunts[3].blocks)
    @test shunts[1].blocks[1].terminals == 1:3
    @test shunts[1].blocks[2].terminals == 4:6
    @test shunts[1].blocks[1].C === shunts[1].blocks[2].C
    resolved = LineCableModels.Resolved(problem, (; line = E.lineinput(problem, blueprints),
        shunt = shunts[2], physical_inputs = E.completed_inputs(problem)))
    @test observe(compute(resolved, other_formula), Y) == observe(results[2], Y)
    # The same physical API supports independent translated assemblies inside
    # one design. Their local operators must be scattered into distinct ports.
    left = internal_shunt_test_design(count = 4, suffix = "_left")
    right = internal_shunt_test_design(count = 4, suffix = "_right")
    pair = build(CableDesign, "two-local-domains",
        assembly(at(left.tree, -0.01, 0.0), at(right.tree, 0.01, 0.0)))
    pair_blueprint = LineCableModels.DataModel.flatten(pair)
    groups = E.ShuntModel.dielectric_groups(first_formula.formulas.shunt_model, pair,
        pair_blueprint, first_formula.formulas, 1).groups
    @test length(groups) == 2
    @test [group.terminals for group in groups] == [1:3, 4:6]
    @test CableConstants(pair; formulation = CableConstantsFormulation(shunt_model = boundary)).C ≈
          fill(constants.C[1], 2) rtol=1e-8
end
