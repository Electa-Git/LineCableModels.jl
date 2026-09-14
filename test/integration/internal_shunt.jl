@testitem "Engine / internal shunt / default workflow and reuse" tags=[:integration] begin
    using LinearAlgebra
    E = LineCableModels.Engine
    include(joinpath(pkgdir(LineCableModels),"test","fixtures","internal_shunt.jl"))
    design = internal_shunt_test_design(count=4)
    system = build(LineCableSystem,[design,design],[(-0.02,-1.0),(0.02,-1.0)];
        connections=[Dict(:inner=>1,:middle=>0,:reference=>0),
            Dict(:inner=>2,:middle=>0,:reference=>0)])
    problem = LineParametersProblem(system;earth_props=homogeneous(rho=100.0),frequencies=[0.1,50.0,1e7])
    physical = (reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    first_formula = Formulation(earth_impedance=:Pollaczek1926,
        earth_admittance=:Pollaczek1926;options=physical)
    other_formula = Formulation(earth_impedance=:Saad1996,
        earth_admittance=:Pollaczek1926;options=physical)
    reduced_formula = Formulation(earth_impedance=:Pollaczek1926,
        earth_admittance=:Pollaczek1926;options=(reduce_bundle=true,
            kron_reduction=true,ideal_transposition=false))
    results = compute(problem,[first_formula,other_formula,reduced_formula];options=(trace=true,))
    first_result = first(results)
    shunt = details(first_result).internal_shunt
    @test shunt.treatment === :resolved_local
    @test shunt.solves == 1
    @test length(shunt.domains) == 2
    @test details(results[2]).trace.Pin == details(first_result).trace.Pin
    @test typeof(results[1]) == typeof(results[2])
    trace = details(first_result).trace
    reduced = E.reduce_primitive_matrices(trace.Z,trace.P,
        problem.system.connection_order,reduced_formula.options)
    @test observe(results[3],Z) ≈ reduced.Z
    for k in eachindex(problem.frequencies)
        @test observe(results[3],Y)[:,:,k] ≈
            (2pi*im*problem.frequencies[k]).*inv(reduced.P[:,:,k]) rtol=1e-10
    end
    @test details(results[3]).internal_shunt.solves == 1
    bp = E.flatten.(Ref(LineCableModelsCoaxial()),problem.system.designs)
    input = E.lineinput(problem,bp)
    # Explicit test of the old internal arithmetic, without a public toggle or
    # a material-law change that would alter the physical reference itself.
    legacy_input = merge(input,(prepared_shunt=nothing,))
    execution = E.computation_options(LineCableModelsCoaxial,(trace=true,))
    legacy = E._compute(LineCableModelsCoaxial(),problem,first_formula,execution,legacy_input)
    @test observe(first_result,Z) == observe(legacy,Z)
    @test details(first_result).trace.Pg == details(legacy).trace.Pg
    @test details(first_result).trace.Pin != details(legacy).trace.Pin
    @test details(first_result).trace.P[3:3:6,:,:] == details(legacy).trace.P[3:3:6,:,:]
    @test details(first_result).trace.P[:,3:3:6,:] == details(legacy).trace.P[:,3:3:6,:]
    @test details(first_result).trace.P[1:3,4:6,:] == details(legacy).trace.P[1:3,4:6,:]
    constants = @inferred CableConstants(design;frequency=50.0)
    local_domain,blueprint = internal_shunt_test_domain(design)
    prepared = E.prepare_internal_shunt([local_domain],3,Formulation().methods,50.0,20.0)
    @test constants.C[1] ≈ prepared.blocks[1].C[1,1] rtol=1e-10
    @test constants.G[1] == 0
    lossy = Formulation(insulation_admittance=:Ametani2004,
        earth_impedance=:Pollaczek1926,earth_admittance=:Pollaczek1926;options=physical)
    mixed = compute(problem,[first_formula,lossy])
    @test length(mixed) == 2
    @test typeof(mixed[1]) == typeof(mixed[2])
    @test details(mixed[2]).internal_shunt.reason === :equivalent_material_law
    # The same physical API supports independent translated assemblies inside
    # one design; their local operators must be scattered into distinct ports.
    left = internal_shunt_test_design(count=4,suffix="_left")
    right = internal_shunt_test_design(count=4,suffix="_right")
    pair = build(CableDesign,"two-local-domains",
        assembly(at(left.origin,-0.01,0.0),at(right.origin,0.01,0.0)))
    pair_blueprint = E.flatten(LineCableModelsCoaxial(),pair)
    domains = E.internal_shunt_domains([pair],[pair_blueprint])
    @test length(domains) == 2
    @test getproperty.(domains,:terminals) == [1:3,4:6]
    @test E._shunt_domain_equal(domains...)
    @test CableConstants(pair).C ≈ fill(constants.C[1],2) rtol=1e-8
end
