@testitem "Engine / solver / multicable reciprocity and modal transformation" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    using LinearAlgebra
    using Logging: NullLogger, with_logger

    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    formulation=Formulation(;
        options = (
        reduce_bundle = true,
        kron_reduction = true,
        ideal_transposition = false
    )
    )
    phase_parameters=compute(problem, formulation; options = (trace = true,))
    trace=details(phase_parameters).trace
    parameters=compute(
        ModalTransformationProblem(phase_parameters),
        ModalTransformationFormulation(:default)
    )
    @test details(parameters).formulations === details(phase_parameters).formulations

    @test domain(parameters) === ModalDomain
    @test size(parameters.Z) == (3, 3, 2)
    @test all(isfinite, parameters.Z)
    @test all(isfinite, parameters.Y)
    @test trace.phase_map == [1, 0, 0, 2, 0, 0, 3, 0, 0]
    @test trace.cable_map == repeat(1:3; inner = 3)
    @test size(trace.Zin) == (9, 9, 2)
    @test size(trace.Pin) == (9, 9, 2)
    @test size(trace.Z) == (9, 9, 2)
    @test size(trace.P) == (9, 9, 2)
    @test size(trace.Zg) == (3, 3, 2)
    @test size(trace.Pg) == (3, 3, 2)
    for identifier in (:default,)
        tracked=@inferred compute(
            ModalTransformationProblem(phase_parameters),
            ModalTransformationFormulation(identifier)
        )
        rebuilt=@inferred compute(ModalTransformationProblem(tracked))
        @test TestNumerics.isapprox_scaled(Z(rebuilt), Z(phase_parameters))
        @test TestNumerics.isapprox_scaled(Y(rebuilt), Y(phase_parameters))
        for matrix in (tracked.Z.values, tracked.Y.values), frequency in 1:2

            slice=@view matrix[:, :, frequency]
            @test norm(slice-Diagonal(diag(slice))) <=
                  1e-6*max(norm(slice), eps(Float64))
        end
    end
    for frequency_index in eachindex(problem.frequencies)
        @test TestNumerics.isapprox_scaled(
            trace.Zg[:, :, frequency_index],
            transpose(trace.Zg[:, :, frequency_index])
        )
        @test TestNumerics.isapprox_scaled(
            trace.Pg[:, :, frequency_index],
            transpose(trace.Pg[:, :, frequency_index])
        )
        @test any(!iszero,
            trace.Zg[:, :, frequency_index] .-
            Diagonal(diag(trace.Zg[:, :, frequency_index])))
        @test any(!iszero,
            trace.Pg[:, :, frequency_index] .-
            Diagonal(diag(trace.Pg[:, :, frequency_index])))
    end

    execution=computation_options(LineCableModelsCoaxial, (;))
    blueprints=LineCableModels.Engine.CableBlueprint{eltype(problem)}[LineCableModels.Engine.flatten(
                                                                          LineCableModelsCoaxial(),
                                                                          design,
                                                                          eltype(problem))
                                                                      for design in problem.system.designs]
    workspace=LineParametersWorkspace(
        problem, Formulation(), execution, blueprints)
    @test workspace.input.phase_map == problem.system.connection_order
    @test workspace.input.cable_map ==
          [entry.cable for entry in problem.system.terminal_order]
    @test workspace.invariants.cable_indices ==
          [findall(entry -> entry.cable == cable, problem.system.terminal_order)
           for cable in 1:ncables(problem.system)]
    capture_allocations(input)=@allocated LineCableModels.Engine._capture_buffers(
        Float64, input, Val(false))
    capture_allocations(workspace.input)
    @test capture_allocations(workspace.input) <= 1024
    @test LineCableModels.Engine._capture_buffers(Float64, workspace.input, Val(false)) ===
          nothing

    allocation_formulation=Formulation(
        earth_impedance = :Pollaczek1926,
        earth_admittance = :IdealGround,
        options = (
            reduce_bundle = true,
            kron_reduction = true,
            ideal_transposition = false
        )
    )
    allocation_workspace=LineParametersWorkspace(
        problem,
        allocation_formulation,
        execution,
        blueprints
    )
    function solve_without_logging(workspace, formulation)
        return with_logger(NullLogger()) do
            LineCableModels.Engine._solve!(workspace, formulation)
        end
    end
    @test @inferred(LineCableModels.Engine._solve!(
        allocation_workspace,
        allocation_formulation
    )) isa LineParameters
    solve_without_logging(allocation_workspace, allocation_formulation)
    allocations=@allocated solve_without_logging(
        allocation_workspace,
        allocation_formulation
    )
    @test allocations <= 4_096
end

@testitem "Engine / Gridpoint / selected line problem reaches scalar compute" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    import LineCableModels.ParametricBuilder as PB

    system=TestFixtures.three_phase_system()
    earth=Grid((EarthModel(10.0), EarthModel(100.0)))
    problem_space=LineParametersProblem(
        system,
        earth;
        frequencies = [50.0]
    )
    @test problem_space isa Gridspace{LineParametersProblem}
    @test length(problem_space) == 2

    point=first(PB.points(problem_space))
    @test point isa LineCableModels.Gridpoint{LineParametersProblem}
    result=compute(point, Formulation())
    @test result isa LineParameters
    @test size(result.Z) == (3, 3, 1)
    @test size(result.Y) == (3, 3, 1)
end

@testitem "Engine / coaxial choreography / local formulas precede earth formulas" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    const EN=LineCableModels.Engine
    const II=EN.InsulationImpedance
    const IA=EN.InsulationAdmittance
    const SA=EN.SemiconAdmittance
    const EZ=EN.EarthImpedance
    const EY=EN.EarthAdmittance

    events=Symbol[]
    insulation_impedance=II.Formula(:default).binding
    local_z=(
        r_in, r_ex,
        mu_r,
        s,
        values,
        options,
        workspace)->begin
        push!(events, :local_z)
        insulation_impedance(r_in, r_ex, mu_r, s, values, options, workspace)
    end
    insulation=IA.Formula(:Ametani2004).binding
    local_insulation_y=(material, frequency, temperature,
        values, options,
        workspace)->begin
        push!(events, :local_y)
        insulation(material, frequency, temperature, values, options, workspace)
    end
    semicon=SA.Formula(:Ametani2004).binding
    local_semicon_y=(material, frequency, temperature,
        values, options,
        workspace)->begin
        push!(events, :local_y)
        semicon(material, frequency, temperature, values, options, workspace)
    end

    earth_z=(
        functor, pair, workspace)->begin
        push!(events, :earth_z)
        functor.binding.equation(functor, pair, workspace)
    end
    earth_y=(
        functor, pair, workspace)->begin
        push!(events, :earth_y)
        functor.binding.equation(functor, pair, workspace)
    end
    for (operation,
        identifier,
        callback) in (
        (II.insulation_impedance, :default, local_z),
        (IA.insulation_material, :Ametani2004, local_insulation_y),
        (SA.semicon_material, :Ametani2004, local_semicon_y))
        @eval LineCableModels.computation_options(
            ::LineCableModels.FormulaMethod{$(QuoteNode(identifier)), typeof($operation)},
            ::$(typeof(callback))) = (;)
    end
    for (operation,
        callback) in (
        (EZ.earth_impedance, earth_z), (EY.earth_potential_coefficient, earth_y))
        @eval LineCableModels.computation_options(
            ::LineCableModels.FormulaMethod{:default, typeof($operation)},
            ::$(typeof(callback))) = (integration = (method = :quad, options = (;)),)
    end
    formulation=Formulation(
        insulation_impedance = formula(:default; hooks = (contribution = local_z,)),
        insulation_admittance = formula(:Ametani2004; hooks = (contribution = local_insulation_y,)),
        semicon_admittance = formula(:Ametani2004; hooks = (contribution = local_semicon_y,)),
        earth_impedance = formula(:default; hooks = (contribution = earth_z,)),
        earth_admittance = formula(:default; hooks = (contribution = earth_y,)),
        options = (ideal_transposition = false,))
    result=compute(
        TestFixtures.line_parameters_problem(frequencies = [50.0]),
        formulation
    )

    local_z_calls=findall(==(:local_z), events)
    earth_z_calls=findall(==(:earth_z), events)
    local_y_calls=findall(==(:local_y), events)
    earth_y_calls=findall(==(:earth_y), events)
    @test all(!isempty, (local_z_calls, earth_z_calls, local_y_calls, earth_y_calls))
    @test maximum(local_z_calls) < minimum(earth_z_calls)
    @test maximum(earth_z_calls) < minimum(local_y_calls)
    @test maximum(local_y_calls) < minimum(earth_y_calls)
    @test all(isfinite, result.Z)
    @test all(isfinite, result.Y)
    single_frequency_events=copy(events)
    empty!(events)
    sweep=compute(
        TestFixtures.line_parameters_problem(frequencies = [1.0, 50.0, 1000.0]),
        formulation
    )
    @test events == repeat(single_frequency_events, 3)
    @test domain(sweep) === PhaseDomain
    @test sweep.Z.values[:, :, 2] == result.Z.values[:, :, 1]
    @test sweep.Y.values[:, :, 2] == result.Y.values[:, :, 1]
end

@testitem "Engine / transform / symmetric two-cable system retains two modes" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    using LinearAlgebra: Diagonal, diag, norm

    design=TestFixtures.mv_cable_design()
    connections=[Dict(terminal=>(terminal===:core ? phase : 0)
                 for terminal in design.terminal_order)
                 for phase in 1:2]
    system=build(
        LineCableSystem,
        [design, design],
        [Pose2(-0.5, -1.0), Pose2(0.5, -1.0)];
        connections,
        system_id = "symmetric-pair"
    )
    problem=LineParametersProblem(
        system;
        earth_props = homogeneous(rho = 100.0),
        frequencies = [50.0, 500.0]
    )
    parameters=compute(
        problem,
        Formulation(options = (ideal_transposition = false,))
    )
    @test size(Z(parameters)) == (2, 2, 2)
    @test size(Y(parameters)) == (2, 2, 2)
    for matrix in (Z(parameters), Y(parameters)), frequency in 1:2

        slice=@view matrix[:, :, frequency]
        @test slice ≈ transpose(slice)
        @test slice[1, 1] ≈ slice[2, 2]
        @test slice[1, 2] ≈ slice[2, 1] atol=100eps(Float64)*norm(slice)
    end

    modal=compute(
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:default); options = (offdiagonal_tolerance = 1e-10,)
    )
    @test domain(modal) === ModalDomain
    for matrix in (Z(modal), Y(modal)), frequency in 1:2

        slice=@view matrix[:, :, frequency]
        @test norm(slice - Diagonal(diag(slice))) <= 1e-10 * max(norm(slice), 1)
    end
    @test size(@observe(modal, (R, diag)[:, :])) == (2, 2)
    @test size(@observe(modal, (L, diag)[:, :])) == (2, 2)
    @test size(@observe(modal, (G, diag)[:, :])) == (2, 2)
    @test size(@observe(modal, (C, diag)[:, :])) == (2, 2)
end

@testitem "Engine / formulation boundary / physical geometry precedes backend support" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    conductor=Material(kind = :conductor, rho = 1.7241e-8)
    design=build(
        CableDesign,
        "elliptical",
        Group(:phase, Region(:elliptical_core, Ellipse(0.01, 0.006), conductor))
    )
    @test design.geometry.regions[1].primitive isa LineCableModels.DataModel.Ellipse

    system=build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0);
        connections = (phase = 1,)
    )
    problem=LineParametersProblem(
        system;
        earth_props = EarthModel(100.0),
        frequencies = [50.0]
    )
    @test problem.system === system
    execution=computation_options(LineCableModelsCoaxial, (;))
    @test_throws ArgumentError LineParametersWorkspace(
        problem,
        Formulation(),
        execution,
        LineCableModels.Engine.CableBlueprint{eltype(problem)}[LineCableModels.Engine.flatten(
                                                                   LineCableModelsCoaxial(), source,
                                                                   eltype(problem)
                                                               )
                                                               for source in
                                                                   problem.system.designs]
    )
    @test_throws ArgumentError compute(problem, Formulation())
    @test_throws ArgumentError CableConstants(design)

    neutral=build(
        LineCableSystem,
        design,
        Pose2(0.0, 0.0);
        connections = (phase = 1,),
        environment = nothing
    )
    @test_throws DomainError LineParametersProblem(
        neutral;
        earth_props = EarthModel(100.0),
        frequencies = [50.0]
    )
end

@testitem "Engine / coaxial profile / explicit radial support boundary" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    conductor=Material(kind = :conductor, rho = 1.7241e-8)
    dielectric=Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    earth=EarthModel(100.0)

    function problem_for(design)
        system=build(
            LineCableSystem,
            design,
            Pose2(0.0, -1.0);
            connections = Dict(first(design.terminal_order)=>1)
        )
        return LineParametersProblem(
            system;
            earth_props = earth,
            frequencies = [50.0]
        )
    end

    bare=build(CableDesign, "bare", Group(
        :phase, Region(:bare, Disk(0.01), conductor)
    ))
    execution=computation_options(LineCableModelsCoaxial, (;))
    workspace(problem,
        formulation = Formulation())=LineParametersWorkspace(
        problem,
        formulation,
        execution,
        LineCableModels.Engine.CableBlueprint{eltype(problem)}[LineCableModels.Engine.flatten(
                                                                   LineCableModelsCoaxial(), source,
                                                                   eltype(problem)
                                                               )
                                                               for source in problem.system.designs]
    )

    bare_input=workspace(problem_for(bare)).input
    @test bare_input.cable.r_ins_in == bare_input.cable.r_ins_ext == [0.01]
    @test bare_input.cable.dielectric_ranges == [1:0]
    @test isempty(bare_input.cable.r_layer_in)

    layered=build(CableDesign,
        "layered",
        Stack(
            Group(:phase,
                Stack(
                    Region(:core, Disk(0.01), conductor),
                    Region(:outer_conductor, Annulus(0.01, 0.011), conductor)
                )),
            Region(:semicon, Shell(0.0005), dielectric),
            Region(:insulation, Shell(0.002), dielectric),
            Group(
                :screen,
                Region(:screen_metal, Annulus(0.0135, 0.014), conductor)
            )
        ))
    layered_input=workspace(problem_for(layered)).input
    @test layered_input.n_phases == 2
    @test layered_input.cable.dielectric_ranges == [1:2, 3:2]
    @test layered_input.cable.r_ins_in[2] ==
          layered_input.cable.r_ins_ext[2] == 0.014

    filled=build(CableDesign,
        "filled",
        Enclosure(
            :pipe,
            Group(:phase, Region(:core, Disk(0.01), conductor));
            primitive = Disk(0.03),
            fill = dielectric
        ))
    filled_input=workspace(problem_for(filled)).input
    @test filled_input.cable.r_layer_in == [0.01]
    @test filled_input.cable.r_layer_ext == [0.03]
    @test getproperty.(filled_input.cable.dielectric_materials, :rho) ==
          [dielectric.rho]

    reappearing=build(CableDesign,
        "reappearing",
        Stack(
            Group(:a, Region(:a_inner, Disk(0.01), conductor)),
            Region(:a_gap, Shell(0.001), dielectric),
            Group(:b, Region(:b_metal, Annulus(0.011, 0.012), conductor)),
            Region(:b_gap, Shell(0.001), dielectric),
            Group(:a, Region(:a_outer, Annulus(0.013, 0.014), conductor))
        ))
    @test_throws ArgumentError workspace(problem_for(reappearing))

    conductor_after_dielectric=build(CableDesign,
        "conductor-after-dielectric",
        Stack(
            Group(:phase, Region(:inner, Disk(0.01), conductor)),
            Region(:gap, Shell(0.001), dielectric),
            Group(:phase, Region(:outer, Annulus(0.011, 0.012), conductor))
        ))
    @test_throws ArgumentError workspace(problem_for(conductor_after_dielectric))

    undeclared_gap=build(CableDesign,
        "undeclared-gap",
        Stack(
            Group(:inner, Region(:inner, Disk(0.01), conductor)),
            Group(:outer, Region(:outer, Annulus(0.012, 0.013), conductor))
        ))
    @test_throws ArgumentError workspace(problem_for(undeclared_gap))

    @test_throws DomainError build(CableDesign,
        "overlap",
        Stack(
            Group(:inner, Region(:inner, Disk(0.01), conductor)),
            Group(:outer, Region(:outer, Annulus(0.009, 0.012), conductor))
        ))
end

@testitem "Engine / solver / bundle-only and singleton reduction policies" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    design=TestFixtures.mv_cable_design()
    duplicate_mapping=Dict("core"=>1, "sheath"=>1, "jacket"=>0)
    duplicate_system=build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0, 0.0);
        connections = duplicate_mapping,
        system_id = "duplicate-bundle",
        line_length = 500.0
    )
    frequencies=[50.0]
    earth=EarthModel(100.0, 10.0, 1.0)
    duplicate_problem=LineParametersProblem(
        duplicate_system;
        earth_props = earth,
        frequencies
    )
    bundle_only=Formulation(;
        temperature_dependence=nothing,
        options = (
        reduce_bundle = true,
        kron_reduction = false,
        ideal_transposition = true
    )
    )
    duplicate_result=@inferred compute(duplicate_problem, bundle_only)
    @test size(duplicate_result.Z) == (2, 2, 1)
    @test all(isfinite, duplicate_result.Z)
    @test all(isfinite, duplicate_result.Y)

    singleton_design=build(
        CableDesign,
        "single-component",
        Stack(deepcopy(design.root.items[1:10]))
    )
    singleton_system=build(
        LineCableSystem,
        singleton_design,
        Pose2(0.0, -1.0, 0.0);
        connections = Dict(:core=>1),
        system_id = "singleton",
        line_length = 100.0
    )
    singleton_problem=LineParametersProblem(
        singleton_system;
        earth_props = earth,
        frequencies
    )
    unreduced=Formulation(;
        temperature_dependence=nothing,
        options = (
        reduce_bundle = false,
        kron_reduction = false,
        ideal_transposition = true
    )
    )
    singleton_result=compute(singleton_problem, unreduced)
    @test size(singleton_result.Z) == (1, 1, 1)
    @test real(singleton_result.Z[1, 1, 1]) > 0
    @test imag(singleton_result.Y[1, 1, 1]) > 0
end

@testitem "Engine / indexed restrictions and Γ overrides reach public compute" tags=[:integration] setup=[
    EngineTestSupport, UseEngineSupport, TestFixtures
] begin
    base=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    explicit=LineParametersProblem(base.system; earth_props = base.earth_props,
        frequencies = base.frequencies, Γ = [1e-5im, 2e-5im])
    @test_throws ArgumentError compute(explicit, Formulation())
    calls=Tuple{Int, Int}[]
    prescription=(s, m, l)->(push!(calls, l); zero(s))
    selected=Formulation(earth_impedance = formula(:default; hooks = (Γ = prescription,)))
    result=compute(base, selected)
    ordinary=compute(base)
    @test result.Z.values == ordinary.Z.values
    @test !isempty(calls) && all(==((2, 2)), calls)
    @test details(result).formulations.modified.earth_impedance
    @test !details(ordinary).formulations.modified.earth_impedance
    zero_problem=LineParametersProblem(base.system; earth_props = base.earth_props,
        frequencies = base.frequencies, Γ = zeros(ComplexF64, 2))
    @test_throws ArgumentError compute(zero_problem, selected)
    @test_throws DimensionMismatch LineParametersProblem(base.system;
        earth_props = base.earth_props, frequencies = base.frequencies, Γ = [0.0])
    design=TestFixtures.mv_cable_design()
    connections(phase)=Dict("core"=>phase, "sheath"=>0, "jacket"=>0)
    mixed=build(LineCableSystem, [design, design], [Pose2(0.0, 1.0), Pose2(1.0, -1.0)];
        connections = [connections(1), connections(2)])
    problem=LineParametersProblem(mixed; earth_props = EarthModel(100.0), frequencies = [50.0])
    @test_throws ArgumentError compute(problem)
    @test_throws ArgumentError compute(problem,
        Formulation(
            earth_impedance = formula(:default; hooks = (contribution = (
            f, p, w)->zero(f.state.jω),))))
end

@testitem "Engine / frequency-dependent earth relation reaches coaxial solve" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    problem=TestFixtures.line_parameters_problem(frequencies = [1.0e6])
    static=compute(problem, Formulation())
    law=(m, f, p, o,
        w)->LineCableModels.Earth.EarthMaterial(m.rho/(1+f/1e5), m.eps_r, m.mu_r)
    @eval LineCableModels.computation_options(
        ::LineCableModels.FormulaMethod{
            :default, typeof(LineCableModels.Earth.FrequencyDependent.earth_material)},
        ::$(typeof(law))) = (;)
    dispersive=compute(problem,
        Formulation(earth_properties = formula(:default; hooks = (contribution = law,))))

    @test all(isfinite, dispersive.Z)
    @test all(isfinite, dispersive.Y)
    @test dispersive.Z.values != static.Z.values
    @test dispersive.Y.values != static.Y.values
end
