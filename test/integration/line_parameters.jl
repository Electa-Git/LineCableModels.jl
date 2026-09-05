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
        ideal_transposition = false,
        temperature_correction = true
    )
    )
    phase_parameters=compute(problem, formulation; options = (trace = true,))
    trace=details(phase_parameters).trace
    parameters=compute(
        ModalTransformationProblem(phase_parameters),
        ModalTransformationFormulation(:Fortescue)
    )
    @test details(parameters) === details(phase_parameters)

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
    for identifier in (:Chrysochos2014, :Fan2009, :Wedepohl1996)
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
    LineCableModels.Engine._capture_buffers(
        Float64, workspace.input, Val(false))
    @test @allocated(LineCableModels.Engine._capture_buffers(
        Float64, workspace.input, Val(false))) <= 1024
    @test LineCableModels.Engine._capture_buffers(
        Float64, workspace.input, Val(false)) === nothing

    allocation_formulation=Formulation(
        earth_impedance = :Pollaczek1926,
        earth_admittance = :IdealGround,
        options = (
            reduce_bundle = true,
            kron_reduction = true,
            ideal_transposition = false,
            temperature_correction = true
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
    insulation_impedance=II.Formula(:Ametani1980).route
    local_z=(
        r_in, r_ex, mu_r, s, values)->begin
        push!(events, :local_z)
        insulation_impedance(r_in, r_ex, mu_r, s, values)
    end
    insulation=IA.Formula(:Ametani2004).route
    local_insulation_y=(material, frequency, temperature,
        values)->begin
        push!(events, :local_y)
        insulation(material, frequency, temperature, values)
    end
    semicon=SA.Formula(:Ametani2004).route
    local_semicon_y=(material, frequency, temperature,
        values)->begin
        push!(events, :local_y)
        semicon(material, frequency, temperature, values)
    end

    earth_z_routes=EZ.routes(Val(:Papadopoulos2010))
    earth_z_self=(functor, pair)->begin
        push!(events, :earth_z)
        earth_z_routes.self(functor, pair)
    end
    earth_z_mutual=(functor, pair)->begin
        push!(events, :earth_z)
        earth_z_routes.mutual(functor, pair)
    end
    earth_y_routes=EY.routes(Val(:Papadopoulos2010))
    earth_y_self=(functor, pair)->begin
        push!(events, :earth_y)
        earth_y_routes.self(functor, pair)
    end
    earth_y_mutual=(functor, pair)->begin
        push!(events, :earth_y)
        earth_y_routes.mutual(functor, pair)
    end

    formulation=Formulation(
        insulation_impedance = formula(:Ametani1980; route = local_z),
        earth_impedance = formula(
            :Papadopoulos2010;
            self = earth_z_self,
            mutual = earth_z_mutual
        ),
        insulation_admittance = formula(
            :Ametani2004;
            route = local_insulation_y
        ),
        semicon_admittance = formula(
            :Ametani2004;
            route = local_semicon_y
        ),
        earth_admittance = formula(
            :Papadopoulos2010;
            self = earth_y_self,
            mutual = earth_y_mutual
        ),
        options = (ideal_transposition = false,)
    )
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
        TestFixtures.line_parameters_problem(frequencies=[1.0, 50.0, 1000.0]),
        formulation
    )
    @test events == repeat(single_frequency_events, 3)
    @test domain(sweep) === PhaseDomain
    @test sweep.Z.values[:, :, 2] == result.Z.values[:, :, 1]
    @test sweep.Y.values[:, :, 2] == result.Y.values[:, :, 1]
end

@testitem "Engine / preservation / two underground wires remain bit exact" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    conductor=Material(
        kind = :conductor,
        rho = eps(Float64),
        eps_r = 1.0,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.0
    )
    insulation=Material(
        kind = :insulator,
        rho = 1.97e14,
        eps_r = 2.3,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.0
    )
    design=build(
        CableDesign,
        "two-bare-wires-preservation",
        Stack(
            Group(:core, Region(:core_metal, Disk(0.0425), conductor)),
            Region(:core_insulation, Shell(1.0e-3), insulation)
        )
    )
    system=build(
        LineCableSystem,
        [design, design],
        [Pose2(0.0, -1.0), Pose2(1.0, -1.0)];
        connections = [Dict(:core=>1), Dict(:core=>2)]
    )
    problem=LineParametersProblem(
        system;
        earth_props = homogeneous(rho = 0.1, eps_r = 1.0, mu_r = 1.0),
        frequencies = [1.0, 50.0, 1000.0]
    )
    @test isconcretetype(fieldtype(typeof(problem), :earth_props))
    @test fieldtype(typeof(problem), :earth_props) === typeof(problem.earth_props)
    formulation=Formulation(
        internal_impedance = :Schelkunoff1934,
        insulation_impedance = :Ametani1980,
        earth_impedance = :Papadopoulos2010,
        insulation_admittance = :Ametani2004,
        semicon_admittance = :Ametani2004,
        earth_admittance = :Papadopoulos2010,
        options = (ideal_transposition = false,)
    )
    parameters=compute(problem, formulation)
    expected_Z=cat(
        ComplexF64[9.972567772765543e-7+1.0667880124003238e-5im 9.97043431689472e-7+6.699003625368725e-6im;
                   9.97043431689472e-7+6.699003625368725e-6im 9.972567772765543e-7+1.0667880124003238e-5im],
        ComplexF64[5.2463186340709196e-5+0.00040739430692732464im 5.226929347335312e-5+0.0002089818015430304im;
                   5.226929347335312e-5+0.0002089818015430304im 5.2463186340709196e-5+0.00040739430692732464im],
        ComplexF64[0.0011582715276183772+0.006050197037293153im 0.0011034983201066505+0.0020950707524030328im;
                   0.0011034983201066505+0.0020950707524030328im 0.0011582715276183772+0.006050197037293153im];
        dims = 3
    )
    expected_Y=cat(
        ComplexF64[1.3725717973458168e-12+3.456886991977605e-8im 1.1151765606538079e-15-1.0103448536100132e-16im;
                   1.1151765606538079e-15-1.0103448536100132e-16im 1.3725717973458172e-12+3.4568869919776055e-8im],
        ComplexF64[3.5935586686360366e-12+1.7284431242009702e-6im 2.0678134023375278e-12-3.7669677894046595e-13im;
                   2.0678134023375278e-12-3.7669677894046595e-13im 3.593558668636037e-12+1.7284431242009702e-6im],
        ComplexF64[5.141397717625919e-10+3.456862375802595e-5im 4.511959006831256e-10-2.452790870106892e-10im;
                   4.511959006831256e-10-2.452790870106892e-10im 5.141397717625921e-10+3.456862375802596e-5im];
        dims = 3
    )
    @test parameters.Z.values == expected_Z
    @test parameters.Y.values == expected_Y
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
        @test slice[1, 2] ≈ slice[2, 1]
    end

    modal=compute(
        ModalTransformationProblem(parameters),
        ModalTransformationFormulation(:Fortescue; tolerance = 1e-10)
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
                                                               for source in problem.system.designs]
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
        formulation = Formulation()) = LineParametersWorkspace(
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
        options = (
        reduce_bundle = true,
        kron_reduction = false,
        ideal_transposition = true,
        temperature_correction = false
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
        options = (
        reduce_bundle = false,
        kron_reduction = false,
        ideal_transposition = true,
        temperature_correction = false
    )
    )
    singleton_result=compute(singleton_problem, unreduced)
    @test size(singleton_result.Z) == (1, 1, 1)
    @test real(singleton_result.Z[1, 1, 1]) > 0
    @test imag(singleton_result.Y[1, 1, 1]) > 0
end

@testitem "Engine / propagation and layer-pair routes reach the solve loop" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    base=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    explicit=LineParametersProblem(
        base.system;
        temperature = base.temperature,
        earth_props = base.earth_props,
        frequencies = base.frequencies,
        Γ = ComplexF64[1.0e-5im, 2.0e-5im]
    )
    @test explicit.Γ == ComplexF64[1.0e-5im, 2.0e-5im]
    @test all(isfinite, compute(explicit, Formulation()).Z)
    @test_throws ArgumentError compute(
        explicit,
        Formulation(earth_impedance = :Pollaczek1926)
    )
    @test_throws DimensionMismatch LineParametersProblem(
        base.system;
        temperature = base.temperature,
        earth_props = base.earth_props,
        frequencies = base.frequencies,
        Γ = [0.0]
    )

    design=TestFixtures.mv_cable_design()
    connections(phase) = Dict("core"=>phase, "sheath"=>0, "jacket"=>0)
    mixed_system=build(
        LineCableSystem,
        [design, design],
        [Pose2(0.0, 1.0), Pose2(1.0, -1.0)];
        connections = [connections(1), connections(2)],
        system_id = "mixed-placement"
    )
    mixed_problem=LineParametersProblem(
        mixed_system;
        earth_props = EarthModel(100.0),
        frequencies = [50.0]
    )
    execution=computation_options(LineCableModelsCoaxial, (;))
    ordinary=Formulation(options = (ideal_transposition = false,))
    blueprints=LineCableModels.Engine.CableBlueprint{eltype(mixed_problem)}[LineCableModels.Engine.flatten(
                                                                                LineCableModelsCoaxial(), design,
                                                                                eltype(mixed_problem)
                                                                            )
                                                                            for design in mixed_problem.system.designs]
    workspace=LineParametersWorkspace(
        mixed_problem, ordinary, execution, blueprints)
    @test getfield.(workspace.invariants.earth_pairs, :layers) ==
          [(1, 1), (1, 2), (2, 2)]
    @test_throws ArgumentError compute(mixed_problem, ordinary)

    impedance=EarthImpedance.Formula(:Papadopoulos2010)
    impedance_default=EarthImpedance.routes(impedance).mutual
    impedance_route=(functor,
        pair)->pair.layers==(2, 2) ?
               impedance_default(functor, pair) :
               zero(complex(pair.separation))
    admittance=EarthAdmittance.Formula(:Papadopoulos2010)
    admittance_default=EarthAdmittance.routes(admittance).mutual
    admittance_route=(functor,
        pair)->pair.layers==(2, 2) ?
               admittance_default(functor, pair) :
               zero(complex(pair.separation))
    experiment=Formulation(
        earth_impedance = EarthImpedance.Formula(
            :Papadopoulos2010;
            self = impedance_route,
            mutual = impedance_route
        ),
        earth_admittance = EarthAdmittance.Formula(
            :Papadopoulos2010;
            self = admittance_route,
            mutual = admittance_route
        ),
        options = (ideal_transposition = false,)
    )
    result=compute(mixed_problem, experiment)
    @test all(isfinite, result.Z)
    @test all(isfinite, result.Y)
end

@testitem "Engine / frequency-dependent earth relation reaches coaxial solve" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    problem=TestFixtures.line_parameters_problem(frequencies = [1.0e6])
    static=compute(problem, Formulation())
    dispersive=compute(
        problem,
        Formulation(earth_properties = :CIGRE2019)
    )

    @test all(isfinite, dispersive.Z)
    @test all(isfinite, dispersive.Y)
    @test dispersive.Z.values != static.Z.values
    @test dispersive.Y.values != static.Y.values
end

@testitem "Engine / pair-local EHEM reaches overhead Coaxial solve" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    const EP=LineCableModels.EarthProps
    const EH=EP.EHEM

    design=TestFixtures.mv_cable_design()
    connections=Dict(
        terminal=>(terminal===:core ? 1 : 0)
    for terminal in design.terminal_order)
    system=build(
        LineCableSystem,
        design,
        Pose2(0.0, 10.0);
        connections,
        system_id = "overhead-ehem"
    )
    earth=build(EarthModel,
        (
            EP.EarthLayer(100.0, 10.0, 1.0, 5.0),
            EP.EarthLayer(500.0, 20.0, 1.0, 10.0),
            EP.EarthLayer(50.0, 5.0, 1.0)
        ))
    problem=LineParametersProblem(
        system;
        earth_props = earth,
        frequencies = [50.0]
    )
    common=(
        earth_impedance = :Pollaczek1926,
        earth_admittance = :Ametani2021,
        options = (ideal_transposition = false,)
    )
    selected=compute(
        problem,
        Formulation(;
            common...,
            equivalent_earth = EH.AfterFD(EH.Layer(-1))
        )
    )
    martins=compute(
        problem,
        Formulation(;
            common...,
            equivalent_earth = EH.AfterFD(:MartinsBritto2020)
        )
    )
    xue=compute(
        problem,
        Formulation(;
            common...,
            equivalent_earth = EH.AfterFD(:Xue2021)
        )
    )

    for result in (selected, martins, xue)
        @test all(isfinite, result.Z)
        @test all(isfinite, result.Y)
    end
    @test !isapprox(martins.Z.values, selected.Z.values; rtol = 1.0e-10)
    @test !isapprox(xue.Z.values, selected.Z.values; rtol = 1.0e-10)
end
