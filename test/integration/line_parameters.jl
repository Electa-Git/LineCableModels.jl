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

    execution=computation_options(Val(LineCableModelsCoaxial), (;))
    workspace=LineParametersWorkspace(
        LineCableModelsCoaxial(), problem, Formulation(), execution)
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
        LineCableModelsCoaxial(),
        problem,
        allocation_formulation,
        execution
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
        earth_props = Earth(rho = 100.0),
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
    execution=computation_options(Val(LineCableModelsCoaxial), (;))
    @test_throws ArgumentError LineParametersWorkspace(
        LineCableModelsCoaxial(), problem, Formulation(), execution)
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
    execution=computation_options(Val(LineCableModelsCoaxial), (;))
    workspace(problem,
        formulation = Formulation()) = LineParametersWorkspace(
        LineCableModelsCoaxial(), problem, formulation, execution)

    bare_input=workspace(problem_for(bare)).input
    @test bare_input.r_ins_in == bare_input.r_ins_ext == [0.01]
    @test bare_input.insulator_layer_ranges == [1:0]
    @test isempty(bare_input.r_ins_layer_in)

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
    @test layered_input.insulator_layer_ranges == [1:2, 3:2]
    @test layered_input.r_ins_in[2] == layered_input.r_ins_ext[2] == 0.014

    filled=build(CableDesign,
        "filled",
        Enclosure(
            :pipe,
            Group(:phase, Region(:core, Disk(0.01), conductor));
            primitive = Disk(0.03),
            fill = dielectric
        ))
    filled_input=workspace(problem_for(filled)).input
    @test filled_input.r_ins_layer_in == [0.01]
    @test filled_input.r_ins_layer_ext == [0.03]
    @test filled_input.rho_ins_layer == [dielectric.rho]

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
    execution=computation_options(Val(LineCableModelsCoaxial), (;))
    ordinary=Formulation(options = (ideal_transposition = false,))
    workspace=LineParametersWorkspace(
        LineCableModelsCoaxial(), mixed_problem, ordinary, execution)
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
    earth=EarthModel(100.0, 10.0, 1.0; thickness = 5.0)
    add!(earth, EP.EarthLayer(500.0, 20.0, 1.0, 10.0))
    add!(earth, EP.EarthLayer(50.0, 5.0, 1.0))
    problem=LineParametersProblem(
        system;
        earth_props = earth,
        frequencies = [1.0e6]
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
    @test martins.Z.values != selected.Z.values
    @test xue.Z.values != selected.Z.values
end
