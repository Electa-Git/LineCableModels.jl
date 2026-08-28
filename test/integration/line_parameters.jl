@testitem "Engine / solver / multicable reciprocity and modal transformation" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    using LinearAlgebra

    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    formulation=Formulation(;
        modal_transform = Fortescue(),
        options = (
            reduce_bundle = true,
            kron_reduction = true,
            ideal_transposition = false,
            temperature_correction = true
        )
    )
    parameters=compute(problem, formulation; options = (trace = true,))
    trace=details(parameters).trace

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

    execution=computation_options(Val(LineCableModelsEngine), (;))
    workspace=LineParametersWorkspace(
        LineCableModelsEngine(), problem, Formulation(), execution)
    LineCableModels.Engine._capture_buffers(
        Float64, workspace.normalized, Val(false))
    @test @allocated(LineCableModels.Engine._capture_buffers(
        Float64, workspace.normalized, Val(false))) <= 1024
    @test LineCableModels.Engine._capture_buffers(
        Float64, workspace.normalized, Val(false)) === nothing
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

    _, modal=Fortescue(tol = 1e-10)(parameters)
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
    execution=computation_options(Val(LineCableModelsEngine), (;))
    @test_throws ArgumentError LineParametersWorkspace(
        LineCableModelsEngine(), problem, Formulation(), execution)
    @test_throws ArgumentError compute(problem, Formulation())
    @test_throws ArgumentError CableConstants(design)
end

@testitem "Engine / native profile / explicit radial support boundary" tags=[:integration] setup=[
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
    execution=computation_options(Val(LineCableModelsEngine), (;))
    workspace(problem, formulation = Formulation()) = LineParametersWorkspace(
        LineCableModelsEngine(), problem, formulation, execution)

    bare_input=workspace(problem_for(bare)).normalized
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
    layered_input=workspace(problem_for(layered)).normalized
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
    filled_input=workspace(problem_for(filled)).normalized
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
