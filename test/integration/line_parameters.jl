@testitem "Engine / solver / multicable reciprocity and modal transformation" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures,
    TestNumerics
] begin
    using LinearAlgebra

    problem=TestFixtures.line_parameters_problem(frequencies = [50.0, 500.0])
    formulation=Formulation(
        Val(:analytical);
        modal_transform = Fortescue(),
        options = (
            reduce_bundle = true,
            kron_reduction = true,
            ideal_transposition = false,
            temperature_correction = true,
            output = :trace
        )
    )
    trace=compute(problem, formulation)
    parameters=trace.result

    @test domain(parameters) === ModalDomain
    @test size(parameters.Z) == (3, 3, 2)
    @test all(isfinite, parameters.Z)
    @test all(isfinite, parameters.Y)
    @test trace isa LineParametersTrace
    @test eltype(trace) === Float64
    @test sprint(show, trace) ==
          "LineParametersTrace(2 frequencies, 9 primitive conductors)"
    @test fieldtype(typeof(trace), :result) <: LineParameters
    @test observe(trace, Z) == observe(parameters, Z)
    @test observe(trace, Y) == observe(parameters, Y)
    @test trace.frequencies == problem.frequencies
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

    input=LineCableModels.Engine.AnalyticalInput(problem)
    LineCableModels.Engine._trace_buffers(Float64, input, Val(:parameters))
    @test @allocated(LineCableModels.Engine._trace_buffers(
        Float64, input, Val(:parameters))) == 0
    @test LineCableModels.Engine._trace_buffers(
        Float64, input, Val(:parameters)) === nothing
end

@testitem "Engine / formulation boundary / physical geometry precedes backend support" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    conductor=Material(kind = :conductor, rho = 1.7241e-8)
    design=CableDesign(
        Group(:phase, Region(:elliptical_core, Ellipse(0.01, 0.006), conductor));
        cable_id = "elliptical"
    )
    @test design.effective === nothing
    @test design.geometry.regions[1].shape.shape isa EllipseShape

    system=LineCableSystem(
        design;
        position = Pose2(0.0, -1.0),
        connections = (phase = 1,)
    )
    problem=LineParametersProblem(
        system;
        earth_props = EarthModel(100.0),
        frequencies = [50.0]
    )
    @test problem.system === system
    @test_throws ArgumentError Engine.AnalyticalInput(problem)
    @test_throws ArgumentError compute(problem, Formulation())

    constants_problem=CableConstantsProblem(design)
    @test constants_problem.design === design
    @test_throws ArgumentError compute(constants_problem, Formulation())
end

@testitem "Engine / analytical profile / explicit radial support boundary" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    conductor=Material(kind = :conductor, rho = 1.7241e-8)
    dielectric=Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    earth=EarthModel(100.0)

    function problem_for(design)
        system=LineCableSystem(
            design;
            position = Pose2(0.0, -1.0),
            connections = Dict(first(design.terminal_order)=>1)
        )
        return LineParametersProblem(
            system;
            earth_props = earth,
            frequencies = [50.0]
        )
    end

    bare=CableDesign(Group(:phase, Region(:bare, Disk(0.01), conductor)))
    bare_dielectric=only(bare.effective).dielectric
    @test bare_dielectric.r_in == bare_dielectric.r_ex == 0.01
    @test isempty(bare_dielectric.layers)
    bare_input=Engine.AnalyticalInput(problem_for(bare))
    @test bare_input.insulator_layer_ranges == [1:0]
    @test isempty(bare_input.r_ins_layer_in)

    layered=CableDesign(Stack(
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
    @test length(layered.effective) == 2
    @test length(layered.effective[1].dielectric.layers) == 2
    @test isempty(layered.effective[2].dielectric.layers)
    @test layered.effective[2].dielectric.r_in ==
          layered.effective[2].dielectric.r_ex == 0.014

    filled=CableDesign(Enclosure(
        :pipe,
        Group(:phase, Region(:core, Disk(0.01), conductor));
        shape = Disk(0.03),
        fill = dielectric
    ))
    fill_layer=only(only(filled.effective).dielectric.layers)
    @test fill_layer.r_in == 0.01
    @test fill_layer.r_ex == 0.03
    @test fill_layer.material === dielectric

    reappearing=CableDesign(Stack(
        Group(:a, Region(:a_inner, Disk(0.01), conductor)),
        Region(:a_gap, Shell(0.001), dielectric),
        Group(:b, Region(:b_metal, Annulus(0.011, 0.012), conductor)),
        Region(:b_gap, Shell(0.001), dielectric),
        Group(:a, Region(:a_outer, Annulus(0.013, 0.014), conductor))
    ))
    @test reappearing.effective === nothing
    @test_throws ArgumentError Engine.AnalyticalInput(problem_for(reappearing))

    conductor_after_dielectric=CableDesign(Stack(
        Group(:phase, Region(:inner, Disk(0.01), conductor)),
        Region(:gap, Shell(0.001), dielectric),
        Group(:phase, Region(:outer, Annulus(0.011, 0.012), conductor))
    ))
    @test conductor_after_dielectric.effective === nothing
    @test_throws ArgumentError Engine.AnalyticalInput(
        problem_for(conductor_after_dielectric)
    )

    undeclared_gap=CableDesign(Stack(
        Group(:inner, Region(:inner, Disk(0.01), conductor)),
        Group(:outer, Region(:outer, Annulus(0.012, 0.013), conductor))
    ))
    @test undeclared_gap.effective === nothing
    @test_throws ArgumentError Engine.AnalyticalInput(problem_for(undeclared_gap))

    @test_throws DomainError CableDesign(Stack(
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
    duplicate_system=LineCableSystem(
        design;
        position = Pose2(0.0, -1.0, 0.0),
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
    bundle_only=Formulation(
        Val(:analytical);
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

    singleton_design=CableDesign(
        Stack(deepcopy(design.root.items[1:10]));
        cable_id = "single-component",
        reference_frequency = design.reference_frequency
    )
    singleton_system=LineCableSystem(
        singleton_design;
        position = Pose2(0.0, -1.0, 0.0),
        connections = Dict(:core=>1),
        system_id = "singleton",
        line_length = 100.0
    )
    singleton_problem=LineParametersProblem(
        singleton_system;
        earth_props = earth,
        frequencies
    )
    unreduced=Formulation(
        Val(:analytical);
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
