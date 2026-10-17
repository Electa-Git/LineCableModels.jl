@testitem "Engine / problems / a gridded keyword lifts the problem" tags=[:unit, :engine] begin
    copper = Material(:conductor, 1.7241e-8)
    xlpe = Material(:insulator, 1.97e14, 2.5)
    design = build(CableDesign, "lifted", Stack(
        Group(:core, Region(:core, Disk(0.01), copper)),
        Region(:insulation, Shell(0.005), xlpe)))
    system = build(LineCableSystem, design, Pose2(0.0, -1.0); connections = Dict(:core => 1))

    problems = LineParametersProblem(system;
        earth_props = homogeneous(rho = Grid((10.0, 100.0))), frequencies = [50.0])
    @test problems isa Gridspace{LineParametersProblem}
    @test [problem.earth_props.layers[2].rho for problem in problems] == [10.0, 100.0]

    constants_problems = CableConstantsProblem(design; temperature = Grid((20.0, 90.0)))
    @test constants_problems isa Gridspace{CableConstantsProblem}
    @test [problem.temperature for problem in constants_problems] == [20.0, 90.0]
end

@testitem "Engine / cable constants / a gridded input returns compute's result space" tags=[:unit, :parametric] begin
    copper = Material(:conductor, 1.7241e-8)
    xlpe = Material(:insulator, 1.97e14, 2.5)
    design = build(CableDesign, "lifted", Stack(
        Group(:core, Region(:core, Disk(0.01), copper)),
        Region(:insulation, Shell(0.005), xlpe),
        Group(:sheath, Region(:sheath, Shell(0.001), copper))))
    temperature = Grid((20.0, 90.0))
    frequency = Grid((50.0, 60.0))
    formulation = CableConstantsFormulation(insulation_admittance = Grid((:lossless, :lossy)))
    # The problem inputs combine by `combine`, and the two formulations multiply them.
    for (combine, count) in ((:product, 8), (:zip, 4))
        results = CableConstants(design; temperature, frequency, formulation, combine)
        @test results isa ParametricResult
        @test length(results) == count
        @test all(value -> value isa CableConstants, results)
        problems = CableConstantsProblem(design; temperature, frequency, combine)
        @test collect(results) == collect(compute(problems, formulation))
    end
end

@testitem "Engine / problems / a one-point Gridspace gives the scalar problem" tags=[:unit, :engine] begin
    copper = Material(:conductor, 1.7241f-8)
    xlpe = Material(:insulator, 1.97f14, 2.5f0)
    design = build(CableDesign, "float32", Stack(
        Group(:core, Region(:core, Disk(0.01f0), copper); at = Pose2(0.0f0, 0.0f0, 0.0f0)),
        Region(:insulation, Shell(0.005f0), xlpe)))
    scalar = CableConstantsProblem(design)
    gridded = only(collect(CableConstantsProblem(Grid((design,)))))
    @test typeof(gridded) === typeof(scalar)
    @test (gridded.temperature, gridded.frequency) == (scalar.temperature, scalar.frequency)
    @test gridded.design === scalar.design

    system = build(LineCableSystem, design, Pose2(0.0f0, -1.0f0);
        connections = Dict(:core => 1), line_length = 1.0f0)
    earth = homogeneous(rho = 100.0)
    scalar = LineParametersProblem(system; earth_props = earth)
    gridded = only(collect(LineParametersProblem(Grid((system,)); earth_props = earth)))
    @test typeof(gridded) === typeof(scalar)
    @test (gridded.temperature, gridded.frequencies) == (scalar.temperature, scalar.frequencies)
    @test gridded.earth_props == scalar.earth_props
    @test gridded.system.positions == scalar.system.positions
end

@testitem "Engine / problems / mixed Float32 and Float64 inputs compute" tags=[:unit, :engine] begin
    copper = Material(:conductor, 1.7241f-8)
    xlpe = Material(:insulator, 1.97f14, 2.5f0)
    # Float32 inputs, with a Float32 pose and length, meet one Float64 input at a time.
    design = build(CableDesign, "float32", Stack(
        Group(:core, Region(:core, Disk(0.01f0), copper); at = Pose2(0.0f0, 0.0f0, 0.0f0)),
        Region(:insulation, Shell(0.005f0), xlpe)))
    system = build(LineCableSystem, design, Pose2(0.0f0, -1.0f0);
        connections = Dict(:core => 1), line_length = 1.0f0)
    earth32 = homogeneous(rho = 100.0f0)
    for problem in (
            LineParametersProblem(system; earth_props = homogeneous(rho = 100.0)),
            LineParametersProblem(system; earth_props = earth32, temperature = 20.0),
            LineParametersProblem(system; earth_props = earth32, frequencies = [50.0]))
        @test @inferred(compute(problem, Formulation())) isa LineParameters
    end
    constants = CableConstantsProblem(design; temperature = 20.0)
    @test @inferred(compute(constants, CableConstantsFormulation())) isa CableConstants
end

@testitem "Engine / problems / a Grid inside an options or connections named tuple lifts" tags=[:unit, :engine] begin
    copper = Material(:conductor, 1.7241e-8)
    xlpe = Material(:insulator, 1.97e14, 2.5)
    design = build(CableDesign, "lifted", Stack(
        Group(:core, Region(:core, Disk(0.01), copper)),
        Region(:insulation, Shell(0.005), xlpe)))

    # A named tuple is varied at its depth, like the whole-options form.
    nested = LineParametersFormulation(; options = (kron_reduction = Grid((true, false)),))
    whole = LineParametersFormulation(;
        options = Grid(((kron_reduction = true,), (kron_reduction = false,))))
    @test nested isa Gridspace{LineParametersFormulation}
    @test [f.options for f in nested] == [f.options for f in whole]
    @test LineParametersFormulation(; options = (kron_reduction = false,)) isa
          LineParametersFormulation

    systems = build(LineCableSystem, design, Pose2(0.0, -1.0);
        connections = (core = Grid((1, 2)),))
    @test systems isa Gridspace{LineCableSystem}
    @test length(systems) == 2

    # The positional values of a problem are the direct parameters of one call, so they
    # zip together at any position.
    earth = homogeneous(rho = 100.0)
    problems = LineParametersProblem(design, Grid((Pose2(0.0, -1.0), Pose2(0.0, -2.0)));
        connections = Dict(:core => 1), temperature = Grid((20.0, 90.0)),
        earth_props = earth, frequencies = [50.0], combine = :zip)
    @test problems isa Gridspace{LineParametersProblem}
    @test [(problem.system.positions[1].y, problem.temperature) for problem in problems] ==
          [(-1.0, 20.0), (-2.0, 90.0)]
end
