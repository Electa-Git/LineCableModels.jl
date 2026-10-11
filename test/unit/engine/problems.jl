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

@testitem "Engine / problems / integer inputs keep the precision of the system" tags=[:unit, :engine] begin
    copper = Material(:conductor, 1.7241f-8)
    xlpe = Material(:insulator, 1.97f14, 2.5f0)
    # The default group pose and system length are Float64, so they are given here.
    design = build(CableDesign, "float32", Stack(
        Group(:core, Region(:core, Disk(0.01f0), copper); at = Pose2(0.0f0, 0.0f0, 0.0f0)),
        Region(:insulation, Shell(0.005f0), xlpe)))
    @test eltype(design) === Float32
    # A scalar input and a one-point Gridspace give problems of one scalar type.
    @test eltype(CableConstantsProblem(design)) === Float32
    @test eltype(only(collect(CableConstantsProblem(Grid((design,)))))) === Float32
    @test eltype(CableConstantsProblem(design; temperature = 90, frequency = 60)) === Float32
    system = build(LineCableSystem, design, Pose2(0.0f0, -1.0f0);
        connections = Dict(:core => 1), line_length = 1.0f0)
    earth = homogeneous(rho = 100.0f0)
    @test eltype(LineParametersProblem(system; earth_props = earth)) === Float32
    @test eltype(only(collect(LineParametersProblem(Grid((system,)); earth_props = earth)))) ===
          Float32
end
