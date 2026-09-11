@testitem "Engine / pipe selections / backend applicability before lowering" tags=[:unit] begin
    using LineCableModels
    const E = LineCableModels.Engine
    copper = Material(:conductor, 1.72e-8)
    dielectric = Material(:insulator, 1e14, 2.3)
    air = Material(:insulator, Inf, 1.0)
    first_core = terminal(:a, solid(copper, Disk(0.005)), insulation(dielectric; t = 0.002))
    second_core = terminal(:b, solid(copper, Disk(0.005)), insulation(dielectric; t = 0.002))
    wall = terminal(:pipe, sheath(copper; t = 0.001), insulation(dielectric; t = 0.002))
    design = build(CableDesign,
        "pipe-preflight",
        pipe(
            at(first_core, -0.01, 0), at(second_core, 0.01, 0);
            shape = Disk(0.025), fill = air, wall))
    system = build(LineCableSystem, design, Pose2(0.0, -1.0);
        connections = Dict(:a=>1, :b=>2, :pipe=>0))
    problem = LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies = [50.0])
    selected = @inferred E.PipeImpedance.Formula(Val(:default))
    expected = "Pipe-type cable formulation is not yet implemented for the coaxial backend. No default formulation is available."
    for execute in (() -> compute(problem), () -> compute(CableConstantsProblem(design)))
        failure = try
            execute()
            nothing
        catch exception
            exception
        end
        @test failure isa ArgumentError
        @test sprint(showerror, failure) == "ArgumentError: $expected"
    end
    # FEM enclosure support is exercised by its geometry and numerical suites.
    @test !haskey(LineCableModelsFEM().methods, :pipe_impedance)
    @test_throws MethodError LineCableModelsFEM(pipe_impedance=formula(:default))
    @test_throws ArgumentError E.PipeImpedance.Formula(:UnimplementedPipe2026)
    @test_throws ArgumentError E.PipeImpedance.Formula(:default; parameters = (invented = true,))
    @test_throws ArgumentError Formulation(pipe_impedance = formula(:default; order = :before))

    concentric = build(
        CableDesign, "concentric-enclosure", pipe(first_core;
            shape = Disk(0.01), fill = air, wall))
    @test E.Formulation(LineCableModelsCoaxial(), selected, concentric) === nothing
    @test all(isfinite, compute(CableConstantsProblem(concentric)).R)

    # A dielectric duct is not a conducting pipe-return formulation.
    ducted = build(CableDesign,
        "dielectric-duct",
        pipe(
            at(first_core, -0.01, 0), at(second_core, 0.01, 0);
            shape = Disk(0.025), fill = air, wall = insulation(dielectric; t = 0.002)))
    @test E.Formulation(LineCableModelsCoaxial(), selected, ducted) === nothing
    for constructor in (Formulation, CableConstantsFormulation)
        selections = constructor(pipe_impedance = Grid((:default, formula(:default))))
        @test length(selections) == 2
        @test all(value -> value.methods.pipe_impedance isa E.PipeImpedanceFormulation, selections)
    end
end
