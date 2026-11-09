using LinearAlgebra: norm

# A touching trefoil of single-core cables buried at 1 m, of core radius, insulation
# thickness and earth resistivity p, at 50 Hz.
function trefoil_problem(p)
    rc, t, ρ = p
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    design = build(CableDesign, "single_core",
        Terminal(:core, Region(:core_metal, Circle(rc), copper)),
        Region(:core_insulation, Offset(t), polymer),
        Terminal(:sheath, Region(:sheath_copper, Annulus(rc + t, rc + t + 1.0e-3), copper)),
        Region(:sheath_jacket, Offset(3.0e-3), polymer))
    places = trefoil(design; spacing = 2r_ex(design.outer), at = Placement(0.0, -1.0))
    system = build(CableSystem, [design, design, design], places;
        connections = [
            (core = 1, sheath = 4), (core = 2, sheath = 5), (core = 3, sheath = 6)],
        id = "trefoil", line_length = 1000.0)
    return LineParametersProblem(system, homogeneous(rho = ρ); frequencies = [50.0])
end

# An overhead line of three phases and an earth wire, of phase radius, earth-wire radius
# and earth resistivity p, at 50 Hz.
function overhead_problem(p)
    rw, re, ρ = p
    aluminium = Material(:conductor, 2.8264e-8, 1, 1, 20, 0.00403)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    phase = build(CableDesign, "phase",
        Terminal(:core, Region(:core_metal, Circle(rw), aluminium)))
    shield = build(CableDesign, "earth_wire",
        Terminal(:core, Region(:core_metal, Circle(re), steel)))
    places = [Placement(-4.0, 12.0), Placement(0.0, 12.0), Placement(4.0, 12.0),
        Placement(0.0, 17.0)]
    system = build(CableSystem, [phase, phase, phase, shield], places;
        connections = [(core = 1,), (core = 2,), (core = 3,), (core = 0,)],
        id = "line", line_length = 1000.0)
    return LineParametersProblem(system, homogeneous(rho = ρ); frequencies = [50.0])
end

# The real and imaginary parts of Z, then of Y.
function line_parameters(problem, formulation, p)
    result = compute(problem(p), formulation)
    return vcat(vec(real(result.Z)), vec(imag(result.Z)), vec(real(result.Y)),
        vec(imag(result.Y)))
end

@testset "ForwardDiff through compute" begin
    trefoil_inputs = [10.0e-3, 8.0e-3, 100.0]
    @testset "$label" for (label, problem, p0, formulation) in (
        ("trefoil, default", trefoil_problem, trefoil_inputs, CoaxialFormulation()),
        ("trefoil, :pollaczek1926", trefoil_problem, trefoil_inputs,
            CoaxialFormulation(earth_impedance = :pollaczek1926)),
        ("trefoil, :papadopoulos2010", trefoil_problem, trefoil_inputs,
            CoaxialFormulation(earth_impedance = :papadopoulos2010,
                earth_admittance = :papadopoulos2010)),
        ("trefoil, :saad1996", trefoil_problem, trefoil_inputs,
            CoaxialFormulation(earth_impedance = :saad1996)),
        ("trefoil, :xue2018", trefoil_problem, trefoil_inputs,
            CoaxialFormulation(earth_impedance = :xue2018, earth_admittance = :xue2018)),
        ("line, default", overhead_problem, [12.0e-3, 5.0e-3, 100.0], CoaxialFormulation()))
        J = ForwardDiff.jacobian(p -> line_parameters(problem, formulation, p), p0)
        @test !any(isnan, J)
        Z = 1:(size(J, 1) ÷ 2)
        Y = (size(J, 1) ÷ 2 + 1):size(J, 1)
        # The loop compares Y for the two radii alone. Its derivative with respect to the
        # earth resistivity is below the error of the central difference.
        for (i, rows) in ((1, Z), (2, Z), (3, Z), (1, Y), (2, Y))
            h = 1.0e-6 * p0[i]
            up = [k == i ? p0[k] + h : p0[k] for k in 1:3]
            down = [k == i ? p0[k] - h : p0[k] for k in 1:3]
            central = (line_parameters(problem, formulation, up) -
                       line_parameters(problem, formulation, down)) / 2h
            @test norm(J[rows, i] - central[rows]) <= 1.0e-6 * norm(central[rows])
        end
    end
end

@testset "the check reads the nominal values of a Dual problem" begin
    Engine = LineCableModels.Engine
    problem = trefoil_problem(ForwardDiff.Dual.([10.0e-3, 8.0e-3, 100.0], 1.0))
    formulation = CoaxialFormulation()
    input = Engine.CoaxialInput(problem, formulation.options)
    given = @inferred Engine.domain(problem, formulation.options)
    tubes = @inferred Engine.domain(input)
    presented = merge(given, tubes)
    @test all(interval -> eltype(interval) === Float64, presented)
    @test (@inferred LineCableModels.validate(presented, formulation)) === presented
    JET.@test_opt target_modules=(LineCableModels,) Engine.domain(problem, formulation.options)
    JET.@test_opt target_modules=(LineCableModels,) Engine.domain(input)
    JET.@test_opt target_modules=(LineCableModels,) LineCableModels.validate(presented, formulation)
end
