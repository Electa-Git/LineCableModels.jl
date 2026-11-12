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

# The trefoil over a layered earth whose upper resistivity is p, at 1 MHz, where the
# equivalent resistivity of :xue2021 is negative.
function layered_problem(p)
    ρ = only(p)
    T = typeof(ρ)
    earth = EarthModel((EarthLayer(T(Inf), 1.0, 1.0), EarthLayer(ρ, 4.0, 1.0, 1.0),
        EarthLayer(T(2000), 30.0, 1.0, 20.0), EarthLayer(T(10), 60.0, 1.0, 3.0),
        EarthLayer(T(300), 8.0, 1.0)))
    return LineParametersProblem(trefoil_problem([10.0e-3, 8.0e-3, 100.0]).system, earth;
        frequencies = [1.0e6])
end

@testset "a negative equivalent resistivity of a Dual problem warns on its nominal value" begin
    formulation = CoaxialFormulation(equivalent_earth = :xue2021)
    J = @test_logs (:warn, r"the equivalent earth of :xue2021") (:warn,
        r"the series impedance Z is not passive at 1 of 1") ForwardDiff.jacobian(
        p -> line_parameters(layered_problem, formulation, p), [50.0])
    @test all(isfinite, J)
end

# The order counts the terms of a sum and keeps its integer type.
@testset "the check reads the nominal values of a Dual problem" begin
    Engine = LineCableModels.Engine
    problem = trefoil_problem(ForwardDiff.Dual.([10.0e-3, 8.0e-3, 100.0], 1.0))
    formulation = CoaxialFormulation()
    input = Engine.CoaxialInput(problem, formulation.options)
    given = @inferred Engine.domain(problem, formulation.options)
    tubes = @inferred Engine.domain(input)
    presented = merge(given, tubes)
    @test all(eltype(presented[key]) === Float64
    for key in keys(presented) if key !== :order)
    @test eltype(presented.order) === Int
    @test (@inferred LineCableModels.validate(presented, formulation)) === presented
    JET.@test_opt target_modules=(LineCableModels,) Engine.domain(problem, formulation.options)
    JET.@test_opt target_modules=(LineCableModels,) Engine.domain(input)
    JET.@test_opt target_modules=(LineCableModels,) LineCableModels.validate(presented, formulation)
end

# The cores of a pipe-type cable, each 3 mm off the origin of its own turned frame. Each member
# lies on the centroid of its core, which the core radius does not move, and its tube ends at
# that radius.
@testset "a member off the origin of its frame has the derivatives of its core radius" begin
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.5)
    steel = Material(:conductor, 7.0e-8, 1.0, 10.0, 20.0, 0.0045)
    function blueprint(r)
        members = [Placement(30.0e-3 * cos(θ),
                       30.0e-3 * sin(θ),
                       θ) => Stack(
                       Terminal(n, Region(n, Circle(r, Placement(3.0e-3, 0.0)), copper)),
                       Region(Symbol(n, :_insulation), Offset(4.0e-3), polymer))
                   for (n, θ) in zip((:a, :b, :c), (0.0, 2π / 3, 4π / 3))]
        design = build(CableDesign, "pipe",
            Enclosure(Region(:fill, Circle(62.0e-3), polymer), Assembly(Tuple(members))),
            Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
            Region(:coating, Offset(5.0e-3), polymer))
        return LineCableModels.flatten(design)
    end
    axis_partials = ForwardDiff.derivative(
        r -> reduce(vcat, [collect(axis) for axis in blueprint(r).positions[2:end]]), 10.0e-3)
    @test all(iszero, axis_partials)
    tube_partials = ForwardDiff.derivative(
        r -> [tube.r_ex for tube in blueprint(r).conductors[1:3]], 10.0e-3)
    @test all(==(1), tube_partials)
end

# A core of seven wires in a steel pipe, at the three positions of a three-core cable. At one
# of them the centre wire lies on the computed axis to the bit. Its distance from the axis then
# takes the derivative zero.
@testset "a member of seven wires has the derivatives of its wire radius" begin
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.5)
    steel = Material(:conductor, 7.0e-8, 1.0, 10.0, 20.0, 0.0045)
    function stranded(r, θ)
        member = Stack(
            Terminal(:a,
                Enclosure(Region(:a_fill, Circle(3r), polymer),
                    Stack(Region(:a_centre, Circle(r), copper),
                        Repeat(CircularPattern(6; r = 2r), Region(:a_wires, Circle(r), copper))))),
            Region(:a_insulation, Offset(3.0e-3), polymer))
        design = build(CableDesign, "stranded",
            Enclosure(Region(:fill, Circle(45.0e-3), polymer),
                Assembly((Placement(22.0e-3 * cos(θ), 22.0e-3 * sin(θ)) => member,))),
            Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
            Region(:coating, Offset(5.0e-3), polymer))
        system = build(CableSystem, [design], [Placement(0.0, -1.0)];
            connections = [(a = 1, pipe = 2)], id = "stranded", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
            frequencies = [50.0])
        return compute(problem, CoaxialFormulation())
    end
    for θ in (0.0, 2π / 3, 4π / 3)
        result = stranded(ForwardDiff.Dual(3.0e-3, 1.0), θ)
        @test all(isfinite(ForwardDiff.partials(part(q))[1])
        for Q in (result.Z, result.Y) for q in Q for part in (real, imag))
        # The reference is a central difference of two Float64 runs 1e-8 m apart. On this build
        # the largest relative gap of the three positions is 3.4e-11.
        up, down = stranded(3.0e-3 + 1.0e-8, θ), stranded(3.0e-3 - 1.0e-8, θ)
        slope = (real(up.Z[1, 1, 1]) - real(down.Z[1, 1, 1])) / 2.0e-8
        partial = ForwardDiff.partials(real(result.Z[1, 1, 1]))[1]
        @test isapprox(partial, slope; rtol = 1.0e-10)
    end
end
