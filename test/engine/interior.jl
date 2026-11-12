# The interiors of the coaxial formulation: the members in the bore of a conductor, at the
# harmonic orders that the option order keeps. Each bound comes from a measurement on Julia
# 1.12.7 with one thread.

# A design of members in the fill of a steel pipe, with the pipe's coating.
function piped_design(id, members, fill, steel, coating)
    return build(CableDesign, id,
        Enclosure(Region(:fill, Circle(30.0e-3), fill), Assembly(Tuple(members))),
        Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
        Region(:coating, Offset(2.0e-3), coating))
end

# The result of a system on a homogeneous earth at four frequencies, with its matrices inside
# the cables.
function interior_result(designs, placements, connections; kw...)
    system = build(CableSystem, designs, placements; connections, id = "interior",
        line_length = 1.0)
    problem = LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
        frequencies = [50.0, 1.0e3, 1.0e5, 1.0e6])
    return compute(problem, CoaxialFormulation(; kw...); trace = true)
end

# The matrix of the loops of the conductors `rows` with the conductor `k` as their return.
looped(M, k, rows) = [M[i, j] - M[i, k] - M[k, j] + M[k, k] for i in rows, j in rows]

@testset "Engine interiors" begin
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    air = material(kind = :insulator, rho = Inf)
    paper = Material(:insulator, 1.0e14, 3.5)
    oil = Material(:insulator, 1.0e14, 2.2)
    single = (a = 1, pipe = 2)

    # A sub-cable declared before a core gives the same pipe matrices as its wires declared directly.
    # Measured: 0 for Z and 1.66e-17 for Y at order 0.
    # Measured: 1.44e-20 for Z and 3.13e-17 for Y at order 20.
    @testset "a sub-cable before a core gives the pipe matrices of its wires declared directly" begin
        wires = [Stack(Terminal(name, Region(name, Circle(1.5e-3), copper)),
                     Region(Symbol(name, :_insulation), Offset(0.5e-3), polymer))
                 for name in (:w1, :w2)]
        cable = Stack(
            Enclosure(Region(:t_fill, Circle(6.0e-3), air),
                Assembly((Placement(3.0e-3, 0.0) => wires[1],
                    Placement(-3.0e-3, 0.0) => wires[2]))),
            Region(:t_jacket, Offset(1.0e-3), air))
        core = Stack(Terminal(:core, Region(:core, Circle(6.0e-3), copper)),
            Region(:core_insulation, Offset(2.0e-3), polymer))
        members = ([Placement(-12.0e-3, 0.0) => cable,
                Placement(14.0e-3, 0.0) => core],
            [Placement(-9.0e-3, 0.0) => wires[1],
                Placement(-15.0e-3, 0.0) => wires[2],
                Placement(14.0e-3, 0.0) => core])
        for order in (0, 20)
            results = [interior_result(
                           [piped_design("pipe", rows, air, steel, polymer)],
                           [Placement(0, -1)], [(w1 = 1, w2 = 2, core = 3, pipe = 4)]; order)
                       for rows in members]
            @test isapprox(results[1].Z, results[2].Z; rtol = 1.0e-15)
            @test isapprox(results[1].Y, results[2].Y; rtol = 1.0e-15)
        end
    end

    # In a jacket without metal, the sub-cable's wires and the core compute as three cables.
    # Measured: 1.22e-20 for Z and 0 for Y at order 0.
    @testset "a sub-cable before a core in a jacket gives the matrices of its wires declared directly" begin
        wires = [Stack(Terminal(name, Region(name, Circle(1.5e-3), copper)),
                     Region(Symbol(name, :_insulation), Offset(0.5e-3), polymer))
                 for name in (:w1, :w2)]
        cable = Stack(
            Enclosure(Region(:t_fill, Circle(6.0e-3), air),
                Assembly((Placement(3.0e-3, 0.0) => wires[1],
                    Placement(-3.0e-3, 0.0) => wires[2]))),
            Region(:t_jacket, Offset(1.0e-3), air))
        core = Stack(Terminal(:core, Region(:core, Circle(6.0e-3), copper)),
            Region(:core_insulation, Offset(2.0e-3), polymer))
        members = ([Placement(-12.0e-3, 0.0) => cable,
                Placement(14.0e-3, 0.0) => core],
            [Placement(-9.0e-3, 0.0) => wires[1],
                Placement(-15.0e-3, 0.0) => wires[2],
                Placement(14.0e-3, 0.0) => core])
        designs = [build(CableDesign, "jacketed",
                       Stack(
                           Enclosure(Region(:fill, Circle(30.0e-3), air), Assembly(Tuple(rows))),
                           Region(:jacket, Offset(2.0e-3), polymer))) for rows in members]
        results = [interior_result([design], [Placement(0, -1)],
                       [(w1 = 1, w2 = 2, core = 3)]) for design in designs]
        @test size(results[1].trace.Zg) == (3, 3, 4)
        @test isapprox(results[1].Z, results[2].Z; rtol = 1.0e-15)
        @test isapprox(results[1].Y, results[2].Y; rtol = 1.0e-15)
    end

    # At order 0 the coaxial formulation refuses a screen with a member and no region
    # around it.
    # Measured: the refusal names cable:1:screen and its missing region.
    @testset "a conductor with members and no region around them is refused at order 0" begin
        core = Stack(Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            Region(:insulation, Offset(2.0e-3), polymer))
        tight = build(CableDesign,
            "tight",
            Stack(Assembly((Placement(0.0, 0.0) => core,)),
                Terminal(:screen, Region(:screen, Offset(1.0e-3), copper)), Region(
                    :jacket, Offset(2.0e-3), polymer)))
        @test_throws ArgumentError("no region lies around the members of cable:1:screen: the coaxial formulation takes the medium of the gap between the members and that conductor from it") interior_result(
            [tight], [Placement(0, -1)], [(core = 1, screen = 2)])
    end

    # A core of radius `a` lies at `b` from the axis of a bore of radius `c`, in one dielectric.
    # Its coefficient is `acosh((a^2 + c^2 - b^2) / (2a * c)) / (2π * ε)`, and its insulation
    # has the filler's medium. Measured: 4.4e-16 at 12 mm and order 20, 6.0e-15 at 17 mm and
    # order 40.
    @testset "the potential of a member off its bore's axis is exact" begin
        core = Stack(Terminal(:a, Region(:a, Circle(10.0e-3), copper)),
            Region(:a_insulation, Offset(1.0e-3), polymer))
        axis = interior_result(
            [piped_design("axis", [Placement(0.0, 0.0) => core], polymer, steel, polymer)],
            [Placement(0, -1)], [single])
        P0 = looped(axis.trace.Pin[:, :, 1], 2, 1:1)[1, 1]
        for (b, order, bound) in ((12.0e-3, 20, 1.0e-15), (17.0e-3, 40, 1.0e-14))
            design = piped_design(
                "offset", [Placement(b, 0.0) => core], polymer, steel, polymer)
            result = interior_result([design], [Placement(0, -1)], [single]; order)
            P = looped(result.trace.Pin[:, :, 1], 2, 1:1)[1, 1]
            @test isapprox(P / P0, acosh((1.0e-3 - b^2) / 6.0e-4) / log(3); rtol = bound)
        end
    end

    # On its bore's axis a member returns no field above order 0, and its cable is the same
    # cable declared concentric. Measured with a coaxial member: 1.3e-17 for Z, 4.5e-16 for Y.
    @testset "a member on its bore's axis is the concentric cable" begin
        concentric = build(CableDesign, "axis",
            Terminal(:a, Region(:a, Circle(10.0e-3), copper)),
            Region(:a_insulation, Offset(2.0e-3), polymer), Region(:fill, Offset(18.0e-3), air),
            Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
            Region(:coating, Offset(2.0e-3), polymer))
        expected = interior_result([concentric], [Placement(0, -1)], [single])
        core = Stack(Terminal(:a, Region(:a, Circle(10.0e-3), copper)),
            Region(:a_insulation, Offset(2.0e-3), polymer))
        design = piped_design("axis", [Placement(0.0, 0.0) => core], air, steel, polymer)
        for order in (0, 1, 20)
            result = interior_result([design], [Placement(0, -1)], [single]; order)
            @test isapprox(result.Z, expected.Z; rtol = 1.0e-14)
            @test isapprox(result.Y, expected.Y; rtol = 1.0e-14)
        end
    end

    # A jacket of the filler's medium leaves the field of every order unchanged. Measured:
    # 3.6e-18 for Z and 3.8e-17 for Y at order 20.
    @testset "a layer of the filler's medium is transparent" begin
        results = map((true, false)) do jacketed
            members = [Placement(x,
                           0.0) => Stack(
                           Terminal(name, Region(name, Circle(8.0e-3), copper)),
                           Region(Symbol(name, :_insulation), Offset(3.0e-3), paper),
                           (jacketed ?
                            (Region(Symbol(name, :_jacket), Offset(2.0e-3), oil),) :
                            ())...)
                       for (x, name) in ((14.0e-3, :a), (-14.0e-3, :b))]
            interior_result([piped_design("jacket", members, oil, steel, polymer)],
                [Placement(0, -1)], [(a = 1, b = 2, pipe = 3)]; order = 20)
        end
        @test isapprox(results[1].Z, results[2].Z; rtol = 1.0e-15)
        @test isapprox(results[1].Y, results[2].Y; rtol = 1.0e-15)
    end

    # A jacket of relative permittivity 1e12 shields like a metal of its outer radius. Its own
    # coefficient differs, by 1.21e-12 of the loop's, measured.
    @testset "a jacket of a perfect dielectric is the equipotential at its outer radius" begin
        perfect = Material(:insulator, 1.0e14, 1.0e12)
        layered = [Placement(x,
                       0.0) => Stack(
                       Terminal(name, Region(name, Circle(8.0e-3), copper)),
                       Region(Symbol(name, :_jacket), Offset(4.0e-3), perfect))
                   for (x, name) in ((13.0e-3, :a), (-13.0e-3, :b))]
        solid = [Placement(x,
                     0.0) => Stack(
                     Terminal(name, Region(name, Circle(12.0e-3), copper)),
                     Region(Symbol(name, :_film), Offset(1.0e-6), oil))
                 for (x, name) in ((13.0e-3, :a), (-13.0e-3, :b))]
        connections = (a = 1, b = 2, pipe = 3)
        results = [interior_result([piped_design("limit", members, oil, steel, polymer)],
                       [Placement(0, -1)], [connections]; order = 20)
                   for members in (layered, solid)]
        for k in 1:4
            @test isapprox(looped(results[1].trace.Pin[:, :, k], 3, 1:2),
                looped(results[2].trace.Pin[:, :, k], 3, 1:2); rtol = 2.0e-12)
        end
    end

    # At order 0 each member's loop with the pipe takes the logarithms of the axes. Measured:
    # 6.1e-16 over the 101 frequencies of the reference case.
    @testset "at order 0 the members take the logarithms of their axes" begin
        problem = Dict(REFERENCES)["three_cores_steel_pipe"]
        schelkunoff = LineCableModels.Engine.InternalImpedance.Formula(Val(:schelkunoff1934))
        result = compute(problem, CoaxialFormulation(); trace = true)
        blueprint = LineCableModels.flatten(only(problem.system.designs))
        tubes = blueprint.conductors
        centers = [complex(position...) for position in blueprint.positions[2:4]]
        c = tubes[4].r_in
        for (k, f) in enumerate(problem.frequencies)
            s = im * (2π * f)
            bore = schelkunoff(c, tubes[4].r_ex, tubes[4].material.rho,
                tubes[4].material.mu_r, s).inner
            expected = [bore + (i == j ?
                         schelkunoff(0.0, tubes[i].r_ex, tubes[i].material.rho,
                            tubes[i].material.mu_r, s).outer +
                         s * LineCableModels.μ₀ / 2π * log(c / tubes[i].r_ex) :
                         s * LineCableModels.μ₀ / 2π *
                         log(c / abs(centers[i] - centers[j])))
                        for i in 1:3, j in 1:3]
            @test isapprox(looped(result.trace.Zin[:, :, k], 4, 1:3), expected;
                rtol = 1.0e-15)
        end
    end

    # The equivalent tubes of the sectors overlap. Their expansions diverge above order 0,
    # and order 0 computes them.
    @testset "the sector cable computes at order 0 and is refused above it" begin
        problem = Dict(REFERENCES)["three_insulated_sectors"]
        @test size(compute(problem, CoaxialFormulation()).Z) == (4, 4, 101)
        @test_throws "the expansions about their axes diverge above order 0; order = 0 computes them" compute(
            problem, CoaxialFormulation(order = 1))
    end

    # Paper insulations in oil touch each other: their bodies touch, and the cable computes.
    @testset "members whose bodies touch compute" begin
        offset = 2 * 28.0e-3 / sqrt(3)
        members = [Placement(offset * cos(θ),
                       offset * sin(θ)) => Stack(
                       Terminal(name, Region(name, Circle(19.5e-3), copper)),
                       Region(Symbol(name, :_insulation), Offset(8.5e-3), paper))
                   for (name, θ) in zip((:a, :b, :c), (0.0, 2π / 3, 4π / 3))]
        design = build(CableDesign, "touching",
            Enclosure(Region(:fill, Circle(62.0e-3), oil), Assembly(Tuple(members))),
            Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
            Region(:coating, Offset(5.0e-3), polymer))
        result = interior_result([design], [Placement(0, -1)],
            [(a = 1, b = 2, c = 3, pipe = 4)]; order = 20)
        @test all(isapprox(result.Y[:, :, k], transpose(result.Y[:, :, k]); rtol = 1.0e-14)
        for k in 1:4)
    end

    # Insulated cores in a common jacket, with no conductor around them, are cables of the
    # system. Each has its own earth return.
    @testset "members that no conductor surrounds are cables of their own" begin
        angles, names = (0, 2π / 3, 4π / 3), (:a, :b, :c)
        cores = [Stack(Terminal(name, Region(name, Circle(3.0e-3), copper)),
                     Region(Symbol(name, :_insulation), Offset(1.0e-3), polymer))
                 for name in names]
        design = build(CableDesign, "three_cores",
            Stack(
                Enclosure(Region(:fill, Circle(13.0e-3), air),
                    Assembly(Tuple(Placement(8.0e-3 * cos(φ), 8.0e-3 * sin(φ)) => core
                    for (φ, core) in zip(angles, cores)))),
                Region(:jacket, Offset(2.0e-3), polymer)))
        result = interior_result([design], [Placement(0, -1)], [(a = 1, b = 2, c = 3)])
        apart = interior_result(
            [build(CableDesign, string(name), core) for (name, core) in zip(names, cores)],
            [Placement(8.0e-3 * cos(φ), -1 + 8.0e-3 * sin(φ)) for φ in angles],
            [(; name => k) for (k, name) in enumerate(names)])
        @test size(result.trace.Zg) == (3, 3, 4)
        @test result.trace.Zg == apart.trace.Zg && result.trace.Pg == apart.trace.Pg
        @test result.Z == apart.Z && result.Y == apart.Y
    end

    # Construction 10 of the blueprints: sector cores with no conductor around them.
    @testset "the sector cores of a cable without metal compute as cables" begin
        sector = LineCableModels.Sector(
            0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))
        c = first(LineCableModels.centroid(sector))
        core = LineCableModels.Sector(
            0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(-c, 0))
        design = build(CableDesign, "three_sector_cores",
            Stack(
                Enclosure(Region(:fill, Circle(11.0e-3), air),
                    Assembly(Tuple(Placement(c * cos(φ),
                                       c * sin(φ),
                                       φ) => Stack(
                                       Terminal(name, Region(:core, core, copper)),
                                       Region(:insulation, Offset(1.0e-3), polymer))
                    for (φ, name) in zip((0, 2π / 3, 4π / 3),
                        (:phase_a, :phase_b, :phase_c))))),
                Region(:jacket, Offset(3.0e-3), polymer)))
        result = interior_result([design], [Placement(0, -1)],
            [(phase_a = 1, phase_b = 2, phase_c = 3)])
        @test size(result.trace.Zg) == (3, 3, 4)
        @test all(iszero, result.trace.Zin[1, 2:3, :]) &&
              all(iszero, result.trace.Pin[1, 2:3, :])
    end

    # A cable of two insulated wires in a sheath is a member of a pipe, beside a coaxial cable.
    # The sheath shields the wires completely from the potentials outside it. Measured: 4.4e-16
    # against the cable alone at order 20.
    @testset "a member with members of its own computes at every order" begin
        wires = [Placement(x,
                     0.0) => Stack(Terminal(name, Region(name, Circle(2.0e-3), copper)),
                     Region(Symbol(name, :_insulation), Offset(1.0e-3), polymer))
                 for (x, name) in ((3.5e-3, :w1), (-3.5e-3, :w2))]
        inner = Stack(
            Enclosure(Region(:inner_fill, Circle(8.0e-3), polymer),
                Assembly(Tuple(wires))),
            Terminal(:sheath, Region(:sheath, Offset(1.0e-3), copper)),
            Region(:inner_jacket, Offset(1.0e-3), polymer))
        b = Stack(Terminal(:b, Region(:b, Circle(6.0e-3), copper)),
            Region(:b_insulation, Offset(3.0e-3), polymer),
            Terminal(:b_screen, Region(:b_screen, Offset(1.0e-3), copper)),
            Region(:b_jacket, Offset(2.0e-3), polymer))
        nested = piped_design("nested",
            [Placement(-14.0e-3, 0.0) => inner, Placement(15.0e-3, 0.0) => b], air, steel,
            polymer)
        connections = (w1 = 1, w2 = 2, sheath = 3, b = 4, b_screen = 5, pipe = 6)
        @test_logs interior_result([nested], [Placement(0, -1)], [connections])
        result = @test_logs (:warn,
            "the coaxial formulation computes the impedances above order 0 of the interior of cable:1:sheath as shielded by that conductor from the interior around it: the field that crosses that conductor is left out") interior_result(
            [nested], [Placement(0, -1)], [connections]; order = 20)
        alone = interior_result([build(CableDesign, "alone", inner)], [Placement(0, -1)],
            [(w1 = 1, w2 = 2, sheath = 3)]; order = 20)
        for k in 1:4
            @test isapprox(looped(result.trace.Pin[:, :, k], 3, 1:2),
                looped(alone.trace.Pin[:, :, k], 3, 1:2); rtol = 1.0e-15)
        end
        problem = LineParametersProblem(
            build(CableSystem, [nested], [Placement(0, -1)]; connections = [connections],
                id = "nested", line_length = 1.0), homogeneous(rho = 100.0, eps_r = 10.0);
            frequencies = [50.0])
        input = LineCableModels.Engine.CoaxialInput(problem, CoaxialFormulation().options)
        buffers = LineCableModels.Engine.initialize_buffers(
            CoaxialFormulation(order = 1), problem, input)
        @test [interior.ranges for interior in buffers.interiors] ==
              [[1:3, 4:5], [1:1, 2:2], UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[]]
        @test input.assemblies == [6:6, 3:3, 1:1, 2:2, 4:5]
    end

    # A formula of order 0 alone serves an interior at order 0, and the formulation refuses
    # it above order 0.
    @testset "the members above order 0 need surfaces at every order" begin
        problem = Dict(REFERENCES)["three_cores_steel_pipe"]
        @test size(compute(problem, CoaxialFormulation(internal_impedance = :wedepohl1973)).Z) ==
              (4, 4, 101)
        @test_throws ArgumentError(
            "the formula :wedepohl1973 gives the surface impedances at order 0 alone; the coaxial formulation computes the members of an interior above order 0 with a formula that gives every order, or at order = 0") compute(
            problem, CoaxialFormulation(internal_impedance = :wedepohl1973, order = 1))
        @test_throws DomainError CoaxialFormulation(order = -1)
    end

    # On a concentric cable only the domains of the formulas read the option order.
    @testset "a formula of order 0 alone warns at a higher order" begin
        problem = first(REFERENCES)[2]
        @test_logs compute(problem, CoaxialFormulation(internal_impedance = :wedepohl1973))
        @test_logs (:warn,
            "the formula :wedepohl1973 holds for order in 0 .. 0; the computation has order in 1 .. 1") compute(
            problem, CoaxialFormulation(internal_impedance = :wedepohl1973, order = 1))
        @test_logs compute(problem, CoaxialFormulation(order = 1))
        @test_logs compute(problem, CoaxialFormulation(order = 20))
    end
end
