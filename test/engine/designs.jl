@testset "Coaxial input of designs" begin
    flatten = LineCableModels.flatten
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1.0e14, 2.3)
    semicon = Material(:semicon, 1000.0, 1000.0)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    air = material(kind = :insulator, rho = Inf)

    coaxial = coaxial_input

    @testset "a wire is its own tube" begin
        wire = coaxial(build(CableDesign, "insulated_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)),
            Region(:core_insulation, Offset(1.0e-3), polymer)))
        tube = only(wire.material)
        @test [only(wire.r_in), only(wire.r_ex)] ≈ [0.0, 0.0425]
        @test tube.rho ≈ copper.rho && tube.mu_r ≈ 1
        @test [only(wire.dielectric_r_in), only(wire.dielectric_r_ex)] ≈ [0.0425, 0.0435]
        @test wire.mu_r_ins == [1.0]
    end

    @testset "a thin annulus keeps its permeability" begin
        # Bounds measured on Julia 1.12.7, x86_64 Linux, with b = 30 mm.
        b = 0.03
        for (metal, ε) in ((copper, 1e-3), (copper, 1e-5), (copper, 1e-8), (steel, 1e-4))
            a = b * (1 - ε)
            design = build(CableDesign, "thin_tube",
                Region(:bore, Circle(a), polymer),
                Terminal(:core, Region(:tube, Annulus(a, b), metal)),
                Region(:jacket, Offset(1e-3), polymer))
            tube = only(coaxial(design).material)
            @test abs(tube.mu_r - metal.mu_r) <= 2eps(metal.mu_r)
        end
        a = b * (1 - 1e-5)
        design = build(CableDesign, "thin_tube",
            Region(:bore, Circle(a), polymer),
            Terminal(:core,
                Region(:tube, Annulus(a, ForwardDiff.Dual(b, 1.0)), copper)),
            Region(:jacket, Offset(1e-3), polymer))
        tube = only(coaxial(design).material)
        @test abs(only(ForwardDiff.partials(tube.mu_r))) < 1e-20
    end

    @testset "a wire inside a tube of its copper is a solid conductor" begin
        core = coaxial(build(CableDesign, "wire_and_tube",
            Terminal(:core,
                Stack(Region(:core_wire, Circle(1.0e-3), copper),
                    Region(:core_tube, Offset(1.0e-3), copper)))))
        tube = only(core.material)
        @test [only(core.r_in), only(core.r_ex)] ≈ [0.0, 2.0e-3]
        @test tube.rho ≈ copper.rho && tube.mu_r ≈ 1
    end

    # On Julia 1.12.7, x86_64 Linux, the residue outer^2 - sum(area, shapes) / π
    # was negative at 52 of these 400 radii. It was positive at 58 and zero at 290.
    @testset "a solid core off the origin has zero inner radius" begin
        for r in range(1e-3, 60e-3; length = 400)
            design = build(CableDesign, "solid_core",
                Terminal(:core,
                    Region(:core_metal, Circle(r, Placement(20e-3 + r, 0.0)), copper)),
                Region(:ins, Offset(8e-3), polymer))
            core = coaxial(design)
            tube = only(core.material)
            @test only(core.r_in) == 0
            @test tube.rho ≈ copper.rho && tube.mu_r ≈ 1
        end
    end

    @testset "a lay over a strand at the radius 0 keeps it straight" begin
        strand = Repeat(CircularPattern(1; r = 0),
            Region(:core_metal, Circle(0.0425), copper))
        laid = build(CableDesign, "laid_wire",
            Terminal(:core, Lay(Helix(LayRatio(12)), strand)),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        straight = build(CableDesign, "straight_wire", Terminal(:core, strand),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        @test only(LineCableModels.resolve(
            NamedTuple[], Placement(0, 0), 1.0, 0.0, first(laid.tree.items))).N == 0
        input, expected = coaxial(laid), coaxial(straight)
        @test only(input.material).rho == only(expected.material).rho
        @test input.mu_r_ins == expected.mu_r_ins == [1.0]
    end

    # On Julia 1.12.7, x86_64 Linux, gmd at τ = 1e-6 gives
    # the rectangular core permeability below.
    @testset "the coaxial cables of the references" begin
        for (strand, core) in (
            (Circle(1.5e-3),
                (2.3580693965940444e-8, 1.0588802203492307, 1.416791518882097)),
            (Rectangle(3.0e-3, 3.0e-3),
                (1.8631839650067386e-8, 1.0034147138418703, 1.4073302667101928)))
            courses = [Lay(Helix(LayRatio(12)),
                           Repeat(CircularPattern(n; r = 3.0e-3k),
                               Region(Symbol(:core_strands_, k), strand, copper)))
                       for (k, n) in enumerate((6, 12, 18))]
            armour = Lay(Helix(LayRatio(10)),
                Repeat(CircularPattern(40; r = 26.5e-3),
                    Region(:armour_wires, Circle(2.0e-3), steel)))
            input = coaxial(build(CableDesign, "coaxial",
                Terminal(:core,
                    Group(Region(:core_fill, Circle(10.5e-3), air),
                        Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
                Region(:core_semicon_inner, Offset(1.0e-3), semicon),
                Region(:core_insulation, Offset(8.0e-3), polymer),
                Region(:core_semicon_outer, Offset(1.0e-3), semicon),
                Terminal(:sheath, Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
                Region(:sheath_bedding, Offset(3.0e-3), polymer),
                Terminal(:armour,
                    Group(Region(:armour_fill, Offset(4.0e-3), air), armour)),
                Region(:armour_jacket, Offset(4.0e-3), polymer)))
            @test input.label == ["cable:1:core", "cable:1:sheath", "cable:1:armour"]
            @test input.r_in ≈ [0.0, 20.5e-3, 24.5e-3]
            @test input.r_ex ≈ [10.5e-3, 21.5e-3, 28.5e-3]
            @test input.dielectric_material ==
                  [air, semicon, polymer, semicon, polymer, air, polymer]
            @test input.dielectric_ranges == [1:4, 5:5, 6:7]
            @test input.dielectric_r_in ≈
                  [10.5e-3, 10.5e-3, 11.5e-3, 19.5e-3, 21.5e-3, 28.5e-3, 28.5e-3]
            @test input.dielectric_r_ex ≈
                  [10.5e-3, 11.5e-3, 19.5e-3, 20.5e-3, 24.5e-3, 28.5e-3, 32.5e-3]
            values = [(tube.rho, tube.mu_r, mu)
                      for (tube, mu) in zip(input.material, input.mu_r_ins)]
            expected = [core, (copper.rho, 1.0, 1.0),
                (2.4270115498150639e-7, 5.8715848278702119, 1.1305521523424971)]
            @test all(isapprox(value, target; rtol = 1.0e-12)
            for (row, targets) in zip(values, expected)
            for (value, target) in zip(row, targets))
            @test [tube.alpha for tube in input.material] ≈ [0.00393, 0.00393, 0.0045]
            @test input.assemblies == [1:3] && input.positions == [(0.0, 0.0)]
            @test isempty(only(input.interior_ranges))
        end
        # 37 round wires in parallel, the 36 laid ones longer by their overlength.
        wire = copper.rho / (π * 1.5e-3^2)
        resistance = wire / (1 + 36 / sqrt(1 + (π / 12)^2))
        @test resistance * π * 10.5e-3^2 ≈ 2.3580693965940444e-8
    end

    @testset "wires smaller than their boundary leave a layer of air" begin
        courses = [Lay(Helix(LayRatio(12)),
                       Repeat(CircularPattern(n; r = 3.0e-3k),
                           Region(Symbol(:core_strands_, k), Circle(1.4e-3), copper)))
                   for (k, n) in enumerate((6, 12, 18))]
        armour = Lay(Helix(LayRatio(10)),
            Repeat(CircularPattern(40; r = 26.5e-3),
                Region(:armour_wires, Circle(2.0e-3), steel)))
        input = coaxial(build(CableDesign, "coaxial_circular_strands",
            Terminal(:core,
                Group(Region(:core_fill, Circle(10.5e-3), air),
                    Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
            Region(:core_semicon_inner, Offset(1.0e-3), semicon),
            Region(:core_insulation, Offset(8.0e-3), polymer),
            Region(:core_semicon_outer, Offset(1.0e-3), semicon),
            Terminal(:sheath, Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
            Region(:sheath_bedding, Offset(3.0e-3), polymer),
            Terminal(:armour, Group(Region(:armour_fill, Offset(4.0e-3), air), armour)),
            Region(:armour_jacket, Offset(4.0e-3), polymer)))
        l = first(input.dielectric_ranges[1])
        @test [input.r_in[1], input.r_ex[1]] ≈ [0.0, 10.4e-3]
        @test [input.dielectric_r_in[l], input.dielectric_r_ex[l]] ≈ [10.4e-3, 10.5e-3]
        @test input.dielectric_material[l] == air
    end

    @testset "only the closing fill of two wire courses forms a layer" begin
        local polymer = Material(:insulator, 1.0e14, 2.5)
        local air = Material(:insulator, 1.0e14, 1.0)
        for eps_r in (1.0, 4.0)
            compound = Material(:insulator, 1.0e14, eps_r)
            input = coaxial(build(CableDesign, "two_courses",
                Terminal(:core,
                    Stack(Region(:centre, Circle(2.0e-3), copper),
                        wires(copper; shape = Circle(1.0e-3), n = 6,
                            lay = LayRatio(12), fill = air, tag = :inner),
                        wires(copper; shape = Circle(1.0e-3), n = 12,
                            lay = LayRatio(12), fill = compound, tag = :outer))),
                Region(:ins, Offset(2.0e-3), polymer),
                Terminal(:screen, Region(:screen, Offset(0.5e-3), copper)),
                Region(:jacket, Offset(2.0e-3), polymer)))
            @test [input.r_in input.r_ex] ≈ [0.0 0.006; 0.008 0.0085] rtol = 1e-12
            @test input.dielectric_ranges == [1:2, 3:3]
            @test [input.dielectric_r_in input.dielectric_r_ex] ≈
                  [0.006 0.006; 0.006 0.008; 0.0085 0.0105] rtol = 1e-12
            @test [layer.eps_r for layer in input.dielectric_material] == [eps_r, 2.5, 2.5]
        end
    end

    @testset "a fill under metal of its terminal stays inside the tube" begin
        local polymer = Material(:insulator, 1.0e14, 2.5)
        compound = Material(:insulator, 1.0e14, 4.0)
        input = coaxial(build(CableDesign, "covered",
            Terminal(:core,
                Stack(
                    Group(Region(:fill, Circle(3.0e-3), compound),
                        Stack(Region(:centre, Circle(1.0e-3), copper),
                            Lay(Helix(LayRatio(12)),
                                Repeat(CircularPattern(6; r = 2.0e-3),
                                    Region(:w, Circle(1.0e-3), copper))))),
                    Region(:tube, Offset(1.0e-3), copper))),
            Region(:ins, Offset(2.0e-3), polymer),
            Terminal(:screen, Region(:screen, Offset(0.5e-3), copper)),
            Region(:jacket, Offset(2.0e-3), polymer)))
        @test [input.r_in input.r_ex] ≈
              [0.0 0.004; 0.006 0.006500000000000001] rtol = 1e-12
        @test input.dielectric_ranges == [1:1, 2:2]
        @test [input.dielectric_r_in input.dielectric_r_ex] ≈
              [0.004 0.006; 0.006500000000000001 0.0085] rtol = 1e-12
        @test [layer.eps_r for layer in input.dielectric_material] == [2.5, 2.5]
    end

    # A wire of 5 mm radius lies 5 mm inside a square fill of 20 mm. The square moves
    # inward by 5 mm to the wire. The conductor is that square of 10 mm, which takes
    # Ametani's cylinder of its area and its perimeter. The fill beyond it is a layer of
    # 5 mm.
    @testset "a wire in a rectangular fill is the rectangle moved in to the wire" begin
        local polymer = Material(:insulator, 1.0e14, 2.5)
        input = coaxial(build(CableDesign, "square",
            Terminal(:core,
                Group(Region(:fill, Rectangle(20.0e-3, 20.0e-3), polymer),
                    Region(:wire, Circle(5.0e-3), copper))),
            Region(:ins, Offset(2.0e-3), polymer)))
        r_e = 40.0e-3 / 2π
        r_i = sqrt(r_e^2 - 100.0e-6 / π)
        @test [only(input.r_in), only(input.r_ex)] ≈ [r_i, r_e] rtol = 1e-12
        @test only(input.material).rho / 100.0e-6 ≈ copper.rho / (π * 5.0e-3^2) rtol = 1e-12
        @test [input.dielectric_r_in input.dielectric_r_ex] ≈
              [r_e r_e+5.0e-3; r_e+5.0e-3 r_e+7.0e-3] rtol = 1e-12
        @test [layer.eps_r for layer in input.dielectric_material] == [2.5, 2.5]
    end

    @testset "opposite lays cancel turns" begin
        relative = map((1, -1)) do dir
            courses = [Lay(Helix(LayRatio(12); dir = k == 2 ? dir : 1),
                           Repeat(CircularPattern(n; r = 3.0e-3k),
                               Region(Symbol(:core_strands_, k), Circle(1.5e-3), copper)))
                       for (k, n) in enumerate((6, 12))]
            core = coaxial(build(CableDesign, "core",
                Terminal(:core,
                    Group(Region(:core_fill, Circle(7.5e-3), air),
                        Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
                Region(:core_insulation, Offset(8.0e-3), polymer)))
            return only(core.mu_r_ins)
        end
        @test relative[2] < relative[1]
    end

    # A rectangular fill contains a straight wire and a course of wires in a lay about it.
    # A wire of the course has the overlength of the lay and makes one turn per lay length.
    # The helix of the flattened conductor has the lay length of the turns of its current.
    @testset "laid wires inside a rectangle keep their turns" begin
        a = 1.0e-3
        design = build(CableDesign, "bar",
            Terminal(:core,
                stranded(copper; shape = Circle(a), n = (1, 6), lay = LayRatio(14),
                    boundary = Rectangle(7.0e-3, 6.2e-3), fill = air)),
            Region(:insulation, Offset(1.0e-3), polymer))
        ℓ = hypot(1, π / 14)
        N = 6 / ℓ / (1 + 6 / ℓ) / (2 * 2a * 14)
        helix = first(flatten(design).items).part.helix
        @test helix.dir / pitch(helix, 1.0) ≈ N rtol = 1e-12
        input = coaxial(design)
        r_e, r_d = only(input.r_ex), only(input.r_ins_ex)
        u = log(r_d / r_e)
        @test only(input.mu_r_ins) ≈ 1 + 4π^2 * N^2 * r_d * r_e * sinh(u) / u rtol = 1e-12
    end

    @testset "designs without coaxial tubes" begin
        strands = Lay(Helix(LayRatio(12)),
            Repeat(CircularPattern(6; r = 1.0e-3), Region(:core_wires, Circle(0.5e-3), copper)))
        core = Stack(
            Terminal(:core,
                Group(Region(:core_fill, Circle(1.5e-3), air),
                    Stack(Region(:core_centre, Circle(0.5e-3), copper), strands))),
            Region(:core_insulation, Offset(1.0e-3), polymer),
            Terminal(:screen, Region(:core_screen, Offset(0.5e-3), copper)))
        cores = Repeat(CircularPattern(3; r = 3.5e-3), core)
        pipe = Region(:pipe_fill, Circle(6.5e-3), polymer)
        three = coaxial(build(CableDesign, "three_cores", Enclosure(pipe, cores)))
        @test three.assemblies == [7:6, 1:2, 3:4, 5:6]
        @test length(three.positions) == 4 && three.positions[1] == (0.0, 0.0)
        @test three.dielectric_material[three.interior_ranges[1]] == [polymer]
        @test_throws ArgumentError flatten(build(CableDesign, "bonded",
            Terminal(:phase, Enclosure(pipe, cores))))
        warm = Material(:conductor, 1.7241e-8, 1, 1, 25, 0.00393)
        @test_throws ArgumentError flatten(build(CableDesign, "two_temperatures",
            Terminal(:core,
                Stack(Region(:core_wire, Circle(1.0e-3), copper),
                    Region(:core_tube, Offset(1.0e-3), warm)))))
    end

    # A sector conductor takes Ametani's tube of its area and perimeter. The layer around it
    # has the thickness of its offset.
    @testset "a sector core is the tube of its area and perimeter" begin
        sector = LineCableModels.Sector(
            0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))
        input = coaxial(build(CableDesign, "sector_core",
            Terminal(:core, Region(:core, sector, copper)),
            Region(:insulation, Offset(1.0e-3), polymer)))
        r_i, r_e = only(input.r_in), only(input.r_ex)
        @test [r_i, r_e] ≈ [2.3105861599593243e-3, 5.652135443739411e-3] rtol = 1e-12
        @test π * (r_e^2 - r_i^2) ≈ 83.59095280005376e-6 rtol = 1e-12
        @test 2π * r_e ≈ 35.51341437429244e-3 rtol = 1e-12
        @test only(input.material).rho ≈ copper.rho rtol = 1e-12
        @test [only(input.dielectric_r_in), only(input.dielectric_r_ex)] ≈
              [r_e, r_e + 1.0e-3] rtol = 1e-12
    end

    # A screen that an offset declares around a sector core is the tube of its thickness. The
    # tube has the resistance and the geometric mean radius of the cross-section of the
    # screen. The values are those of the commit 3e15c5a1b.
    @testset "a screen around a sector core is the tube of its thickness" begin
        aluminum = Material(:conductor, 2.8264e-8, 1, 1, 20, 0.00429)
        sector = LineCableModels.Sector(
            0, 10.24e-3, deg2rad(119.0), 1.1e-3 / sind(60), 1.02e-3, Placement(0, 0))
        input = coaxial(build(CableDesign, "screened_sector",
            Terminal(:core, Region(:core, sector, aluminum)),
            Region(:insulation, Offset(1.1e-3), polymer),
            Terminal(:screen, Region(:screen, Offset(0.2e-3), copper)),
            Region(:jacket, Offset(2.0e-3), polymer)))
        @test input.r_in ≈ [2.584596240685967e-3, 6.948789605543471e-3] rtol = 1e-12
        @test input.r_ex ≈ [5.84878960554347e-3, 7.1487896055434695e-3] rtol = 1e-12
        @test [medium.rho for medium in input.material] ≈ [2.8264e-8, 1.7241e-8] rtol = 1e-12
        @test [medium.mu_r for medium in input.material] ≈
              [1.749196663688495, 5.401337489441598] rtol = 1e-12
    end

    # A fill whose container is the outline of a sector core at 2 mm is the layer of 2 mm
    # declared after the core. Measured on Julia 1.12.7, x86_64 Linux: the radii differ by
    # 1.6e-15 and the resistivity by 6.7e-16, relative.
    @testset "a fill around a sector core is the layer declared after it" begin
        sector = LineCableModels.Sector(
            0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))
        container = LineCableModels.resolve(sector, Offset(2.0e-3)).outer
        inside = coaxial(build(CableDesign, "inside",
            Terminal(:core,
                Group(Region(:fill, container, polymer), Region(:core, sector, copper))),
            Region(:insulation, Offset(1.0e-3), polymer)))
        after = coaxial(build(CableDesign, "after",
            Terminal(:core, Region(:core, sector, copper)),
            Region(:fill, Offset(2.0e-3), polymer),
            Region(:insulation, Offset(1.0e-3), polymer)))
        @test [inside.r_in inside.r_ex] ≈ [after.r_in after.r_ex] rtol = 1e-12
        @test [inside.dielectric_r_in inside.dielectric_r_ex] ≈
              [after.dielectric_r_in after.dielectric_r_ex] rtol = 1e-12
        @test only(inside.material).rho ≈ only(after.material).rho rtol = 1e-12
        @test only(inside.material).mu_r ≈ only(after.material).mu_r rtol = 1e-12
        @test inside.dielectric_material == after.dielectric_material
    end

    # The sector core of a 1 kV cable of three cores has 17 round wires on four arcs about
    # the center of the sector. The wires of the outer arc touch the back of the sector. The
    # conductor is then the sector, with the tube of the solid sector, and the fill beyond it
    # is a layer of no width. Wires 0.1 mm smaller at the same places end 0.1 mm inside the
    # sector. The tube ends 0.1 mm lower, and the fill beyond it is a layer of 0.1 mm.
    @testset "the conductor of a stranded sector ends at its wires" begin
        aluminum = Material(:conductor, 2.8264e-8, 1, 1, 20, 0.00429)
        pvc = Material(:insulator, Inf, 8.0)
        wedge = LineCableModels.Sector(
            0, 10.24e-3, deg2rad(119.0), 1.1e-3 / sind(60), 1.02e-3, Placement(0, 0))
        function strands(a)
            wire = Region(:wires, Circle(a), aluminum)
            arcs = [Repeat(CircularPattern(n; r, φ0 = -φ, span = 2φ), wire)
                    for (n, r, φ) in (
                (3, 4.99e-3, 0.594), (5, 7.09e-3, 0.731), (8, 9.19e-3, 0.8029))]
            return Group(Region(:core_fill, wedge, air),
                Stack(Repeat(CircularPattern(1; r = 2.89e-3), wire), arcs...))
        end
        r_i, r_e = 2.584596240685967e-3, 5.84878960554347e-3
        for (a, gap) in ((1.05e-3, 0.0), (0.95e-3, 0.1e-3))
            input = coaxial(build(CableDesign, "sector", Terminal(:core, strands(a)),
                Region(:insulation, Offset(1.1e-3), pvc)))
            @test [only(input.r_in), only(input.r_ex)] ≈ [r_i, r_e - gap] rtol = 1e-9
            @test only(input.material).rho /
                  (π * (only(input.r_ex)^2 - only(input.r_in)^2)) ≈
                  aluminum.rho / (17π * a^2) rtol = 1e-12
            @test [input.dielectric_r_in input.dielectric_r_ex] ≈
                  [r_e-gap r_e; r_e r_e+1.1e-3] rtol = 1e-9
            @test input.dielectric_material == [air, pvc]
        end
        # A cable has such a core for each phase, under a neutral of 30 wires. Each core
        # has the tube of the core alone and lies about the centroid of its sector. The
        # tubes of the cores overlap. Order 0 computes the cable, and a higher order
        # refuses it.
        r_wire, r_center = 0.79e-3, 14.36e-3
        member = Terminal(:core,
            Stack(strands(1.05e-3), Region(:insulation, Offset(1.1e-3), pvc)))
        cable = build(CableDesign, "three_sectors",
            Enclosure(Region(:filler, Circle(r_center - r_wire), pvc),
                cores(member; r = 0, names = (:a, :b, :c))),
            Terminal(:neutral,
                Group(Region(:neutral_fill, Offset(2r_wire), pvc),
                    Repeat(CircularPattern(30; r = Offset(0)),
                        Region(:neutral_wires, Circle(r_wire), copper)))),
            Region(:jacket, Offset(2.1e-3), pvc))
        input = coaxial(cable)
        @test [input.r_in[1:3] input.r_ex[1:3]] ≈ repeat([r_i r_e], 3) rtol = 1e-9
        c = first(LineCableModels.centroid(wedge))
        @test [collect(p) for p in input.positions[2:4]] ≈
              [[c * cos(φ), c * sin(φ)] for φ in (0, 2π / 3, 4π / 3)] rtol = 1e-9
        system = build(CableSystem, [cable], [Placement(0.0, -1.0)];
            connections = [(a = 1, b = 2, c = 3, neutral = 4)], id = "three_sectors",
            line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
            frequencies = [50.0, 1.0e4])
        @test size(compute(problem, CoaxialFormulation()).Z) == (4, 4, 2)
        @test_throws "the expansions about their axes diverge above order 0; order = 0 computes them" compute(
            problem, CoaxialFormulation(order = 1))
    end

    # A circle is its own tube, and a solid ellipse takes Ametani's cylinder. The impedance
    # and the admittance of the ellipse approach those of the circle in proportion to ε,
    # the relative difference of its semi-axes. Measured on Julia 1.12.7, x86_64 Linux, as
    # the largest entry from ε = 1e-1 down to 1e-10: at most 0.065 ε for Z and 0.343 ε
    # for Y.
    @testset "a solid ellipse computes as the circle that it approaches" begin
        function result(shape)
            design = build(CableDesign, "core",
                Terminal(:core, Region(:core, shape, copper)),
                Region(:insulation, Offset(4.0e-3), polymer),
                Terminal(:sheath, Region(:sheath, Offset(1.0e-3), copper)),
                Region(:jacket, Offset(2.0e-3), polymer))
            system = build(CableSystem, [design], [Placement(0.0, -1.0)];
                connections = [(core = 1, sheath = 2)], id = "core", line_length = 1.0)
            return compute(
                LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
                    frequencies = [50.0, 1.0e3, 1.0e5]),
                CoaxialFormulation())
        end
        a = 10.0e-3
        circle = result(Circle(a))
        for ε in (1.0e-2, 1.0e-6, 1.0e-10)
            ellipse = result(LineCableModels.Ellipse(a * (1 + ε), a))
            @test maximum(abs, ellipse.Z - circle.Z) <= 0.07ε * maximum(abs, circle.Z)
            @test maximum(abs, ellipse.Y - circle.Y) <= 0.35ε * maximum(abs, circle.Y)
        end
    end

    # A core has 37 laid wires of 1.4 mm in a round fill of 10.5 mm. A layer of air of
    # 0.1 mm lies around its tube. In an elliptical fill the ellipse moves inward to the
    # wires. The impedance and the admittance approach those of the round fill in proportion
    # to ε, the relative difference of the semi-axes. Measured on Julia 1.12.7, x86_64 Linux,
    # as the largest entry from ε = 1e-2 down to 1e-10: at most 0.048 ε for Z and 0.252 ε
    # for Y.
    @testset "wires in an elliptical fill compute as in the circle that it approaches" begin
        function result(container)
            courses = [Lay(Helix(LayRatio(12)),
                           Repeat(CircularPattern(n; r = 3.0e-3k),
                               Region(Symbol(:wires_, k), Circle(1.4e-3), copper)))
                       for (k, n) in enumerate((6, 12, 18))]
            design = build(CableDesign, "core",
                Terminal(:core,
                    Group(Region(:core_fill, container, air),
                        Stack(Region(:center, Circle(1.4e-3), copper), courses...))),
                Region(:insulation, Offset(8.0e-3), polymer),
                Terminal(:sheath, Region(:sheath, Offset(1.0e-3), copper)),
                Region(:jacket, Offset(3.0e-3), polymer))
            system = build(CableSystem, [design], [Placement(0.0, -1.0)];
                connections = [(core = 1, sheath = 2)], id = "core", line_length = 1.0)
            return compute(
                LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
                    frequencies = [50.0, 1.0e3, 1.0e5]),
                CoaxialFormulation())
        end
        b = 10.5e-3
        circle = result(Circle(b))
        for ε in (1.0e-2, 1.0e-6, 1.0e-10)
            ellipse = result(LineCableModels.Ellipse(b * (1 + ε), b))
            @test maximum(abs, ellipse.Z - circle.Z) <= 0.05ε * maximum(abs, circle.Z)
            @test maximum(abs, ellipse.Y - circle.Y) <= 0.26ε * maximum(abs, circle.Y)
        end
    end

    # An open tape and its fill tile one layer. Under a radial field the two carry their
    # currents in parallel, and the admittivity of the layer is the area-weighted sum of theirs.
    @testset "an open dielectric tape is one layer of the composite material" begin
        polyester = material(kind = :insulator, rho = 1.0e15, eps_r = 3.3, tan_delta = 2.0e-3)
        compound = material(kind = :insulator, rho = 1.0e13, eps_r = 2.2, tan_delta = 5.0e-4)
        input = coaxial(build(CableDesign, "open_tape",
            Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            tape(polyester; w = 10.0e-3, t = 0.2e-3, lay = LayRatio(10), fill = compound),
            Region(:insulation, Offset(1.0e-3), polymer)))
        @test length(input.dielectric_material) == 2
        layer = first(input.dielectric_material)
        f = 10.0e-3 * 0.2e-3 / (π * (10.2e-3^2 - 10.0e-3^2))
        @test layer isa Material{Float64, :insulator}
        @test [input.dielectric_r_in[1], input.dielectric_r_ex[1]] ≈ [10.0e-3, 10.2e-3] rtol = 1e-12
        @test 1 / layer.rho ≈ f / 1.0e15 + (1 - f) / 1.0e13 rtol = 1e-12
        @test layer.eps_r ≈ f * 3.3 + (1 - f) * 2.2 rtol = 1e-12
        @test layer.eps_r * layer.tan_delta ≈
              f * 3.3 * 2.0e-3 + (1 - f) * 2.2 * 5.0e-4 rtol = 1e-12
        for ω in 2π .* (50.0, 1.0e3, 1.0e6)
            @test 1 / layer.rho +
                  ω * LineCableModels.ε₀ * layer.eps_r *
                  (layer.tan_delta + im) ≈
                  f * (1 / 1.0e15 + ω * LineCableModels.ε₀ * 3.3 * (2.0e-3 + im)) +
                  (1 - f) * (1 / 1.0e13 + ω * LineCableModels.ε₀ * 2.2 * (5.0e-4 + im)) rtol = 1e-12
        end
    end

    # The same tape and fill give one layer in an interval and outside the filler of an
    # interior. The interior is that of a cable without a conductor, alone and under the fill
    # and the wall of a duct.
    @testset "a tiled tape has the same layer in an interval and an interior" begin
        polyester = material(kind = :insulator, rho = 1.0e15, eps_r = 3.3, tan_delta = 2.0e-3)
        compound = material(kind = :insulator, rho = 1.0e13, eps_r = 2.2, tan_delta = 5.0e-4)
        part = tape(polyester; w = 10.0e-3, t = 0.2e-3, lay = LayRatio(10), fill = compound)
        interval = coaxial(build(CableDesign, "interval_tape",
            Terminal(:core, Region(:core, Circle(7.0e-3), copper)),
            part, Region(:insulation, Offset(1.0e-3), polymer)))
        d = 2 * 3.0e-3 / sqrt(3)
        core = Stack(Terminal(:core, Region(:core, Circle(2.0e-3), copper)),
            Region(:insulation, Offset(1.0e-3), polymer))
        cable = Stack(
            Enclosure(Region(:cable_fill, Circle(7.0e-3), air),
                Assembly(Tuple(Placement(d * cos(θ), d * sin(θ)) => core
                for θ in (π / 2, 7π / 6, 11π / 6)))),
            part, Region(:cable_jacket, Offset(1.0e-3), polymer))
        interior = coaxial(build(CableDesign, "interior_tape", cable))
        hosted = coaxial(build(CableDesign, "hosted_tape", cable,
            Region(:duct_fill, Offset(3.8e-3), air),
            Terminal(:duct, Region(:duct, Offset(2.0e-3), steel)),
            Region(:coating, Offset(1.0e-3), polymer)))
        @test length(only(interval.dielectric_ranges)) == 2
        @test length(interior.interior_ranges[1]) == 3
        @test length(hosted.interior_ranges[1]) == 4
        l = first(only(interval.dielectric_ranges))
        expected = interval.dielectric_material[l]
        f = 10.0e-3 * 0.2e-3 / (π * (7.2e-3^2 - 7.0e-3^2))
        @test [interval.dielectric_r_in[l], interval.dielectric_r_ex[l]] ≈ [7.0e-3, 7.2e-3] rtol = 1e-12
        @test expected.eps_r ≈ f * 3.3 + (1 - f) * 2.2 rtol = 1e-12
        @test 1 / expected.rho ≈ f / 1.0e15 + (1 - f) / 1.0e13 rtol = 1e-12
        for input in (interior, hosted)
            k = input.interior_ranges[1][2]
            @test input.dielectric_r_in[k] == interval.dielectric_r_in[l]
            @test input.dielectric_r_ex[k] == interval.dielectric_r_ex[l]
            @test input.dielectric_material[k].eps_r == expected.eps_r
            @test input.dielectric_material[k].rho == expected.rho
        end
    end

    @testset "two dielectrics combine into the class of their merged material" begin
        values = (
            1.0e13, 2.4, polymer.mu_r, polymer.T0, polymer.alpha, polymer.rho_thermal,
            polymer.theta_max, 1.0e-4, polymer.sigma_solar)
        @test Material(polymer, air, values...) === Material(Val(:insulator), values...)
        @test Material(polymer, semicon, values...) === Material(Val(:semicon), values...)
        @test Material(semicon, semicon, values...) === Material(Val(:semicon), values...)
        @test_throws MethodError Material(copper, polymer, values...)
    end

    # Construction 10. The cores lie at their centroids. The root has no conductor, so its fill
    # and its common jacket form its interior.
    @testset "the cores of a three-core sector cable sit at their centroids" begin
        sector = LineCableModels.Sector(
            0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))
        c = first(LineCableModels.centroid(sector))
        core = LineCableModels.Sector(
            0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(-c, 0))
        angles = (0, 2π / 3, 4π / 3)
        input = coaxial(build(CableDesign, "three_sector_cores",
            Stack(
                Enclosure(Region(:fill, Circle(11.0e-3), air),
                    Assembly(Tuple(Placement(c * cos(φ),
                                       c * sin(φ),
                                       φ) => Stack(
                                       Terminal(name, Region(:core, core, copper)),
                                       Region(:insulation, Offset(1.0e-3), polymer))
                    for (φ, name) in zip(angles, (:phase_a, :phase_b, :phase_c))))),
                Region(:jacket, Offset(3.0e-3), polymer))))
        @test all(isapprox([r_i, r_e],
                      [2.3105861599593243e-3, 5.652135443739411e-3]; rtol = 1e-12)
        for (r_i, r_e) in zip(input.r_in, input.r_ex))
        @test [collect(p) for p in input.positions] ≈
              vcat([[0.0, 0.0]], [[c * cos(φ), c * sin(φ)] for φ in angles]) rtol = 1e-12
        @test isempty(input.assemblies[1])
        @test input.assemblies[2:4] == [1:1, 2:2, 3:3]
        @test input.dielectric_ranges == [1:1, 2:2, 3:3]
        @test input.interior_ranges[1] == 4:5 &&
              all(isempty, input.interior_ranges[2:4])
        @test input.dielectric_material == [polymer, polymer, polymer, air, polymer]
        @test all(isapprox(input.dielectric_r_ex[l], input.dielectric_r_in[l] + 1.0e-3;
                      rtol = 1e-12) for l in 1:3)
    end

    # Construction 15: three coaxials in a trefoil, each on the wall of a steel pipe. The
    # members' conductors precede the pipe's, and the fill is the interior of the root.
    @testset "three coaxials in a pipe flatten to their rows" begin
        d = 2 * 22.0e-3 / sqrt(3)
        θs = (π / 2, 7π / 6, 11π / 6)
        coax = Stack(Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            Region(:insulation, Offset(8.0e-3), polymer),
            Terminal(:screen, Region(:screen, Offset(1.0e-3), copper)),
            Region(:jacket, Offset(3.0e-3), polymer))
        design = build(CableDesign, "pipe",
            Stack(
                Enclosure(Region(:fill, Circle(22.0e-3 * (1 + 2 / sqrt(3))), air),
                    Assembly(Tuple(Placement(d * cos(θ), d * sin(θ)) => coax for θ in θs))),
                Terminal(:pipe, Region(:pipe, Offset(8.0e-3), steel)),
                Region(:coating, Offset(3.0e-3), polymer)))
        input = coaxial(design)
        R = 0.04740341184434353
        @test [[r_i, r_e] for (r_i, r_e) in zip(input.r_in, input.r_ex)] ≈
              [[0, 0.010], [0.018, 0.019], [0, 0.010], [0.018, 0.019], [0, 0.010],
            [0.018, 0.019], [R, R + 0.008]] rtol = 1e-12
        @test [[a, b] for (a, b) in zip(input.dielectric_r_in, input.dielectric_r_ex)] ≈
              [[0.010, 0.018], [0.019, 0.022], [0.010, 0.018], [0.019, 0.022],
            [0.010, 0.018], [0.019, 0.022], [0, R], [R + 0.008, R + 0.011]] rtol = 1e-12
        @test input.dielectric_material ==
              [polymer, polymer, polymer, polymer, polymer, polymer, air, polymer]
        @test input.dielectric_ranges == [1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 8:8]
        @test input.assemblies == [7:7, 1:2, 3:4, 5:6]
        @test input.interior_ranges == [7:7, 1:0, 3:2, 5:4]
        @test [collect(p) for p in input.positions] ≈
              vcat([[0.0, 0.0]], [[d * cos(θ), d * sin(θ)] for θ in θs]) rtol = 1e-12
        @test input.mu_r_ins ≈ ones(7) rtol = 1e-12
        system = build(CableSystem, [design], [Placement(0, -1)];
            connections = [(core = 1, screen = 0, pipe = 0)], id = "pipe", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
            frequencies = [50.0])
        # The jackets touch each other and the pipe, in air. Measured against order 80, as the
        # largest entry, Zin differs by 4.2e-13 at order 20 and Pin by 5.6e-8 at order 40. The
        # asymmetry of Pin is 4.9e-16, and Zin has none.
        results = Dict(order => compute(problem, CoaxialFormulation(; order); trace = true)
        for order in (20, 40, 80))
        (; Zin, Pin) = results[80].trace
        @test maximum(abs, results[20].trace.Zin - Zin) <= 1.0e-12 * maximum(abs, Zin)
        @test maximum(abs, results[40].trace.Pin - Pin) <= 1.0e-7 * maximum(abs, Pin)
        @test maximum(abs, Pin[:, :, 1] - transpose(Pin[:, :, 1])) <=
              1.0e-15 * maximum(abs, Pin)
        @test maximum(abs, Zin[:, :, 1] - transpose(Zin[:, :, 1])) <=
              1.0e-15 * maximum(abs, Zin)
    end

    # Construction 16: three round cores in a steel pipe, the geometry of the reference
    # case `three_cores_steel_pipe`. Each member is at the origin, turned, with its core's
    # circle off its origin. The member lies about the circle's center, the centroid of its
    # terminal. Cores at their own centers give the same input, equal to the bit on
    # Julia 1.12.7.
    @testset "a member lies about the centroid of its first terminal" begin
        offset = 2 * 28.0e-3 / sqrt(3)
        angles = (0.0, 2π / 3, 4π / 3)
        names = (:a, :b, :c)
        off = Circle(19.5e-3, Placement(offset, 0.0))
        turned = Tuple(Placement(0.0,
                           0.0,
                           θ) => Stack(
                           Terminal(name, Region(name, off, copper)),
                           Region(Symbol(name, :_insulation), Offset(8.5e-3), polymer))
        for (name, θ) in zip(names, angles))
        placed = Tuple(Placement(offset * cos(θ),
                           offset * sin(θ),
                           θ) => Stack(
                           Terminal(name, Region(name, Circle(19.5e-3), copper)),
                           Region(Symbol(name, :_insulation), Offset(8.5e-3), polymer))
        for (name, θ) in zip(names, angles))
        pipe,
        own = [coaxial(build(CableDesign, "pipe_type",
                   Enclosure(Region(:filler, Circle(62.0e-3), polymer), Assembly(members)),
                   Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
                   Region(:jacket, Offset(5.0e-3), polymer)))
               for members in (turned, placed)]
        centers = vcat([[0.0, 0.0]], [[offset * cos(θ), offset * sin(θ)] for θ in angles])
        @test [collect(p) for p in pipe.positions] ≈ centers rtol = 1e-12
        @test [pipe.r_in pipe.r_ex] ≈
              [0 19.5e-3; 0 19.5e-3; 0 19.5e-3; 62.0e-3 65.0e-3] rtol = 1e-12
        @test [pipe.dielectric_r_in pipe.dielectric_r_ex] ≈
              [19.5e-3 28.0e-3; 19.5e-3 28.0e-3; 19.5e-3 28.0e-3; 0 62.0e-3;
               65.0e-3 70.0e-3] rtol = 1e-12
        @test pipe.assemblies == [4:4, 1:1, 2:2, 3:3]
        @test pipe.interior_ranges == [4:4, 1:0, 2:1, 3:2]
        @test pipe.dielectric_ranges == [1:1, 2:2, 3:3, 5:5]
        @test pipe.members == [[2, 3, 4], Int[], Int[], Int[]]
        @test own.positions == pipe.positions
        @test [[own.r_in[k], own.r_ex[k], tube.rho, tube.mu_r]
               for (k, tube) in enumerate(own.material)] ==
              [[pipe.r_in[k], pipe.r_ex[k], tube.rho, tube.mu_r]
               for (k, tube) in enumerate(pipe.material)]
        @test [own.dielectric_r_in own.dielectric_r_ex] ==
              [pipe.dielectric_r_in pipe.dielectric_r_ex]
        @test own.mu_r_ins == pipe.mu_r_ins
        @test (own.assemblies, own.interior_ranges, own.dielectric_ranges, own.members) ==
              (pipe.assemblies, pipe.interior_ranges, pipe.dielectric_ranges, pipe.members)
    end

    # Construction 17: three wires around a central rod in a fill, under a binder tape. The
    # rod and the fill form the filler, one layer from the center to the tape. Its material
    # weighs the two by their areas, π mm² and 13.3125π mm², equal to the bit on Julia 1.12.7.
    # The tape is the annulus outside the filler.
    @testset "an interior is its filler and the annuli outside it" begin
        wire = Stack(Terminal(:core, Region(:core, Circle(0.5e-3), copper)),
            Region(:insulation, Offset(0.25e-3), polymer))
        input = coaxial(build(
            CableDesign, "rod", Region(:rod, Circle(1.0e-3), polymer),
            Enclosure(Region(:fill, Offset(3.0e-3), air),
                Assembly(Tuple(Placement(2.5e-3 * cos(θ), 2.5e-3 * sin(θ)) => wire
                for θ in (0.0, 2π / 3, 4π / 3)))),
            Region(:binder, Offset(0.2e-3), polymer)))
        span = input.interior_ranges[1]
        filler, binder = input.dielectric_material[span]
        @test [input.dielectric_r_in[span] input.dielectric_r_ex[span]] ≈
              [0 4.0e-3; 4.0e-3 4.2e-3] rtol = 1e-12
        @test filler.eps_r ≈ (polymer.eps_r + 13.3125 * air.eps_r) / 14.3125 rtol = 1e-14
        @test filler.rho ≈ 14.3125 * polymer.rho rtol = 1e-14
        @test binder == polymer
    end

    # Construction 18: a three-core cable without a metallic layer inside a steel duct. The
    # cable is a member of the duct's assembly, off the duct's axis and then on it. It has no
    # terminal, so its cores are members of the duct's assembly at their own centers. Its fill
    # and jacket join the duct's filler. The filler weighs the air, 102π mm², and the jacket,
    # 15π mm², equal to the bit on Julia 1.12.7.
    #
    # The same cable on the duct's axis, declared as the layers of one stack, has its fill as
    # the filler. Its jacket and the duct's fill are then annuli outside the filler.
    @testset "a sub-cable with no conductor of its own gives its members and regions to its container" begin
        d = 2 * 3.0e-3 / sqrt(3)
        angles = (π / 2, 7π / 6, 11π / 6)
        core = Stack(Terminal(:core, Region(:core, Circle(2.0e-3), copper)),
            Region(:insulation, Offset(1.0e-3), polymer))
        cable = Stack(
            Enclosure(Region(:cable_fill, Circle(7.0e-3), air),
                Assembly(Tuple(Placement(d * cos(θ), d * sin(θ)) => core for θ in angles))),
            Region(:cable_jacket, Offset(1.0e-3), polymer))
        for x in (2.0e-3, 0.0)
            input = coaxial(build(CableDesign, "duct",
                Enclosure(Region(:duct_fill, Circle(12.0e-3), air),
                    Assembly((Placement(x, 0.0) => cable,))),
                Terminal(:duct, Region(:duct, Offset(2.0e-3), steel)),
                Region(:coating, Offset(1.0e-3), polymer)))
            @test input.members == [[3, 4, 5], Int[], Int[], Int[], Int[]]
            @test [collect(p) for p in input.positions[3:5]] ≈
                  [[x + d * cos(θ), d * sin(θ)] for θ in angles] rtol = 1e-12
            @test isempty(input.assemblies[2]) && isempty(input.interior_ranges[2])
            f = only(input.interior_ranges[1])
            @test [input.dielectric_r_in[f], input.dielectric_r_ex[f]] ≈ [0, 12.0e-3] rtol = 1e-12
            @test input.dielectric_material[f].eps_r ≈
                  (102 * air.eps_r + 15 * polymer.eps_r) / 117 rtol = 1e-14
        end
        centered = coaxial(build(CableDesign, "duct", cable,
            Region(:duct_fill, Offset(4.0e-3), air),
            Terminal(:duct, Region(:duct, Offset(2.0e-3), steel)),
            Region(:coating, Offset(1.0e-3), polymer)))
        span = centered.interior_ranges[1]
        @test [centered.dielectric_r_in[span] centered.dielectric_r_ex[span]] ≈
              [0 7.0e-3; 7.0e-3 8.0e-3; 8.0e-3 12.0e-3] rtol = 1e-12
        @test [layer.eps_r for layer in centered.dielectric_material[span]] ≈
              [air.eps_r, polymer.eps_r, air.eps_r] rtol = 1e-14
    end

    # A sub-cable before a core gives the pipe's members in conductor order.
    # Measured: the pipe lists rows 3, 4, and 5, in that order.
    @testset "a sub-cable before a core gives its container members in conductor order" begin
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
        members = [Placement(-12.0e-3, 0.0) => cable,
            Placement(14.0e-3, 0.0) => core]
        input = coaxial(build(CableDesign, "pipe",
            Enclosure(Region(:fill, Circle(30.0e-3), air), Assembly(Tuple(members))),
            Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
            Region(:coating, Offset(2.0e-3), polymer)))
        @test input.members == [[3, 4, 5], Int[], Int[], Int[], Int[]]
    end

    # Construction 19: a terminal of two courses of 6 and 12 wires at the lay ratio 12, around
    # a polymer rod. The rod lies outside the terminal, in a hole of the fill around them. The
    # fill's container makes the block an annulus, and its tube starts at the axis. The tube has
    # the resistance of the 18 wires, each longer by its overlength. It is an equivalent tube,
    # as Ametani's is.
    @testset "a fill around a rod outside its terminal gives a tube from the axis" begin
        courses = [Lay(Helix(LayRatio(12)),
                       Repeat(CircularPattern(n; r = 3.0e-3k),
                           Region(Symbol(:wires_, k), Circle(1.5e-3), copper)))
                   for (k, n) in enumerate((6, 12))]
        input = coaxial(build(CableDesign, "rod",
            Group(Region(:core_fill, Circle(7.5e-3), air),
                Stack(Region(:rod, Circle(1.5e-3), polymer),
                    Terminal(:core, Stack(courses...)))),
            Region(:insulation, Offset(3.0e-3), polymer)))
        @test [only(input.r_in), only(input.r_ex)] ≈ [0, 7.5e-3] rtol = 1e-12
        resistance = copper.rho * sqrt(1 + (π / 12)^2) / (18π * 1.5e-3^2)
        @test only(input.material).rho / (π * 7.5e-3^2) ≈ resistance rtol = 1e-12
    end

    # Construction 21: one insulated core in a steel pipe, on its axis and 6 mm off it. The
    # pipe's fill is the filler of the root, from its axis to the pipe.
    @testset "a single core in a pipe has the pipe's fill as its filler" begin
        core = Stack(Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            Region(:insulation, Offset(2.0e-3), polymer))
        for offset in (0.0, 6.0e-3)
            input = coaxial(build(CableDesign, "pipe",
                Enclosure(Region(:fill, Circle(30.0e-3), air),
                    Assembly((Placement(offset, 0.0) => core,))),
                Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
                Region(:coating, Offset(2.0e-3), polymer)))
            @test [input.r_in input.r_ex] ≈ [0.0 10.0e-3; 30.0e-3 33.0e-3] rtol = 1e-12
            f = only(input.interior_ranges[1])
            @test [input.dielectric_r_in[f], input.dielectric_r_ex[f]] ≈ [0, 30.0e-3] rtol = 1e-12
            @test input.dielectric_material[f] == air
            @test collect(input.positions[2]) ≈ [offset, 0.0] rtol = 1e-12
        end
    end

    # Construction 20: one coaxial in a steel pipe, on its axis and 5 mm off it. The filler is
    # the pipe's fill, from the axis to the pipe.
    @testset "the filler of a pipe reaches from its axis to the pipe" begin
        coax = Stack(Terminal(:core, Region(:core, Circle(10.0e-3), copper)),
            Region(:insulation, Offset(8.0e-3), polymer),
            Terminal(:screen, Region(:screen, Offset(1.0e-3), copper)),
            Region(:jacket, Offset(3.0e-3), polymer))
        for offset in (0.0, 5.0e-3)
            enclosure = Enclosure(Region(:fill, Circle(30.0e-3), air),
                Assembly((Placement(offset, 0.0) => coax,)))
            input = coaxial(build(CableDesign, "pipe",
                Stack(enclosure, Terminal(:pipe, Region(:pipe, Offset(8.0e-3), steel)),
                    Region(:coating, Offset(3.0e-3), polymer))))
            f = only(input.interior_ranges[1])
            @test [input.dielectric_r_in[f], input.dielectric_r_ex[f]] ≈ [0, 30.0e-3] rtol = 1e-12
        end
    end

    # The cores of these cables lie in one filler and one jacket, with no metallic layer
    # around them. The flexible cable has three phase cores that touch each other and three
    # earth cores in the gaps between them. The second cable has six cores around a pilot
    # wire. The control cable has 8 and 16 cores in two layers.
    @testset "the cores of one jacket compute as cables of their own at their places" begin
        pvc = Material(:insulator, 1.0e13, 4.0)
        phase = Stack(core(copper; r = 1.5e-3), insulation(pvc; t = 0.8e-3))
        earth = Stack(core(copper; r = 0.4e-3), insulation(pvc; t = 0.3e-3))
        signal = Stack(core(copper; r = 0.7e-3), insulation(pvc; t = 0.5e-3))
        r, g = 2.3e-3 / sinpi(1 / 3), 3.6e-3
        flexible = [Placement(d * cospi(θ), d * sinpi(θ)) => Terminal(name, part)
                    for (name, part, d, θ) in (
            (:a, phase, r, 0), (:b, phase, r, 2 / 3), (:c, phase, r, 4 / 3),
            (:e1, earth, g, 1 / 3), (:e2, earth, g, 1), (:e3, earth, g, 5 / 3))]
        around = cores(Terminal(:core, signal); r = 2.4e-3,
            names = (:c1, :c2, :c3, :c4, :c5, :c6))
        pilot = [Placement(0.0, 0.0) => Terminal(:pilot, signal); around.members...]
        layers = [cores(Terminal(:core, signal); r = 1.2e-3 / sinpi(1 / n),
                      names = Tuple(Symbol(prefix, k) for k in 1:n))
                  for (prefix, n) in ((:inner, 8), (:outer, 16))]
        control = [first(layers).members..., last(layers).members...]
        for (members, container) in ((flexible, 5.5e-3), (pilot, 4.0e-3), (control, 8.0e-3))
            names = [member.name for (place, member) in members]
            cable = build(CableDesign, "cores",
                Enclosure(Region(:filler, Circle(container), pvc), Assembly(Tuple(members))),
                jacket(pvc; t = 1.5e-3))
            @test [row.name for row in cable.terminals] == names
            at = Placement(0.0, -1.0)
            whole = build(CableSystem, [cable], [at];
                connections = [NamedTuple{Tuple(names)}(Tuple(eachindex(names)))],
                id = "whole", line_length = 1.0)
            apart = build(CableSystem,
                [build(CableDesign, String(name), member.second)
                 for (name, member) in zip(names, members)],
                [at * place for (place, member) in members];
                connections = [NamedTuple{(name,)}((k,)) for (k, name) in enumerate(names)],
                id = "apart", line_length = 1.0)
            results = map((whole, apart)) do system
                return compute(
                    LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
                        frequencies = [50.0, 1.0e4]),
                    CoaxialFormulation())
            end
            @test size(first(results).Z) == (length(names), length(names), 2)
            @test first(results).Z ≈ last(results).Z rtol = 1e-12
            @test first(results).Y ≈ last(results).Y rtol = 1e-12
        end
    end

    # A rectangular duct of 60 mm by 30 mm contains three insulated cores side by side. With
    # a wall of polymer, or with none, no conductor surrounds the cores, and they compute as
    # cables of their own at their places.
    #
    # A wall of aluminum makes the duct a pipe. Its filler ends at the outer radius of
    # Ametani's cylinder of the rectangle, which is 180 mm / 2π. The tube of the wall begins
    # there and has the thickness of the wall. That tube is Ametani's cylinder of the wall,
    # with its area and its outer perimeter. The pipe computes as the round pipe of that
    # tube, with the resistance and the geometric mean radius of the wall. Each pair of
    # results is equal to the bit on Julia 1.12.7.
    @testset "a rectangular duct is the pipe of its perimeter" begin
        aluminum = Material(:conductor, 2.8264e-8, 1, 1, 20, 0.00429)
        part = Stack(core(copper; r = 5.0e-3), insulation(polymer; t = 3.0e-3))
        members = [Placement(x, 0.0) => Terminal(name, part)
                   for (name, x) in ((:a, -18.0e-3), (:b, 0.0), (:c, 18.0e-3))]
        contents = Assembly(Tuple(members))
        w, h, t = 60.0e-3, 30.0e-3, 3.0e-3
        ducts = [build(CableDesign, "duct", duct(air; shape = Rectangle(w, h), contents),
                     layers...) for layers in ((jacket(polymer; t, tag = :wall),), ())]
        trunk = build(CableDesign, "trunk", duct(air; shape = Rectangle(w, h), contents),
            Terminal(:trunk, sheath(aluminum; t, tag = :trunk)), jacket(polymer; t = 2.0e-3))
        input = coaxial(trunk)
        c = (w + h) / π
        f = only(input.interior_ranges[1])
        @test [input.dielectric_r_in[f], input.dielectric_r_ex[f]] ≈ [0, c] rtol = 1e-12
        @test [input.r_in[4], input.r_ex[4]] ≈ [c, c + t] rtol = 1e-12
        @test [input.dielectric_r_in[end], input.dielectric_r_ex[end]] ≈
              [c + t, c + t + 2.0e-3] rtol = 1e-12
        wall = only(placed for placed in trunk.regions if placed.source.tag == :trunk)
        @test π * (input.r_ex[4]^2 - input.r_in[4]^2) ≈ area(wall.shape) rtol = 1e-12
        @test 2π * input.r_ex[4] ≈ perimeter(boundary(wall.shape)) rtol = 1e-12
        tube = build(CableDesign, "tube", duct(air; shape = Circle(c), contents),
            Terminal(:trunk,
                Region(:trunk, Annulus(c, c + t), aluminum.rho / area(wall.shape),
                    LineCableModels.DataModel.gmr(wall.shape, aluminum.mu_r),
                    aluminum.T0, aluminum.alpha)),
            jacket(polymer; t = 2.0e-3))
        at = Placement(0.0, -1.0)
        apart = build(CableSystem,
            [build(CableDesign, String(member.name), member) for (place, member) in members],
            [at * place for (place, member) in members];
            connections = [(a = 1,), (b = 2,), (c = 3,)], id = "apart", line_length = 1.0)
        systems = [apart;
                   [build(CableSystem, [design], [at];
                        connections = [NamedTuple{names}(Tuple(eachindex(names)))],
                        id = design.id, line_length = 1.0)
                    for (design, names) in (
                       (first(ducts), (:a, :b, :c)), (last(ducts), (:a, :b, :c)),
                       (trunk, (:a, :b, :c, :trunk)), (tube, (:a, :b, :c, :trunk)))]]
        cables, walled,
        bare,
        piped,
        expected = map(systems) do system
            return compute(
                LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
                    frequencies = [50.0, 1.0e4]),
                CoaxialFormulation())
        end
        for result in (walled, bare)
            @test result.Z ≈ cables.Z rtol = 1e-12
            @test result.Y ≈ cables.Y rtol = 1e-12
        end
        @test size(piped.Z) == (4, 4, 2)
        @test piped.Z ≈ expected.Z rtol = 1e-12
        @test piped.Y ≈ expected.Y rtol = 1e-12
    end

    # A tube of 1 mm lies on a rectangular filler of 10 mm by 5 mm. No conductor and no
    # layer lies before it, and it takes Ametani's cylinder of its own cross-section. The
    # inner radius of the cylinder is the perimeter of the filler over 2π. A layer that an
    # offset declares on such a filler has no radius to begin at, and the input refuses it.
    #
    # A tube of 2 mm on an elliptical filler approaches the tube on the round filler in
    # proportion to ε, the relative difference of the semi-axes. Measured on Julia 1.12.7,
    # x86_64 Linux, as the largest entry from ε = 1e-1 down to 1e-10: at most 0.044 ε for Z
    # and 0.218 ε for Y.
    @testset "a tube on a filler that is not round is the cylinder of its cross-section" begin
        bar = filler(polymer; shape = Rectangle(10.0e-3, 5.0e-3))
        tube = Terminal(:tube, sheath(copper; t = 1.0e-3))
        design = build(CableDesign, "hollow_bar", bar, tube, jacket(polymer; t = 2.0e-3))
        input = coaxial(design)
        r = 30.0e-3 / 2π
        @test [only(input.r_in), only(input.r_ex)] ≈ [r, r + 1.0e-3] rtol = 1e-12
        @test only(input.material).rho / (π * ((r + 1.0e-3)^2 - r^2)) ≈
              copper.rho / area(design.regions[2].shape) rtol = 1e-12
        @test [only(input.dielectric_r_in), only(input.dielectric_r_ex)] ≈
              [r + 1.0e-3, r + 3.0e-3] rtol = 1e-12
        lined = build(CableDesign, "lined_bar", bar,
            insulation(polymer; t = 0.5e-3, tag = :liner), tube)
        @test_throws "the region :liner declares an Offset from the outline of a filler that is not round" coaxial(lined)
        function result(shape)
            hollow = build(CableDesign, "hollow", filler(polymer; shape),
                Terminal(:core, solid(copper; shape = Offset(2.0e-3))),
                insulation(polymer; t = 4.0e-3),
                Terminal(:sheath, sheath(copper; t = 1.0e-3)), jacket(polymer; t = 2.0e-3))
            system = build(CableSystem, [hollow], [Placement(0.0, -1.0)];
                connections = [(core = 1, sheath = 2)], id = "hollow", line_length = 1.0)
            return compute(
                LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
                    frequencies = [50.0, 1.0e3, 1.0e5]),
                CoaxialFormulation())
        end
        a = 5.0e-3
        circle = result(Circle(a))
        for ε in (1.0e-2, 1.0e-6, 1.0e-10)
            ellipse = result(LineCableModels.Ellipse(a * (1 + ε), a))
            @test maximum(abs, ellipse.Z - circle.Z) <= 0.05ε * maximum(abs, circle.Z)
            @test maximum(abs, ellipse.Y - circle.Y) <= 0.22ε * maximum(abs, circle.Y)
        end
    end
end
